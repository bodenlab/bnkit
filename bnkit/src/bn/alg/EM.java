/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package bn.alg;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.JDF;
import bn.JPT;
import bn.Variable;
import bn.alg.CGVarElim.CGVarElimRuntimeException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.List;

/**
 * Implementation of Expectation-Maximization for learning in Bayesian network.
 * Note that this class is still under development and has undoubtedly bugs.
 * It is loosely based on the superseded bioalg: ml.bayes implementation that I (MB)
 * wrote years ago.
 * @author mikael
 */
public class EM extends LearningAlg {

    private Inference inf;

    /**
     * @param bn Bayesian network
     * @param inference the inference algorithm
     */
    public EM(BNet bn, Inference inference) {
        super(bn);
        this.inf = inference;
    }

    /**
     * @param bn
     */
    public EM(BNet bn) {
        this(bn, new CGVarElim());
    }

    /**
     * EM convergence is true (and training stops) when the improvement is
     * smaller than this value. Note that smaller value leads to longer
     * convergence times but closer to optimal solutions.
     */
    public static double EM_CONVERGENCE_CRITERION = 0.00001;

    /**
     * EM maximum number of rounds (iterations). Training stops after the
     * specified number of rounds independently of convergence.
     */
    public static int EM_MAX_ROUNDS = 1000;
    
    /**
     * EM status print-outs
     */
    public static boolean EM_PRINT_STATUS = true;

    /**
     * EM option: currently two different approaches to determine expectations.
     * 1. one node at a time = one query per node
     * 2. all nodes at a time = one big query
     */
    public static int EM_OPTION = 1;

    /**
     * Train the BN using EM.
     *
     * @see bn.alg.LearningAlg#train(java.lang.Object[][], bn.Variable[], long)
     * @param values the values of the variables [row][variable], if a value is
     * null, it means "unspecified"
     * @param vars the variables that correspond to the values
     * @param seed the seed used to randomly select initial parameters (same
     * seed, same data --> same result)
     */
    @Override
    public void train(Object[][] values, Variable[] vars, long seed) {

        int nSample = values.length; // this is how training samples we have
        // we only need to initialize the relevant nodes.
        //	setRandom(seed); // init network CPTs and CDTs

        // below we create the list of nodes that will be updated during training
        // note that only named nodes (variables) will be included, and their parents or children (the latter as latent variables)
        Set<Variable> updateVars = new HashSet<>();
        for (int i = 0; i < vars.length; i++) {
            updateVars.add(vars[i]);
        }
        Map<BNode, Object[]> update = new HashMap<>(); 			// the set of nodes that may need to be updated:
        for (BNode node : bn.getNodes()) {      		// check all nodes
            Variable var = node.getVariable();			// this is the variable of the node
            if (updateVars.contains(var)) {			// is it in the data set?
                update.put(node, null);				// if yes, add it and move to next node
                continue;
            }
            List<EnumVariable> parents = node.getParents();	// if no, determine which are the parents of the node
            if (parents != null) {
                for (EnumVariable parent : parents) {		// go through the parents
                    if (updateVars.contains(parent)) {		// check if parent is in data set
                        update.put(node, null);			// if yes, add the node to the set of nodes to update
                        break;
                    }
                }
            }
        } // finished constructing the node update set

        Map<BNode, Object[]> node_keys = new HashMap<>();       // hold the evidence for each of the nodes for efficient access

        for (BNode node : update.keySet()) {
            if (node.isTrainable()) {
                node.randomize(System.currentTimeMillis());
            }
        }

        // randomize nodes that will be updated as part of learning ???
        // initialize learning parameters
        double log_prob = -999998;
        double prev_prob = -999999;
        double conv_rate = 0.02;
        int round = 0;
        boolean latentVariablesExist = false;

        
        // start training (keep going until convergence or stop criterion is met)
        while (conv_rate > EM_CONVERGENCE_CRITERION && round < EM_MAX_ROUNDS) {
            round++;
            prev_prob = log_prob;
            log_prob = 0;

            // for each sample with observations...
            for (int i = 0; i < values.length; i++) {
                // set variables and keys according to observations
                for (int j = 0; j < vars.length; j++) {
                    BNode instantiate_me = bn.getNode(vars[j]);
                    if (instantiate_me == null) {
                        throw new EMRuntimeException("Variable \"" + vars[i].getName() + "\" is not part of Bayesian network");
                    }
                    if (values[i][j] != null) { // check so that the observation is not null
                        // the node is instantiated to the value in the data set
                        instantiate_me.setInstance(values[i][j]);
                    } else { // the observation is null
                        // the node is reset, i.e. un-instantiated 
                        instantiate_me.resetInstance();
                    }
                }

                /*
                 * There are two principal ways of performing EM here, both based on the available evidence E:
                 * 
                 * 1. For each node, with variable X0 conditioned on X1, X2, ..., Xm; we refer to the set {X0, ..., Xm} as X.
                 * 	  For each node, infer the probability P(X'|E') where X' is {X} - {E} and E' is the union of {E} and {X}.
                 * 	  For each node, assign each value-permutation of X (X0 = x0, X1 = x1, ..., Xm = xm) the inferred probability P(X' = x'|E)
                 * 
                 * 2. Gather all nodes, to identify Y the total set of variables in the BN
                 * 	  Once, infer the probability P(Y'|E) where Y' is {Y} - {E}. 
                 * 	     Note this table can be large (size increases exponentially with number of variables)
                 * 	  For each node, marginalise over {Y'} - {X'} for P(Y'|E) where {Y'} - {X'} is the set of query variables not in the node.
                 * 	     Assign a probability to each value permutation.
                 * 
                 * Which of the ways is better depends on the number of variables in Y' and the number of nodes N.
                 */
                inf.instantiate(bn);
                Variable.Assignment[] evidence = Variable.Assignment.array(vars, values[i]); // the evidence here
                
                switch (EM_OPTION) {
                    case 1: 

                        // Principal way 1: go through each of the BN nodes... pose of query for, and update each...
                        
                        for (BNode node : update.keySet()) {

                            if (node.isTrainable()) {
                                
                                // identify what variables that we need to infer, to generate expectations
                                List<Variable> query_vars = new ArrayList<>();
                                Object[] evid_key = null; // only applicable if the node has parents
                                if (!node.isRoot()) {
                                    evid_key = EnumTable.getKey(node.getParents(), evidence);
                                    // add to the variables that must be queried during inference
                                    for (int key_index = 0; key_index < evid_key.length; key_index ++) {
                                        if (evid_key[key_index] == null)
                                            query_vars.add(node.getParents().get(key_index));
                                    }
                                }
                                Variable var = node.getVariable();
                                Object ovalue = node.getInstance(); // check the value, if instantiated
                                if (ovalue == null)
                                    query_vars.add(var);
                                
                                // check if inference is required
                                if (query_vars.size() > 0) { // there are unspecified/latent variables for this node
                                    try {
                                        latentVariablesExist = true;
                                        Variable[] query_arr = new Variable[query_vars.size()];
                                        query_vars.toArray(query_arr);
                                        Query q = inf.makeQuery(query_arr);
                                        CGTable qr = (CGTable) inf.infer(q);
                                        int[] indices = qr.getIndices(); 
                                        // for each permutation of the enumerable query variables
                                        for (int qr_index : indices) {
                                            Object[] qr_key = qr.getKey(qr_index);
                                            double p = qr.getFactor(qr_index);
                                            JDF jdf = null;
                                            if (qr.hasNonEnumVariables())
                                                jdf = qr.getJDF(qr_index);

                                            if (!node.isRoot()) { // if node has parents
                                                // we need to construct a key for the update of the node
                                                // first, put in the result from the inference, but in the order of the node's parents
                                                Variable.Assignment[] expected = Variable.Assignment.array(qr.getEnumVariables(), qr_key);
                                                Object[] inf_key = EnumTable.getKey(node.getParents(), expected);
                                                // second, overlay the evidence
                                                EnumTable.overlay(evid_key, inf_key);
                                                if (ovalue != null)
                                                    node.countInstance(evid_key, ovalue, p);
                                                else { // we don't know the value so use expected value
                                                    try {
                                                        EnumVariable evar = (EnumVariable) var;
                                                        for (Variable.Assignment assigned : expected) {
                                                            if (assigned.var.equals(evar)) {
                                                                node.countInstance(evid_key, assigned.val, p);
                                                                break;
                                                            }
                                                        }
                                                    } catch (ClassCastException e) {
                                                        // we think it is a continuous variable, so we should have a distrib for it
                                                        Distrib d = jdf.getDistrib(var);
                                                        node.countInstance(evid_key, d, p);
                                                    }
                                                }
                                            } else { // node IS root
                                                if (ovalue != null)
                                                    node.countInstance(null, ovalue, p);
                                                else { // we don't know the value so use expected value
                                                    try {
                                                        EnumVariable evar = (EnumVariable) var;
                                                        Variable.Assignment[] expected = Variable.Assignment.array(qr.getEnumVariables(), qr_key);
                                                        for (Variable.Assignment assigned : expected) {
                                                            if (assigned.var.equals(evar)) {
                                                                node.countInstance(null, assigned.val, p);
                                                                break;
                                                            }
                                                        }
                                                    } catch (ClassCastException e) {
                                                        // we think it is a continuous variable, but it is a root node!
                                                        throw new EMRuntimeException("Failed query for sample #"+(i+1)+": " + var.getName() + " is a non-enumerable root node");
                                                    }
                                                }
                                            }
                                        } 
                                    } catch (CGVarElimRuntimeException e) {
                                        throw new EMRuntimeException("Failed query for sample #"+(i+1)+" and node " + node.getName()+ ": " + e.getMessage());
                                    }
                                } else { // all variables are instantiated, no need to do inference
                                    node.countInstance(evid_key, ovalue);
                                }
                            }
                        }
                        if ((EM_PRINT_STATUS && round % 10 == 0) || round == 1)
                            log_prob += Math.log(((CGVarElim)inf).likelihood());

                        break; // end EM_OPTION == 1
                        
                    case 2:
                        // Principal way 2: collect query variables from all (trainable) nodes... pose one query and update all (trainable) nodes
                        Set<Variable> query_vars = new HashSet<>(); // all query variables are stored here

                        for (BNode node : update.keySet()) {
                            // check if the node should be updated, and if so collect query variables from it
                            if (node.isTrainable()) {
                                // if the node has parents, we need to check out the variables of its parents too
                                if (!node.isRoot()) {
                                    Object[] evid_key = EnumTable.getKey(node.getParents(), evidence);
                                    update.put(node, evid_key); // associate each node with a (potentially partial) key for evidence, expected values later...
                                    // add to the variables that must be queried during inference
                                    for (int key_index = 0; key_index < evid_key.length; key_index ++) {
                                        if (evid_key[key_index] == null)
                                            query_vars.add(node.getParents().get(key_index));
                                    }
                                }
                                Object ovalue = node.getInstance(); // check the value, if instantiated
                                if (ovalue == null) 
                                    query_vars.add(node.getVariable());
                            }
                        }
                        
                        if (query_vars.size() > 0) { // there are unspecified/latent variables for this node
                            // so we need to perform inference, which can go bad (hence potential for exception)
                            try {
                                latentVariablesExist = true;
                                // perform inference, ALL query variables in one go
                                Variable[] query_arr = new Variable[query_vars.size()];
                                query_vars.toArray(query_arr);
                                Query q = inf.makeQuery(query_arr);
                                CGTable qr = (CGTable) inf.infer(q); 

                                int[] indices = qr.getIndices(); // FIXME: For efficiency, at least initially, consider only looking at indices of events that are more probable...
                                // for each permutation of the enumerable query variables
                                for (int qr_index : indices) {
                                    Object[] qr_key = qr.getKey(qr_index);
                                    double p = qr.getFactor(qr_index);
                                    JDF jdf = null;
                                    if (qr.hasNonEnumVariables())
                                        jdf = qr.getJDF(qr_index);
                                    Variable.Assignment[] assignment = Variable.Assignment.array(qr.getEnumVariables(), qr_key);
                                    for (BNode node : update.keySet()) {
                                        // check if the node should be updated, and if so put together expectations for maximisation...
                                        if (node.isTrainable()) {
                                            // if the node has parents, we need to check out the variables of its parents too
                                            if (!node.isRoot()) {
                                                Object[] evid_key = update.get(node);
                                                Object[] inf_key = EnumTable.getKey(node.getParents(), assignment);
                                                EnumTable.overlay(evid_key, inf_key);
                                                Object value = node.getInstance();
                                                if (value != null)
                                                    node.countInstance(evid_key, value, p);
                                                else { // we don't know the value so use expected value
                                                    Variable var = node.getVariable();
                                                    try {
                                                        EnumVariable evar = (EnumVariable) var;
                                                        for (Variable.Assignment assigned : assignment) {
                                                            if (assigned.var.equals(evar)) {
                                                                node.countInstance(evid_key, assigned.val, p);
                                                                break;
                                                            }
                                                        }
                                                    } catch (ClassCastException e) {
                                                        // we think it is a continuous variable, so we should have a distrib for it
                                                        Distrib d = jdf.getDistrib(var);
                                                        node.countInstance(evid_key, d, p);
                                                    }
                                                }
                                            } else { // node IS root
                                                Object value = node.getInstance();
                                                if (value != null)
                                                    node.countInstance(null, value, p);
                                                else { // we don't know the value so use expected value
                                                    Variable var = node.getVariable();
                                                    try {
                                                        EnumVariable evar = (EnumVariable) var;
                                                        for (Variable.Assignment assigned : assignment) {
                                                            if (assigned.var.equals(evar)) {
                                                                node.countInstance(null, assigned.val, p);
                                                                break;
                                                            }
                                                        }
                                                    } catch (ClassCastException e) {
                                                        // we think it is a continuous variable, but it is a root node!
                                                        throw new EMRuntimeException("Failed query for sample #"+(i+1)+": " + var.getName() + " is a non-enumerable root node");
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            } catch (CGVarElimRuntimeException e) {
                                throw new EMRuntimeException("Failed query for sample #"+(i+1)+": " + e.getMessage());
                            }
                        } else { // nothing needs to be inferred
                            for (BNode node : update.keySet()) {
                                // check if the node should be updated, and if so put together expectations for maximisation...
                                if (node.isTrainable()) {
                                    // if the node has parents, we need to check out the variables of its parents too
                                    if (!node.isRoot()) {
                                        Object[] evid_key = update.get(node);
                                        node.countInstance(evid_key, node.getInstance());
                                    }
                                }
                            }
                        }
                        if ((EM_PRINT_STATUS && round % 10 == 0) || round == 1)
                            log_prob += Math.log(((CGVarElim)inf).likelihood());

                        break; // end EM_OPTION == 2
                }                        
            }
            // finally complete the M-step by transferring counts to probabilities
            for (BNode node : update.keySet()) {
                if (node.isTrainable()) {
                    node.maximizeInstance();
                }
            }

            if ((EM_PRINT_STATUS && round % 10 == 0) || round == 1) {
                conv_rate = Math.abs(log_prob - prev_prob); // use abs because the joint prob may exceed 1 (it is not normalized)
                System.err.println("Completed " + round + " round(s), L=" + log_prob);
            }
        }

        // un-set instances
        for (BNode node : bn.getNodes()) {
            node.resetInstance();
        }
        if (EM_PRINT_STATUS) {
            System.err.println("Completed " + round + " rounds, L=" + log_prob + ". Done.");
        }
    }

    public class EMRuntimeException extends RuntimeException {

        private static final long serialVersionUID = 1L;

        public EMRuntimeException(String message) {
            super(message);
        }
    }

}
