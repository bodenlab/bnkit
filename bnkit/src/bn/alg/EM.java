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
import dat.EnumTable;
import dat.EnumVariable;
import bn.JDF;
import dat.Variable;
import bn.alg.VarElim.VarElimRuntimeException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Implementation of Expectation-Maximization for learning in Bayesian network.
 * Current implementation uses variable elimination for inference.
 * It has options for multi-threaded operation.
 * @author mikael
 * @author alex
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
        this(bn, new VarElim());
    }

    /**
     * EM convergence is true (and training stops) when the improvement is
     * smaller than this value (expressed in percent; based on the average
     * log-likelihood over the last 5 x 10 rounds; only sampled every 10 rounds). 
     * Note that smaller value leads to longer convergence times but closer 
     * to optimal solutions.
     */
    public double EM_CONVERGENCE_CRITERION = 0.00005; //


    /**
     * EM maximum number of rounds (iterations). Training stops after the
     * specified number of rounds independently of convergence.
     */
    public int EM_MAX_ROUNDS = 1000;
    
    /**
     * EM status print-outs
     */
    public boolean EM_PRINT_STATUS = true;
    
    /**
     * Threads to be used in EM - case 3
     */
    public int EM_THREAD_COUNT = 1;

    /**
     * EM option: currently two different approaches to determine expectations.
     * In either case, nodes are then updated one at a time.
     * 1. one node at a time = one query per node
     * 2. all nodes at a time = one big query
     */
    public int EM_OPTION = 1;

    public void setConvergeCrit(double value){
        this.EM_CONVERGENCE_CRITERION = (value > 0.0) ? value : 0.00005;
    }

    /**
     * Set the print status for the training
     * @param status
     */
    public void setPrintStatus(Boolean status){
        this.EM_PRINT_STATUS = status;
    }
    
    /**
     * Set the threads to be used in EM - case 1
     * @param thread
     */
    public void setThreadCount(int thread) {
    	this.EM_THREAD_COUNT = thread;
    }

    /**
     * Set EM training option
     * @param option
     */
    public void setEMOption(int option){
    	this.EM_OPTION = (option == 1) ? 
    			option : 
    				(option == 2 ? 2 : 3);
    }

    /**
     * Set the max number of rounds for training.
     * @param rounds
     */
    public void setMaxRounds(int rounds){
        this.EM_MAX_ROUNDS = (rounds >= 1) ? rounds : 1;
    }

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

        for (BNode node : bn.getNodes()) {
            if (node.isTrainable()) {
                node.randomize(seed);
            }
        }

        // initialize learning parameters
        double[] last_LL = new double[] {-999999, -999999, -999999, -999999, -999999};
        
        int round = 0;
        boolean latentVariablesExist = false;
        
        boolean EM_TERMINATE = false;

        // start training (keep going until convergence or stop criterion is met)
        while (round < EM_MAX_ROUNDS && !EM_TERMINATE) {
            round++;

            double log_likelihood = 0;
            // the set of nodes to be updated:
            Map<BNode, Object[]> updateMap = new HashMap<>();
            Map<BNode, Object[]> update = Collections.synchronizedMap(updateMap); //THREAD SAFE MAP
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
                        //Has evidence so needs to be updated
                        update.put(instantiate_me, null);
                        List<EnumVariable> parents = instantiate_me.getParents();	// if no, determine which are the parents of the node
                        if (parents != null) {
                            for (EnumVariable parent : parents) {		// go through the parents and add them to the update set
                                update.put(bn.getNode(parent), null);
                            }
                        }
                        Set<String> children = bn.getChildrenNames(instantiate_me);	// if no, determine which are the parents of the node
                        if (children != null) {
                            for (String child : children) {		// go through the children and add them to the update set
                                BNode c = bn.getNode(child);
                                update.put(c, null);
                            }
                        }
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
                 * Should benchmark this so way is chose automatically depending on the query complexity.
                 */
                inf.instantiate(bn);
                Variable.Assignment[] evidence = Variable.Assignment.array(vars, values[i]); // the evidence here
                
                switch (EM_OPTION) {
                case 1: 

                    // Principal way 1: go through each of the BN nodes... pose a query for, and update each...
                    if (i == 11)
                        i = 11;
                    
                    for (BNode node : update.keySet()) {

                        if (node.isTrainable()) {

                            // identify what variables that we need to infer, to generate expectations
                            List<Variable> query_vars = new ArrayList<>();
                            Object[] evid_key = null; // only applicable if the node has parents
                            if (!node.isRoot()) {
                                evid_key = EnumTable.getKey(node.getParents(), evidence);
                                // add to the variables that must be queried during inference
                                for (int key_index = 0; key_index < evid_key.length; key_index++) {
                                    if (evid_key[key_index] == null) {
                                        query_vars.add(node.getParents().get(key_index));
                                    }
                                }
                            }
                            Variable var = node.getVariable();
                            Object ovalue = node.getInstance(); // check the value, if instantiated
                            if (ovalue == null) {
                                query_vars.add(var);
                            }

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
                                        if (p == 0 || Double.isNaN(p)) // count is zero (or the log prob was so small that conversion failed)
                                            continue;
                                        JDF jdf = null;
                                        if (qr.hasNonEnumVariables()) {
                                            jdf = qr.getJDF(qr_index);
                                        }

                                        if (!node.isRoot()) { // if node has parents
                                            // we need to construct a key for the update of the node
                                            // first, put in the result from the inference, but in the order of the node's parents
                                            Variable.Assignment[] expected = Variable.Assignment.array(qr.getEnumVariables(), qr_key);
                                            Object[] inf_key = EnumTable.getKey(node.getParents(), expected);
                                            // second, overlay the evidence
                                            EnumTable.overlay(evid_key, inf_key);
                                            if (ovalue != null) {
                                                node.countInstance(evid_key, ovalue, p);
                                            } else { // we don't know the value so use expected value
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
                                            if (ovalue != null) {
                                                node.countInstance(null, ovalue, p);
                                            } else { // we don't know the value so use expected value
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
                                                    throw new EMRuntimeException("Failed query for sample #" + (i + 1) + ": " + var.getName() + " is a non-enumerable root node");
                                                }
                                            }
                                        }
                                    }
                                } catch (RuntimeException e) {
                                    throw new EMRuntimeException("Failed query for sample #" + (i + 1) + " and node " + node.getName() + ": " + e.getLocalizedMessage());
                                }
                            } else { // all variables are instantiated, no need to do inference
                                node.countInstance(evid_key, ovalue);
                            }
                        }
                    }

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
                                    for (int key_index = 0; key_index < evid_key.length; key_index++) {
                                        if (evid_key[key_index] == null) {
                                            query_vars.add(node.getParents().get(key_index));
                                        }
                                    }
                                }
                                Object ovalue = node.getInstance(); // check the value, if instantiated
                                if (ovalue == null) {
                                    query_vars.add(node.getVariable());
                                }
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
                                    if (qr.hasNonEnumVariables()) {
                                        jdf = qr.getJDF(qr_index);
                                    }
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
                                                if (value != null) {
                                                    node.countInstance(evid_key, value, p);
                                                } else { // we don't know the value so use expected value
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
                                                if (value != null) {
                                                    node.countInstance(null, value, p);
                                                } else { // we don't know the value so use expected value
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
                                                        throw new EMRuntimeException("Failed query for sample #" + (i + 1) + ": " + var.getName() + " is a non-enumerable root node");
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            } catch (RuntimeException e) {
                                throw new EMRuntimeException("Failed query for sample #" + (i + 1) + ": " + e.getMessage());
                            }
                        } else { // nothing needs to be inferred
                            for (BNode node : update.keySet()) {
                                // check if the node should be updated, and if so put together expectations for maximisation...
                                if (node.isTrainable()) {
                                    Object[] evid_key = null;
                                    // if the node has parents, we need to check out the variables of its parents too
                                    if (!node.isRoot()) {
                                        evid_key = update.get(node);
                                    }
                                    node.countInstance(evid_key, node.getInstance());
                                }
                            }
                        }
                        break; // end EM_OPTION == 2

                    case 3:
                        // Principal way 1 with parallel implementation: go through each of the BN nodes... pose of query for, and update each...
                        ExecutorService executor = Executors.newFixedThreadPool(EM_THREAD_COUNT);
                        for (BNode node : update.keySet()) {
//                		node.getTable().getSize();
                            EMc1 work = new EMc1(node, evidence, i);
                            executor.execute(work);
                        }
                        executor.shutdown();
                        while (!executor.isTerminated()) {
                        }
                        //                        System.out.println("Finished all threads");

                        break; // end EM_OPTION == 3
                }                        
            }
            
            // finally complete the M-step by transferring counts to probabilities
            for (BNode node : update.keySet()) {
                if (node.isTrainable()) {

                    node.maximizeInstance();
                }
            }

            if (round % 10 == 0) { // || round == 1) {
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
                    double sample_likelihood = ((VarElim) inf).logLikelihood();
                    if (Double.isNaN(sample_likelihood)) {
                        System.err.println("Sample " + i + "/" + values.length + " log-likelihood is " + sample_likelihood);
//                        for (int j = 0; j < vars.length; j++) {
//                            System.err.println("\t" + vars[j].getName() + " = " + values[i][j]);
//                        }
                    }
                    log_likelihood += sample_likelihood;
                    if (Double.isInfinite(log_likelihood)) {
                        System.err.println("Log-likelihood is infinite: " + log_likelihood);
                    }
                }
            }

            if (round % 10 == 0) { // & round %50 == 0
                // summarise progress
                double sd_LL = 0;
                // copy previous LL (log-likelihood of data)
                if (round <= 50) {
                    for (int i = 0; i < last_LL.length - 1; i++) {
                        last_LL[i] = last_LL[i + 1]; // shuffle llhs through list - only ever record 5
                    }
                    last_LL[last_LL.length - 1] = log_likelihood; // add new log_likelihood
//            		System.out.println(Arrays.toString(last_LL));
                } else {
                    double mean_LL = last_LL[0] / last_LL.length;
                    for (int i = 0; i < last_LL.length - 1; i++) {
                        last_LL[i] = last_LL[i + 1];
                        mean_LL += (last_LL[i] / last_LL.length);
                    }
                    double[] sdl_LL = new double[last_LL.length];
                    for (int j = 0; j < last_LL.length; j++) {
                    	sdl_LL[j] = (last_LL[j] - mean_LL)*(last_LL[j] - mean_LL);
                    }
                    sd_LL = sdl_LL[0] / sdl_LL.length;
                    for (int i = 0; i < last_LL.length - 1; i++) {
                        sdl_LL[i] = sdl_LL[i + 1];
                        sd_LL += (sdl_LL[i] / sdl_LL.length);
                    }
                    last_LL[last_LL.length - 1] = log_likelihood;
//                    if ((-mean_LL - -log_likelihood) < (EM_CONVERGENCE_CRITERION * 0.01 * -mean_LL)) // percent improvement < EM_CONVERGENCE_CRITERION
                    if(mean_LL > 0) {
                        if (sd_LL < (EM_CONVERGENCE_CRITERION * 0.01 * mean_LL)) // percent improvement < EM_CONVERGENCE_CRITERION
                        {
                            EM_TERMINATE = true;
                        }
                    } else {
                        if (sd_LL < (EM_CONVERGENCE_CRITERION * 0.01 * -mean_LL)) // percent improvement < EM_CONVERGENCE_CRITERION
                        {
                            EM_TERMINATE = true;
                        }
                    }
                }

                if (EM_PRINT_STATUS || EM_TERMINATE) {
                    System.err.println("Completed " + round + " round(s), L = " + log_likelihood);
                    if (EM_TERMINATE)
                        System.err.println("SD(L) = " + sd_LL);
                }
            }
        }

        // un-set instances
        for (BNode node : bn.getNodes()) {
            node.resetInstance();
        }
        if (EM_PRINT_STATUS) {
            System.err.println("Done.");
        }
    }

    public class EMRuntimeException extends RuntimeException {

        private static final long serialVersionUID = 1L;

        public EMRuntimeException(String message) {
            super(message);
        }
    }
    
    /**
     * Runnable 'worker' for case 3 of EM
     * In order for the code to be threaded, the process of the for loop for case 1 had
     * to be removed to a new class
     * This class can then be handed as work to each thread
     * @author Alex
     *
     */
    public class EMc1 implements Runnable {

        private boolean init = false;
        private BNode node;
        private Variable.Assignment[] evidence;
        private int i;

        public EMc1(BNode node, Variable.Assignment[] evidence, int i) {
            this.init = true;
            this.node = node;
            this.evidence = evidence;
            this.i = i;
        }

        @Override
        public void run() {
//    		System.out.println(Thread.currentThread().getName()+" Start. Node = "+node.toString());
            case3(node, evidence, i);
//    		System.out.println(Thread.currentThread().getName()+" End.");
        }

        public void case3(BNode node, Variable.Assignment[] evidence, int i) {

            boolean latentVariablesExist = false;

            if (node.isTrainable()) {

                // identify what variables that we need to infer, to generate expectations
                List<Variable> query_vars = new ArrayList<>();
                Object[] evid_key = null; // only applicable if the node has parents
                if (!node.isRoot()) {
                    evid_key = EnumTable.getKey(node.getParents(), evidence);
                    // add to the variables that must be queried during inference
                    for (int key_index = 0; key_index < evid_key.length; key_index++) {
                        if (evid_key[key_index] == null) {
                            query_vars.add(node.getParents().get(key_index));
                        }
                    }
                }
                Variable var = node.getVariable();
                Object ovalue = node.getInstance(); // check the value, if instantiated
                if (ovalue == null) {
                    query_vars.add(var);
                }

                // check if inference is required
                if (query_vars.size() > 0) { // there are unspecified/latent variables for this node
                    try {
                        latentVariablesExist = true;
                        Variable[] query_arr = new Variable[query_vars.size()];
                        query_vars.toArray(query_arr); //FIXME NOT THREADSAFE - does inference update these nodes in any way?
                        Query q = inf.makeQuery(query_arr); 
                        CGTable qr = (CGTable) inf.infer(q);
                        int[] indices = qr.getIndices();
                        // for each permutation of the enumerable query variables
                        for (int qr_index : indices) {
                            Object[] qr_key = qr.getKey(qr_index);
                            double p = qr.getFactor(qr_index);
                            JDF jdf = null;
                            if (qr.hasNonEnumVariables()) {
                                jdf = qr.getJDF(qr_index);
                            }

                            if (!node.isRoot()) { // if node has parents
                                // we need to construct a key for the update of the node
                                // first, put in the result from the inference, but in the order of the node's parents
                                Variable.Assignment[] expected = Variable.Assignment.array(qr.getEnumVariables(), qr_key);
                                Object[] inf_key = EnumTable.getKey(node.getParents(), expected);
                                // second, overlay the evidence
                                EnumTable.overlay(evid_key, inf_key); //FIXME
                                if (ovalue != null) {
                                    node.countInstance(evid_key, ovalue, p);
                                } else { // we don't know the value so use expected value
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
                                if (ovalue != null) {
                                    node.countInstance(null, ovalue, p);
                                } else { // we don't know the value so use expected value
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
                                        throw new EMRuntimeException("Failed query for sample #" + (i + 1) + ": " + var.getName() + " is a non-enumerable root node");
                                    }
                                }
                            }
                        }
                    } catch (RuntimeException e) {
                        throw new EMRuntimeException("Failed query for sample #" + (i + 1) + " and node " + node.getName() + ": " + e.getMessage());
                    }
                } else { // all variables are instantiated, no need to do inference
                    node.countInstance(evid_key, ovalue);
                }
            }
        }
    }
}
