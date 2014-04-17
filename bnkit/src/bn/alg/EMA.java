

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
import bn.CPT;
import bn.Distrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.FactorTable;
import bn.JPT;
import bn.Variable;
import bn.file.BNBuf;

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
public class EMA extends LearningAlg {

    private Inference inf;

    /**
     * @param bn Bayesian network
     * @param inference the inference algorithm
     */
    public EMA(BNet bn, Inference inference) {
        super(bn);
        this.inf = inference;
    }

    /**
     * @param bn
     */
    public EMA(BNet bn) {
        this(bn, new ApproxInfer());
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
    public static int EM_MAX_ROUNDS = 10000;

    /**
     * EM status print-outs
     */
    public static boolean EM_PRINT_STATUS = true;

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
        Set<Variable> updateVars = new HashSet<Variable>();
        for (int i = 0; i < vars.length; i++) {
            updateVars.add(vars[i]);
        }
        Set<BNode> update = new HashSet<BNode>(); 			// the set of nodes that may need to be updated:
        for (BNode node : bn.getNodes()) {      			// check all nodes
            Variable var = node.getVariable();				// this is the variable of the node
            if (updateVars.contains(var)) {					// is it in the data set?
                update.add(node);							// if yes, add it and move to next node
                continue;
            }
            List<EnumVariable> parents = node.getParents();	// if no, determine which are the parents of the node
            if (parents != null) {
                for (EnumVariable parent : parents) {		// go through the parents
                    if (updateVars.contains(parent)) {		// check if parent is in data set
                        update.add(node);					// if yes, add the node to the set of nodes to update
                        break;
                    }
                }
            }
        } // finished constructing the node update set

        for (BNode node : update) {
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
                // set variables according to observations
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
//                bn.getNode("Alpha_latent").setInstance(true);
//                bn.getNode("Age_latent").setInstance(false);
//                bn.getNode("Cpl_latent").setInstance(true);
//                BNet mb = bn.getMB(bn.getNode("Class").getVariable());
//                Double test = mb.getMBProb(mb, bn.getNode("Class"));
//                System.out.println("End of test");

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
                // Principal way 1: go through each of the BN nodes...
                inf.instantiate(bn);
                for (BNode node : update) {

                    if (node.isTrainable()) {
//                    	System.out.println("Node");
//                    	System.out.println(node.toString());
//                    	System.out.println(node.getInstance());
                        List<Variable> query_vars = new ArrayList<Variable>();
                        List<EnumVariable> parents = node.getParents();
                        Object[] parent_key = null;
                        int[] parent_map = null;
                        if (parents != null) {
//                        	System.out.println("parents");
//                            for (EnumVariable pa : parents) {
//                            	System.out.println(bn.getNode(pa).toString());
//                            }
                            parent_key = new Object[parents.size()];
                            parent_map = new int[parents.size()];
                            int index_in_query = 0;
                            for (int p = 0; p < parents.size(); p++) {
                                Variable vpar = parents.get(p);
                                BNode npar = bn.getNode(vpar);
                                parent_key[p] = npar.getInstance();
                                if (parent_key[p] == null) {
                                    query_vars.add(vpar);
                                    parent_map[p] = index_in_query++;
                                } else {
                                    parent_map[p] = -1;
                                }
                            }
                        }
                        Variable var = node.getVariable();
                        Object ovalue = node.getInstance(); // check the value, if instantiated
                        if (ovalue == null) {
                            query_vars.add(var);
                        }
                        if (query_vars.size() > 0) { // there are unspecified/latent variables for this node
                            latentVariablesExist = true;
                            Variable[] query_arr = new Variable[query_vars.size()];
                            query_vars.toArray(query_arr);
                            Query q = inf.makeQuery(query_arr);
                            QueryResult qr = inf.infer(q);
                            JPT jpt = qr.getJPT();
                            log_prob += inf.getLogLikelihood();
//                            Set<Map.Entry<Integer, Double>> test = jpt.table.getMapEntries();
                            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                                Object[] jpt_key = jpt.table.getKey(entry.getKey().intValue());
				// jpt_key will contain values for latent variables
                                // other variables are already instantiated
                                if (node.getInstance() == null) {
                                    ovalue = jpt_key[jpt_key.length - 1];
                                }
                                double prob = entry.getValue().doubleValue();
                                //Some sort of bug here - null as key?
                                if (node.isRoot()) { // if prior
                                    node.countInstance(null, ovalue, prob);
                                } else { // if node has parents
                                    for (int p = 0; p < parent_key.length; p++) {
                                        if (parent_map[p] >= 0) {
                                            parent_key[p] = jpt_key[parent_map[p]];
                                        }
                                    }
                                    try {
                                        node.countInstance(parent_key, ovalue, prob);
                                    } catch (java.lang.RuntimeException e) {
                                        throw new EMRuntimeException("Problem with sample #"+(i+1)+": " + e.getMessage());
                                    }
                                }
                            }

                        } else { // all variables are instantiated, no need to do inference
                            node.countInstance(parent_key, ovalue, 1.0);
                        }
                    }
                }
                // Principal way 2: not yet implemented but see old bioalg:ml.bayes which does it this way.
                ;;

            }
            // finally complete the M-step by transferring counts to probabilities
            for (BNode node : update) {
                if (node.isTrainable()) {
                        node.maximizeInstance();
                }
            }

            conv_rate = Math.abs(log_prob - prev_prob); // use abs because the joint prob may exceed 1 (it is not normalized)
            if ((EM_PRINT_STATUS && round % 10 == 0) || round == 1) {
                System.err.println("Completed " + round + " round(s), L=" + log_prob);
//                BNBuf.save(bn, "antonTrain1.new");
            }
        }

        // un-set instances
        for (BNode node : bn.getNodes()) {
            node.resetInstance();
        }
        if (EM_PRINT_STATUS || round == 180) {
            System.err.println("Completed " + round + " rounds, L=" + log_prob + ". Done.");
            BNBuf.save(bn, "antonTrain2.new");
        }
    }

    public class EMRuntimeException extends RuntimeException {

        private static final long serialVersionUID = 1L;

        public EMRuntimeException(String message) {
            super(message);
        }
    }

}


