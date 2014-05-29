package bn.alg;

import bn.BNet;
import bn.BNode;
import bn.CPT;
import bn.CountTable;
import bn.DataSample;
import bn.Distrib;
import bn.EnumDistrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.Enumerable;
import bn.FactorTable;
import bn.GaussianDistrib;
import bn.JPT;
import bn.MixtureDistrib;
import bn.Predef;
import bn.SampleTable;
import bn.SampleTable.Sample;
import bn.Variable;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;

/**
 * Approximate inference in Bayesian network by Gibbs algorithm (MCMC). In
 * accordance with the method described in Russell and Norvig, Artificial
 * Intelligence: A Modern Approach, 3e, 2009.
 *
 * Functional with a hybrid Bayesian network Does not currently conform to
 * CGTable output of query
 *
 * Convergence of algorithm is incomplete
 *
 * @author Alex
 *
 */
public class ApproxInference implements Inference {

    public BNet bn;
    private double logLikelihood = 1;
    private Random randomGenerator = new Random();
    public static int iterations = 500;

    /**
     * Approximate inference in Bayesian network by Gibbs algorithm (MCMC). In
     * accordance with the method described in Russell and Norvig, Artificial
     * Intelligence: A Modern Approach, 3e, 2009.
     *
     * Convergence of algorithm is incomplete
     */
    public void instantiate(BNet bn) {
        this.bn = bn;
        this.bn.compile();
    }

    /**
     * Construct the data structure for the query X - query variable
     * (non-evidence) E - evidence variables in bn Z - all non-evidence
     * variables
     *
     * @param qvars
     */
    @Override
    public Query makeQuery(Variable... qvars) {
        List<Variable> X = new ArrayList<>(); // Query variables
        List<Variable> E = new ArrayList<>(); // Evidence variables
        List<BNode> Z = new ArrayList<>(); // 
        BNet qbn = bn.getRelevant(qvars); // create new BN with variables that are relevant to query, 
        try {
            X.addAll(Arrays.asList(qvars));
            for (BNode node : qbn.getOrdered()) { // topological order: top-down
                Variable var = node.getVariable();
                if (node.getInstance() != null) {
                    E.add(var);
                    //need to keep track of ALL non-evidence variables
                    //this includes the query variable
                } else {
                    Z.add(qbn.getNode(var));
                }
            }
        } catch (RuntimeException e) {
            throw new RuntimeException("makeQuery, ApproxInfer didn't work");
        }

        return new AQuery(X, E, Z, qbn);
    }

    /**
     * Perform Approximate inference using Gibbs sampling algorithm
     */
    @SuppressWarnings("rawtypes")
    public AResult infer(Query query) {
        JPT answer = null;
        AQuery q = (AQuery) query;
        // BN that will be queried
        BNet cbn = q.qbn;
        // First set all non-evidenced nodes including query
        cbn.sampleInstance();
        
        DataSample data = new DataSample(q.X); //Storage class - maintains instance of each query node for each 'state' the chain passes through

        //Iterations of sampling
        int N = iterations;

        //The main loop for the algorithm (see Russell and Norvig 2e p. 517)
        for (int j = 0; j < N; j++) {
            //Iterate over all non-evidenced nodes, including query nodes
            for (BNode node : q.Z) {
                // These variables are in "topological order" (a node is never seen until all its parents have been seen)
                // Get the Markov blanket for the node
                Set<BNode> mbVar = cbn.getMB(node);
                // Sample from the mb distribution
                Object result = bn.getMBProb(mbVar, node); //FIXME - this is the most time consuming element of ApproxInfer
                if (result != null) {
                    cbn.getNode(node.getVariable()).setInstance(result);
                }
            }

            //Record instances of query node for current state of network
            for (Variable cQuery : q.X) {
                data.addValue(cQuery, cbn.getNode(cQuery).getInstance());
            }
        }

        //Reset all unevidenced nodes in network
        for (BNode node : q.Z) 
            node.resetInstance();

        data.createData(); //Take raw counts and sort into appropriate data structures
//		logLikelihood = logLikelihood*0.8;//Random value for training purposes when needed

            //Process the results stored in data
        //FIXME adapt this to CGTable output?
        if (data.allContinuous()) { //Query contains only continuous nodes
            Map<Variable, Distrib> allCont = data.getGaussianDistrib();
            data.map = null;
            return new AResult(allCont);
        } else if (data.getNonEnumTable() == null) { //Query contains only discrete nodes
            EnumTable<Double> disc = data.getNormalizedCounts();
            JPT result = new JPT(disc);
            data.map = null;
            data.counts = null;
            return new AResult(result);
        } else { //Mixed/hybrid query
            EnumTable<Double> disc = data.getNormalizedCounts();
            Map<Variable, EnumTable<Distrib>> cont = data.getMixedGDistrib();
            JPT result = new JPT(disc);
            data.map = null;
            data.counts = null;
            data.nonEnumTables = null;
            return new AResult(result, cont);
        }

//		//FIXME RESET Data structures?  
//		FactorTable result = null;
            //TODO - Convergence of algorithm
    }

//    public BNet instantiateNet(List<Variable> uninstantiated, BNet cbn) {
//    	//Occurs only once during inference
//        //Should be done randomly**
//        for (Variable var : uninstantiated) {
//            String pre = var.getPredef();
//            //getParams() did not return suitable results for boolean nodes
//            //Set manually
//            if (pre.equals("Boolean")) {
//                Object[] params = {true, false};
//                //Generate a pseudo random number to select parameter
//                int index = randomGenerator.nextInt(params.length);
//                //Randomly set the instance for the node
//                cbn.getNode(var).setInstance(params[index]);
//            } else if (pre.equals("Real")) {
//                // Ultimately want to sample from one of the possible distributions
//                // FIXME - better way to do this sampling? Something without odd casts
//                Object[] nodeDistribs = cbn.getNode(var).getTable().getValues().toArray();
//                int index = randomGenerator.nextInt(nodeDistribs.length);
//                //Get individual mu and sigma from random distribution
//                String[] values = nodeDistribs[index].toString().split(";");
//                //Create distribution to sample from
//                Distrib d = new GaussianDistrib(Double.parseDouble(values[0]), Double.parseDouble(values[1]));
//                //If this distribution is randomly generated it will not be a very accurate value here
//                cbn.getNode(var).setInstance(d.sample());
//            } else {
//                String parm = var.getParams();
//                String[] params;
//                if (parm != null) {
//                    params = parm.split(";");
//                    //Generate a pseudo random number to select parameter
//                    int index = randomGenerator.nextInt(params.length);
//                    //Randomly set the instance for the node
//                    cbn.getNode(var).setInstance(params[index]);
//                } else {
//                    throw new ApproxInferRuntimeException("Node must contain parameters");
//                }
//            }
//        }
//        return cbn;
//    }

    /**
     * Get the number of iterations sampling will complete
     *
     * @return iterations
     */
    public static int getIterations() {
        return iterations;
    }

    /**
     * Set the number of iterations sampling will complete
     *
     * @return iterations
     */
    public void setIterations(int iter) {
        iterations = iter;
    }

    public double getLogLikelihood() {
        return logLikelihood;
    }

    public class ApproxInferRuntimeException extends RuntimeException {

        private static final long serialVersionUID = 1L;

        public ApproxInferRuntimeException(String message) {
            super(message);
        }
    }

    public class AQuery implements Query {

        final List<Variable> X;
        final List<Variable> E;
        final List<BNode> Z;
        final BNet qbn;

        AQuery(List<Variable> X, List<Variable> E, List<BNode> Z, BNet qbn) {
            this.X = X;
            this.E = E;
            this.Z = Z;
            this.qbn = qbn;
        }
    }

    public class AResult implements QueryResult {

        final private JPT jpt;
        final private Map<Variable, EnumTable<Distrib>> nonEnumTables;
        final private Map<Variable, Distrib> nonEnumDistribs;

        public AResult(JPT jpt) {
            this.jpt = jpt;
            this.nonEnumTables = null;
            this.nonEnumDistribs = null;
        }

        public AResult(JPT jpt, Map<Variable, EnumTable<Distrib>> nonEnum) {
            this.jpt = jpt;
            this.nonEnumTables = nonEnum;
            this.nonEnumDistribs = null;
        }

        public AResult(Map<Variable, Distrib> nonEnum) {
            this.jpt = null;
            this.nonEnumTables = null;
            this.nonEnumDistribs = nonEnum;
        }

        @Override
        public JPT getJPT() {
            return this.jpt;
        }

        public Map<Variable, EnumTable<Distrib>> getNonEnum() {
            return this.nonEnumTables;
        }

        public Map<Variable, Distrib> getNonEnumDistrib() {
            return this.nonEnumDistribs;
        }

        /**
         * Method to retrieve a single distribution with all other variables in
         * original query unspecified.
         *
         * @param query the query variable
         * @return
         */
        public Distrib getDistrib(Variable query) {
            return getDistrib(query, null);
        }

        /**
         * Method to retrieve a single distribution GIVEN some optional
         * evidence.
         *
         * @param query the query variable
         * @return
         */
        public Distrib getDistrib(Variable query, Evidence... evid) {

            return null;
        }

        public class Evidence {

            public Variable var;
            Object val;

            public Evidence(Variable var, Object val) {
                this.var = var;
                this.val = val;
            }
        }
    }

}