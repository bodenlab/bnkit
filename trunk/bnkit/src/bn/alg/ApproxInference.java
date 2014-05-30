package bn.alg;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.EnumTable;
import bn.JPT;
import bn.SampleTrace;
import bn.Variable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

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
        List<BNode> X = new ArrayList<>(); // Query variables
        List<BNode> E = new ArrayList<>(); // Evidence variables
        List<BNode> Z = new ArrayList<>(); // 
        BNet qbn = bn.getRelevant(qvars); // create new BN with variables that are relevant to query, 
        try {
            for (Variable x : qvars) 
                X.add(qbn.getNode(x));
            for (BNode node : qbn.getOrdered()) { // topological order: top-down
                Variable var = node.getVariable();
                if (node.getInstance() != null) {
                    E.add(qbn.getNode(var));
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
     * @param query
     */
    @SuppressWarnings("rawtypes")
    @Override
    public AResult infer(Query query) {
        AQuery q = (AQuery) query;
        // BN that will be queried
        BNet cbn = q.qbn;
        // First set all non-evidenced nodes including query
        cbn.sampleInstance();
        
        SampleTrace data = new SampleTrace(q.X); //Storage class - maintains instance of each query node for each 'state' the chain passes through

        Map<BNode, Object> nextInstance = new HashMap<>();
        //Iterations of sampling
        int N = iterations;

        //The main loop for the algorithm (see Russell and Norvig 2e p. 517)
        for (int j = 0; j < N; j++) {
            //Iterate over all non-evidenced nodes, including query nodes
            for (BNode node : q.Z) {
                // These variables are in "topological order" (a node is never seen until all its parents have been seen)
                // Get the Markov blanket for the node
                //Set<BNode> mbVar = cbn.getMB(node); // don't need it with the getMBProb below
                // Sample from the mb distribution
                //Object result = bn.getMBProb(mbVar, node); //FIXME - this is the most time consuming element of ApproxInfer
                Object result = bn.getMBProb(node);
                if (result != null) {
                    nextInstance.put(node, result);
                    //node.setInstance(result);// wait until after the full round
                }
            }

            // instantiate the nodes according to the MB-based sampling
            for (Map.Entry<BNode, Object> inst : nextInstance.entrySet()) 
                inst.getKey().setInstance(inst.getValue());
            
            //Record instances of query node for current state of network
            for (BNode cQuery : q.X) {
                data.addValue(cQuery, cQuery.getInstance());
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

    }


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

        final List<BNode> X;
        final List<BNode> E;
        final List<BNode> Z;
        final BNet qbn;

        AQuery(List<BNode> X, List<BNode> E, List<BNode> Z, BNet qbn) {
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
