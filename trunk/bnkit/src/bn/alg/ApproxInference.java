package bn.alg;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import dat.EnumTable;
import bn.factor.Factor;
import bn.JPT;
import bn.SampleTrace;
import dat.Variable;
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
//        BNet qbn = bn.getRelevant(qvars); // create new BN with variables that are relevant to query, 
        List<BNode> rnl = bn.getDconnected(qvars); //List of relevant nodes to be used based on FULL network
        try {
            for (Variable x : qvars) 
                X.add(bn.getNode(x));
            for (BNode node : rnl) { // topological order: top-down
                Variable var = node.getVariable();
                if (node.getInstance() != null) {
                    E.add(bn.getNode(var));
                    //need to keep track of ALL non-evidence variables
                    //this includes the query variable
                } else {
                    Z.add(bn.getNode(var));
                }
            }
        } catch (RuntimeException e) {
            throw new RuntimeException("makeQuery, ApproxInfer didn't work");
        }
        return new AQuery(X, E, Z, rnl);
    }

    /**
     * Perform approximate inference using Gibbs sampling algorithm
     * @param query
     */
    @SuppressWarnings("rawtypes")
    @Override
    public CGTable infer(Query query) {
        AQuery q = (AQuery) query;
        // BN that will be queried
//        BNet cbn = q.qbn;
        List<BNode> rnl = q.rnl;
        // First set all non-evidenced nodes including query
        bn.sampleInstance();
//        bn.sampleInstance(rnl); // will instantiate all nodes
        SampleTrace data = new SampleTrace(q.X, iterations);    // Storage class - maintains instance of each query node for each 'state' the chain passes through
        data.count();                            // Observe the current instantiation, starting 'state' of chain
        
        Map<BNode, Object> nextInstance = new HashMap<>();
        // Iterations of sampling
        int N = iterations;
        //The main loop for the algorithm (see Russell and Norvig 2e p. 517)
        for (int j = 0; j < N - 1; j++) {
            //Iterate over all non-evidenced nodes, including query nodes
            for (BNode node : q.Z) {
                // These variables are in "topological order" (a node is never seen until all its parents have been seen)
                // Get the Markov blanket for the node 
            	//The minimal set of nodes which d-separates node A from all other nodes is A's Markov blanket (MB)
                // Sample from the Markov blanket distribution of the node
                Object result = bn.getMBProb(node);
                if (result != null) {
                    nextInstance.put(node, result); // save the sample for later
                }
            }
            // instantiate the nodes according to the MB-based sampling
            for (Map.Entry<BNode, Object> inst : nextInstance.entrySet()) 
                inst.getKey().setInstance(inst.getValue());
            //Record instances of query node for current state of network
            data.count();
        }
        // Reset all unevidenced nodes in network
        for (BNode node : q.Z) 
            node.resetInstance();
        // put everything in place
        Factor f = data.getFactor();
        return new CGTable(f);
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
        // FIXME: Not implemented
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
        final List<BNode> rnl;

        AQuery(List<BNode> X, List<BNode> E, List<BNode> Z, List<BNode> rnl) {
            this.X = X;
            this.E = E;
            this.Z = Z;
            this.rnl = rnl;
        }
    }

}
