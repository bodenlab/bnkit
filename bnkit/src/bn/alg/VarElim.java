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

import asr.GRASP;
import bn.BNet;
import bn.BNode;
import bn.ctmc.SubstNode;
import dat.EnumVariable;
import bn.factor.Factor;
import bn.JPT;
import dat.Variable;
import bn.factor.AbstractFactor;
import bn.factor.Factorize;
import bn.factor.Factorize.FactorProductTree;
import util.MilliTimer;

import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Exact inference in Bayesian network by variable elimination, more
 * specifically a variant of "Bucket Elimination" that involves conditional non-enumerable 
 * variables, including Gaussians. 
 * 
 * Implementation of belief updating in accordance with the method described in 
 * Dechter, R. Bucket Elimination: A Unifying Framework for Probabilistic Inference, 
 * in Uncertainty in Artificial Intelligence, 1998.
 * 
 * The idea of conditional Gaussians is discussed in Lauritzen SL and Jensen F, 
 * Stable local computation with conditional Gaussian distributions, Statistics and
 * Computing 2001 11:191-203.
 * 
 * @author mikael
 */
public class VarElim implements Inference {

    public BNet bn;
    
    // Query status settings:
    static int STATUS_BEL = 0; // Belief (probability of)
    static int STATUS_MPE = 1; // Most probable explanation

    @Override
    public void instantiate(BNet bn) {
        this.bn = bn;
        this.bn.compile();
    }

    /**
     * Construct the data structure for the specified variables in preparation
     * of inference of belief. There are three types of variables (given the BN): 
     *  1. Assignment variables E--which have been assigned values via the BN 
     *  2. Query variables Q--for which probabilities are sought P(Q|E) 
     *  3. Other variables X--which will be summed out during inference P(Q|E) = SUM_X P(Q|E,X).
     *  Note that inference should return a JPT with variables in the *same* order 
     *  as that specified by Q.
     * @param qvars variables to include in query
     */
    @SuppressWarnings("rawtypes")
    @Override
    public Query makeQuery(Variable... qvars) {
	// Find out which variables in the BN that will be summed out and organise them into "buckets".
        // They will be listed in "topological order" (parents before children) as per heuristics given in Dechter.
        // Each margin variable will be assigned a bucket.
        List<Variable> Q = new ArrayList<>(); // Query, all nodes identified by user of this function
        List<Variable.Assignment> E = new ArrayList<>(); // Assignment, all nodes that are instantiated with values AND relevant (not d-separated from any query node)
        List<Variable> X = new ArrayList<>(); // Un-instantiated but relevant nodes (not d-separated from any query node), to-be summed out
        Q.addAll(Arrays.asList(qvars));
        List<BNode> rnl = null;
        rnl = bn.getDconnected(qvars); //relevant *ordered* node list, based on the concept of D-separation, and topological ordering
        for (BNode node : rnl) {
            Variable var = node.getVariable();
            Object val = node.getInstance();
            if (val != null) {
                E.add(new Variable.Assignment(var, val));
            } else if (!Q.contains(var)) {
                X.add(var);
            }
        }
        return new CGQuery(Q, E, X);
    }
    
    /**
     * Construct the data structure for the specified variables in preparation
     * of inference of the most probable explanation (also known as MAP query).
     * There are three types of variables (given the BN): 
     *  1. Assignment variables E--which have been assigned values via the BN 
     *  2. Query variables Q--typically empty set, but if not, these are not maxed-out but inferred
     *  3. Other variables X--which will be maxed out during inference P(Q|E) = MAX_X P(Q|E,X).
     * Most probable assignments to X can be accessed via CGTable.getMPE.
     * TODO: Investigate if "hybrid" queries, known as marginal MAP (see Koller and Friedman, p. 27) can fit
     * as a variant query under this or if a new query form needs to be implemented. 
     * Marginal MAP means that variables are either "summed out" or 
     * "maxed out" (if assigned values to form part of a diagnosis).
     * @param qvars variables to include in query, can be empty
     */
    @SuppressWarnings("rawtypes")
    public Query makeMPE(Variable... qvars) {
	// Find out which variables in the BN that will be max:ed out and organise them into "buckets".
        // They will be listed in "topological order" (parents before children) as per heuristics given in Dechter.
        // Each margin variable will be assigned a bucket.
        List<Variable> Q = new ArrayList<>(); // Query
        List<Variable.Assignment> E = new ArrayList<>(); // Assignment, all nodes that are instantiated with values AND relevant (not d-separated from any query node)
        List<Variable> X = new ArrayList<>(); // Unspecified, to-be summed out
        Q.addAll(Arrays.asList(qvars));
//        BNet qbn = bn.getRelevant(qvars);
        List<BNode> rnl = bn.getOrdered(); // with MPE all nodes are relevant by definition // getDconnected(qvars); //relevant ordered node list
        //BNet qbn = bn;
        for (BNode node : rnl) {
            Variable var = node.getVariable();
            Object val = node.getInstance();
            if (val != null) {
                E.add(new Variable.Assignment(var, val));
            } else if (!Q.contains(var)) {
                X.add(var);
            }
        }
        CGQuery q = new CGQuery(Q, E, X);
        q.setStatus(STATUS_MPE);
        return q;
    }
    
    /**
     * Construct the data structure for the specified variables in preparation
     * of inference of the most probable explanation (also known as MAP query)
     * but only over nominated variables.
     * There are three types of variables (given the BN): 
     *  1. Assignment variables E--which have been assigned values via the BN 
     *  2. Query variables Q--empty set in this variant of MPE
     *  3. Other variables X--which, IF nominated, will be maxed out during inference
     * Most probable assignments to X can be accessed via CGTable.getMPE.
     * @param nomvars all variables to include
     * @return the query data structure
     */
    @SuppressWarnings("rawtypes")
    public Query makeNominatedMPEQuery(Variable... nomvars) {
	// Find out which variables in the BN that will be max:ed out and organise them into "buckets".
        List<Variable> Q = new ArrayList<>(); // Query, empty in this variant
        List<Variable.Assignment> E = new ArrayList<>(); // Assignment, all nodes that are instantiated with values AND relevant (not d-separated from any query node)
        List<Variable> X = new ArrayList<>(); // Unspecified, to-be summed out
        List<BNode> rnl = bn.getDconnected(nomvars); //relevant *ordered* node list, based on the concept of D-separation, and topological ordering
        for (BNode node : rnl) {
            Variable var = node.getVariable();
            Object val = node.getInstance();
            if (val != null) {
                E.add(new Variable.Assignment(var, val));
            } else {
                X.add(var);
            }
        }
        CGQuery q = new CGQuery(Q, E, X);
        q.setStatus(STATUS_MPE);
        return q;
    }
    
    

    /**
     * Perform exact inference for the specified variables, with support for continuous variables 
     * in the form of conditional Gaussians ({@see bn.FactorTable}).
     * @param query the query (indicating query variable(s))
     */
    @SuppressWarnings("rawtypes")
    @Override
    public QueryResult infer(Query query) {
//        MilliTimer timer = new MilliTimer();
        int[] nprods = new int[3];
	    // All CPTs will be converted to "factors", and put in the bucket which is the first to sum-out any of the variables in the factor.
        // Assignment will be incorporated into the factor when it is constructed.
        // Create list of buckets: first one has query variables, then all sum-outs in topological ordering (as sorted by the constructor)
        CGQuery q = (CGQuery) query;
        Variable[] qarr = new Variable[q.Q.size()];
        q.Q.toArray(qarr);
        //BNet rel_bn = bn.getRelevant(qarr);
        List<Bucket> buckets = new ArrayList<>();
        Bucket first_bucket = new Bucket(q.Q);
        buckets.add(first_bucket);
        for (Variable x : q.X) { // the list of unspecified variables appear in a topological order
            // only create buckets for enumerable variables
            try {
                buckets.add(new Bucket((EnumVariable)x));
            } catch (ClassCastException e) {
                ;
            }
        }
        int nBuckets = buckets.size();

        // Fill buckets backwards with appropriate factor tables (instantiated when "made")
//        timer.start("factors");
        Map<Variable, Object> relmap = q.getRelevant();
        for (Variable var : relmap.keySet()) {
            BNode node = bn.getNode(var);
            // next call is causing delays with threading
            AbstractFactor ft = node.makeDenseFactor(relmap); // forces new makeFactor method to be used on only relevant nodes
            boolean added = false;
            if (!ft.hasEnumVars()) { // // the FT is empty of enumerable variables, hence will only "scale" factors
                buckets.get(0).put(ft); // we will need to keep non-enumerable variables for later though
                added = true;
                continue;
            }
            // go through buckets in reverse order, choosing the first (from end) which "matches" the variables of the FT
            for (int i = nBuckets - 1; i >= 0 && !added; i--) {
                Bucket b = buckets.get(i);
                if (b.match(ft)) {
                    b.put(ft);
                    added = true;
                    break;
                }
            }
            if (!added) { // if not added as per sum-out variable
                // FT is somehow corrupt, e.g. no variables
                throw new VarElimRuntimeException("Node can not be eliminated in inference: " + node.getName());
            }
        }
//        timer.stop("factors");
        // Purge buckets, merge sum-out variables
//        System.out.println(buckets.size());
        List<Integer> remove = new ArrayList<>();
//        timer.start("merge");
        for (int i = 1; i < nBuckets; i++) { // ignore query bucket
//        	try {
//        		Bucket b = buckets.get(i);
//        	} catch (IndexOutOfBoundsException e) {
//        		System.out.println("i = "+i+", nBuckets = "+nBuckets);
//        		System.out.println(buckets.size());
//        	}
            Bucket b = buckets.get(i);
            if (b.factors.isEmpty()) { // no factors, put sum-out variables in other bucket(s)
                for (Variable sumout : b.vars) { // check each sum-out variable
                    for (int jj = i + 1; jj < nBuckets; jj++) { // search suitable bucket for sum-out
                        try {
                            Bucket b2 = buckets.get(jj);
                            if (b2.hasFactorWith(sumout) || jj == nBuckets - 1) { // we've found a bucket with a factor with sum-out variable
                                b2.vars.add(sumout);
                            }
                        } catch (IndexOutOfBoundsException e) {
                            throw new VarElimRuntimeException("Bucket elimination failed during purging and merging: Variable is " + sumout.getName());
                        }
                    }
                }
                remove.add(i);
//                buckets.remove(i); // we can safely remove bucket since variables have been assigned to others
                //FIXME if you remove a bucket the loop fails because nBuckets has not been updated
            }
        }
        buckets.removeAll(remove);
//        timer.stop("merge");
        nBuckets = buckets.size(); // update bucket number
        // Create a factor of each bucket, by performing factor products and marginalisation as appropriate
//        timer.start("products");
        for (int i = nBuckets - 1; i >= 0; i--) {
            Bucket b = buckets.get(i);
//            boolean ignore = true; //  It is safe to ignore buckets with no evidence and no newly computed factor (function of other factors)
//            for (AbstractFactor ft : b.factors) {
//                if (ft.function || ft.evidenced || i == 0) {
//                    ignore = false;
//                    break;
//                }
//            }
            int nFactors = b.factors.size();
//            if (nFactors > 0 && !ignore) {
            if (nFactors > 0) {
                // Perform product of all factors in bucket
                AbstractFactor result = null;
//                timer.start("product");
                if (nFactors == 1) {
                    nprods[0] += 1;
                    result = b.factors.get(0);
                } else if (nFactors == 2) {
                    nprods[1] += 1;
                    result = Factorize.getProduct(b.factors.get(0), b.factors.get(1));
                } else {
                    nprods[2] += 1;
                    AbstractFactor[] fs = new AbstractFactor[b.factors.size()];
                    b.factors.toArray(fs);
                    result = Factorize.getProduct(fs);
                }
//                timer.stop("product");
                if (i > 0) { // not the last bucket, so normal operation
//                    timer.start("margin");
                    // The code below assumes that all buckets except the first have only enumerable variables to be summed out
                    // If continuous variables are unspecified (X) they should have been placed in the first bucket.
                    // If they need summing out, it needs to be done later.
                    try {
                        Variable[] margin = new Variable[b.vars.size()];
                        for (int j = 0; j < margin.length; j ++) 
                            margin[j] = b.vars.get(j);
                        if (q.getStatus() == STATUS_MPE)
                            result = Factorize.getMaxMargin(result, margin); // max-out variables of bucket
                        else
                            result = Factorize.getMargin(result, margin);    // sum-out variables of bucket
                        
                        if (!result.hasEnumVars())      // if no enumerable variables, we may still have non-enumerables
                            buckets.get(0).put(result); // so we put the factor in the first bucket
                        else {                          // there are enumerables so...
                            for (int jj = i - 1; jj >= 0; jj--) { // find a new bucket for the result factor
                                Bucket b2 = buckets.get(jj);      
                                if (b2.match(result)) {
                                    b2.put(result);
                                    break;
                                }
                            }
                        }
                    } catch (ClassCastException e) {
                        throw new VarElimRuntimeException("Cannot marginalize or maximize-out continuous variables");
                    }
//                    timer.stop("margin");
                } else {    
                    // This is the final (first) bucket so we should not marginalize/maximize out query variables (and non-enumerables), 
                    // instead we should extract query results from the final factor, including a JPT.
                    // If Q is only a non-enumerable variable or list there-of, we will not be able to create a JPT.
                    // The first section below is just making sure that the variables are presented in the same order as that in the query
//                    timer.stop("products");
//                    System.out.println("Products 1x " + nprods[0] + " 2x " + nprods[1] + " nx " + nprods[2]);
//                    timer.report(true);
                    return new CGTable(result, q.Q);
                }
            }
        }
        throw new VarElimRuntimeException("Variable elimination failed");
    }
    
    @SuppressWarnings("rawtypes")
    public QueryResult infer(Variable[] query_vars) {
        return infer(makeQuery(query_vars));
    }

    public QueryResult infer(Variable query_var) {
        return infer(new Variable[]{query_var});
    }

    public QueryResult infer(BNode[] query_nodes) {
        Variable[] vars = new Variable[query_nodes.length];
        for (int i = 0; i < vars.length; i++) {
            vars[i] = query_nodes[i].getVariable();
        }
        return infer(makeQuery(vars));
    }

    public QueryResult infer(BNode query_node) {
        return infer(new BNode[]{query_node});
    }

    /**
     * Determine the probability of the instantiated variables. 
     * This function does not integrate over/marginalize continuous (and other density-based) variables. 
     * Hence, such variables may have influence the result if they are not instantiated.
     * @return the log likelihood of the evidence (instantiated nodes)
     */
    public double logLikelihood() {
	// All CPTs will be converted to "factors", and put in the bucket which is the first to sum-out any of the variables in the factor.
        // Assignment will be incorporated into the factor when it is constructed.
        List<Variable> X = new ArrayList<>(); // unspecified variables, to-be summed-out
        Map<Variable, Object> R = new HashMap<>(); // all variables that are relevant with corresponding instantiations
        for (BNode node : bn.getOrdered()) {
            Variable var = node.getVariable();
            Object instance = node.getInstance();
            R.put(var, instance); // currently we consider all variables are relevant, even when not specified
            if (instance == null) 
                X.add(var);
        }
        List<Bucket> buckets = new ArrayList<>();
        for (Variable x : X) { 
            // only create buckets for enumerable, unspecified variables to-be marginalised-out
            try {
                buckets.add(new Bucket((EnumVariable)x));
            } catch (ClassCastException e) {
            }
        }
        int nBuckets = buckets.size();
        if (nBuckets == 0) // if no variables are unknown, then we need to 
            buckets.add(new Bucket((EnumVariable)null));
        // Fill buckets backwards with appropriate factor tables (instantiated when "made")
        for (BNode node : bn.getNodes()) {
            // node is converted into a factor, all nodes are considered relevant
            AbstractFactor ft = node.makeDenseFactor(R); // Factor ft = node.makeFactor(bn, true);
            boolean added = false;
            if (!ft.hasEnumVars()) { // this happen if the FT has no variables
                buckets.get(0).put(ft);
                added = true;
                continue;
            }
            // go through buckets in reverse order, choosing the first (from end) which "matches" the variables of the FT
            for (int i = nBuckets - 1; i >= 0 && !added; i--) { 
                Bucket b = buckets.get(i);
                if (b.match(ft)) {
                    b.put(ft);
                    added = true;
                }
            }
            if (!added) { // if not added as per sum-out variable, so put in penultimate bucket
                buckets.get(0).put(ft);
            }
        }
        // Purge buckets, merge sum-out variables
        //FIXME error same as above
        for (int i = 1; i < nBuckets; i++) { // do not re-organise penultimate bucket
            Bucket b = buckets.get(i);
            if (b.factors.isEmpty()) { // no factors, put sum-out variables in other bucket(s)
                for (Variable sumout : b.vars) { // check each sum-out variable
                    for (int jj = i + 1; jj < nBuckets; jj++) { // search suitable bucket for sum-out
                        Bucket b2 = buckets.get(jj);
                        if (b2.hasFactorWith(sumout) || jj == nBuckets - 1) { // we've found a bucket with a factor with sum-out variable
                            b2.vars.add(sumout);
                        }
                    }
                }
                buckets.remove(i); // we can safely remove bucket since variables have been assigned to others
            }
        }
        nBuckets = buckets.size(); // update bucket number
        // Create a factor of each bucket, by performing factor products and marginalisation as appropriate
        for (int i = nBuckets - 1; i >= 0; i--) {
            Bucket b = buckets.get(i);
//            boolean ignore = true; //  It is safe to ignore buckets with no evidence and no newly computed factor (function of other factors)
//            for (Factor ft : b.factors) {
//                if (ft.function || ft.evidenced || i == 0) {
//                    ignore = false;
//                    break;
//                }
//            }
            int nFactors = b.factors.size();
//            if (nFactors > 0 && !ignore) {
            if (nFactors > 0) {
                // Perform product of all factors in bucket
                AbstractFactor result = null;
                if (nFactors == 1) {
                    result = b.factors.get(0);
                } else if (nFactors == 2) { 
                    result = Factorize.getProduct(b.factors.get(0), b.factors.get(1));
                } else {
                    AbstractFactor[] fs = new AbstractFactor[b.factors.size()];
                    b.factors.toArray(fs);
                    result = Factorize.getProduct(fs);
                }
                if (i > 0) { // not the last bucket, so normal operation 
                    // The code below assumes that all buckets except the first have only enumerable variables to be summed out
                    // If continuous variables are unspecified (X) they should have been placed in the first bucket.
                    // If they need summing out, it needs to be done later.
                    try {
                        Variable[] margin = new Variable[b.vars.size()];
                        for (int j = 0; j < margin.length; j ++) 
                            margin[j] = b.vars.get(j);
                        result = Factorize.getMargin(result, margin);   // sum-out variables of bucket
                        for (int jj = i - 1; jj >= 0; jj--) { // find a new bucket for the result factor
                            Bucket b2 = buckets.get(jj);      // FIXME: problem if FT is atomic (ie no variables)
                            if (b2.match(result)) {
                                b2.put(result);
                                break;
                            }
                        }
                    } catch (ClassCastException e) {
                        throw new VarElimRuntimeException("Cannot marginalize continuous variables");
                    }
                } else {    
                    // This is the final (first) bucket 
                    double sum = result.getLogSum();
                    return sum;
                }
            }
        }
        throw new VarElimRuntimeException("Exited variable elimination prematurely");
    }
    
    public class FTCompare implements Comparator<Factor> {
        @Override
        public int compare(Factor o1, Factor o2) {
            return o1.getEnumVariables().size() - o2.getEnumVariables().size();
        }
    }

    class CGQuery implements Query {
        // Note only *relevant* evidence and unspecified variables are included here
        final List<Variable> Q; // query variables, for which values are sought
        final Map<Variable, Object> E; // evidence variables; their nodes are instantiated to values
        final List<Variable> X; // unspecified variables, incl all not listed as query or evidence
        private int status = STATUS_BEL;
        
        CGQuery(List<Variable> Q, List<Variable.Assignment> E, List<Variable> X) {
            this.Q = Q;
            this.E = new HashMap<>();
            // evidence variables passed to this method are relevant
            for (Variable.Assignment assign : E) 
                this.E.put(assign.var, assign.val);
            // all query variables are relevant by default, but no values are assigned
            for (Variable q : Q)
                this.E.put(q, null);
            // unspecified variables passed to this method are relevant, but no values are assigned
            for (Variable x : X)
                this.E.put(x, null);
            this.X = X;
        }
        
        void setStatus(int status) {
            this.status = status;
        }
        
        int getStatus() {
            return this.status;
        }
        
        /**
         * Get relevant variables for the query as a map, with evidence if available.
         * @return a map of all relevant variables with assigned values
         */
        Map<Variable, Object> getRelevant() {
            return E;
        }
        
        public String toString() {
            StringBuilder sbuf = new StringBuilder();
            sbuf.append("Q:");
            for (Variable v:Q)
                sbuf.append(v.getName()).append(",");
            sbuf.append("|E:");
            for (Map.Entry<Variable, Object> v:E.entrySet()) {
                if (v.getValue() != null)
                    sbuf.append(v.getKey().toString()).append("=").append(v.getValue().toString()).append(",");
            }
            sbuf.append("|X:");
            for (Variable v:X)
                sbuf.append(v.getName()).append(",");
            sbuf.append("|Status:").append(status);
            return sbuf.toString();
        }
    }

    
    /**
     * A bucket holds factor tables and one variable
     *
     * @author mikael
     */
    public class Bucket {

        List<AbstractFactor> factors;
        List<Variable> vars = new ArrayList<>();

        /**
         * Create a bucket by identifying the variable it processes
         *
         * @param var the variable
         */
        Bucket(Variable var) {
            this.vars.add(var);
            this.factors = new ArrayList<>();
        }

        /**
         * Create a bucket by identifying the variables it processes
         *
         * @param vars the variables
         */
        Bucket(List<Variable> vars) {
            this.vars = vars;
            this.factors = new ArrayList<>();
        }

        /**
         * Decides if the factor table belongs to this bucket (not exclusively
         * though).
         *
         * @param f factor table
         * @return true if the bucket can be used to process the factor table,
         * false otherwise
         */
        boolean match(AbstractFactor f) {
            for (EnumVariable fvar : f.getEnumVars()) {
                if (vars.contains(fvar)) 
                    return true;
            }
            // currently the matching does NOT include continuous variables since they are not summed out anyway
            // for (Variable nvar : nvars) {
            //    if (vars.contains(nvar)) 
            //        return true;
            // }
            return false;
        }

        boolean hasFactorWith(Variable var) {
            for (AbstractFactor ft : factors) {
                if (ft.hasVariable(var)) 
                    return true;
            }
            return false;
        }

        void put(AbstractFactor f) {
            factors.add(f);
        }

        List<AbstractFactor> get() {
            return factors;
        }
    }
    
    public class VarElimRuntimeException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        public VarElimRuntimeException(String message) {
            super(message);
        }
    }

}
