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
import bn.EnumVariable;
import bn.Factor;
import bn.JPT;
import bn.Variable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Exact inference in Bayesian network by variable elimination, more
 * specifically a variant of "Bucket Elimination" that involves conditional Gaussians. 
 * Implementation of belief updating in accordance with the method described in 
 * Dechter, R. Bucket Elimination: A Unifying Framework for Probabilistic Inference, 
 * in Uncertainty in Artificial Intelligence, 1998.
 * The idea of conditional Gaussians is discussed in Lauritzen SL and Jensen F, 
 * Stable local computation with conditional Gaussian distributions, Statistics and
 * Computing 2001 11:191-203.
 * 
 * @author mikael
 */
public class CGVarElim implements Inference {

    public BNet bn;

    @Override
    public void instantiate(BNet bn) {
        this.bn = bn;
        this.bn.compile();
    }

    /**
     * Construct the data structure for the specified variables in preparation
     * of inference. There are three types of variables (given the BN): 
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
        // Each sumout variable will be assigned a bucket.
        List<Variable> Q = new ArrayList<>(); // Query
        List<Variable> E = new ArrayList<>(); // Assignment
        List<Variable> X = new ArrayList<>(); // Unspecified, to-be summed out
        Q.addAll(Arrays.asList(qvars));
        BNet qbn = bn.getRelevant(qvars);
        //BNet qbn = bn;
        for (BNode node : qbn.getOrdered()) {
            Variable var = node.getVariable();
            if (node.getInstance() != null) {
                E.add(var);
            } else if (!Q.contains(var)) {
                X.add(var);
            }
        }
        return new CGQuery(Q, E, X, qbn);
    }

    /**
     * Perform exact inference for the specified variables, with support for continuous variables 
     * in the form of conditional Gaussians ({@see bn.FactorTable}).
     * @param query the query (indicating query variable(s))
     */
    @SuppressWarnings("rawtypes")
    @Override
    public QueryResult infer(Query query) {
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
        for (Variable x : q.X) {
            // only create buckets for enumerable variables
            try {
                buckets.add(new Bucket((EnumVariable)x));
            } catch (ClassCastException e) {
                ;
            }
        }        
        int nBuckets = buckets.size();
        // Fill buckets backwards with appropriate factor tables (instantiated when "made")
        for (BNode node : q.qbn.getNodes()) {
            Factor ft = node.makeFactor(q.qbn);
            boolean added = false;
            if (ft.isAtomic()) { // // the FT is empty of enumerable variables, hence will only "scale" factors
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
                throw new CGVarElimRuntimeException("Node can not be eliminated in inference: " + node.getName());
            }
        }
        // Purge buckets, merge sum-out variables
        for (int i = 1; i < nBuckets; i++) { // ignore query bucket
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
                            throw new CGVarElimRuntimeException("Bucket elimination failed during purging and merging: Variable is " + sumout.getName());
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
            boolean ignore = true; //  It is safe to ignore buckets with no evidence and no newly computed factor (function of other factors)
            for (Factor ft : b.factors) {
                if (ft.function || ft.evidenced || i == 0) {
                    ignore = false;
                    break;
                }
            }
            int nFactors = b.factors.size();
            if (nFactors > 0 && !ignore) {
                List<Factor> fts = new ArrayList<>(b.factors);
                // The order in which FTs are multiplied can have big influence on efficiency, the strategy below is naive
                Collections.sort(fts, new FTCompare()); // sort the factors in order of variable-count (smaller-to-greater)				
                Factor result = fts.get(0); // perform products in that order
                for (int j = 1; j < nFactors; j++) {
                    result = Factor.product(result, fts.get(j));
                }
                if (i > 0) { // not the last bucket, so normal operation 
                    // The code below assumes that all buckets except the first have only enumerable variables to be summed out
                    // If continuous variables are unspecified (X) they should have been placed in the first bucket.
                    // If they need summing out, it needs to be done later.
                    try {
                        List<EnumVariable> evars = new ArrayList<>(b.vars.size());
                        for (Variable bvar : b.vars) 
                            evars.add((EnumVariable)bvar);
                        result = result.marginalize(evars);   // sum-out variables of bucket
                        if (result.isAtomic())          // if no enumerable variables, we may still have non-enumerables
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
                        throw new CGVarElimRuntimeException("Cannot marginalize continuous variables");
                    }
                } else {    
                    // This is the final (first) bucket so we should not marginalize out query variables (and non-enumerables), 
                    // instead we should extract query results from the final factor, including a JPT.
                    // If Q is only a non-enumerable variable or list there-of, we will not be able to create a JPT.
                    // The first section below is just making sure that the variables are presented in the same order as that in the query
                    
                    //return extractResult(result, q.Q);
                    return new CGTable(result.rehash(q.Q));
                }
            }
        }
        throw new CGVarElimRuntimeException("Variable elimination failed");
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
     * Determine the log-transformed probability of the instantiated variables. 
     * This function does not integrate over/marginalize continuous variables in case they are not instantiated.
     * @return the likelihood of the evidence (instantiated nodes)
     */
    public double likelihood() {
        JPT answer = null;
	// All CPTs will be converted to "factors", and put in the bucket which is the first to sum-out any of the variables in the factor.
        // Assignment will be incorporated into the factor when it is constructed.
        List<Variable> E = new ArrayList<>(); // evidence variables
        List<Variable> X = new ArrayList<>(); // unspecified variables, to-be summed-out
        for (BNode node : bn.getOrdered()) {
            Variable var = node.getVariable();
            if (node.getInstance() != null) 
                E.add(var);
            else
                X.add(var);
        }
        List<Bucket> buckets = new ArrayList<>();
        for (Variable x : X) {
            // only create buckets for enumerable, unspecified variables to-be summed-out
            try {
                buckets.add(new Bucket((EnumVariable)x));
            } catch (ClassCastException e) {
                ;
            }
        }
        int nBuckets = buckets.size();
        if (nBuckets == 0) // if no variables are unknown, then we need to 
            buckets.add(new Bucket((EnumVariable)null));
        // Fill buckets backwards with appropriate factor tables (instantiated when "made")
        for (BNode node : bn.getNodes()) {
            Factor ft = node.makeFactor(bn);
            boolean added = false;
            if (ft.isAtomic()) { // this happen if the FT has no variables
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
            boolean ignore = true; //  It is safe to ignore buckets with no evidence and no newly computed factor (function of other factors)
            for (Factor ft : b.factors) {
                if (ft.function || ft.evidenced || i == 0) {
                    ignore = false;
                    break;
                }
            }
            int nFactors = b.factors.size();
            if (nFactors > 0 && !ignore) {
                List<Factor> fts = new ArrayList<>(b.factors);
                // The order in which FTs are multiplied can have big influence on efficiency, the strategy below is naive
                Collections.sort(fts, new FTCompare()); // sort the factors in order of variable-count (smaller-to-greater)				
                Factor result = fts.get(0); // perform products in that order
                for (int j = 1; j < nFactors; j++) {
                    result = Factor.product(result, fts.get(j));
                }
                if (i > 0) { // not the last bucket, so normal operation 
                    // The code below assumes that all buckets except the first have only enumerable variables to be summed out
                    // If continuous variables are unspecified (X) they should have been placed in the first bucket.
                    // If they need summing out, it needs to be done later.
                    try {
                        List<EnumVariable> evars = new ArrayList<>(b.vars.size());
                        for (Variable bvar : b.vars) 
                            evars.add((EnumVariable)bvar);
                        result = result.marginalize(evars);   // sum-out variables of bucket
                        for (int jj = i - 1; jj >= 0; jj--) { // find a new bucket for the result factor
                            Bucket b2 = buckets.get(jj);      // FIXME: problem if FT is atomic (ie no variables)
                            if (b2.match(result)) {
                                b2.put(result);
                                break;
                            }
                        }
                    } catch (ClassCastException e) {
                        throw new CGVarElimRuntimeException("Cannot marginalize continuous variables");
                    }
                } else {    
                    // This is the final (first) bucket 
                    return result.getSum();
                }
            }
        }
        throw new CGVarElimRuntimeException("Exited variable elimination prematurely");
    }
    
    public class FTCompare implements Comparator<Factor> {
        @Override
        public int compare(Factor o1, Factor o2) {
            return o1.getEnumVariables().size() - o2.getEnumVariables().size();
        }
    }

    public class CGQuery implements Query {
        final List<Variable> Q;
        final List<Variable> E;
        final List<Variable> X;
        final BNet qbn;
        CGQuery(List<Variable> Q, List<Variable> E, List<Variable> X, BNet qbn) {
            this.Q = Q;
            this.E = E;
            this.X = X;
            this.qbn = qbn;
        }
    }

    
    /**
     * A bucket holds factor tables and one variable
     *
     * @author mikael
     */
    public class Bucket {

        List<Factor> factors;
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
        boolean match(Factor f) {
            List<EnumVariable> fvars = f.getEnumVariables();
            List<Variable> nvars = f.getNonEnumVariables();
            for (EnumVariable fvar : fvars) {
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
            for (Factor ft : factors) {
                if (ft.getEnumVariables().contains(var)) 
                    return true;
                if (ft.getNonEnumVariables().contains(var))
                    return true;
            }
            return false;
        }

        void put(Factor f) {
            factors.add(f);
        }

        List<Factor> get() {
            return factors;
        }
    }
    
    public class CGVarElimRuntimeException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        public CGVarElimRuntimeException(String message) {
            super(message);
        }
    }

}
