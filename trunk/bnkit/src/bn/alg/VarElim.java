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
import bn.CountTable;
import bn.EnumTable;
import bn.EnumVariable;
import bn.FactorTable;
import bn.JPT;
import bn.Variable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Exact inference in Bayesian network by variable elimination, more
 * specifically "Bucket Elimination". Implementation of belief updating in
 * accordance with the method described in Dechter, R. Bucket Elimination: A
 * Unifying Framework for Probabilistic Inference, in Uncertainty in Artificial
 * Intelligence, 1998.
 *
 * @author mikael
 */
public class VarElim implements Inference {

    public BNet bn;
    private double logLikelihood = 0;

    @Override
    public void instantiate(BNet bn) {
        this.bn = bn;
        this.bn.compile();
    }

    /**
     * Construct the data structure for the specified variables in preparation
     * of inference. There are three types of variables (given the BN): 1.
     * Evidence variables E--which have been assigned values via the BN 2. Query
     * variables Q--for which probabilities are sought P(Q|E) 3. Other variables
     * X--which will be summed out during inference P(Q|E) = SUM_X P(Q|E,X).
     * Note that inference should return a JPT with variables in the *same* order 
     * as that specified by Q.
     */
    @SuppressWarnings("rawtypes")
    @Override
    public Query makeQuery(Variable... qvars) {
	// Find out which variables in the BN that will be summed out and organise them into "buckets".
        // They will be listed in "topological order" (parents before children) as per heuristics given in Dechter.
        // Each sumout variable will be assigned a bucket.
        List<EnumVariable> Q = new ArrayList<EnumVariable>();
        List<Variable> E = new ArrayList<Variable>();
        List<EnumVariable> X = new ArrayList<EnumVariable>();
        List<Variable> ordered = new ArrayList<Variable>();
        try {
            for (int i = 0; i < qvars.length; i++) {
                Q.add((EnumVariable) qvars[i]);
            }
            for (BNode node : bn.getOrdered()) {
                Variable var = node.getVariable();
                ordered.add(var);
                if (node.getInstance() != null) {
                    E.add(var);
                } else if (!Q.contains(var)) {
                    X.add((EnumVariable) var);
                }
            }
        } catch (RuntimeException e) {
            throw new RuntimeException("Cannot perform variable elimination on non-enumerable nodes");
        }
        return new VEQuery(Q, E, X);
    }

    /**
     * Perform exact inference for the specified variables.
     */
    @SuppressWarnings("rawtypes")
    @Override
    public JPT infer(Query query) {
        JPT answer = null;
	// All CPTs will be converted to "factors", and put in the bucket which is the first to sum-out any of the variables in the factor.
        // Evidence will be incorporated into the factor when it is constructed.
        // Create list of buckets: first one has query variables, then all sum-outs in topological ordering
        VEQuery q = (VEQuery) query;
        List<Bucket> buckets = new ArrayList<Bucket>();
        buckets.add(new Bucket(q.Q));
        for (EnumVariable x : q.X) {
            buckets.add(new Bucket(x));
        }
        int nBuckets = buckets.size();
        // Fill buckets backwards with appropriate factor tables (instantiated when "made")
        for (BNode node : bn.getNodes()) {
            FactorTable ft = node.makeFactor(bn);
            if (ft == null) // the node is eliminated
            {
                continue;
            }
            boolean added = false;
            for (int i = nBuckets - 1; i >= 0; i--) {
                Bucket b = buckets.get(i);
                if (b.match(ft)) {
                    b.put(ft);
                    added = true;
                    break;
                }
            }
            if (!added) // if not added as per sum-out variable 
            {
                if (ft.getParents().size() > 0) {
                    throw new RuntimeException("Node can not be eliminated in inference: " + node.getName());
                } // put it in the last bucket?
                //buckets.get(buckets.size() - 1).put(ft);
                else
                    ; // all variables had been factored out so ignore
            }
        }
        // Purge buckets, merge sum-out variables
        for (int i = 1; i < nBuckets; i++) { // ignore query bucket
            Bucket b = buckets.get(i);
            if (b.factors.isEmpty()) { // no factors, put sum-out variables in other bucket(s)
                for (EnumVariable sumout : b.vars) { // check each sum-out variable
                    for (int jj = i + 1; jj < nBuckets; jj++) { // search suitable bucket for sum-out
                        Bucket b2 = buckets.get(jj);
                        if (b2.hasFactorWith(sumout) || jj == nBuckets - 1) // we've found a bucket with a factor with sum-out variable
                        {
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
            for (FactorTable ft : b.factors) {
                if (ft.function || ft.evidenced || i == 0) {
                    ignore = false;
                    break;
                }
            }
            int nFactors = b.factors.size();
            if (nFactors > 0 && !ignore) {
                List<FactorTable> fts = new ArrayList<>(b.factors);
                Collections.sort(fts, new FTCompare()); // sort the factors in order of variable-count (smaller-to-greater)				
                FactorTable result = fts.get(0); // perform products in that order
                for (int j = 1; j < nFactors; j++) {
                    result = FactorTable.product(result, fts.get(j));
                }
                if (i > 0) { // not the last bucket, so normal operation 
                    result = result.marginalize(b.vars); // sum-out variables of bucket
                    for (int jj = i - 1; jj >= 0; jj--) { // find a new bucket for the result factor
                        Bucket b2 = buckets.get(jj);
                        if (b2.match(result)) {
                            b2.put(result);
                            break;
                        }
                    }
                } else { // last bucket so we should not marginalize out query variables, instead return the JPT from the final factor
                    List<EnumVariable> f_parents = result.getParents();
                    List<EnumVariable> q_parents = q.Q;
                    int[] map2q = new int[q_parents.size()];
                    for (int jf = 0; jf < map2q.length; jf ++) {
                        map2q[jf] = -1;
                        for (int jq = 0; jq < map2q.length; jq ++) {
                            if (f_parents.get(jf).getName().equals(q_parents.get(jq).getName()))  {
                                map2q[jf] = jq;
                                break;
                            }
                        }
                    }
                    for (int j = 0; j < map2q.length; j ++) {
                        if (map2q[j] == -1)
                            throw new RuntimeException("Invalid inference result");
                    }
                    EnumTable<Double> et = new EnumTable<Double>(q_parents);
                    for (Map.Entry<Integer, Double> entry : result.getMapEntries()) {
                        Object[] fkey = result.getKey(entry.getKey().intValue());
                        Object[] qkey = new Object[fkey.length];
                        for (int j = 0; j < fkey.length; j ++)
                            qkey[map2q[j]] = fkey[j];
                        et.setValue(qkey, entry.getValue());
                    }
                    answer = new JPT(et);
                    logLikelihood = result.getLogLikelihood();
                    return answer;
                }
            }
        }
        logLikelihood = 0.0;
        return answer;
    }

    @SuppressWarnings("rawtypes")
    public JPT infer(Variable[] query_vars) {
        return infer(makeQuery(query_vars));
    }

    public JPT infer(Variable query_var) {
        return infer(new Variable[]{query_var});
    }

    public JPT infer(BNode[] query_nodes) {
        Variable[] vars = new Variable[query_nodes.length];
        for (int i = 0; i < vars.length; i++) {
            vars[i] = query_nodes[i].getVariable();
        }
        return infer(makeQuery(vars));
    }

    public JPT infer(BNode query_node) {
        return infer(new BNode[]{query_node});
    }

    public class FTCompare implements Comparator<FactorTable> {

        @Override
        public int compare(FactorTable o1, FactorTable o2) {
            return o1.getParents().size() - o2.getParents().size();
        }
    }

    public class VEQuery implements Query {

        final List<EnumVariable> Q;
        final List<Variable> E;
        final List<EnumVariable> X;

        VEQuery(List<EnumVariable> Q, List<Variable> E, List<EnumVariable> X) {
            this.Q = Q;
            this.E = E;
            this.X = X;
        }
    }

    /**
     * A bucket holds factor tables and one variable
     *
     * @author mikael
     */
    public class Bucket {

        List<FactorTable> factors;
        List<EnumVariable> vars = new ArrayList<EnumVariable>();

        /**
         * Create a bucket by identifying the variable it processes
         *
         * @param var the variable
         */
        Bucket(EnumVariable var) {
            this.vars.add(var);
            this.factors = new ArrayList<FactorTable>();
        }

        /**
         * Create a bucket by identifying the variable it processes
         *
         * @param var the variable
         */
        Bucket(List<EnumVariable> vars) {
            this.vars = vars;
            this.factors = new ArrayList<FactorTable>();
        }

        /**
         * Decides if the factor table belongs to this bucket (not exclusively
         * though).
         *
         * @param f factor table
         * @return true if the bucket can be used to process the factor table,
         * false otherwise
         */
        boolean match(FactorTable f) {
            List<EnumVariable> fvars = f.getParents();
            for (EnumVariable fvar : fvars) {
                if (vars.contains(fvar)) {
                    return true;
                }
            }
            return false;
        }

        boolean hasFactorWith(EnumVariable var) {
            for (FactorTable ft : factors) {
                if (ft.getParents().contains(var)) {
                    return true;
                }
            }
            return false;
        }

        void put(FactorTable f) {
            factors.add(f);
        }

        List<FactorTable> get() {
            return factors;
        }
    }

    @Override
    public double getLogLikelihood() {
        return logLikelihood;
    }
}
