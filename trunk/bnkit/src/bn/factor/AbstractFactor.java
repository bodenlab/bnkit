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
package bn.factor;

import bn.Distrib;
import bn.EnumVariable;
import bn.GaussianDistrib;
import bn.JDF;
import bn.Predef;
import bn.Variable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Table for storing and retrieving doubles based on Enumerable keys.
 * 
 * The intended use is for "factors", representing the factorisation of (conditional) probabilities.
 * They form the backbone data structure in algorithms such as variable elimination and 
 * message-passing in belief propagation. These operations need to be efficient. 
 * 
 * AbstractFactor is an abstract class that defines methods for specific implementations but 
 * implements basic operations that are used externally on factors, including product, and 
 * variable marginalisation (summing-out or maxing-out).
 * 
 * The class is sensitive to the order of variables, and exposes the user to some
 * details so caution must be exercised.
 * 
 * TODO: Consider improving efficiency further to exploit the fact that now all variables are sorted
 * in the constructor. (Currently, some code does not assume order.)
 * TODO: Choose type of implementation when new factors are constructed as a result of operations,
 * informed of their likely requirements.
 *
 * @author mikael
 */
public abstract class AbstractFactor {

    protected static int PRODUCT_OPTION = -1; // choose strategy for complex cases by timing ("-1") or by fixed option (currently "0" and "1")
    protected final int nEVars; // number of enumerable variables
    protected static boolean VERBOSE = false;
    protected final int nNVars; // number of non-enumerable variables

    protected final EnumVariable[] evars; // enumerable variables, which form the "keys" to entries in the map
    protected final Variable[] nvars; // non-enumerable variables, which densities are attached to the entries via JDF
    protected final int[] period;
    protected final int[] step;
    protected final int[] domsize; // size of domain
    
    public boolean evidenced = false;

    /**
     * Construct a new table without any variables. 
     * This type of factor is used when all variables are summed out. 
     * They can appear as part of products to "scale" its opposite.
     */
    protected AbstractFactor() {
        this.evars = null;
        this.nEVars = 0;
        this.nvars = null;
        this.nNVars = 0;
        this.step = null;
        this.period = null;
        this.domsize = null;
    }

    /**
     * Construct a new table with the specified variables.
     * Enumerable variables will form keys to index the entries, 
     * non-enumerables will be associated in the form of their densities with
     * each entry. The order in which the variables are given is significant.
     * The internal order between enumerable variables is maintained, and can
     * impact on the efficiency of various operations including product.
     * As a rule, expect that presenting variables in the same order between
     * different tables will be beneficial.
     * The order of non-enumerable variables plays no role, however.
     * 
     * @param useVariables variables, either enumerable or non-enumerable, 
     * potentially unsorted and redundant
     */
    protected AbstractFactor(Variable... useVariables) {
        // sort and weed out duplicates
        Variable[] uniqueVars = sortRemoveDuplicates(useVariables);
        // count the number of enumerable and non-enumerable variables
        int cnt_evars = 0, cnt_nvars = 0; 
        for (Variable var : uniqueVars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                cnt_evars ++;
            } catch (ClassCastException e) {
                cnt_nvars ++;
            }
        }
        // hybrid table: both enumerable and non-enumerable variables
        if (cnt_nvars > 0 && cnt_evars > 0) { 
            this.evars = new EnumVariable[cnt_evars];
            this.nvars = new Variable[cnt_nvars];
            cnt_evars = 0; cnt_nvars = 0;
            for (Variable uvar : uniqueVars) {
                try {
                    EnumVariable evar = (EnumVariable) uvar;
                    this.evars[cnt_evars ++] = evar;
                } catch (ClassCastException e) {
                    this.nvars[cnt_nvars ++] = uvar;
                }
            }
            this.nEVars = evars.length;
            this.nNVars = nvars.length;
        } else if (cnt_evars > 0) { // only enumerable variables
            this.nNVars = 0;
            this.nvars = null;
            this.evars = new EnumVariable[cnt_evars];
            System.arraycopy(uniqueVars, 0, this.evars, 0, cnt_evars);
            this.nEVars = cnt_evars;
        } else if (cnt_nvars > 0) { // only non-enumerable variables
            this.nEVars = 0;
            this.evars = null;
            this.nvars = new Variable[cnt_nvars];
            System.arraycopy(uniqueVars, 0, this.nvars, 0, cnt_nvars);
            this.nNVars = cnt_nvars;
        } else { // no variables
            this.evars = null;
            this.nEVars = 0;
            this.nvars = null;
            this.nNVars = 0;
        }
        // now determine the indexing for quick access to entries based on enumerable variables
        this.step = new int[this.nEVars];
        this.period = new int[this.nEVars];
        this.domsize = new int[this.nEVars];
        int prod = 1;
        for (int i = 0; i < nEVars; i++) {
            int parent = nEVars - i - 1;
            this.domsize[parent] = evars[parent].size();
            this.step[parent] = prod;
            prod *= this.domsize[parent];
            this.period[parent] = prod;
        }
    }

    /**
     * Construct an array of variables sorted according to their canonical index,
     * that contains all of the unique variables in the provided array.
     * This method is used by the constructor of this class.
     * @param A potentially unsorted array of variables, potentially with duplicates
     * @return array of the same variables, but sorted and unique
     */
    protected static Variable[] sortRemoveDuplicates(Variable[] A) {
        if (A.length < 2) {
            return A;
        }
        int j = 0;
        int i = 1;
        try { 
            Arrays.sort(A);
        } catch (java.lang.NullPointerException e) {
            System.out.println();
        }
        while (i < A.length) {
            if (A[i] == A[j]) {
                i++;
            } else {
                j++;
                A[j] = A[i];
                i++;
            }
        }
        return Arrays.copyOf(A, j + 1);
    }

    /**
     * Gauges the alignment of enumerable variables in two tables to enable
     * optimization of products, in particular.
     * Assumes that variables are ordered per their canonical index.
     * @param Xvars enumerable variables of a factor X
     * @param Yvars enumerable variables of a factor Y
     * @return number of enumerable variables that overlap
     */
    public static int getOverlap(EnumVariable[] Xvars, EnumVariable[] Yvars) {
        if (Xvars == null) Xvars = new EnumVariable[0];
        if (Yvars == null) Yvars = new EnumVariable[0];
        int x = 0;
        int y = 0;
        int overlap = 0;
        while (Xvars.length > x && Yvars.length > y) {
            if (Xvars[x] == Yvars[y]) {
                x++;
                y++;
                overlap++;
            } else if (Xvars[x].getCanonicalIndex() < Yvars[y].getCanonicalIndex()) {
                x++;
            } else {
                y++;
            }
        }
        return overlap;
    }

    /**
     * Gauges the alignment of enumerable variables in two tables to enable
     * optimization of products, in particular.
     * Assumes that variables are ordered per their canonical index.
     * @param X factor X
     * @param Y factor Y
     * @return number of enumerable variables that overlap
     */
    public static int getOverlap(AbstractFactor X, AbstractFactor Y) {
        int x = 0;
        int y = 0;
        int overlap = 0;
        while (X.nEVars > x && Y.nEVars > y) {
            if (X.evars[x] == Y.evars[y]) {
                x++;
                y++;
                overlap++;
            } else if (X.evars[x].getCanonicalIndex() < Y.evars[y].getCanonicalIndex()) {
                x++;
            } else {
                y++;
            }
        }
        return overlap;
    }

    /**
     * Gauges the computational cost of multiplying the two tables, to enable
     * optimization of the order in which to products. 
     * This "gauge" does not work well at the moment for structuring a binary tree.
     * 
     * Assumes that variables are ordered per their canonical index.
     * Note that the cost of joining shared variables is optional.
     *
     * @param Xvars enumerable variables of a factor X
     * @param Yvars enumerable variables of a factor Y
     * @param includeJoin set to true if operations that involved share variables should be counted, false otherwise
     * @return complexity of the product X * Y
     */
    protected static int getComplexity(EnumVariable[] Xvars, EnumVariable[] Yvars, boolean includeJoin) {
        if (Xvars == null) Xvars = new EnumVariable[0];
        if (Yvars == null) Yvars = new EnumVariable[0];
        int x = 0;
        int y = 0;
        int multiplier = 1;
        while (Xvars.length > x && Yvars.length > y) {
            if (Xvars[x] == Yvars[y]) {
                if (includeJoin)
                    multiplier *= Xvars[x].size();
                x++;
                y++;
            } else if (Xvars[x].getCanonicalIndex() < Yvars[y].getCanonicalIndex()) {
                multiplier *= Xvars[x].size();
                x++;
            } else {
                multiplier *= Yvars[y].size();
                y++;
            }
        }
        while (Xvars.length > x) {
            multiplier *= Xvars[x++].size();
        }
        while (Yvars.length > y) {
            multiplier *= Yvars[y++].size();
        }
        return multiplier;
    }

    /**
     * Gauges the computational cost of multiplying the two tables, to enable
     * optimization of the order in which to products. 
     * Assumes that variables are ordered per their canonical index.
     * Note that the cost of joining shared variables is optional.
     *
     * @param X factor X
     * @param Y factor Y
     * @param includeJoin set to true if operations that involved share variables should be counted, false otherwise
     * @return complexity of the product X * Y
     */
    public static int getComplexity(AbstractFactor X, AbstractFactor Y, boolean includeJoin) {
        return getComplexity(X.getEnumVars(), Y.getEnumVars(), includeJoin);
    }

    /**
     * Create a binary tree with all the pairwise products that are required to complete the full product.
     * The order of the products are determined so to minimise the computational cost over the full product.
     * @param factors all the factors that are to be multiplied
     * @return a binary tree that defines a good order
     */
    public static FactorProductTree getProductTree(AbstractFactor[] factors) {
        int N = factors.length;
        // deal with special cases
        if (N == 0)
            return null;
        if (N == 1)
            return new FactorProductTree(factors[0]);
        if (N == 2) 
            return new FactorProductTree(new FactorProductTree(factors[0]), new FactorProductTree(factors[1]));
        // more than two so optimisation is performed...
        Map<FactorProduct, Integer> cmplx = new HashMap<>();
        // construct the initial working set, of:
        // nodes that MUST be used in the final tree
        List<FactorProductTree> fpool = new ArrayList<>();
        // first leaves for all factors that must take part of the tree
        for (int i = 0; i < N; i ++) 
            fpool.add(new FactorProductTree(factors[i]));
        // now, the main part of repeatedly identifying the best node and integrating that into the final tree
        // this will always be N - 1 products/internal nodes (where N is the number of factors)
        FactorProductTree node = null;
        for (int rank = 0; rank < N - 1; rank ++) {
            int lowest = Integer.MAX_VALUE;
            int a = -1, b = -1;
            for (int i = 0; i < fpool.size(); i ++) {
                EnumVariable[] evars_i = fpool.get(i).getEnumVars();
                for (int j = i + 1; j < fpool.size(); j ++) {
                    EnumVariable[] evars_j = fpool.get(j).getEnumVars();
                    Integer cost = getComplexity(evars_i, evars_j, true); // this is quick so no real need to cache these numbers
                    if (cost < lowest) {
                        a = i; 
                        b = j;
                        lowest = cost;
                    }
                }
            }
            FactorProductTree child_a = fpool.get(a);
            FactorProductTree child_b = fpool.get(b);
            node = new FactorProductTree(child_a, child_b);
            fpool.remove(b); // have to remove in opposite order, to not change the indices where "a" is
            fpool.remove(a);
            fpool.add(node);
        }
        return node;
    }
    
    static private class FactorProduct {
        final EnumVariable[] Xvars, Yvars;
        FactorProduct(EnumVariable[] Xvars, EnumVariable[] Yvars) {
            this.Xvars = Xvars; 
            this.Yvars = Yvars;
        }
        @Override
	public boolean equals(Object obj) {
            if (this == obj) return true;
            if (obj == null || (obj.getClass() != this.getClass())) return false;
            FactorProduct fp = (FactorProduct)obj;
            return ((Xvars == fp.Xvars && Yvars == fp.Yvars) || ((Xvars == fp.Yvars) && (Yvars == fp.Xvars)));
        }
        @Override
        public int hashCode() {
            int hash = 7;
            hash = 31 * hash + (Xvars == null ? 0 : Arrays.hashCode(Xvars)) + (Yvars == null ? 0 : Arrays.hashCode(Yvars));
            return hash;
        }
    }
    
    static public class FactorProductTree {
        private final FactorProductTree x, y;
        private AbstractFactor f = null;
        private EnumVariable[] evars = null;
        /**
         * Construct a node in the tree
         * @param x child
         * @param y child
         */
        FactorProductTree(FactorProductTree x, FactorProductTree y) {
            this.x = x;
            this.y = y;
            Variable[] vars = sortRemoveDuplicates(concat(x.getEnumVars(), y.getEnumVars()));
            this.evars = new EnumVariable[vars.length];
            for (int i = 0; i < vars.length; i ++)
                this.evars[i] = (EnumVariable) vars[i];
        }
        
        /**
         * Construct a leaf node in the tree
         * @param f the factor
         */
        FactorProductTree(AbstractFactor f) {
            this.x = null;
            this.y = null;
            this.f = f;
        }
        EnumVariable[] getEnumVars() {
            if (f == null)
                return this.evars;
            else
                return f.getEnumVars();
        }
        /**
         * Assign the factor to a node in the tree
         * @param f the factor
         */
        void setFactor(AbstractFactor f) {
            this.f = f;
        }
        AbstractFactor getFactor() {
            return f;
        }
    }
    
    /**
     * Concatenate two arrays of variables into one.
     * Does not consider order or duplicates.
     * @param one
     * @param two
     * @return
     */
    protected static Variable[] concat(Variable[] one, Variable[] two) {
        if (one == null && two == null) {
            return new Variable[0];
        } else if (one == null) {
            Variable[] arr = new Variable[two.length];
            System.arraycopy(two, 0, arr, 0, two.length);
            return arr;
        } else if (two == null) {
            Variable[] arr = new Variable[one.length];
            System.arraycopy(one, 0, arr, 0, one.length);
            return arr;
        } else {
            Variable[] arr = new Variable[one.length + two.length];
            System.arraycopy(one, 0, arr, 0, one.length);
            System.arraycopy(two, 0, arr, one.length, two.length);
            return arr;
        }
    }

    /**
     * Create indices that allow quick cross-referencing between two lists of variables.
     * @param xvars variables in a table X
     * @param xcross indices to go from X to Y, i.e. at [x] you find y (or -1 if no mapping)
     * @param yvars variables in a table Y
     * @param ycross indices to go from Y to X, i.e. at [y] you find x (or -1 if no mapping)
     * @return number of positions that overlap
     */
    protected static int crossReference(Variable[] xvars, int[] xcross, Variable[] yvars, int[] ycross) {
        if (xcross != null && ycross != null) {
            if (xvars.length != xcross.length || yvars.length != ycross.length) {
                throw new RuntimeException("Lists must be equal");
            }
            int noverlap = 0; // number of overlapping variables
            for (int j = 0; j < yvars.length; j++) {
                ycross[j] = -1; // Assume Yj does not exist in X
            }
            for (int i = 0; i < xvars.length; i++) {
                xcross[i] = -1; // Assume Xi does not exist in Y
                for (int j = 0; j < yvars.length; j++) {
                    if (xvars[i].equals(yvars[j])) {
                        // this variable Xi is at Yj
                        xcross[i] = j;
                        ycross[j] = i;
                        noverlap++;
                        break;
                    }
                }
            }
            return noverlap;
        } else if (xcross != null) {
            // ycross == null
            if (xvars.length != xcross.length) {
                throw new RuntimeException("X Lists must be equal");
            }
            int noverlap = 0; // number of overlapping variables
            for (int i = 0; i < xvars.length; i++) {
                xcross[i] = -1; // Assume Xi does not exist in Y
                for (int j = 0; j < yvars.length; j++) {
                    if (xvars[i].equals(yvars[j])) {
                        // this variable Xi is at Yj
                        xcross[i] = j;
                        noverlap++;
                        break;
                    }
                }
            }
            return noverlap;
        } else if (ycross != null) {
            // xcross == null
            if (yvars.length != ycross.length) {
                throw new RuntimeException("Y Lists must be equal");
            }
            int noverlap = 0; // number of overlapping variables
            for (int i = 0; i < yvars.length; i++) {
                ycross[i] = -1; // Assume Yi does not exist in X
                for (int j = 0; j < xvars.length; j++) {
                    if (yvars[i].equals(xvars[j])) {
                        // this variable Yi is at Xj
                        ycross[i] = j;
                        noverlap++;
                        break;
                    }
                }
            }
            return noverlap;
        }
        return 0;
    }

    /**
     * Construct a new table that is the result of a factor product of the two specified tables.
     * Works in the general case but there are significant efficiency gains if joint variables are ordered.
     * The implementation currently identifies different cases some of which can be run efficiently,
     * whilst others require tricky indexing. The method is thus quite long but most of the time only a fraction
     * of the code is actually executed.
     * 
     * TODO: Make informed choices as to what implementation of AbstractFactor should be used. 
     * Currently only DenseFactor is used.
     *
     * @param X one table
     * @param Y other table
     * @return the product of one and the other table
     */
    public static AbstractFactor getProduct(AbstractFactor X, AbstractFactor Y) {
        // First resolve cases with tables without enumerable variables
        if (X.nEVars == 0 && Y.nEVars == 0) {
            // both tables lack enumerable variables
            AbstractFactor dt = new DenseFactor(concat(X.nvars, Y.nvars));
            dt.setValue(X.getValue() * Y.getValue());
            if (X.isJDF() && Y.isJDF()) {
                dt.setJDF(JDF.combine(X.getJDF(), Y.getJDF()));
            } else if (X.isJDF()) {
                dt.setJDF(X.getJDF());
            } else if (Y.isJDF()) {
                dt.setJDF(Y.getJDF());
            }
            if (X.isTraced()) {
                dt.setAssign(X.getAssign());
            }
            if (Y.isTraced()) {
                dt.setAssign(Y.getAssign());
            }
            return dt;
        } else if (X.nEVars == 0) {
            // only X is without enumerables
            AbstractFactor dt = new DenseFactor(concat(Y.evars, concat(X.nvars, Y.nvars)));
            for (int j = 0; j < Y.getSize(); j++) {
                dt.setValue(j, Y.getValue(j) * X.getValue());
                if (X.isJDF() && Y.isJDF()) {
                    dt.setJDF(j, JDF.combine(X.getJDF(), Y.getJDF(j)));
                } else if (X.isJDF()) {
                    dt.setJDF(j, X.getJDF());
                } else if (Y.isJDF()) {
                    dt.setJDF(j, Y.getJDF(j));
                }
                if (X.isTraced()) {
                    dt.setAssign(j, X.getAssign());
                }
                if (Y.isTraced()) {
                    dt.setAssign(j, Y.getAssign(j));
                }
            }
            return dt;
        } else if (Y.nEVars == 0) {
            // only Y is without enumerables
            AbstractFactor dt = new DenseFactor(concat(X.evars, concat(X.nvars, Y.nvars)));
            for (int j = 0; j < X.getSize(); j++) {
                dt.setValue(j, X.getValue(j) * Y.getValue());
                if (X.isJDF() && Y.isJDF()) {
                    dt.setJDF(j, JDF.combine(X.getJDF(j), Y.getJDF()));
                } else if (X.isJDF()) {
                    dt.setJDF(j, X.getJDF(j));
                } else if (Y.isJDF()) {
                    dt.setJDF(j, Y.getJDF());
                }
                if (X.isTraced()) {
                    dt.setAssign(j, X.getAssign(j));
                }
                if (Y.isTraced()) {
                    dt.setAssign(j, Y.getAssign());
                }
            }
            return dt;
        }
        // Both tables have enumerable variables, so let's work out how they relate
        int[] xcross2y = new int[X.nEVars]; // map from X to Y indices [x] = y
        int[] ycross2x = new int[Y.nEVars]; // map from Y to X indices [y] = x
        int noverlap = 0; // number of overlapping variables
        noverlap = crossReference(X.evars, xcross2y, Y.evars, ycross2x);
        if (VERBOSE) {
            System.err.println("#overlap = " + noverlap);
        }
        // check if ordered
        boolean ordered = true; // true if *all* overlapping variables are in the same order
        int prev = -1;
        for (int i = 0; i < xcross2y.length && ordered; i++) {
            if (xcross2y[i] != -1) {
                if (xcross2y[i] <= prev) {
                    ordered = false;
                    break;
                }
                prev = xcross2y[i];
            }
        }
        if (VERBOSE) {
            System.err.println("ordered = " + ordered);
        }
        // Before handling complex scenarios, consider some special, simpler cases
        // 1. two tables with the same variables in the same order
        if (noverlap == Math.min(X.nEVars, Y.nEVars)) {
            // at least one table is "contained" by the other
            if (X.nEVars == Y.nEVars) {
                Object[] ykey = new Object[Y.nEVars];
                AbstractFactor dt = new DenseFactor(concat(X.evars, concat(X.nvars, Y.nvars)));
                for (int x = 0; x < X.getSize(); x++) {
                    int y = x; // if ordered the inices in the tables are identical
                    if (!ordered) {
                        // re-index since the variables are listed in a different order
                        Object[] xkey = X.getKey(x);
                        for (int i = 0; i < xkey.length; i++) {
                            ykey[xcross2y[i]] = xkey[i];
                        }
                        y = Y.getIndex(ykey);
                    }
                    dt.setValue(x, X.getValue(x) * Y.getValue(y));
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(x, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(x, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(x, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.setAssign(x, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.setAssign(x, Y.getAssign(y));
                    }
                }
                if (VERBOSE) {
                    System.err.println("DT: Complete overlap, ordered = " + ordered);
                }
                return dt;
            } else if (X.nEVars > Y.nEVars) {
                // Y is more compact than X
                // some variables in X are not in Y
                Set<EnumVariable> notInY = new HashSet<>();
                for (int i = 0; i < xcross2y.length; i++) {
                    if (xcross2y[i] == -1) {
                        notInY.add(X.evars[i]);
                    }
                }
                Object[] ykey = new Object[Y.nEVars];
                AbstractFactor dt = new DenseFactor(concat(X.evars, concat(X.nvars, Y.nvars)));
                for (int x = 0; x < X.getSize(); x++) {
                    double xval = X.getValue(x);
                    if (xval == 0) {
                        continue; // no point in continuing since product will always be zero for entries with this x-value
                    }
                    int y; // there can only be one index in Y (if it is contained in X)
                    if (!ordered) {
                        // re-index considering that variables are listed in a different order
                        Object[] xkey = X.getKey(x);
                        for (int i = 0; i < xkey.length; i++) {
                            if (xcross2y[i] != -1) {
                                ykey[xcross2y[i]] = xkey[i];
                            }
                        }
                        y = Y.getIndex(ykey);
                    } else {
                        // re-index but ordered, not sure if this is quicker than above... TODO: test
                        y = X.maskIndex(x, notInY);
                    }
                    double yval = Y.getValue(y);
                    if (yval == 0) {
                        continue; // the product will be zero no what
                    }
                    int idx = x;
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        if (!dt.isTraced())
                            dt.setTraced(true);
                        dt.setAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        if (!dt.isTraced())
                            dt.setTraced(true);
                        dt.setAssign(idx, Y.getAssign(y));
                    }
                }
                if (VERBOSE) {
                    System.err.println("DT: Partial overlap (X>Y), ordered = " + ordered);
                }
                return dt;
            } else if (Y.nEVars > X.nEVars) {
                // X is more compact than Y
                // some variables in Y are not in X
                Set<EnumVariable> notInX = new HashSet<>();
                for (int i = 0; i < ycross2x.length; i++) {
                    if (ycross2x[i] == -1) {
                        notInX.add(Y.evars[i]);
                    }
                }
                Object[] xkey = new Object[X.nEVars];
                AbstractFactor dt = new DenseFactor(concat(Y.evars, concat(X.nvars, Y.nvars)));
                for (int y = 0; y < Y.getSize(); y++) {
                    double yval = Y.getValue(y);
                    if (yval == 0) {
                        continue; // no point in continuing since product will always be zero for entries with this x-value
                    }
                    int x; // there can only be one index in X (if it is contained in Y)
                    if (!ordered) {
                        // re-index considering that variables are listed in a different order
                        Object[] ykey = Y.getKey(y);
                        for (int i = 0; i < ykey.length; i++) {
                            if (ycross2x[i] != -1) {
                                xkey[ycross2x[i]] = ykey[i];
                            }
                        }
                        x = X.getIndex(xkey);
                    } else {
                        // re-index but ordered, not sure if this is quicker than above... TODO: test
                        x = Y.maskIndex(y, notInX);
                    }
                    double xval = X.getValue(x);
                    if (xval == 0) {
                        continue; // the product will be zero no what
                    }
                    int idx = y;
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.setAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.setAssign(idx, Y.getAssign(y));
                    }
                }
                if (VERBOSE) {
                    System.err.println("DT: Partial overlap (X<Y), ordered = " + ordered);
                }
                return dt;
            }
        }
        // Failing the above, we must construct a table which is an amalgamate of the two.
        AbstractFactor dt = new DenseFactor(concat(concat(X.evars, Y.evars), concat(X.nvars, Y.nvars)));
        Object[] reskey = new Object[dt.nEVars]; // this will initialise all elements to null
        int[] xcross2dt = new int[X.nEVars]; // map from X index to dt index
        int[] ycross2dt = new int[Y.nEVars]; // map from Y index to dt index
        crossReference(X.evars, xcross2dt, dt.evars, null);
        crossReference(Y.evars, ycross2dt, dt.evars, null);
        // 2. two tables have nothing in common
        if (noverlap == 0) {
            for (int x = 0; x < X.getSize(); x++) {
                double xval = X.getValue(x);
                if (xval == 0) {
                    continue; // no point in continuing since product will always be zero for entries with this x-value
                }
                Object[] xkey = X.getKey(x);
                for (int y = 0; y < Y.getSize(); y++) {
                    double yval = Y.getValue(y);
                    if (yval == 0) {
                        continue; // the product will be zero no what
                    }
                    Object[] ykey = Y.getKey(y);
                    for (int i = 0; i < xkey.length; i ++)
                        reskey[xcross2dt[i]] = xkey[i];
                    for (int i = 0; i < ykey.length; i ++)
                        reskey[ycross2dt[i]] = ykey[i];
                    int idx = dt.getIndex(reskey);
                    // **************************
                    //     This does NOT work
                    // **************************
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.setAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.setAssign(idx, Y.getAssign(y));
                    }
                }
            }
            if (VERBOSE) {
                System.err.println("DT: No overlap.");
            }
            return dt;
        }
        // 3. General case, if none of those implemented above has worked
        Object[] searchkey = new Object[Y.nEVars]; // this will initialise all elements to null
        int option = PRODUCT_OPTION;
        long start0 = 0;
        long start1 = 0;
        long total0 = 0;
        long total1 = 0;
        for (int x = 0; x < X.getSize(); x++) {
            double xval = X.getValue(x);
            if (xval == 0) {
                continue; // no point in continuing since product will always be zero for entries with this x-value
            }
            Object[] xkey = X.getKey(x);
            for (int i = 0; i < xcross2y.length; i++) {
                int idx = xcross2y[i];
                if (idx > -1) {
                    searchkey[idx] = xkey[i];
                }
                reskey[xcross2dt[i]] = xkey[i];
            }
            if (start0 == 0 && option == -1) {
                start0 = System.nanoTime();
                option = 0;
            } else if (start1 == 0 && option == -1) {
                start1 = System.nanoTime();
                option = 1;
            }
            // before computing all *matching* indices in Y, weigh the cost of traversing all indices instead, to subsequently check for match
            if (option == 0) {
                int[] yindices = Y.getIndices(searchkey);
                for (int y : yindices) {
                    double yval = Y.getValue(y);
                    if (yval == 0) {
                        continue; // the product will be zero no what
                    }
                    Object[] ykey = Y.getKey(y);
                    for (int i = 0; i < ykey.length; i++) {
                        if (ycross2x[i] == -1) {
                            reskey[ycross2dt[i]] = ykey[i];
                        }
                    }
                    int idx = dt.getIndex(reskey);
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.setAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.setAssign(idx, Y.getAssign(y));
                    }
                }
                if (start0 != 0) {
                    total0 = System.nanoTime() - start0;
                    option = -1;
                }
            } else if (option == 1) {
                for (int y = 0; y < Y.getSize(); y++) {
                    if (Y.isMatch(searchkey, y)) {
                        double yval = Y.getValue(y);
                        if (yval == 0) {
                            continue; // the product will be zero no what
                        }
                        Object[] ykey = Y.getKey(y);
                        for (int i = 0; i < ykey.length; i++) {
                            if (ycross2x[i] == -1) {
                                reskey[ycross2dt[i]] = ykey[i];
                            }
                        }
                        int idx = dt.getIndex(reskey);
                        dt.setValue(idx, xval * yval);
                        if (X.isJDF() && Y.isJDF()) {
                            dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                        } else if (X.isJDF()) {
                            dt.setJDF(idx, X.getJDF(x));
                        } else if (Y.isJDF()) {
                            dt.setJDF(idx, Y.getJDF(y));
                        }
                        if (X.isTraced()) {
                            dt.setAssign(idx, X.getAssign(x));
                        }
                        if (Y.isTraced()) {
                            dt.setAssign(idx, Y.getAssign(y));
                        }
                    }
                }
                if (start1 != 0) {
                    total1 = System.nanoTime() - start1;
                    option = -1;
                }
            }
            if (start1 != 0) {
                // done with timing
                if (total0 > total1) {
                    option = 1;
                } else {
                    option = 0;
                }
            }
        }
        if (VERBOSE) {
            System.err.println("DT: Generic case. Option = " + option + " (Option 0 took " + total0 + "ns. Option 1 took " + total1 + "ns.)");
        }
        return dt;
    }

    /**
     * Construct a new table from an existing, by summing-out specified variable/s.
     * 
     * TODO: Make informed choices as to what implementation of AbstractFactor should be used. 
     * Currently only DenseFactor is used.
     *
     * @param X existing table
     * @param anyvars variables to sum-out
     * @return the resulting "margin" of the table
     */
    public static AbstractFactor getMargin(AbstractFactor X, Variable... anyvars) {
        // sort and weed out duplicates
        Variable[] uniqueVars = sortRemoveDuplicates(anyvars);
        // first check that X has enumerable variables, because that would be required
        if (X.nEVars == 0) {
            return X;
        }
        // if X is ok, get rid of non-enumerables
        int cnt_evars = 0;
        for (Variable var : uniqueVars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                cnt_evars++;
            } catch (ClassCastException e) {
            }
        }
        EnumVariable[] evars = new EnumVariable[cnt_evars];
        if (cnt_evars == 0) {
            return X;
        } else {
            cnt_evars = 0;
            for (Variable var : uniqueVars) {
                try {
                    EnumVariable evar = (EnumVariable) var;
                    evars[cnt_evars++] = evar;
                } catch (ClassCastException e) {
                }
            }
        }
        // now, resolve the relationships between X's variables and those to be summed-out
        Variable[] yvars = new Variable[X.nEVars + X.nNVars - evars.length];
        int[] xkeyidx = new int[X.nEVars - evars.length]; // map from evars to X indices [evar] = x
        int cnt = 0;
        for (int i = 0; i < X.nEVars; i++) {
            boolean keep = true;
            for (EnumVariable evar : evars) {
                if (X.evars[i].equals(evar)) {
                    keep = false;
                    break;
                }
            }
            if (keep) {
                xkeyidx[cnt] = i;
                yvars[cnt++] = X.evars[i];
            }
        }
        for (int i = 0; i < X.nNVars; i++) {
            yvars[cnt++] = X.nvars[i];
        }
        if (cnt != X.nEVars + X.nNVars - evars.length) {
            throw new AbstractFactorRuntimeException("Invalid variable list");
        }
        AbstractFactor Y = new DenseFactor(yvars);
        Object[] xkey_search = new Object[X.nEVars];
        for (int y = 0; y < Y.getSize(); y++) { // FIXME: not effective for sparse factors
            Object[] ykey = Y.getKey(y);
            for (int i = 0; i < xkeyidx.length; i++) {
                xkey_search[xkeyidx[i]] = ykey[i];
            }
            double sum = 0;
            int[] indices = X.getIndices(xkey_search);
            JDF first = null, subsequent = null;
            JDF mixture = null;
            double prev_weight = 0;
            for (int x : indices) {
                double xval = X.getValue(x);
                sum += xval;
            }
            Y.setValue(y, sum);
            if (X.isJDF()) {
                for (int x : indices) {
                    double nxval = X.getValue(x) / sum; // normalized
                    if (nxval != 0) {
                        if (first == null) {
                            // this will be true only once, for the first entry that will be collapsed
                            first = X.getJDF(x); // save the JDF for later
                            prev_weight = nxval; // save the value associated with this entry to weight the distributions
                        } else if (mixture == null) {
                            // for the subsequent entry, the mixture is created, both JDF weighted
                            subsequent = X.getJDF(x);
                            if (subsequent != null)
                                mixture = JDF.mix(first, prev_weight, subsequent, nxval);
                        } else {
                            // for all subsequent entries, the mixture is mixed unchanged with the new entry's JDF weighted
                            subsequent = X.getJDF(x);
                            if (subsequent != null)
                                mixture = JDF.mix(mixture, subsequent, nxval);
                        }
                    }
                }
                if (first != null) {
                    if (mixture != null)
                        Y.setJDF(y, mixture);
                    else
                        Y.setJDF(y, first);
                }
            }
        }
        return Y;
    }

    /**
     * Construct a new table from an existing, by maxing-out specified variable/s, and tracing
     * the assignment that provided the maximum value in the resulting table.
     * This code is based on that of getMargin.
     *
     * @param X existing table
     * @param anyvars variables to sum-out
     * @return the resulting "margin" of the table
     */
    public static AbstractFactor getMaxMargin(AbstractFactor X, Variable... anyvars) {
        // sort and weed out duplicates
        Variable[] uniqueVars = sortRemoveDuplicates(anyvars);
        // first check that X has enumerable variables, because that would be required
        if (X.nEVars == 0) {
            return X;
        }
        // if X is ok, get rid of non-enumerables
        int cnt_evars = 0;
        for (Variable var : uniqueVars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                cnt_evars++;
            } catch (ClassCastException e) {
            }
        }
        EnumVariable[] evars = new EnumVariable[cnt_evars];
        if (cnt_evars == 0) {
            return X;
        } else {
            cnt_evars = 0;
            for (Variable var : uniqueVars) {
                try {
                    EnumVariable evar = (EnumVariable) var;
                    evars[cnt_evars++] = evar;
                } catch (ClassCastException e) {
                }
            }
        }
        // now find the relationships between the enumerables
        Variable[] yvars = new Variable[X.nEVars + X.nNVars - evars.length];
        int[] xkeyidx = new int[X.nEVars - evars.length]; // map from Y to X indices [y] = x
        int[] outidx = new int[evars.length]; // index of enumerable variables in X that will be summed out
        int cnt_keep = 0;
        int cnt_dontkeep = 0;
        for (int i = 0; i < X.nEVars; i++) {
            boolean keep = true;
            for (EnumVariable evar : evars) {
                if (X.evars[i].equals(evar)) {
                    keep = false;
                    break;
                }
            }
            if (keep) {
                xkeyidx[cnt_keep] = i;
                yvars[cnt_keep++] = X.evars[i];
            } else {
                outidx[cnt_dontkeep++] = i;
            }
        }
        for (int i = 0; i < X.nNVars; i++) {
            yvars[cnt_keep++] = X.nvars[i];
        }
        if (cnt_keep != X.nEVars + X.nNVars - evars.length) {
            throw new AbstractFactorRuntimeException("Invalid variable list");
        }
        AbstractFactor Y = new DenseFactor(yvars);
        Y.setTraced(true);
        Object[] xkey_search = new Object[X.nEVars];
        for (int y = 0; y < Y.getSize(); y++) {
            Object[] ykey = Y.getKey(y);
            for (int i = 0; i < xkeyidx.length; i++) {
                xkey_search[xkeyidx[i]] = ykey[i];
            }
            double max = Double.NEGATIVE_INFINITY;
            int maxidx = 0;
            int[] indices = X.getIndices(xkey_search);
            for (int x : indices) {
                double xval = X.getValue(x);
                if (xval > max) {
                    max = xval;
                    maxidx = x;
                }
            }
            Y.setValue(y, max);
            if (Y.isJDF()) {
                Y.setJDF(y, X.getJDF(maxidx));
            }
            // Trace implied assignments
            Object[] xkey = X.getKey(maxidx);
            for (int i = 0; i < evars.length; i++) {
                Variable.Assignment a = new Variable.Assignment(evars[i], xkey[outidx[i]]);
                Y.addAssign(y, a);
            }
        }
        return Y;
    }

    public static AbstractFactor getNormal(AbstractFactor X) {
        DenseFactor Y = new DenseFactor(concat(X.evars, X.nvars));
        if (X.hasEnumVars()) {
            double sum = X.getSum();
            for (int i = 0; i < X.getSize(); i ++) {
                Y.setValue(i, X.getValue(i) / sum);
                if (X.isJDF())
                    Y.setJDF(i, X.getJDF(i));
            }
        } else {
            Y.setValue(X.getValue());
            if (X.isJDF())
                Y.setJDF(X.getJDF());
        }
        return Y;
    }

    /**
     * Get a key from a partial assignment of variables defined for the table.
     * @param vararr variables in order
     * @param evid the evidence
     * @return the key that encodes the values for the provided variables
     */
    public static Object[] getKey(EnumVariable[] vararr, Variable.Assignment[] evid) {
        Object[] key = new Object[vararr.length]; // key is initialised to nulls by default
        if (key.length <= evid.length) {
            // for efficiency, we check what to iterate over
            for (Variable.Assignment e : evid) {
                try {
                    EnumVariable evar = (EnumVariable) e.var;
                    int var_index = Arrays.binarySearch(vararr, evar);
                    if (var_index >= 0) {
                        key[var_index] = e.val;
                    }
                } catch (ClassCastException exception) {
                    ; // ignore non-enumerable variables
                }
            }
        } else {
            // evidence is longer than key, so let's iterate over key
            for (int i = 0; i < key.length; i++) {
                Variable var = vararr[i];
                for (Variable.Assignment evid1 : evid) {
                    if (evid1.var.equals(var)) {
                        key[i] = evid1.val;
                        break;
                    }
                }
            }
        }
        return key;
    }

    /**
     * Calculate the sum of factors
     * @return the sum
     */
    public double getSum() {
        double sum =0;
        if (!hasEnumVars())
            return getValue();
        for (int i = 0; i < getSize(); i ++) {
            sum += getValue(i);
        }
        return sum;
    }
    
    /**
     * Copy over all non-null values from source to target key.
     *
     * @param target
     * @param source
     * @return target
     **/
    public static Object[] overlay(Object[] target, Object[] source) {
        if (target.length != source.length) {
            throw new AbstractFactorRuntimeException("Invalid operation since keys are of difference lengths (" + target.length + " vs " + source.length + ")");
        }
        for (int i = 0; i < target.length; i++) {
            if (source[i] != null) {
                target[i] = source[i];
            }
        }
        return target;
    }

    /**
     * Find out if this table has JDFs or not.
     * Equivalent to asking if the table has non-enumerable variables.
     * @return true if non-enumerable variables are included in this table, false otherwise
     */
    public boolean isJDF() {
        return this.nNVars > 0;
    }

    /**
     * Check if factor is defined by the specified variable.
     * @param var variable
     * @return true if in factor (enumerable or non-enumerable), false otherwise
     */
    public boolean hasVariable(Variable var) {
        for (int i = 0; i < nEVars; i ++)
            if (evars[i].equals(var))
                return true;
        for (int i = 0; i < nNVars; i ++)
            if (nvars[i].equals(var))
                return true;
        return false;
    }
    
    /**
     * Get the theoretical number of entries in this table. Note this number is
     * always equal to the actual number of entries, but some may never have been explicitly
     * set to a value.
     *
     * @return the size (number of entries)
     */
    public int getSize() {
        if (period != null)
            if (period.length > 0)
                return period[0];
        return 1;
    }

    /**
     * Find out how full the table is.
     * @return the percentage of table that is occupied with non-zero values
     */
    public double getCapacity() {
        return getOccupied() / getSize();
    }
    
    /**
     * Get the enumerable variables of the table.
     * @return the variables in original order.
     */
    public EnumVariable[] getEnumVars() {
        if (nEVars > 0)
            return evars;
        return new EnumVariable[0];
    }
    
    /**
     * Find out if there are enumerable variables defining the table.
     * @return true if enumerable variables define the table, false otherwise
     */
    public boolean hasEnumVars() {
        return nEVars > 0;
    }

    /**
     * Get the non-enumerable variables of the table.
     * @return the variables in original order.
     */
    public Variable[] getNonEnumVars() {
        if (nNVars > 0)
            return nvars;
        return new Variable[0];
    }

    /**
     * Find out if there are non-enumerable variables defining the table.
     * @return true if non-enumerable variables define the table, false if not
     */
    public boolean hasNonEnumVars() {
        return nNVars > 0;
    }

    /**
     * Get the canonical names of the parent variables (names + "." + index)
     * @return names of variables (in order)
     */
    public String[] getLabels() {
        String[] labels = new String[nEVars + nNVars];
        Variable[] all = concat(this.evars, this.nvars);
        for (int i = 0; i < labels.length; i++) {
            labels[i] = all[i].toString();
        }
        return labels;
    }

    /**
     * Retrieve the index for the specified key
     *
     * @param key the values by which the index is identified (order the same as
     * when constructing the factor table)
     * @return the index for the instantiated key
     */
    public int getIndex(Object[] key) {
        if (getSize() == 1) {
            throw new AbstractFactorRuntimeException("Invalid key: no variables");
        }
        if (key.length != nEVars) {
            throw new AbstractFactorRuntimeException("Invalid key: length is " + key.length + " not " + nEVars);
        }
        int sum = 0;
        for (int i = 0; i < nEVars; i++) {
            if (key[i] == null) {
                throw new AbstractFactorRuntimeException("Null in key");
            }
            sum += (evars[i].getIndex(key[i]) * step[i]);
        }
        return sum;
    }

    /**
     * Instantiate the key for the specified index. The variables in the factor
     * table must be Enumerable.
     *
     * @param index the index of an entry in the factor table
     * @return the values of the key corresponding to the entry with the index
     */
    public Object[] getKey(int index) {
        if (index >= getSize() || index < 0 || this.getSize() == 1) {
            throw new AbstractFactorRuntimeException("Invalid index");
        }
        int remain = index;
        Object[] key = new Object[nEVars];
        for (int i = 0; i < nEVars; i++) {
            int keyindex = remain / step[i];
            key[i] = evars[i].getDomain().get(keyindex);
            remain -= keyindex * step[i];
        }
        return key;
    }

    /**
     * Method to check if a key instance has the specified index.
     * Must be fast since it may be called frequently.
     *
     * @param key
     * @param index
     * @return true if the key instance maps to the index
     */
    public boolean isMatch(Object[] key, int index) {
        if (key.length != nEVars || index < 0 || index >= getSize()) {
            throw new AbstractFactorRuntimeException("Invalid index or key");
        }
        int remain = index;
        for (int i = 0; i < nEVars; i++) {
            if (key[i] != null) {
                int keyindex = evars[i].getIndex(key[i]);
                if (keyindex != remain / step[i]) {
                    return false;
                }
                remain -= keyindex * step[i];
            } else {
                // key[i] == null
                int missing = remain / step[i];
                remain -= missing * step[i];
            }
        }
        return true;
    }

    /**
     * Takes an entry index of the current table and "masks" out a subset of
     * parents, to determine the index in a table with only the parents that are
     * not masked.
     * Note that the order of variables is assumed to be the same.
     * Time complexity is O(3n) where n is the number of parents in
     * the current table. (This computation could be done marginally more
     * efficiently.)
     *
     * @param origindex
     * @param maskMe
     * @return index in other, more compact table
     */
    public int maskIndex(int origindex, Set<EnumVariable> maskMe) {
        int origremain = origindex;
        int sum = 0;
        int jn = 0;
        int[] newstep = new int[nEVars - maskMe.size()];
        int[] newvale = new int[nEVars - maskMe.size()];
        for (int i = 0; i < nEVars; i++) {
            if (!maskMe.contains(evars[i])) {
                newvale[jn++] = domsize[i];
            }
        }
        jn = newstep.length - 1;
        int prod = 1;
        for (int i = nEVars - 1; i >= 0; i--) {
            if (!maskMe.contains(evars[i])) {
                // if NOT masked-out
                newstep[jn] = prod;
                prod *= newvale[jn--];
            }
        }
        jn = 0;
        for (int i = 0; i < nEVars; i++) {
            int key = origremain / step[i];
            origremain -= key * step[i];
            if (!maskMe.contains(evars[i])) {
                sum += (key * newstep[jn++]);
            }
        }
        return sum;
    }

    /**
     * Get a key from a partial assignment of variables defined for the table.
     * @param evid
     * @return
     */
    public Object[] getKey(Variable.Assignment[] evid) {
        return AbstractFactor.getKey(this.evars, evid);
    }

    /**
     * Print the table.
     */
    public void display() {
        System.out.print("Idx ");
        for (int j = 0; j < this.nEVars; j++) {
            System.out.print(String.format("[%8s]", this.evars[j].getName()));
        }
        System.out.print(" F     ");
        for (int i = 0; i < nNVars; i++) {
            System.out.print(String.format("[%8s]", nvars[i].getName()));
        }
        System.out.println();
        for (int i = 0; i < this.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            if (this.getSize() == 1) {
                System.out.print(String.format(" %5.3f ", this.getValue()));
                for (int j = 0; j < nNVars; j++) {
                    Distrib d = this.getDistrib(nvars[j]);
                    System.out.print(String.format(" %s ", d == null ? "-" : d.toString()));
                }
                if (this.isTraced()) {
                    for (Variable.Assignment a : this.getAssign()) {
                        System.out.print(a + ";");
                    }
                }
            } else {
                Object[] key = this.getKey(i);
                for (Object key1 : key) {
                    System.out.print(String.format(" %-8s ", key1.toString()));
                }
                System.out.print(String.format(" %5.3f ", this.getValue(i)));
                for (int j = 0; j < nNVars; j++) {
                    Distrib d = this.getDistrib(i, nvars[j]);
                    System.out.print(String.format(" %s ", d == null ? "-" : d.toString()));
                }
                if (this.isTraced()) {
                    for (Variable.Assignment a : this.getAssign(i)) {
                        System.out.print(a + ";");
                    }
                }
            }
            System.out.println();
        }
    }
    
    /**
     * Retrieve the value of the factor, if table is without enumerable variables.
     * @return the only value of the factor
     */
    public abstract double getValue();
    
    /**
     * Retrieve the JDF of the factor without enumerable variables.
     * @return probability distribution for variable
     */
    public abstract JDF getJDF();
    
    /**
     * Retrieve the value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    public abstract double getValue(int index);

    /**
     * Retrieve the JDF
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @return probability distribution for variable
     */
    public abstract JDF getJDF(int index);
    
    /**
     * Retrieve the value of the entry identified by the instantiated key
     *
     * @param key the entry
     * @return the value of the entry
     */
    public double getValue(Object[] key) {
        int index = getIndex(key);
        return getValue(index);
    }


    /**
     * Retrieve the value of the entry identified by the instantiated key
     *
     * @param key the entry
     * @return the value of the entry
     */
    public double getJDF(Object[] key) {
        int index = getIndex(key);
        return getValue(index);
    }


    /**
     * Set the only value associated with a table without enumerable variables.
     * @param value
     * @return 
     */
    public abstract int setValue(double value);
    
    /**
     * Associate the specified key-index with the given value. Note that using
     * getValue and setValue with index is quicker than with key, if more than
     * one operation is done.
     *
     * @param key_index
     * @param value
     * @return the index at which the value was stored
     */
    public abstract int setValue(int key_index, double value);


    /**
     * Associate the specified key with the given value.
     * @param key
     * @param value
     * @return the index at which the value was stored
     */
    public int setValue(Object[] key, double value) {
        if (getSize() == 1)
            throw new AbstractFactorRuntimeException("Invalid key: no variables");
        int index = getIndex(key);
        return setValue(index, value);
    }


    /**
     * Set the JDF for a table without enumerable variables.
     * @param value JDF
     * @return the index (always 0)
     */
    public abstract int setJDF(JDF value);
    
    /**
     * Associate the specified key-index with the given JDF. 
     *
     * @param key_index index for entry
     * @param value
     * @return the index at which the value was stored
     */
    public abstract int setJDF(int key_index, JDF value);


    /**
     * Associate the specified key with the given JDF.
     * @param key
     * @param value JDF
     * @return the index at which the JDF was stored
     */
    public int setJDF(Object[] key, JDF value) {
        if (getSize() == 1)
            throw new AbstractFactorRuntimeException("Invalid key: no variables");
        int index = getIndex(key);
        return setJDF(index, value);
    }

    /**
     * Activate tracing of implied assignments.
     * @param status true if activated, false otherwise
     */
    public abstract void setTraced(boolean status);

    /**
     * Find out if tracing of implied assignments is active.
     * @return true if activated, false otherwise
     */
    public abstract boolean isTraced();

        /**
     * Associate the only entry with a collection of assignments.
     * @param assign assignment
     * @return index of entry
     */
    public abstract int setAssign(Collection<Variable.Assignment> assign);
    /**
     * Associate the entry with a collection of assignments.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    public abstract int setAssign(int key_index, Collection<Variable.Assignment> assign);
    
    /**
     * Associate the entry with a collection of assignments.
     * @param key key for entry
     * @param assign assignment
     * @return index of entry
     */
    public int setAssign(Object[] key, Collection<Variable.Assignment> assign) {
        return setAssign(getIndex(key), assign);
    }
    
    
    /**
     * Retrieve a collection of assignments from the single entry.
     * @return set of assignments
     */
    public abstract Set<Variable.Assignment> getAssign();

    /**
     * Retrieve a collection of assignments from entry.
     * @param key_index index of entry
     * @return set of assignments
     */
    public abstract Set<Variable.Assignment> getAssign(int key_index);
    
    /**
     * Retrieve a collection of assignments from the specified entry.
     * @param key the entry
     * @return set of assignments
     */
    public Set<Variable.Assignment> getAssign(Object[] key) {
        return getAssign(getIndex(key));
    }
    
    /**
     * Tag the entry with an assignment.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    public abstract int addAssign(int key_index, Variable.Assignment assign);
    
    /**
     * Tag the entry with an assignment.
     * @param key key for entry
     * @param assign assignment
     * @return index of entry
     */
    public int addAssign(Object[] key, Variable.Assignment assign) {
        return addAssign(getIndex(key), assign);
    }
    
    /**
     * Retrieve the distribution of the specified non-enumerable variable.
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public abstract Distrib getDistrib(Variable nvar);
    
    /**
     * Retrieve the distribution of the specified non-enumerable variable, 
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public abstract Distrib getDistrib(int index, Variable nvar);
    
    /**
     * Retrieve the distribution of the specified non-enumerable variable, 
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param key the key of the table, containing the instantiation
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public Distrib getDistrib(Object[] key, Variable nvar) {
        int index = getIndex(key);
        return this.getDistrib(index, nvar);
    }
    
   /**
     * Set the conditional distribution of a non-enumerable variable. 
     * Use with care as the method may not do all integrity checks.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public abstract int setDistrib(int key_index, Variable nonenum, Distrib d);

    /**
     * Set the conditional distribution of a non-enumerable variable. 
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(Object[] key, Variable nonenum, Distrib d) {
        if (this.getSize() == 1)
            throw new AbstractFactorRuntimeException("Invalid key: no variables");
        int index = getIndex(key);
        return setDistrib(index, nonenum, d);
    }


    /**
     * Set the distribution of a non-enumerable variable for a factor without enumerable variables. 
     * @param nvar the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry, always 0
     */
    public abstract int setDistrib(Variable nvar, Distrib d);
    
    /**
     * Find out how many entries that occupied.
     */
    public abstract int getOccupied();
    
    /**
     * Set all entries to 0.
     */
    public abstract void setEmpty();

    /**
     * Identify each index that is linked to the specified key (which may
     * include "wildcards", indicated by null values).
     * This function finds all indices that match the key.
     *
     * @param key
     * @return an array with all matching indices
     */
    public abstract int[] getIndices(Object[] key);
    
    /*************************** Methods to test AbstractFactor **************************/
    
    protected static Variable[] getVariablePool(long seed, int n) {
        Random random = new Random(seed);
        Variable[] vars = new Variable[n];
        for (int i = 0; i < vars.length; i++) {
            int type = random.nextInt(5);
            Variable var = null;
            switch (type) {
                case 0:
                    var = Predef.Boolean();
                    break;
                case 1:
                    var = Predef.Nominal("a", "b", "c");
                    break;
                case 2:
                    var = Predef.NucleicAcid();
                    break;
                case 3:
                    var = Predef.Real();
                    break;
                case 4:
                    var = Predef.Number(random.nextInt(8) + 2);
                    break;
            }
            vars[i] = var;
        }
        return vars;
    }

    protected static Variable[] getSubset(long seed, Variable[] vars, int n) {
        Random random = new Random(seed);
        Set<Variable> unique = new HashSet<>();
        n = Math.min(vars.length, n);
        Variable[] subset = new Variable[n];
        while (unique.size() < n) {
            unique.add(vars[random.nextInt(vars.length)]);
        }
        unique.toArray(subset);
        return subset;
    }

    protected static AbstractFactor[] getFactorPool(long seed, Variable[] vars, int n) {
        Random random = new Random(seed);
        int M = Math.abs((int) (random.nextGaussian() * n) + 1);
        int N = Math.abs((int) (random.nextGaussian() * n) + 1);
        AbstractFactor[] dfs = new DenseFactor[N];
        for (int i = 0; i < dfs.length; i++) {
            int nvars = random.nextInt(Math.max(1, M));
            if (nvars > 0) {
                Variable[] myvars = getSubset(random.nextInt(), vars, nvars);
                dfs[i] = new DenseFactor(myvars);
            } else {
                dfs[i] = new DenseFactor();
            }
            int npop = dfs[i].getSize(); //random.nextInt(dfs[i].getSize()); // number of entries to populate
            for (int j = 0; j < npop; j++) {
                if (dfs[i].getSize() == 1) {
                    dfs[i].setValue(Math.abs(random.nextGaussian()) / npop);    
                    if (dfs[i].isJDF()) {
                        for (Variable nvar : dfs[i].getNonEnumVars()) {
                            dfs[i].setDistrib(nvar, new GaussianDistrib(random.nextGaussian() * random.nextInt(100), Math.abs(random.nextGaussian() * (random.nextInt(10) + 1))));
                        }
                    }
                } else {
                    int index = random.nextInt(npop);
                    dfs[i].setValue(index, Math.abs(random.nextGaussian()) / npop);
                    if (dfs[i].isJDF()) {
                        for (Variable nvar : dfs[i].getNonEnumVars()) {
                            dfs[i].setDistrib(index, nvar, new GaussianDistrib(random.nextGaussian() * random.nextInt(100), Math.abs(random.nextGaussian() * (random.nextInt(10) + 1))));
                        }
                    }
                }
                
            }
        }
        return dfs;
    }

    protected static AbstractFactor getProductBenchmarked(FactorProductTree node) {
        if (node.getFactor() != null)
            return node.getFactor();
        AbstractFactor X = getProductBenchmarked(node.x);
        AbstractFactor Y = getProductBenchmarked(node.y);
        long startTime = System.nanoTime();
        AbstractFactor f = AbstractFactor.getProduct(X, Y);
        long endTime = System.nanoTime();
        for (int j = 0; j < 20; j ++) {
            if (testProductIntegrity(j, X, Y, f) == false)
                System.err.println("Test failed");
        }
        int overlap = getOverlap(X, Y);
        int minevars = Math.min(X.nEVars, Y.nEVars);
        int maxevars = Math.max(Y.nEVars, Y.nEVars);
        System.out.println(maxevars + "\t" + minevars + "\t" + overlap + "\t" + (minevars == 0 ? 0.0 : overlap / (float)minevars) + "\t" + getComplexity(X, Y, false) + "\t" + getComplexity(X, Y, true) + "\t" + (endTime - startTime) / 100000.0);
        node.setFactor(f);
        return f;
    }
    
    protected static AbstractFactor getProductBenchmarked(AbstractFactor[] factors) {
        if (factors.length == 0)
            return null;
        AbstractFactor R = factors[0];
        for (int i = 1; i < factors.length; i ++) {
            AbstractFactor X = R;
            AbstractFactor Y = factors[i];
            long startTime = System.nanoTime();
            R = AbstractFactor.getProduct(X, Y);
            long endTime = System.nanoTime();
            for (int j = 0; j < 20; j ++) {
                if (testProductIntegrity(j, X, Y, R) == false)
                    System.err.println("Test failed");
            }
            int overlap = getOverlap(X, Y);
            int minevars = Math.min(X.nEVars, Y.nEVars);
            int maxevars = Math.max(Y.nEVars, Y.nEVars);
            System.out.println(maxevars + "\t" + minevars + "\t" + overlap + "\t" + (minevars == 0 ? 0.0 : overlap / (float)minevars) + "\t" + getComplexity(X, Y, false) + "\t" + getComplexity(X, Y, true) + "\t" + (endTime - startTime) / 100000.0);
        }
        return R;
    }
    
    public static AbstractFactor getProduct(AbstractFactor[] factors) {
        if (factors.length == 0)
            return null;
        AbstractFactor R = factors[0];
        for (int i = 1; i < factors.length; i ++) {
            AbstractFactor X = R;
            AbstractFactor Y = factors[i];
            R = AbstractFactor.getProduct(X, Y);
        }
        return R;
    }
    
    
    /** Calculate products linearly, don't consider order */
    static final int POOL_OPTION_LINEAR = 0; 
    /** Calculate products according to binary tree, optimised to minimise computational cost */
    static final int POOL_OPTION_TREE = 1; 
    
    protected static AbstractFactor productPool(AbstractFactor[] dfs, int option) {
        if (dfs.length == 0)
            return null;
        if (dfs.length == 1)
            return dfs[0];
        AbstractFactor f = null;
        long startTime = System.nanoTime();
        switch (option) {
            case POOL_OPTION_LINEAR:
                f = AbstractFactor.getProductBenchmarked(dfs);
                break;
            case POOL_OPTION_TREE:
                FactorProductTree tree = AbstractFactor.getProductTree(dfs);
                f = AbstractFactor.getProductBenchmarked(tree);
                break;
        }
        long endTime = System.nanoTime();
        System.out.println("\t\t\t\t\t\t\t" + (endTime - startTime) / 100000.0);
        if (f.nEVars > 0) {
            Random rand = new Random(endTime);
            int n = rand.nextInt(f.nEVars);
            Variable[] sumout = new Variable[n];
            for (int i = 0; i < n; i ++)
                sumout[i] = f.evars[rand.nextInt(n)];
            if (rand.nextBoolean())
                AbstractFactor.getMargin(f, sumout);
            else
                AbstractFactor.getMaxMargin(f, sumout);
            return f;
        }
        return f;
    }

    protected static boolean testProductIntegrity(long seed, AbstractFactor X, AbstractFactor Y, AbstractFactor PROD) {
        Random rand = new Random(seed);
        int p = rand.nextInt(PROD.getSize());
        Object[] pkey;
        if (PROD.getSize() == 1)
            pkey = new Object[0];
        else 
            pkey = PROD.getKey(p);
        EnumVariable[] pevars = PROD.getEnumVars();
        EnumVariable[] xevars = X.getEnumVars();
        EnumVariable[] yevars = Y.getEnumVars();
        Object[] xkey = new Object[X.nEVars];
        Object[] ykey = new Object[Y.nEVars];
        for (int i = 0; i < pkey.length; i ++) {
            if (!pevars[i].getDomain().isValid(pkey[i]))
                return false;
            int xi = -1;
            for (int j = 0; j < X.nEVars; j ++) {
                if (xevars[j].equals(pevars[i]))
                    xi = j;
            }
            int yi = -1;
            for (int j = 0; j < Y.nEVars; j ++) {
                if (yevars[j].equals(pevars[i]))
                    yi = j;
            }
            if (yi != -1)
                ykey[yi] = pkey[i];
            if (xi != -1)
                xkey[xi] = pkey[i];
        }
        double xval ;
        if (xkey.length != 0)
            xval = X.getValue(X.getIndex(xkey));
        else
            xval = X.getValue();
        double yval;
        if (ykey.length != 0)
            yval = Y.getValue(Y.getIndex(ykey));
        else
            yval = Y.getValue();
        double pval;
        if (pkey.length == 0) 
            pval = PROD.getValue();
        else
            pval = PROD.getValue(p);
        //System.out.print(p + "\t");
        return xval * yval == pval;
    }
    
    public static void main(String[] args) {
        System.out.println("maxEV\tminEV\tOverlap\tContain\tProduct\tPJoin\tTime (ms)");
        for (long seed = 0; seed < 200; seed++) {
            Variable[] vars = getVariablePool(seed, 10);
            AbstractFactor[] dfs = getFactorPool(seed, vars, 8);
            AbstractFactor f1 = productPool(dfs, POOL_OPTION_LINEAR);
            AbstractFactor f2 = productPool(dfs, POOL_OPTION_TREE);
            if (f1 == null && f2 == null)
                continue;
            if (f1.getSize() != f2.getSize())
                System.err.println("Invalid product size");
            if (f1.getSize() == 1) {
                if (f1.getValue() < f2.getValue() * 0.999 || f1.getValue() > f2.getValue() *1.001) {
                    System.err.println("Invalid atomic product: " + f1.getValue() + " v " + f2.getValue());
                    System.exit(1);
                }
            } else {
                for (int i = 0; i < f1.getSize(); i ++) {
                    if (f1.getValue(i) < f2.getValue(i) * 0.999 || f1.getValue(i) > f2.getValue(i) *1.001) {
                        System.err.println("Invalid product: " + f1.getValue(i) + " v " + f2.getValue(i));
                        System.exit(1);
                    }
                }
            }                    
        }
    }
    
}

class AbstractFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public AbstractFactorRuntimeException(String string) {
        message = string;
    }
}
