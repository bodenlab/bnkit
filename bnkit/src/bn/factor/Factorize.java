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

import bn.JDF;
import bn.Predef;
import bn.prob.GaussianDistrib;
import dat.EnumVariable;
import dat.Variable;

import java.util.*;

/**
 * This class is intended to coordinate and implement high-level operations on factors.
 * It contains routines for main types of operations, e.g. product and marginalisation.
 * This class also contains (mostly static) utility routines for transforming lists of variables etc.
 * Central coordination of factor operations should enable multi-threading of processes, 
 * e.g. factor products organised into binary trees, but not yet implemented.
 * @author mikael
 */
public class Factorize {
    
    protected static final double LOG0 = Double.NEGATIVE_INFINITY;
    protected static boolean isLOG0(double x) { return Double.isInfinite(x); }
    
    protected static boolean VERBOSE = false;

    protected static int PRODUCT_OPTION = -1; // choose strategy for complex cases by timing ("-1") or by fixed option (currently "0" and "1")

    /**
     * Empty constructor.
     */
    public Factorize() {
    }

    /**
     * Determine the difference between two arrays.
     * @param A the primary array
     * @param B the secondary array
     * @return the array which contain all elements of A except those in B
     */
    public static Variable[] getDifference(Variable[] A, Variable[] B) {
        List<Variable> accept = new ArrayList<>();
        for (int i = 0; i < A.length; i++) {
            boolean present_in_B = false;
            for (int j = 0; j < B.length; j++) {
                if (A[i].equals(B[j])) {
                    present_in_B = true;
                    break;
                }
            }
            if (!present_in_B) {
                accept.add(A[i]);
            }
        }
        Variable[] C = new Variable[accept.size()];
        accept.toArray(C);
        return C;
    }

    /**
     * Construct a non-redundant array of variables 
     * that contains all of the unique variables in the provided array.
     * @param A potentially unsorted array of variables, potentially with duplicates
     * @return array of the same variables, but with one copy of each unique variable
     */
    protected static Variable[] getNonredundant(Variable[] A) {
        List<Variable> accept = new ArrayList<>();
        for (int i = 0; i < A.length; i++) {
            if (!accept.contains(A[i])) {
                accept.add(A[i]);
            }
        }
        Variable[] B = new Variable[accept.size()];
        accept.toArray(B);
        return B;
    }

    /**
     * Construct an array of variables sorted according to their canonical index,
     * that contains all of the unique variables in the provided array.
     * @param A potentially unsorted array of variables, potentially with duplicates
     * @return array of the same variables, but sorted and unique
     */
    protected static Variable[] getNonredundantSorted(Variable[] A) {
        if (A.length < 2) {
            return A;
        }
        int j = 0;
        int i = 1;
        try {
            Arrays.sort(A);
        } catch (NullPointerException e) {
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

    public static EnumVariable[] getEnumVars(Variable[] A) {
        int cnt_evars = 0;
        for (Variable var : A) {
            try {
                EnumVariable evar = (EnumVariable) var;
                cnt_evars++;
            } catch (ClassCastException e) {
            }
        }
        EnumVariable[] evars = new EnumVariable[cnt_evars];
        if (cnt_evars == 0) {
            return evars;
        } else {
            cnt_evars = 0;
            for (Variable var : A) {
                try {
                    EnumVariable evar = (EnumVariable) var;
                    evars[cnt_evars++] = evar;
                } catch (ClassCastException e) {
                }
            }
        }
        return evars;
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
        if (Xvars == null) {
            Xvars = new EnumVariable[0];
        }
        if (Yvars == null) {
            Yvars = new EnumVariable[0];
        }
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
        if (Xvars == null) {
            Xvars = new EnumVariable[0];
        }
        if (Yvars == null) {
            Yvars = new EnumVariable[0];
        }
        int x = 0;
        int y = 0;
        int multiplier = 1;
        while (Xvars.length > x && Yvars.length > y) {
            if (Xvars[x] == Yvars[y]) {
                if (includeJoin) {
                    multiplier *= Xvars[x].size();
                }
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
    protected static FactorProductTree getProductTree(AbstractFactor[] factors) {
        int N = factors.length;
        // deal with special cases
        if (N == 0) 
            return null;
        if (N == 1) 
            return new FactorProductTree(factors[0]);
        if (N == 2) 
            return new FactorProductTree(new FactorProductTree(factors[0]), new FactorProductTree(factors[1]));
        // more than two so optimisation is performed...
        // construct the initial working set, of:
        // nodes that MUST be used in the final tree
        List<FactorProductTree> fpool = new ArrayList<>();
        // first leaves for all factors that must take part of the tree
        for (int i = 0; i < N; i ++) {
            fpool.add(new FactorProductTree(factors[i]));
        }
        // now, the main part of repeatedly identifying the best node and integrating that into the final tree
        // this will always be N - 1 products/internal nodes (where N is the number of factors)
        FactorProductTree node = null;
        for (int rank = 0; rank < N - 1; rank++) {
            int lowest = Integer.MAX_VALUE;
            int a = -1;
            int b = -1;
            for (int i = 0; i < fpool.size(); i++) {
                EnumVariable[] evars_i = fpool.get(i).getEnumVars();
                for (int j = i + 1; j < fpool.size(); j++) {
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

    public static boolean isValid(AbstractFactor f) {
        int n = f.getSize();
        boolean allNegInf = true;
        if (n == 1) {
            double y = f.getLogValue();
            allNegInf = Double.isInfinite(y);
        } else {
            for (int i = 0; i < n; i ++) {
                double y = f.getLogValue(i);
                allNegInf = allNegInf && Double.isInfinite(y);
            }
        }
        return !allNegInf;
    }
    
    public static void exitIfInvalid(AbstractFactor f, String msg) {
        if (!isValid(f)) {
            System.err.println("Invalid factor: " + f + " (" + msg + ")");
            f.display();
            System.exit(111);
        }
    }
    /**
     * Recursive function to perform a single product of factors, potentially resulting 
     * from products themselves.
     * @return the product of the factors represented by this node in a binary tree
     */
    protected static AbstractFactor getProduct(FactorProductTree node) {
        if (node.getFactor() != null) {
            return node.getFactor();
        }
        AbstractFactor X = getProduct(node.x);
        AbstractFactor Y = getProduct(node.y);
        AbstractFactor f = getProduct(X, Y);
        //exitIfInvalid(f, "getProduct");
        node.setFactor(f);
        return f;
    }

    /**
     * Perform products of factors, by organising them into a binary tree of single products,
     * which minimises the number of elements in the resulting, intermittent factors.
     * @param factors the factors that are multiplied
     * @return the product of the factors
     */
    public static AbstractFactor getProduct(AbstractFactor[] factors) {
        if (factors.length == 0) {
            return null;
        }
        FactorProductTree tree = getProductTree(factors);
        return getProduct(tree);
    }


    /**
     * Concatenate two arrays of variables into one.
     * Does not consider order or duplicates.
     * @param one
     * @param two
     * @return
     */
    protected static Variable[] getConcat(Variable[] one, Variable[] two) {
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
    public static int getCrossref(Variable[] xvars, int[] xcross, Variable[] yvars, int[] ycross) {
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
            if (xcross.length > 0) {
                if (xvars.length != xcross.length) {
                    throw new RuntimeException("Lists must be equal");
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
            }
        } else if (ycross != null) {
            // xcross == null
            if (ycross.length > 0) {
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
     * The product is performed entirely in log space. Underflow issues will only arise in extreme 
     * situations.
     *
     * TODO: Make informed choices as to what implementation of AbstractFactor should be used.
     * Currently only DenseFactor is used.
     *
     * @param X one table
     * @param Y other table
     * @return the product of one and the other table
     */
    public static AbstractFactor getProduct(AbstractFactor X, AbstractFactor Y) {
        // System.err.println("Factor product " + X.toString() + " x " + Y.toString());
        // First resolve cases with tables without enumerable variables
        if (X.nEVars == 0 && Y.nEVars == 0) {
            // both tables lack enumerable variables
            AbstractFactor dt = new DenseFactor(getConcat(X.nvars, Y.nvars));
            if (X.isTraced() || Y.isTraced()) {
                dt.setTraced(true);
            }
            dt.setLogValue(X.getLogValue() + Y.getLogValue());
            if (X.isJDF() && Y.isJDF()) {
                dt.setJDF(JDF.combine(X.getJDF(), Y.getJDF()));
            } else if (X.isJDF()) {
                dt.setJDF(X.getJDF());
            } else if (Y.isJDF()) {
                dt.setJDF(Y.getJDF());
            }
            if (X.isTraced()) {
                dt.addAssign(X.getAssign());
            }
            if (Y.isTraced()) {
                dt.addAssign(Y.getAssign());
            }
            return dt;
        } else if (X.nEVars == 0) {
            // only X is without enumerables
            AbstractFactor dt = new DenseFactor(getConcat(Y.evars, getConcat(X.nvars, Y.nvars)));
            if (X.isTraced() || Y.isTraced()) {
                dt.setTraced(true);
            }
            for (int j = 0; j < Y.getSize(); j++) {
                dt.setLogValue(j, Y.getLogValue(j) + X.getLogValue());
                if (X.isJDF() && Y.isJDF()) {
                    dt.setJDF(j, JDF.combine(X.getJDF(), Y.getJDF(j)));
                } else if (X.isJDF()) {
                    dt.setJDF(j, X.getJDF());
                } else if (Y.isJDF()) {
                    dt.setJDF(j, Y.getJDF(j));
                }
                if (X.isTraced()) {
                    dt.addAssign(j, X.getAssign());
                }
                if (Y.isTraced()) {
                    dt.addAssign(j, Y.getAssign(j));
                }
            }
            return dt;
        } else if (Y.nEVars == 0) {
            // only Y is without enumerables
            AbstractFactor dt = new DenseFactor(getConcat(X.evars, getConcat(X.nvars, Y.nvars)));
            if (X.isTraced() || Y.isTraced()) {
                dt.setTraced(true);
            }
            for (int j = 0; j < X.getSize(); j++) {
                dt.setLogValue(j, X.getLogValue(j) + Y.getLogValue());
                if (X.isJDF() && Y.isJDF()) {
                    dt.setJDF(j, JDF.combine(X.getJDF(j), Y.getJDF()));
                } else if (X.isJDF()) {
                    dt.setJDF(j, X.getJDF(j));
                } else if (Y.isJDF()) {
                    dt.setJDF(j, Y.getJDF());
                }
                if (X.isTraced()) {
                    dt.addAssign(j, X.getAssign(j));
                }
                if (Y.isTraced()) {
                    dt.addAssign(j, Y.getAssign());
                }
            }
            return dt;
        }
        // Both tables have enumerable variables, so let's work out how they relate
        int[] xcross2y = new int[X.nEVars]; // map from X to Y indices [x] = y
        int[] ycross2x = new int[Y.nEVars]; // map from Y to X indices [y] = x
        int noverlap = 0; // number of overlapping variables
        noverlap = getCrossref(X.evars, xcross2y, Y.evars, ycross2x);
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
        // Before handling complex scenarios, consider some special, simpler cases
        // 1. two tables with the same variables in the same order
        if (noverlap == Math.min(X.nEVars, Y.nEVars)) {
            // at least one table is "contained" by the other
            if (X.nEVars == Y.nEVars) {
                Object[] ykey = new Object[Y.nEVars];
                AbstractFactor dt = new DenseFactor(getConcat(X.evars, getConcat(X.nvars, Y.nvars)));
                if (X.isTraced() || Y.isTraced()) {
                    dt.setTraced(true);
                }
                for (int x = 0; x < X.getSize(); x++) {
                    int y = x; // if ordered the indices in the tables are identical
                    if (!ordered) {
                        // re-index since the variables are listed in a different order
                        Object[] xkey = X.getKey(x);
                        for (int i = 0; i < xkey.length; i++) {
                            ykey[xcross2y[i]] = xkey[i];
                        }
                        y = Y.getIndex(ykey);
                    }
                    double xval = X.getLogValue(x);
                    double yval = Y.getLogValue(y);
                    double XY = xval + yval;  // <=== Factor product in log space
                    dt.setLogValue(x, XY);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(x, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(x, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(x, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.addAssign(x, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.addAssign(x, Y.getAssign(y));
                    }
                }
                if (VERBOSE) {
                    System.err.println("DT: Complete overlap, ordered = " + ordered + " : " + X.toString() + " x " + Y.toString());
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
                AbstractFactor dt = new DenseFactor(getConcat(X.evars, getConcat(X.nvars, Y.nvars)));
                if (X.isTraced() || Y.isTraced()) {
                    dt.setTraced(true);
                }
                for (int x = 0; x < X.getSize(); x++) {
                    double xval = X.getLogValue(x);
                    if (isLOG0(xval)) {
                        continue; // no point in continuing since product will always be log-zero for entries with this x-value
                        // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
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
                    double yval = Y.getLogValue(y);
                    if (isLOG0(yval)) { // <=== Factor check
                        continue; // the product will be zero no what
                        // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
                    }
                    int idx = x;
                    double XY = xval + yval; // <=== Factor product in log space
                    dt.setLogValue(idx, XY);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.addAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.addAssign(idx, Y.getAssign(y));
                    }
                }
                if (VERBOSE) {
                    System.err.println("DT: Partial overlap (X>Y), ordered = " + ordered + " : " + X.toString() + " x " + Y.toString());
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
                AbstractFactor dt = new DenseFactor(getConcat(Y.evars, getConcat(X.nvars, Y.nvars)));
                if (X.isTraced() || Y.isTraced()) {
                    dt.setTraced(true);
                }
                for (int y = 0; y < Y.getSize(); y++) {
                    double yval = Y.getLogValue(y);
                    if (isLOG0(yval)) { // <=== Factor check
                        continue; // no point in continuing since product will always be zero for entries with this x-value
                        // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
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
                    double xval = X.getLogValue(x);
                    if (isLOG0(xval)) { // <=== Factor check
                        continue; // the product will be zero no what
                        // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
                    }
                    int idx = y;
                    double XY = xval + yval; // <=== Factor product in log space
                    dt.setLogValue(idx, XY);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.addAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.addAssign(idx, Y.getAssign(y));
                    }
                }
                if (VERBOSE) {
                    System.err.println("DT: Partial overlap (X<Y), ordered = " + ordered + " : " + X.toString() + " x " + Y.toString());
                }
                return dt;
            }
        }
        // Failing the above, we must construct a table which is an amalgamate of the two.
        AbstractFactor dt = new DenseFactor(getConcat(getConcat(X.evars, Y.evars), getConcat(X.nvars, Y.nvars)));
        if (X.isTraced() || Y.isTraced()) {
            dt.setTraced(true);
        }
        Object[] reskey = new Object[dt.nEVars]; // this will initialise all elements to null
        int[] xcross2dt = new int[X.nEVars]; // map from X index to dt index
        int[] ycross2dt = new int[Y.nEVars]; // map from Y index to dt index
        getCrossref(X.evars, xcross2dt, dt.evars, null);
        getCrossref(Y.evars, ycross2dt, dt.evars, null);
        // 2. two tables have nothing in common
        if (noverlap == 0) {
            for (int x = 0; x < X.getSize(); x++) {
                double xval = X.getLogValue(x);
                if (isLOG0(xval)) { // <=== Factor check
                    continue; // no point in continuing since product will always be zero for entries with this x-value
                    // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
                }
                Object[] xkey = X.getKey(x);
                for (int y = 0; y < Y.getSize(); y++) {
                    double yval = Y.getLogValue(y);
                    if (isLOG0(yval)) { // <=== Factor check
                        continue; // the product will be zero no what
                        // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
                    }
                    Object[] ykey = Y.getKey(y);
                    for (int i = 0; i < xkey.length; i++) {
                        reskey[xcross2dt[i]] = xkey[i];
                    }
                    for (int i = 0; i < ykey.length; i++) {
                        reskey[ycross2dt[i]] = ykey[i];
                    }
                    int idx = dt.getIndex(reskey);
                    // **************************
                    //     This does NOT work??
                    // **************************
                    double XY = xval + yval; // <=== Factor product in log space
                    dt.setLogValue(idx, XY);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.addAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.addAssign(idx, Y.getAssign(y));
                    }
                }
            }
            if (VERBOSE) {
                System.err.println("DT: No overlap : " + X.toString() + " x " + Y.toString());
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
            double xval = X.getLogValue(x);
            if (isLOG0(xval)) { // <=== Factor check
                continue; // no point in continuing since product will always be zero for entries with this x-value
                // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
            }
            Object[] xkey = X.getKey(x);
            for (int i = 0; i < xcross2y.length; i++) {
                if (xcross2y[i] > - 1) {
                    searchkey[xcross2y[i]] = xkey[i];
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
                    double yval = Y.getLogValue(y);
                    if (isLOG0(yval)) { // <=== Factor check
                        continue; // the product will be zero no what
                        // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
                    }
                    Object[] ykey = Y.getKey(y);
                    for (int i = 0; i < ykey.length; i++) {
                        reskey[ycross2dt[i]] = ykey[i];
                    }
                    int idx = dt.getIndex(reskey);
                    double XY = xval + yval; // <=== Factor product in log space
                    dt.setLogValue(idx, XY);
                    if (X.isJDF() && Y.isJDF()) {
                        dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                    } else if (X.isJDF()) {
                        dt.setJDF(idx, X.getJDF(x));
                    } else if (Y.isJDF()) {
                        dt.setJDF(idx, Y.getJDF(y));
                    }
                    if (X.isTraced()) {
                        dt.addAssign(idx, X.getAssign(x));
                    }
                    if (Y.isTraced()) {
                        dt.addAssign(idx, Y.getAssign(y));
                    }
                }
                if (start0 != 0) {
                    total0 = System.nanoTime() - start0;
                    option = -1;
                }
            } else if (option == 1) {
                for (int y = 0; y < Y.getSize(); y++) {
                    if (Y.isMatch(searchkey, y)) {
                        double yval = Y.getLogValue(y);
                        if (isLOG0(yval)) { // <=== Factor check
                            continue; // the product will be zero no what
                            // Note: it is OK to leave the cell in the factortable un-set since it was initialised to LOG0
                        }
                        Object[] ykey = Y.getKey(y);
                        for (int i = 0; i < ykey.length; i++) {
                            reskey[ycross2dt[i]] = ykey[i];
                        }
                        int idx = dt.getIndex(reskey);
                        double XY = xval + yval; // <=== Factor product in log space
                        dt.setLogValue(idx, XY);
                        if (X.isJDF() && Y.isJDF()) {
                            dt.setJDF(idx, JDF.combine(X.getJDF(x), Y.getJDF(y)));
                        } else if (X.isJDF()) {
                            dt.setJDF(idx, X.getJDF(x));
                        } else if (Y.isJDF()) {
                            dt.setJDF(idx, Y.getJDF(y));
                        }
                        if (X.isTraced()) {
                            dt.addAssign(idx, X.getAssign(x));
                        }
                        if (Y.isTraced()) {
                            dt.addAssign(idx, Y.getAssign(y));
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
            System.err.println("DT: Generic case: " + X.toString() + " x " + Y.toString() + " Option = " + option + " (Option 0 took " + total0 + "ns. Option 1 took " + total1 + "ns.)");
        }
        return dt;
    }

    /**
     * Determine log(x + y) from log(x) and log(y).
     * While this method does perform exponentiation, it does so by minimising risk of underflow.
     * @param logx
     * @param logy
     * @return the log of a sum
     */
    public static double logSumOfLogs(double logx, double logy) {
        if (isLOG0(logx) && !isLOG0(logy))
            return logy;
        else if (!isLOG0(logx) && isLOG0(logy))
            return logx;
        else if (isLOG0(logx) && isLOG0(logy)) // added to overcome NaNs in the calculation below (ask Alex)
            return LOG0;
        if (logx > logy)
            return logx + Math.log(1.0 + Math.exp(logy - logx));
        else
            return logy + Math.log(1.0 + Math.exp(logx - logy));
    }
    
    /**
     * Construct a new table from an existing, by summing-out specified variable/s.
     * This operation is performed primarily in log space, and tries to avoid numerical underflow.
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
        Variable[] uniqueVars = getNonredundantSorted(anyvars);
        // first check that X has enumerable variables, because that would be required
        if (X.nEVars == 0) {
            return X;
        }
        // if X is ok, get rid of non-enumerables
        EnumVariable[] evars = getEnumVars(uniqueVars);
        // now, resolve the relationships between X's variables and those to be summed-out
        Variable[] yvars = getDifference(getConcat(X.evars, X.nvars), evars);
        // construct new table
        AbstractFactor Y = new DenseFactor(yvars);
        if (Y.getSize() == 1) {
            // we are creating a factor with no enumerable variables)
            double logsum = 1; // cumulative in log space
            JDF first = null;
            JDF subsequent = null;
            JDF mixture = null;
            double prev_weight = 0;
            for (int x = 0; x < X.getSize(); x++) {
                double xval = X.getLogValue(x);
                if (logsum > 0)
                    logsum = xval;
                else
                    logsum = logSumOfLogs(logsum, xval); // <=== Factor addition in log space
            }
            Y.setLogValue(logsum);  // <=== still in log space
            if (X.isJDF()) {
                for (int x = 0; x < X.getSize(); x++) {
                    double lognxval = X.getLogValue(x) - logsum;
                    if (!isLOG0(lognxval)) {
                        double nxval = Math.exp(lognxval); // normalization in log space
                        if (first == null) {
                            // this will be true only once, for the first entry that will be collapsed
                            first = X.getJDF(x); // save the JDF for later
                            prev_weight = nxval; // save the value associated with this entry to weight the distributions
                        } else if (mixture == null) {
                            // for the subsequent entry, the mixture is created, both JDF weighted
                            subsequent = X.getJDF(x);
                            if (subsequent != null) {
                                mixture = JDF.mix(first, prev_weight, subsequent, nxval);
                            }
                        } else {
                            // for all subsequent entries, the mixture is mixed unchanged with the new entry's JDF weighted
                            subsequent = X.getJDF(x);
                            if (subsequent != null) {
                                mixture = JDF.mix(mixture, subsequent, nxval);
                            }
                        }
                    }
                }
                if (first != null) {
                    if (mixture != null) {
                        Y.setJDF(mixture);
                    } else {
                        Y.setJDF(first);
                    }
                }
            }
        } else {
            // work out how variables in Y map to X so that we can search X from each index in Y
            int[] ycross2x = new int[Y.nEVars];
            getCrossref(Y.evars, ycross2x, X.evars, null);
            Object[] xkey_search = new Object[X.nEVars];
            for (int y = 0; y < Y.getSize(); y++) {
                // FIXME: loop's not efficient for sparse factors, start with X, be selective about Y?
                Object[] ykey = Y.getKey(y);
                for (int i = 0; i < ykey.length; i++) {
                    xkey_search[ycross2x[i]] = ykey[i];
                }
                double logsum = 1; // cumulative in log space
                int[] indices = X.getIndices(xkey_search);
                JDF first = null;
                JDF subsequent = null;
                JDF mixture = null;
                double prev_weight = 0;
                for (int x : indices) {
                    double xval = X.getLogValue(x);
                    if (logsum > 0)
                        logsum = xval;
                    else
                        logsum = logSumOfLogs(logsum, xval); // <=== Factor addition in log space
                }
                Y.setLogValue(y, logsum);
                if (X.isJDF()) {
                    for (int x : indices) {
                        double lognxval = X.getLogValue(x) - logsum;
                        if (!isLOG0(lognxval)) {
                            double nxval = Math.exp(lognxval); // normalization in log space
                            if (first == null) {
                                // this will be true only once, for the first entry that will be collapsed
                                first = X.getJDF(x); // save the JDF for later
                                prev_weight = nxval; // save the value associated with this entry to weight the distributions
                            } else if (mixture == null) {
                                // for the subsequent entry, the mixture is created, both JDF weighted
                                subsequent = X.getJDF(x);
                                if (subsequent != null) {
                                    mixture = JDF.mix(first, prev_weight, subsequent, nxval);
                                }
                            } else {
                                // for all subsequent entries, the mixture is mixed unchanged with the new entry's JDF weighted
                                subsequent = X.getJDF(x);
                                if (subsequent != null) {
                                    mixture = JDF.mix(mixture, subsequent, nxval);
                                }
                            }
                        }
                    }
                    if (first != null) {
                        if (mixture != null) {
                            Y.setJDF(y, mixture);
                        } else {
                            Y.setJDF(y, first);
                        }
                    }
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
     * This operation is performed entirely in log space.
     *
     * @param X existing table
     * @param anyvars variables to max-out
     * @return the resulting factor
     */
    public static AbstractFactor getMaxMargin(AbstractFactor X, Variable... anyvars) {
        // sort and weed out duplicates
        Variable[] uniqueVars = getNonredundantSorted(anyvars);
        // first check that X has enumerable variables, because that would be required
        if (X.nEVars == 0) {
            return X;
        }
        // if X is ok, get rid of non-enumerables
        EnumVariable[] evars = getEnumVars(uniqueVars);
        int[] ecross2x = new int[evars.length];
        getCrossref(evars, ecross2x, X.evars, null); // so we can map maxed-out vars to the original table
        // now, resolve the relationships between X's variables and those to be summed-out
        Variable[] yvars = getDifference(getConcat(X.evars, X.nvars), evars);
        // construct new table
        AbstractFactor Y = new DenseFactor(yvars);
        Y.setTraced(true); // make sure we save assignments made in max-outs
        if (!Y.hasEnumVars()) {
            double max = Double.NEGATIVE_INFINITY;
            int maxidx = 0;
            for (int x = 0; x < X.getSize(); x++) {
                double xval = X.getLogValue(x);
                if (xval > max) {
                    max = xval;
                    maxidx = x;
                }
            }
            Y.setLogValue(max);
            if (Y.isJDF()) {
                Y.setJDF(X.getJDF(maxidx));
            }
            // transfer old assignments
            if (X.isTraced())
                Y.addAssign(X.getAssign(maxidx));
            // Trace implied assignments
            Object[] xkey = X.getKey(maxidx);
            for (int i = 0; i < evars.length; i++) {
                Variable.Assignment a = new Variable.Assignment(evars[i], xkey[ecross2x[i]]);
                Y.addAssign(a);
            }
        } else {
            // work out how variables in Y map to X so that we can search X from each index in Y
            int[] ycross2x = new int[Y.nEVars];
            int[] xcross2y = new int[X.nEVars];
            getCrossref(Y.evars, ycross2x, X.evars, xcross2y);
            Object[] xkey_search = new Object[X.nEVars];
            for (int y = 0; y < Y.getSize(); y++) {
                // FIXME: loop's not efficient for sparse factors
                Object[] ykey = Y.getKey(y);
                for (int i = 0; i < ykey.length; i++) {
                    xkey_search[ycross2x[i]] = ykey[i];
                }
                double max = Double.NEGATIVE_INFINITY;
                int maxidx = 0;
                int[] indices = X.getIndices(xkey_search);
                for (int x : indices) {
                    double xval = X.getLogValue(x);
                    if (xval > max) {
                        max = xval;
                        maxidx = x;
                    }
                }
                Y.setLogValue(y, max);
                if (Y.isJDF()) {
                    Y.setJDF(y, X.getJDF(maxidx));
                }
                // transfer old assignments
                if (X.isTraced())
                    Y.addAssign(y, X.getAssign(maxidx));
                // Trace implied assignments
                Object[] xkey = X.getKey(maxidx);
                for (int i = 0; i < evars.length; i++) {
                    if (ecross2x[i] >= 0) { // check so that the enumerable variable can be assigned by X
                        Variable.Assignment a = new Variable.Assignment(evars[i], xkey[ecross2x[i]]);
                        Y.addAssign(y, a);
                    }
                }
            }
        }
        return Y;
    }

    private static class NominateTree {
        final Variable[] vars;
        int N;
        NominateTree(Variable... vars) {
            this.vars = vars;
        }
        Set<Variable[]> nominate(int N) {
            if (N > vars.length)
                throw new RuntimeException("Invalid N=" + N + " exceeds number of variables in factor");
            this.N = N;
            Set<Integer> set = getBinary(0, 0);
            Set<Variable[]> all = new HashSet<>();
            for (Integer b : set) {
                Variable[] list = new Variable[vars.length - N];
                int idx = 0;
                for (int i = 1; i <= vars.length; i ++) {
                    if ((b & (1 << i)) == 0) // not set, so sum-out
                        list[idx ++] = vars[i - 1];
                }
                all.add(list);
            }
            return all;
        }
        Set<Integer> getBinary(int n, int b) {
            Set<Integer> ret = new HashSet<>();
            if (n == N) {
                ret.add(b);
            } else {
                for (int i = 1; i <= vars.length; i ++) {
                    if ((b & (1 << i)) != 0) // already set
                        continue;
                    int bb = b;
                    bb |= (1 << i);
                    ret.addAll(getBinary(n + 1, bb));
                }
            }
            return ret;
        }
    }

    /**
     * Decompose a factor into N-term sub-factors.
     * @param f the factor
     * @param N the number of variables in the decomposed factors
     * @return an array of sub-factors
     */
    public static AbstractFactor[] decompose(AbstractFactor f, int N) {
        NominateTree ntree = new NominateTree(f.evars);
        Set<Variable[]> permut = ntree.nominate(N);
        if (permut == null)
            return null;
        AbstractFactor[] ff = new AbstractFactor[permut.size()];
        int cnt = 0;
        for (Variable[] one : permut) {
            ff[cnt] = Factorize.getMargin(f, one);
            cnt += 1;
        }
        return ff;
    }

    private static double log2(double x) {
        return Math.log(x) / Math.log(2);
    }

    /**
     * Determine the Mutual Information (en.wikipedia.org/wiki/Mutual_information)
     * between all pairs of variables in the factor.
     * The matrix is ordered by variables as ordered in the factor.
     * @param f the factor
     * @return the standard mutual information score
     */
    public static double[][] getMIMatrix(AbstractFactor f) {
        double[][] ret = new double[f.nEVars][f.nEVars];
        AbstractFactor[] f1 = decompose(f, 1);
        Map<Variable, Integer> f1map = new HashMap<>();
        for (int i = 0; i < f.nEVars; i ++) {
            Variable v1 = f.evars[i];
            for (int k = 0; k < f1.length; k ++) {
                if (f1[k].hasVariable(v1)) {
                    f1map.put(v1, k);
                    break;
                }
            }
        }
        AbstractFactor[] f2 = decompose(f, 2);
        for (int i = 0; i < f.nEVars; i ++) {
            Variable v1 = f.evars[i];
            for (int j = 0; j < i; j ++) {
                Variable v2 = f.evars[j];
                AbstractFactor myf2 = null;
                for (int k = 0; k < f2.length; k++) {
                    if (f2[k].hasVariable(v1) && f2[k].hasVariable(v2)) {
                        myf2 = f2[k];
                        break;
                    }
                }
                if (myf2 == null)
                    throw new RuntimeException("Factor decomposition error in getMIMatrix");
                // loop through all possible assignments of v1 and v2, to calculate MI
                double sum = 0;
                for (int index : myf2) {
                    Object[] values = myf2.getKey(index);
                    double p_joint = myf2.getValue(index);
                    int v1idx = f1map.get(myf2.evars[0]);
                    double p_1 = f1[v1idx].getValue(new Object[] {values[0]});
                    int v2idx = f1map.get(myf2.evars[1]);
                    double p_2 = f1[v2idx].getValue(new Object[] {values[1]});
                    sum += p_joint * log2(p_joint / (p_1 * p_2));
                }
                ret[i][j] = ret[j][i] = sum;
            }
        }
        return ret;
    }

    /**
     * Determine a normalised version of the factor.
     * This operation is performed primarily in log space, and is very 
     * unlikely to cause numerical underflow.
     * 
     * @param X factor
     * @return a normalised copy of the provided factor
     */
    public static AbstractFactor getNormal(AbstractFactor X) {
        AbstractFactor Y = new DenseFactor(getConcat(X.evars, X.nvars));
        if (X.isTraced())
            Y.setTraced(true);
        if (X.hasEnumVars()) {
            double logmax = Double.NEGATIVE_INFINITY;
            double logmin = Double.POSITIVE_INFINITY;
            for (int i = 0; i < X.getSize(); i++) {
                double xval = X.getLogValue(i);
                if (xval > logmax)
                    logmax = xval;
                if (xval < logmin)
                    logmin = xval;
            }
            double sum = 0;
            double[] yvals = new double[X.getSize()];
            for (int i = 0; i < X.getSize(); i++) {
                yvals[i] = X.getLogValue(i) - logmax;
                if (i == 0)
                    sum = yvals[i];
                else
                    sum = logSumOfLogs(sum, yvals[i]);
            }
            if (isLOG0(sum))
                System.err.println("getNormal fails");

            for (int i = 0; i < X.getSize(); i++) {
                Y.setLogValue(i, yvals[i] - sum);
                if (X.isJDF()) {
                    Y.setJDF(i, X.getJDF(i));
                }
                if (X.isTraced()) {
                    Y.addAssign(i, X.getAssign(i));
                }
            }
        } else {
            Y.setLogValue(X.getLogValue());
            if (X.isJDF()) {
                Y.setJDF(X.getJDF());
            }
            if (X.isTraced()) {
                Y.addAssign(X.getAssign());
            }
        }
        return Y;
    }


    /**
     * To represent a pair of factors. Originally intended to cache products, currently not used.
     */
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
    
    /**
     * Class to structure the binary tree of factor products.
     */
    static public class FactorProductTree {
        public final FactorProductTree x, y;
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
            Variable[] vars = Factorize.getNonredundantSorted(Factorize.getConcat(x.getEnumVars(), y.getEnumVars()));
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
    
}

