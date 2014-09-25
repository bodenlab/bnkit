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
package bn;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Table for storing and retrieving double based on Enumerable keys.
 * It can be used efficiently when data are dense, i.e. when almost all possible key 
 * combinations are used, though storage requirements grow exponentially with number of 
 * variables.
 * The class is sensitive to the order of variables, and exposes the user to other
 * details so caution must be exercised.
 *
 * @author mikael
 */
public class DenseFactor {

    protected final double[] map; // the factors
    protected Set<Variable.Assignment>[] assigned = null; // the latent assignments associated with each factor, disabled by default
    protected JDF[] jdf = null; // the Joint Density Function associated with each factor, disabled by default
    protected final int nEVars; // number of enumerable variables
    protected final int nNVars; // number of non-enumerable variables
    protected final EnumVariable[] evars; // enumerable variables, which form the "keys" to entries in the map
    protected final Variable[] nvars; // non-enumerable variables, which densities are attached to the entries via JDF
    protected final int[] period;
    protected final int[] step;
    protected final int[] domsize; // size of domain

    protected static int PRODUCT_OPTION = -1; // choose strategy for complex cases by timing ("-1") or by fixed option (currently "0" and "1")
    protected static boolean VERBOSE = false;
    
    /**
     * Construct a new table without any variables.
     */
    public DenseFactor() {
        this.evars = null;
        this.nEVars = 0;
        this.nvars = null;
        this.nNVars = 0;
        this.step = null;
        this.period = null;
        this.domsize = null;
        this.map = new double[1];
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
     * @param useVariables variables, either enumerable or non-enumerable
     */
    public DenseFactor(Variable... useVariables) {
        int cnt_evars = 0, cnt_nvars = 0;
        for (Variable useVariable : useVariables) {
            try {
                EnumVariable evar = (EnumVariable) useVariable;
                cnt_evars ++;
            } catch (ClassCastException e) {
                cnt_nvars ++;
            }
        }
        if (cnt_nvars > 0 && cnt_evars > 0) { 
            this.evars = new EnumVariable[cnt_evars];
            this.nvars = new Variable[cnt_nvars];
            cnt_evars = 0; cnt_nvars = 0;
            for (Variable useVariable : useVariables) {
                try {
                    EnumVariable evar = (EnumVariable) useVariable;
                    this.evars[cnt_evars ++] = evar;
                } catch (ClassCastException e) {
                    this.nvars[cnt_nvars ++] = useVariable;
                }
            }
            this.nEVars = evars.length;
            this.nNVars = nvars.length;
        } else if (cnt_evars > 0) { // only enumerable
            this.nNVars = 0;
            this.nvars = null;
            this.evars = new EnumVariable[cnt_evars];
            System.arraycopy(useVariables, 0, this.evars, 0, cnt_evars);
            this.nEVars = cnt_evars;
        } else { // only non-enumerable
            this.nEVars = 0;
            this.evars = null;
            this.nvars = new Variable[cnt_nvars];
            System.arraycopy(useVariables, 0, this.nvars, 0, cnt_nvars);
            this.nNVars = cnt_nvars;
        }
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
        if (this.nEVars > 0)
            this.map = new double[period[0]];
        else // no enumerable variables, so need only one entry
            this.map = new double[1];
        if (this.nNVars > 0) {
            this.jdf = new JDF[this.map.length];
            for (int i = 0; i < this.jdf.length; i ++) 
                this.jdf[i] = new JDF(nvars);
        }
    }

    /**
     * Activate tracing of implied assignments.
     * @param status true if activated, false otherwise
     */
    public void setTraced(boolean status) {
        if (status) 
            this.assigned = new Set[this.getSize()];
        else
            this.assigned = null;
    }
    
    /**
     * Find out if tracing of implied assignments is active.
     * @return true if activated, false otherwise
     */
    public boolean isTraced() {
        return (this.assigned != null);
    }
    
    public boolean isJDF() {
        return (this.jdf != null);
    }
    
    /**
     * Retrieve the value of the factor without enumerable variables.
     * @return the only value of the factor
     */
    public double getValue() {
        if (this.getSize() == 1)
            return this.map[0];
        throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }
    
    /**
     * Retrieve the distribution of the specified non-enumerable variable.
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public Distrib getDistrib(Variable nvar) {
        if (this.getSize() == 1)
            return this.jdf[0].getDistrib(nvar);
        throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }
    
    /**
     * Retrieve the value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    public double getValue(int index) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index");
        double value = map[index];
        return value;
    }

    /**
     * Retrieve the distribution of the specified non-enumerable variable, 
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public Distrib getDistrib(int index, Variable nvar) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index");
        return this.jdf[index].getDistrib(nvar);
    }
    
    /**
     * Retrieve the value of the entry identified by the instantiated key
     *
     * @param key the entry
     * @return the value of the entry
     */
    public double getValue(Object[] key) {
        if (getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key: no variables");
        if (key.length != nEVars)
            throw new DenseFactorRuntimeException("Invalid key: length is " + key.length + " not " + nEVars);
        int index = getIndex(key);
        return getValue(index);
    }

    /**
     * Retrieve the distribution of the specified non-enumerable variable, 
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param key the key of the table, containing the instantiation
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public Distrib getDistrib(Object[] key, Variable nvar) {
        if (getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key: no variables");
        if (key.length != nEVars)
            throw new DenseFactorRuntimeException("Invalid key: length is " + key.length + " not " + nEVars);
        int index = this.getIndex(key);
        return this.getDistrib(index, nvar);
    }
    
    /**
     * Set the only value associated with a table without enumerable variables.
     * @param value
     * @return 
     */
    public int setValue(double value) {
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
        map[0] = value;
        return 0;
    }
    
    /**
     * Set the distribution of a non-enumerable variable for a factor without enumerable variables. 
     * @param nvar the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry, always 0
     */
    public int setDistrib(Variable nvar, Distrib d) {
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
        this.jdf[0].setDistrib(d, nvar);
        return 0;
    }
    
    /**
     * Add a non-enumerable distribution as required by marginalization (removal) of enumerable variables.
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @param weight the weight assigned to the distribution
     * @return the index for the key
     */
    public int addDistrib(Variable nonenum, Distrib d, double weight) {
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
        this.jdf[0].mixDistrib(nonenum, d, weight);
        return 0; // no index since there are no enum variables
    }
    
    /**
     * Associate the specified key-index with the given value. Note that using
     * getValue and setValue with index is quicker than with key, if more than
     * one operation is done.
     *
     * @param key_index
     * @param value
     * @return the index at which the value was stored
     */
    public int setValue(int key_index, double value) {
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        map[key_index] = value;
        return key_index;
    }

    /**
     * Set the conditional distribution of a non-enumerable variable. 
     * Use with care as the method may not do all integrity checks.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(int key_index, Variable nonenum, Distrib d) {
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        this.jdf[key_index].setDistrib(d, nonenum);
        return key_index; // 
    }

    /**
     * Associate the specified key with the given value.
     * @param key
     * @param value
     * @return the index at which the value was stored
     */
    public int setValue(Object[] key, double value) {
        if (getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key: no variables");
        int index = this.getIndex(key);
        map[index] = value;
        return index;
    }

    /**
     * Set the conditional distribution of a non-enumerable variable. 
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(Object[] key, Variable nonenum, Distrib d) {
        if (this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key: no variables");
        int index = this.getIndex(key);
        return setDistrib(index, nonenum, d);
    }

    public int addAssign(int key_index, Variable.Assignment assign) {
        if (!isTraced()) 
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        if (assigned[key_index] == null)
            assigned[key_index] = new HashSet<>();
        assigned[key_index].add(assign);
        return key_index;
    }
    
    public int addAssign(Object[] key, Variable.Assignment assign) {
        return addAssign(this.getIndex(key), assign);
    }
    
    
    /**
     * Determine the mixture of distributions for all entries that match key
     * @param key the specific entries as identified by a (partial) key
     * @return the joint mixture distribution
     */
//    public JDF getJDFSum(Object[] key) {
//        if (densityTable != null) {
//            return getJDFSum(densityTable.getIndices(key));
//        } else {
//            return atomicDensity;
//        }
//    }
    
    /**
     * Determine the mixture of distributions for all entries that match key
     * @param indices the specific entries as identified by their indices
     * @return the joint mixture distribution
     */
//    public JDF getJDFSum(int[] indices) {
//        double sum = getSum(indices);
//        if (densityTable != null) {
//            JDF jdf = new JDF(getNonEnumVariables());
//            for (int index : indices) {
//                JDF current = densityTable.getValue(index);
//                double weight = factorTable.getValue(index) / sum;
//                jdf = JDF.mix(jdf, current, weight); 
//            }
//            return jdf;
//        } else {
//            return atomicDensity;
//        }
//    }
    
    private static Variable[] concat(Variable[] one, Variable[] two) {
        if (one == null && two == null)
            return new Variable[0];
        else if (one == null) {
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
     * Construct a new table that is the result of a factor product of the two specified tables.
     * Works in the general case but there are significant efficiency gains if joint variables are ordered.
     * @param X one table
     * @param Y other table
     * @return the product of one and the other table
     */
    public static DenseFactor getProduct(DenseFactor X, DenseFactor Y) {
        // First resolve cases with tables without enumerable variables
        if (X.nEVars == 0 && Y.nEVars == 0) { // both tables lack enumerable variables
            DenseFactor dt = new DenseFactor();
            dt.setValue(X.getValue() * Y.getValue());
            if (X.isJDF() && Y.isJDF())
                dt.jdf[0] = JDF.combine(X.jdf[0], Y.jdf[0]);
            else if (X.isJDF())
                dt.jdf[0] = X.jdf[0];
            else if (Y.isJDF()) 
                dt.jdf[0] = Y.jdf[0];
            if (X.isTraced())
                dt.assigned[0].addAll(X.assigned[0]);
            if (Y.isTraced())
                dt.assigned[0].addAll(Y.assigned[0]);
            return dt;
        } else if (X.nEVars == 0) { // only X is without enumerables
            DenseFactor dt = new DenseFactor(concat(Y.evars, Y.nvars));
            for (int j = 0; j < Y.getSize(); j ++) {
                dt.map[j] = Y.getValue(j) * X.getValue();
                if (X.isJDF() && Y.isJDF())
                    dt.jdf[j] = JDF.combine(X.jdf[0], Y.jdf[j]);
                else if (X.isJDF())
                    dt.jdf[j] = X.jdf[0];
                else if (Y.isJDF()) 
                    dt.jdf[j] = Y.jdf[j];
                if (X.isTraced())
                    dt.assigned[j].addAll(X.assigned[0]);
                if (Y.isTraced())
                    dt.assigned[j].addAll(Y.assigned[j]);
            }
            return dt;
        } else if (Y.nEVars == 0) { // only Y is without enumerables
            DenseFactor dt = new DenseFactor(concat(X.evars, X.nvars));
            for (int j = 0; j < X.getSize(); j ++) {
                dt.setValue(j, X.getValue(j) * Y.getValue());
                if (X.isJDF() && Y.isJDF())
                    dt.jdf[j] = JDF.combine(X.jdf[j], Y.jdf[0]);
                else if (X.isJDF())
                    dt.jdf[j] = X.jdf[j];
                else if (Y.isJDF()) 
                    dt.jdf[j] = Y.jdf[0];
                if (X.isTraced())
                    dt.assigned[j].addAll(X.assigned[j]);
                if (Y.isTraced())
                    dt.assigned[j].addAll(Y.assigned[0]);
            }
            return dt;
        }
        // Both tables have enumerable variables, so let's work out how they relate
        int[] ykeyidx = new int[X.nEVars]; // map from X to Y indices [x] = y
        int[] xkeyidx = new int[Y.nEVars]; // map from Y to X indices [y] = x
        int noverlap = 0; // number of overlapping variables
        for (int j = 0; j < Y.nEVars; j ++)
            xkeyidx[j] = -1; // Assume Yj does not exist in X
        for (int i = 0; i < X.nEVars; i ++) {
            ykeyidx[i] = -1; // Assume Xi does not exist in Y
            for (int j = 0; j < Y.nEVars; j ++) {
                if (X.evars[i].equals(Y.evars[j])) { // this variable Xi is at Yj
                    ykeyidx[i] = j;
                    xkeyidx[j] = i;
                    noverlap ++;
                    break;
                }
            }
        }
        if (VERBOSE)
            System.err.println("#overlap = " + noverlap);
        // check if ordered
        boolean ordered = true; // true if *all* overlapping variables are in the same order
        int prev = -1;
        for (int i = 0; i < ykeyidx.length && ordered; i ++) {
            if (ykeyidx[i] != -1) {
                if (ykeyidx[i] <= prev) {
                    ordered = false;
                    break;
                }
                prev = ykeyidx[i];
            }
        }
        if (VERBOSE)
            System.err.println("ordered = " + ordered);
        // Before handling complex scenarios, consider some special, simpler cases
        // 1. two tables with the same variables in the same order
        if (noverlap == Math.min(X.nEVars, Y.nEVars)) { // at least one table is "contained" by the other
            if (X.nEVars == Y.nEVars) {
                Object[] ykey = new Object[Y.nEVars];
                DenseFactor dt = new DenseFactor(X.evars);
                for (int x = 0; x < X.map.length; x ++) {
                    int y = x; // if ordered the inices in the tables are identical
                    if (!ordered) { // re-index since the variables are listed in a different order
                        Object[] xkey = X.getKey(x);
                        for (int i = 0; i < xkey.length; i ++)
                            ykey[ykeyidx[i]] = xkey[i];
                        y = Y.getIndex(ykey);
                    }
                    dt.map[x] = X.getValue(x) * Y.getValue(y);
                    if (X.isJDF() && Y.isJDF())
                        dt.jdf[x] = JDF.combine(X.jdf[x], Y.jdf[y]);
                    else if (X.isJDF())
                        dt.jdf[x] = X.jdf[x];
                    else if (Y.isJDF()) 
                        dt.jdf[x] = Y.jdf[y];
                    if (X.isTraced())
                        dt.assigned[x].addAll(X.assigned[x]);
                    if (Y.isTraced())
                        dt.assigned[x].addAll(Y.assigned[y]);
                }
                if (VERBOSE)
                    System.err.println("DT: Complete overlap, ordered = " + ordered);
                return dt;
            } else if (X.nEVars > Y.nEVars) { // Y is more compact than X
                // some variables in X are not in Y
                Set<EnumVariable> notInY = new HashSet<>();
                for (int i = 0; i < ykeyidx.length; i ++)
                    if (ykeyidx[i] == -1)
                        notInY.add(X.evars[i]);
                Object[] ykey = new Object[Y.nEVars];
                DenseFactor dt = new DenseFactor(X.evars);
                for (int x = 0; x < X.getSize(); x ++) {
                    double xval = X.getValue(x);
                    if (xval == 0)
                        continue; // no point in continuing since product will always be zero for entries with this x-value
                    int y; // there can only be one index in Y (if it is contained in X)
                    if (!ordered) { // re-index considering that variables are listed in a different order
                        Object[] xkey = X.getKey(x);
                        for (int i = 0; i < xkey.length; i ++) {
                            if (ykeyidx[i] != -1)
                                ykey[ykeyidx[i]] = xkey[i];
                        }
                        y = Y.getIndex(ykey);
                    } else { // re-index but ordered, not sure if this is quicker than above... TODO: test
                        y = X.maskIndex(x, notInY); 
                    }
                    double yval = Y.getValue(y);
                    if (yval == 0)
                        continue; // the product will be zero no what
                    int idx = x; 
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF())
                        dt.jdf[idx] = JDF.combine(X.jdf[x], Y.jdf[y]);
                    else if (X.isJDF())
                        dt.jdf[idx] = X.jdf[x];
                    else if (Y.isJDF()) 
                        dt.jdf[idx] = Y.jdf[y];
                    if (X.isTraced())
                        dt.assigned[idx].addAll(X.assigned[x]);
                    if (Y.isTraced())
                        dt.assigned[idx].addAll(Y.assigned[y]);
                }
                if (VERBOSE)
                    System.err.println("DT: Partial overlap (X>Y), ordered = " + ordered);
                return dt;
            } else if (Y.nEVars > X.nEVars) { // X is more compact than Y
                // some variables in Y are not in X
                Set<EnumVariable> notInX = new HashSet<>();
                for (int i = 0; i < xkeyidx.length; i ++)
                    if (xkeyidx[i] == -1)
                        notInX.add(Y.evars[i]);
                Object[] xkey = new Object[X.nEVars];
                DenseFactor dt = new DenseFactor(Y.evars);
                for (int y = 0; y < Y.getSize(); y ++) {
                    double yval = Y.getValue(y);
                    if (yval == 0)
                        continue; // no point in continuing since product will always be zero for entries with this x-value
                    int x; // there can only be one index in X (if it is contained in Y)
                    if (!ordered) { // re-index considering that variables are listed in a different order
                        Object[] ykey = Y.getKey(y);
                        for (int i = 0; i < ykey.length; i ++) {
                            if (xkeyidx[i] != -1)
                                xkey[xkeyidx[i]] = ykey[i];
                        }
                        x = X.getIndex(xkey);
                    } else { // re-index but ordered, not sure if this is quicker than above... TODO: test
                        x = Y.maskIndex(y, notInX); 
                    }
                    double xval = X.getValue(x);
                    if (xval == 0)
                        continue; // the product will be zero no what
                    int idx = y; 
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF())
                        dt.jdf[idx] = JDF.combine(X.jdf[x], Y.jdf[y]);
                    else if (X.isJDF())
                        dt.jdf[idx] = X.jdf[x];
                    else if (Y.isJDF()) 
                        dt.jdf[idx] = Y.jdf[y];
                    if (X.isTraced())
                        dt.assigned[idx].addAll(X.assigned[x]);
                    if (Y.isTraced())
                        dt.assigned[idx].addAll(Y.assigned[y]);
                }
                if (VERBOSE)
                    System.err.println("DT: Partial overlap (X<Y), ordered = " + ordered);
                return dt;
            }
        }

        // Failing the above, we must construct a table which is an amalgamate of the two.
        // Now, we prepare the variable "key" for the new, aggregate table.
        EnumVariable[] rkey = new EnumVariable[X.nEVars + Y.nEVars - noverlap];
        System.arraycopy(X.evars, 0, rkey, 0, X.nEVars);
        int ycnt = 0;
        for (int j = 0; j < Y.nEVars; j ++) {
            if (xkeyidx[j] == -1)
                rkey[X.nEVars + ycnt ++] = Y.evars[j];
        }
        DenseFactor dt = new DenseFactor(rkey);
        
        // 2. two tables have nothing in common
        if (noverlap == 0) {
            for (int x = 0; x < X.getSize(); x ++) {
                double xval = X.getValue(x);
                if (xval == 0)
                    continue; // no point in continuing since product will always be zero for entries with this x-value
                for (int y = 0; y < Y.getSize(); y ++) {
                    double yval = Y.getValue(y);
                    if (yval == 0)
                        continue; // the product will be zero no what
                    int idx = x * Y.getSize() + y; 
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF())
                        dt.jdf[idx] = JDF.combine(X.jdf[x], Y.jdf[y]);
                    else if (X.isJDF())
                        dt.jdf[idx] = X.jdf[x];
                    else if (Y.isJDF()) 
                        dt.jdf[idx] = Y.jdf[y];
                    if (X.isTraced())
                        dt.assigned[idx].addAll(X.assigned[x]);
                    if (Y.isTraced())
                        dt.assigned[idx].addAll(Y.assigned[y]);
                }
            }  
            if (VERBOSE)
                System.err.println("DT: No overlap.");
            return dt;
        }
        
        // 3. General case, if none of those implemented above has worked
        Object[] searchkey = new Object[Y.nEVars]; // this will initialise all elements to null
        Object[] reskey = new Object[X.nEVars + Y.nEVars - noverlap]; // this will initialise all elements to null
        int option = PRODUCT_OPTION;
        long start0 = 0, start1 = 0, total0 = 0, total1 = 0;
        for (int x = 0; x < X.getSize(); x ++) {
            double xval = X.getValue(x);
            if (xval == 0)
               continue; // no point in continuing since product will always be zero for entries with this x-value
            Object[] xkey = X.getKey(x);
            for (int i = 0; i < ykeyidx.length; i ++) {
                int idx = ykeyidx[i];
                if (idx > -1) 
                    searchkey[idx] = xkey[i];
                reskey[i] = xkey[i];
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
                    if (yval == 0)
                        continue; // the product will be zero no what
                    Object[] ykey = Y.getKey(y);
                    int matches = 0;
                    for (int i = 0; i < ykey.length; i ++) {
                        if (xkeyidx[i] == -1)
                            reskey[xkey.length + (matches ++)] = ykey[i];
                    }
                    int idx = dt.getIndex(reskey);
                    dt.setValue(idx, xval * yval);
                    if (X.isJDF() && Y.isJDF())
                        dt.jdf[idx] = JDF.combine(X.jdf[x], Y.jdf[y]);
                    else if (X.isJDF())
                        dt.jdf[idx] = X.jdf[x];
                    else if (Y.isJDF()) 
                        dt.jdf[idx] = Y.jdf[y];
                    if (X.isTraced())
                        dt.assigned[idx].addAll(X.assigned[x]);
                    if (Y.isTraced())
                        dt.assigned[idx].addAll(Y.assigned[y]);
                }
                if (start0 != 0) {
                    total0 = System.nanoTime() - start0;
                    option = -1;
                }
            } else if (option == 1) {
                for (int y = 0; y < Y.getSize(); y ++) {
                    if (Y.isMatch(searchkey, y)) {
                        double yval = Y.getValue(y);
                        if (yval == 0)
                            continue; // the product will be zero no what
                        Object[] ykey = Y.getKey(y);
                        int matches = 0;
                        for (int i = 0; i < ykey.length; i ++) {
                            if (xkeyidx[i] == -1)
                                reskey[xkey.length + (matches ++)] = ykey[i];
                        }
                        int idx = dt.getIndex(reskey);
                        dt.setValue(idx, xval * yval);
                        if (X.isJDF() && Y.isJDF())
                            dt.jdf[idx] = JDF.combine(X.jdf[x], Y.jdf[y]);
                        else if (X.isJDF())
                            dt.jdf[idx] = X.jdf[x];
                        else if (Y.isJDF()) 
                            dt.jdf[idx] = Y.jdf[y];
                        if (X.isTraced())
                            dt.assigned[idx].addAll(X.assigned[x]);
                        if (Y.isTraced())
                            dt.assigned[idx].addAll(Y.assigned[y]);
                    }
                }
                if (start1 != 0) {
                    total1 = System.nanoTime() - start1;
                    option = -1;
                }
            }
            if (start1 != 0) { // done with timing
                if (total0 > total1)
                    option = 1;
                else
                    option = 0;
            }
        }
        if (VERBOSE)
            System.err.println("DT: Generic case. Option = " + option + " (Option 0 took " + total0 + "ns. Option 1 took " + total1 + "ns.)");
        return dt;
    }
    
    /**
     * Construct a new table from an existing, by summing-out specified variable/s.
     * 
     * @param X existing table
     * @param anyvars variables to sum-out
     * @return the resulting "margin" of the table
     */
    public static DenseFactor getMargin(DenseFactor X, Variable ... anyvars) {
        // first check that X has enumerable variables, because that would be required
        if (X.nEVars == 0) // if no enumerables, nothing can be done, so we'll return the same factor
            return X;
        // if X is ok, get rid of non-enumerables
        int cnt_evars = 0;
        for (Variable var : anyvars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                cnt_evars ++;
            } catch (ClassCastException e) {
            }
        }
        EnumVariable[] evars = new EnumVariable[cnt_evars];
        if (cnt_evars == 0) 
            return X;
        else { 
            cnt_evars = 0; 
            for (Variable var : anyvars) {
                try {
                    EnumVariable evar = (EnumVariable) var;
                    evars[cnt_evars ++] = evar;
                } catch (ClassCastException e) {
                }
            }
        }
        // now, resolve the relationships between X's variables and those to be summed-out
        Variable[] yvars = new Variable[X.nEVars + X.nNVars - evars.length];
        int[] xkeyidx = new int[X.nEVars - evars.length]; // map from evars to X indices [evar] = x
        int cnt = 0;
        for (int i = 0; i < X.nEVars; i ++) {
            boolean keep = true;
            for (EnumVariable evar : evars) {
                if (X.evars[i].equals(evar)) {
                    keep = false;
                    break;
                }
            }
            if (keep) {
                xkeyidx[cnt] = i;
                yvars[cnt ++] = X.evars[i];
            } 
        }
        for (int i = 0; i < X.nNVars; i ++)
            yvars[cnt ++] = X.nvars[i];
        if (cnt != X.nEVars + X.nNVars - evars.length)
            throw new DenseFactorRuntimeException("Invalid variable list");
        DenseFactor Y = new DenseFactor(yvars);
        Object[] xkey_search = new Object[X.nEVars];
        for (int y = 0; y < Y.getSize(); y ++) {
            Object[] ykey = Y.getKey(y);
            for (int i = 0; i < xkeyidx.length; i ++)
                xkey_search[xkeyidx[i]] = ykey[i];
            double sum = 0;
            int[] indices = X.getIndices(xkey_search);
            JDF first = null, mixture = null;
            double prev_weight = 0;
            for (int x : indices) {
                double xval = X.map[x];
                sum += xval;
                if (X.isJDF()) {
                    if (first == null) {  // this will be true only once, for the first entry that will be collapsed
                        first = X.jdf[x]; // save the JDF for later
                        prev_weight = xval; // save the value associated with this entry to weight the distributions
                    } else if (mixture == null) { // for the second entry, the mixture is created, both JDF weighted
                        mixture = JDF.mix(mixture, prev_weight, X.jdf[x], xval);
                    } else { // for all subsequent entries, the mixture is mixed unchanged with the new entry's JDF weighted
                        mixture = JDF.mix(mixture, X.jdf[x], xval);
                    }
                }
            }
            Y.map[y] = sum;
            if (Y.isJDF())
                Y.jdf[y] = mixture;
        }
        return Y;
    }
    
    /**
     * Construct a new table from an existing, by maxing-out specified variable/s, and tracing
     * the assignment that provided the maximum value in the resulting table.
     * This code is based on that of getMargin.
     * 
     * @param X existing table
     * @param evars variables to sum-out
     * @return the resulting "margin" of the table
     */
    public static DenseFactor getMaxMargin(DenseFactor X, Variable ... anyvars) {
        // first check that X has enumerable variables, because that would be required
        if (X.nEVars == 0) // if no enumerables, nothing can be done, so we'll return the same factor
            return X;
        // if X is ok, get rid of non-enumerables
        int cnt_evars = 0;
        for (Variable var : anyvars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                cnt_evars ++;
            } catch (ClassCastException e) {
            }
        }
        EnumVariable[] evars = new EnumVariable[cnt_evars];
        if (cnt_evars == 0) 
            return X;
        else { 
            cnt_evars = 0; 
            for (Variable var : anyvars) {
                try {
                    EnumVariable evar = (EnumVariable) var;
                    evars[cnt_evars ++] = evar;
                } catch (ClassCastException e) {
                }
            }
        }
        // now find the relationships between the enumerables
        Variable[] yvars = new Variable[X.nEVars + X.nNVars - evars.length];
        int[] xkeyidx = new int[X.nEVars - evars.length]; // map from Y to X indices [y] = x
        int[] outidx = new int[evars.length]; // index of enumerable variables in X that will be summed out
        int cnt_keep = 0, cnt_dontkeep = 0;
        for (int i = 0; i < X.nEVars; i ++) {
            boolean keep = true;
            for (EnumVariable evar : evars) {
                if (X.evars[i].equals(evar)) {
                    keep = false;
                    break;
                }
            }
            if (keep) {
                xkeyidx[cnt_keep] = i;
                yvars[cnt_keep ++] = X.evars[i];
            } else {
                outidx[cnt_dontkeep ++] = i;
            }
        }
        for (int i = 0; i < X.nNVars; i ++)
            yvars[cnt_keep ++] = X.nvars[i];
        if (cnt_keep != X.nEVars + X.nNVars - evars.length)
            throw new DenseFactorRuntimeException("Invalid variable list");
        DenseFactor Y = new DenseFactor(yvars);
        Y.setTraced(true);
        Object[] xkey_search = new Object[X.nEVars];
        for (int y = 0; y < Y.getSize(); y ++) {
            Object[] ykey = Y.getKey(y);
            for (int i = 0; i < xkeyidx.length; i ++)
                xkey_search[xkeyidx[i]] = ykey[i];
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
            Y.map[y] = max;
            if (Y.isJDF())
                Y.jdf[y] = X.jdf[maxidx];
            // Trace implied assignments
            Object[] xkey = X.getKey(maxidx);
            for (int i = 0; i < evars.length; i ++) {
                Variable.Assignment a = new Variable.Assignment(evars[i], xkey[outidx[i]]);
                Y.addAssign(y, a);
            }
        }
        return Y;
    }
    

    /**
     * Set all entries to 0.
     */
    public void setEmpty() {
        for (int i = 0; i < map.length; i ++)
            map[i] = 0;
    }
    
    /**
     * Get the theoretical number of entries in this table. Note this number is
     * always equal to the actual number of entries, but some may never have been explicitly 
     * set to a value.
     *
     * @return the size (number of entries)
     */
    public int getSize() {
        return map.length;
    }

    /**
     * Get the enumerable variables of the table.
     * @return the variables in original order.
     */
    public EnumVariable[] getEnumVars() {
        return evars;
    }

    /**
     * Get the non-enumerable variables of the table.
     * @return the variables in original order.
     */
    public Variable[] getNonEnumVars() {
        return nvars;
    }

    /**
     * Get the canonical names of the parent variables (names + "." + index)
     * @return names of variables (in order)
     */
    public String[] getLabels() {
        String[] labels = new String[nEVars + nNVars];
        Variable[] all = concat(this.evars, this.nvars);
        for (int i = 0; i < labels.length; i++) 
            labels[i] = all[i].toString();
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
        if (getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key: no variables");
        if (key.length != nEVars)
            throw new DenseFactorRuntimeException("Invalid key: length is " + key.length + " not " + nEVars);
        int sum = 0;
        for (int i = 0; i < nEVars; i++) {
            if (key[i] == null)
                throw new DenseFactorRuntimeException("Null in key");
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
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index");
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
            throw new DenseFactorRuntimeException("Invalid index or key");
        }
        int remain = index;
        for (int i = 0; i < nEVars; i++) {
            if (key[i] != null) {
                int keyindex = evars[i].getIndex(key[i]);
                if (keyindex != remain / step[i]) {
                    return false;
                }
                remain -= keyindex * step[i];
            } else { // key[i] == null
                int missing = remain / step[i];
                remain -= missing * step[i];
            }
        }
        return true;
    }

    /**
     * Identify each "theoretical" index that is linked to the specified key (which may
     * include "wildcards", indicated by null values).
     * In contrast to only identifying populated entries like getIndices does, this function
     * finds all indices that in theory match the key.
     *
     * @param key
     * @return an array with all matching indices
     */
    public int[] getIndices(Object[] key) {
        if (key.length != nEVars)
            throw new DenseFactorRuntimeException("Invalid key for EnumTable: key should be " + nEVars + " but is " + key.length + " values");
        int[] idx; // size a function of domain sizes of null columns
        int startentry = 0; // determined from non-null entries, where counting will start
        int tot = 1;
        for (int i = 0; i < key.length; i ++) {
            if (key[i] == null) {
                tot *= domsize[i];
            } else {
                int keyidx = evars[i].getIndex(key[i]);
                startentry += keyidx * step[i];
            }
        }
        idx = new int[tot];
        if (tot == 1) // no null values
            idx[0] = startentry;
        else
            getIndicesRecursive(key, idx, 0, startentry, 0);
        return idx;
    }

    private synchronized int getIndicesRecursive(Object[] key, int[] idx, int my_idx, int my_tab, int parent) {
        if (key.length <= parent) 
            return -1;
        // parent is real
        if (key[parent] == null) {
            for (int i = 0; i < domsize[parent]; i ++) {
                int start = getIndicesRecursive(key, idx, my_idx, my_tab, parent + 1);
                if (start != -1) {
                    my_idx = start;
                    my_tab += step[parent];
                } else { // start == null meaning that we're at leaf
                    idx[my_idx ++] = my_tab;
                    my_tab += step[parent];
                }
            }
        } else {
            return getIndicesRecursive(key, idx, my_idx, my_tab, parent + 1);
        }
        return my_idx;
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
            if (!maskMe.contains(evars[i])) // if NOT masked-out
                newvale[jn++] = domsize[i];
        }
        jn = newstep.length - 1;
        int prod = 1;
        for (int i = nEVars - 1; i >= 0; i--) {
            if (!maskMe.contains(evars[i])) { // if NOT masked-out
                newstep[jn] = prod;
                prod *= newvale[jn--];
            }
        }
        jn = 0;
        for (int i = 0; i < nEVars; i++) {
            int key = origremain / step[i];
            origremain -= key * step[i];
            if (!maskMe.contains(evars[i])) // if NOT masked-out
                sum += (key * newstep[jn++]);
        }
        return sum;
    }

    /**
     * Get a key from a partial assignment of variables defined for the table.
     * @param evid
     * @return 
     */
    public Object[] getKey(Variable.Assignment[] evid) {
        return DenseFactor.getKey(this.evars, evid);
    }

    /**
     * Print the table.
     */
    public void display() {
        System.out.print("Idx ");
        for (int j = 0; j < this.nEVars; j++) {
            System.out.print(String.format("[%8s]", this.evars[j].getName()));
        }
        System.out.println(" P");
        for (int i = 0; i < this.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = this.getKey(i);
            for (int j = 0; j < key.length; j++) {
                System.out.print(String.format(" %-8s ", key[j].toString()));
            }
            System.out.println(String.format(" %5.3f", this.getValue(i)));
        }
    }

    /**
     * Get a key from a partial assignment of variables defined for the table.
     * @param vararr variables in order
     * @param evid the evidence
     * @return the key that encodes the values for the provided variables
     */
    public static Object[] getKey(EnumVariable[] vararr, Variable.Assignment[] evid) {
        Object[] key = new Object[vararr.length]; // key is initialised to nulls by default
        if (key.length <= evid.length) { // for efficiency, we check what to iterate over
            for (Variable.Assignment e : evid) {
                try {
                    EnumVariable evar = (EnumVariable) e.var;
                    int var_index = Arrays.binarySearch(vararr, evar);
                    if (var_index >= 0)
                        key[var_index] = e.val;
                } catch (ClassCastException exception) {
                    ; // ignore non-enumerable variables
                }
            }
        } else { // evidence is longer than key, so let's iterate over key
            for (int i = 0; i < key.length; i ++) {
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
     * Copy over all non-null values from source to target key.
     * 
     * @param target 
     * @param source
     * @return target
     **/
    public static Object[] overlay(Object[] target, Object[] source) {
        if (target.length != source.length)
            throw new DenseFactorRuntimeException("Invalid operation since keys are of difference lengths (" + target.length + " vs " + source.length + ")");
        for (int i = 0; i < target.length; i ++)
            if (source[i] != null)
                target[i] = source[i];
        return target;
    }
    
    public static void main(String[] args) {
        java.util.Random random = new java.util.Random();
        EnumVariable x1 = Predef.Boolean("X1");
        EnumVariable x2 = Predef.AminoAcid("AA1");
        EnumVariable y1 = Predef.Number(2);
        EnumVariable y2 = Predef.Nominal(new String[] {"a", "b", "c"});
        EnumVariable z1 = Predef.NucleicAcid("NA1");
        EnumVariable z2 = Predef.AminoAcid("AA2");
        DenseFactor dt0 = new DenseFactor(x1,y1,y2);
        for (int key_index = 0; key_index < dt0.getSize(); key_index ++)
            dt0.setValue(key_index, random.nextDouble());
        dt0.display();
    
        DenseFactor mt0 = DenseFactor.getMargin(dt0, y1, x1);
        mt0.display();
        
        DenseFactor dt1 = new DenseFactor(y1,x1);
        for (int key_index = 0; key_index < dt1.getSize(); key_index ++)
            dt1.setValue(key_index, random.nextDouble());
        dt1.display();
        
        DenseFactor mt1 = DenseFactor.getMargin(dt1, x1);
        mt1.display();
        
        DenseFactor dt2 = new DenseFactor(x2,y2,z1,x1);
        for (int key_index = 0; key_index < dt2.getSize(); key_index ++)
            dt2.setValue(key_index, random.nextDouble());
        dt2.display();
        DenseFactor mt2 = DenseFactor.getMargin(dt2, x2, y2);
        mt2.display();
        
        DenseFactor dt3 = new DenseFactor(z1,x2,z2);
        for (int key_index = 0; key_index < dt3.getSize(); key_index ++)
            dt3.setValue(key_index, random.nextDouble());
//        dt3.display();
        DenseFactor dt4 = DenseFactor.getProduct(dt2, dt1);
//        dt4.display();
        DenseFactor dt5 = DenseFactor.getProduct(dt3, dt4);
//        dt5.display();
        DenseFactor dt6 = DenseFactor.getProduct(dt3, dt4);
//        dt6.display();
    }
}

class DenseFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public DenseFactorRuntimeException(String string) {
        message = string;
    }
}

