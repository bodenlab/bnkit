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
import bn.JDF;
import dat.EnumTable;
import dat.EnumVariable;
import dat.Variable;

import java.util.*;

/**
 * Table for storing and retrieving doubles based on Enumerable keys.
 * 
 * The intended use is for "factors", representing the factorisation of (conditional) probabilities.
 * They form the backbone data structure in algorithms such as variable elimination and 
 * message-passing in belief propagation. These operations need to be efficient. 
 * 
 * AbstractFactor is an abstract class that defines methods for specific implementations but 
 * implements basic operations on the data structure that are used externally on factors, incl.
 * finding indices, or keys. See {@see bn.factor.Factorize} for numerical operations on factors.
 * 
 * The class is sensitive to the order of variables, and exposes the user to some
 * details so caution must be exercised. Of specific importance is that the table stores the
 * natural logarithm of each factor, so that operations can be performed entirely in log space,
 * to avoid numerical issues, e.g. underflow.
 * 
 * TODO: Consider improving efficiency further to exploit the fact that now all variables are sorted
 * in the constructor. (Currently, some code does not assume order.)
 * TODO: Choose type of implementation when new factors are constructed as a result of operations,
 * informed of their likely requirements.
 *
 * @author mikael
 */
public abstract class AbstractFactor implements Iterable<Integer> {
    @Override
    abstract public Iterator<Integer> iterator();

    protected final static double LOG0 = Double.NEGATIVE_INFINITY;
    protected static boolean isLOG0(double x) { return Double.isInfinite(x); }

    protected final int nEVars; // number of enumerable variables
    protected final int nNVars; // number of non-enumerable variables

    protected final EnumVariable[] evars; // enumerable variables, which form the "keys" to entries in the map
    protected final Variable[] nvars; // non-enumerable variables, which densities are attached to the entries via JDF
    protected final int[] period;
    protected final int[] step;
    protected final int[] domsize; // size of domain
    
    public boolean evidenced = false;

    protected long cacheID = 0L; // 0L means not possible to cache

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
        Variable[] uniqueVars = Factorize.getNonredundantSorted(useVariables);
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

    public String toString() {
        StringBuilder sb = new StringBuilder("F(");
        for (int i = 0; i < nEVars; i ++)
            sb.append(evars[i].getName() + ((i == nEVars - 1)?(";"):(",")));
        for (int i = 0; i < nNVars; i ++)
            sb.append(nvars[i].getName() + ((i == nNVars - 1)?(""):(",")));
        sb.append(")");
        if (!hasEnumVars())
            sb.append("=" + this.getValue());
        return sb.toString();
    }

    /**
     * Set the ID used to cache factors, intended to conserve memory and time when factors are essentially the same
     * (and do not require separate storage and calculation when used as part of repeated operations between them)
     * @param id
     */
    protected void setCacheID(long id) {
        cacheID = id;
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
        return Math.exp(getLogSum());
    }
    
    /**
     * Calculate the log sum of factors: log(a + b + ... )
     * @return the sum
     */
    public double getLogSum() {
        double sum = 1;
        if (!hasEnumVars())
            return getLogValue();
        for (int i = 0; i < getSize(); i ++) {
            double logp = getLogValue(i);
            if (Double.isNaN(logp)) // skip entries that are un-instantiated
                continue;
            if (i == 0)
                sum = logp;
            else
                sum = Factorize.logSumOfLogs(sum, logp);
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
     * Find index of the specified variable.
     * @param var variable
     * @return index in factor (positive if enumerable, negative if non-enumerable)
     */
    public int getVariableIndex(Variable var) {
        for (int i = 0; i < nEVars; i ++)
            if (evars[i].equals(var))
                return i;
        for (int i = 0; i < nNVars; i ++)
            if (nvars[i].equals(var))
                return -i;
        throw new RuntimeException("Invalid variable " + var + " for factor");
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
        Variable[] all = Factorize.getConcat(this.evars, this.nvars);
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
     * Takes an entry index of the current table and a revised ordering of variables,
     * to determine the index in the re-ordered table.
     * Note that the number of variables is expected to be the same.
     *
     * @param origindex the index in this table
     * @param newvars new ordering, same variables
     * @return index in other, re-ordered table
     */
    public int reIndex(int origindex, EnumVariable[] newvars) {
        if (newvars.length != this.nEVars)
            throw new AbstractFactorRuntimeException("Invalid variable list");
        EnumVariable[] x = this.evars;
        EnumVariable[] y = newvars;
        int[] xcross2y = new int[x.length];
        int[] ycross2x = new int[y.length];
        Factorize.getCrossref(x, xcross2y, y, ycross2x);
        Object[] xkey = this.getKey(origindex);
        Object[] ykey = new Object[y.length];
        for (int i = 0; i < xkey.length; i ++)
            ykey[xcross2y[i]] = xkey[i];
        return EnumTable.getIndex(ykey, y);
    }
    
    /**
     * Takes an entry index of the current table and a revised ordering of variables,
     * to determine the index in the re-ordered table.
     * Note that the number of variables is expected to be the same.
     *
     * @param origindex the index in this table
     * @param newvars new ordering, same variables
     * @return index in other, re-ordered table
     */
    public int reIndex(int origindex, List<EnumVariable> newvars) {
        EnumVariable[] arr = new EnumVariable[newvars.size()];
        newvars.toArray(arr);
        return reIndex(origindex, arr);
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
    public double getValue() {
        double y = getLogValue();
        if (Double.isInfinite(y))
            return 0;
        return Math.exp(y);
    }
    
    /**
     * Retrieve the log value of the factor, if table is without enumerable variables.
     * @return the only value of the factor
     */
    public abstract double getLogValue();
    
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
    public double getValue(int index) {
        double y = getLogValue(index);
        if (Double.isInfinite(y))
            return 0;
        return Math.exp(y);
    }

    /**
     * Retrieve the log value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    public abstract double getLogValue(int index);

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
    public int setValue(double value) {
        if (value == 0)
            return setLogValue(LOG0);
        return setLogValue(Math.log(value));
    }
    
    /**
     * Set the only value associated with a table without enumerable variables.
     * @param value
     * @return 
     */
    public abstract int setLogValue(double value);
    
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
        if (value == 0)
            return setLogValue(key_index, LOG0);
        return setLogValue(key_index, Math.log(value));
    }

    /**
     * Associate the specified key-index with the given log value. Note that using
     * getValue and setValue with index is quicker than with key, if more than
     * one operation is done.
     *
     * @param key_index
     * @param value
     * @return the index at which the value was stored
     */
    public abstract int setLogValue(int key_index, double value);


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
     * Associate the specified key with the given log value.
     * @param key
     * @param value
     * @return the index at which the value was stored
     */
    public int setLogValue(Object[] key, double value) {
        if (getSize() == 1)
            throw new AbstractFactorRuntimeException("Invalid key: no variables");
        int index = getIndex(key);
        return setLogValue(index, value);
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
    public abstract int addAssign(Collection<Variable.Assignment> assign);
    /**
     * Associate the entry with a collection of assignments.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    public abstract int addAssign(int key_index, Collection<Variable.Assignment> assign);
    
    /**
     * Associate the entry with a collection of assignments.
     * @param key key for entry
     * @param assign assignment
     * @return index of entry
     */
    public int setAssign(Object[] key, Collection<Variable.Assignment> assign) {
        return addAssign(getIndex(key), assign);
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
     * Tag the sole entry with an assignment.
     * @param assign assignment
     * @return index of entry
     */
    public abstract int addAssign(Variable.Assignment assign);

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
        return AbstractFactor.this.addAssign(getIndex(key), assign);
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
     * Find out how many entries that are occupied.
     * @return  number of entries that are instantiated
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
        
}

class AbstractFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public AbstractFactorRuntimeException(String string) {
        message = string;
    }
}
