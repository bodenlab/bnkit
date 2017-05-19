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
import dat.Variable;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Table for storing and retrieving doubles based on Enumerable keys.
 * 
 * The intended use is for "factors", representing the factorisation of (conditional) probabilities.
 * They form the backbone data structure in algorithms such as variable elimination and 
 * message-passing in belief propagation. The main operations that are supported include 
 * product, and variable marginalisation (summing-out or maxing-out).
 * These operations need to be efficient. 
 * 
 * SparseFactor is intended for very large tables that are sparse, i.e have few values that are non-zero.
 *
 * @author mikael
 */

public class SparseFactor extends AbstractFactor {

    private final EnumTable<Double> fac;      // the factors for each permutation of the enumerable variables
    private Double fac_atomic = null;         // the factor when there are no enumerable variables
    private EnumTable<JDF> jdf = null;     // the densities for each permutation of the enumerables
    private JDF jdf_atomic = null;            // the JDF when there are no enumerable variables
    private EnumTable<Set<Variable.Assignment>> ass = null; // the traced variable assignments
    private Set<Variable.Assignment> ass_atomic = null;     // the traced variable assignment when there are no enumerables

    @Override
    public Iterator<Integer> iterator() {
        return fac.iterator();
    }

    /**
     * Construct a new table without any variables.
     * This type of factor is used when all variables are summed out. 
     * They can appear as part of products to "scale" its opposite.
     */
    public SparseFactor() {
        super();
        this.fac = new EnumTable<>();
    }

    
    /**
     * Construct a new table with the specified variables.
     * Enumerable variables will form keys to index the entries, 
     * non-enumerables will be associated in the form of their densities with
     * each entry. 
     * 
     * @param useVariables variables, either enumerable or non-enumerable, 
     * potentially unsorted and redundant
     */
    public SparseFactor(Variable... useVariables) {
        super(useVariables);
        // most of the time, there would be one or more enumerable variables so allocate a map
        // for all factor values
        if (this.nEVars > 0)
            this.fac = new EnumTable<>(super.getEnumVars());
        else { // sometimes there are no enumerable variables
            this.fac = null;
            this.fac_atomic = LOG0;
        }
        // if one or more non-enumerable variables, we allocate a matching map for JDFs
        if (this.nNVars > 0 && this.nEVars > 0) 
            this.jdf = new EnumTable<>();
        else if (this.nNVars > 0) 
            this.jdf_atomic = new JDF(this.nvars);
    }


    @Override
    public double getLogValue() {
        if (this.getSize() == 1) {
            Double value = fac_atomic;
            if (value == null)
                return LOG0;
            return fac_atomic;
        }
        throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }

    @Override
    public JDF getJDF() {
        if (this.getSize() == 1)
            return this.jdf_atomic;
        throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }

    @Override
    public double getLogValue(int index) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid index");
        Double value = fac.getValue(index);
        if (value == null)
            return LOG0;
        return value;
    }

    @Override
    public JDF getJDF(int index) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index");
        JDF y = jdf.getValue(index);
        return y;
    }

    @Override
    public int setLogValue(double value) {
        if (this.getSize() == 1) {
            fac_atomic = value;
            return 0;
        }
        throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }

    @Override
    public int setLogValue(int key_index, double value) {
        if (key_index >= getSize() || key_index < 0 || this.getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid index");
        if (value == LOG0)
            fac.removeValue(key_index);
        else 
            fac.setValue(key_index, value);
        return key_index;
    }

    @Override
    public int setJDF(JDF value) {
        if (this.getSize() == 1) {
            jdf_atomic = value;
            return 0;
        }
        throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }

    @Override
    public int setJDF(int key_index, JDF value) {
        if (key_index >= getSize() || key_index < 0 || this.getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid index");
        return jdf.setValue(key_index, value);
    }

    @Override
    public void setTraced(boolean status) {
        if (status && getSize() != 1) {
            this.ass = new EnumTable<>();
        } else if (status) {
            this.ass_atomic = new HashSet<>();
        } else {
            this.ass = null;
            this.ass_atomic = null;
        }
    }

    @Override
    public boolean isTraced() {
        return (this.ass != null || this.ass_atomic != null);
    }

    @Override
    public int addAssign(Collection<Variable.Assignment> assign) {
        if (!isTraced()) 
            throw new SparseFactorRuntimeException("Tracing is not enabled");
        if (this.getSize() != 1)
            throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
        ass_atomic.addAll(assign);
        return 0;
    }

    @Override
    public int addAssign(int key_index, Collection<Variable.Assignment> assign) {
        if (!isTraced()) 
            throw new SparseFactorRuntimeException("Tracing is not enabled");
        if (this.getSize() != 1)
            throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
        Set<Variable.Assignment> a = ass.getValue(key_index);
        if (a == null)
            a = new HashSet<>();
        a.addAll(assign);
        return key_index;
    }

    @Override
    public Set<Variable.Assignment> getAssign() {
        if (this.getSize() == 1)
            return ass_atomic;
        throw new SparseFactorRuntimeException("Table has variables that must be used to index access");
    }

    @Override
    public Set<Variable.Assignment> getAssign(int key_index) {
        if (key_index >= getSize() || key_index < 0 || getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid key index: outside map");
        return ass.getValue(key_index);
    }

    @Override
    public int addAssign(Variable.Assignment assign) {
        if (!isTraced()) 
            throw new SparseFactorRuntimeException("Invalid key index: outside map");
        if (getSize() == 1)
            throw new SparseFactorRuntimeException("Table has variables that must be used to index access");
        ass_atomic.add(assign);
        return 0;
    }

    @Override
    public int addAssign(int key_index, Variable.Assignment assign) {
        if (!isTraced()) 
            throw new SparseFactorRuntimeException("Invalid key index: outside map");
        if (key_index >= getSize() || key_index < 0 || getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid key index: outside map");
        Set<Variable.Assignment> a = ass.getValue(key_index);
        if (a == null)
            a = new HashSet<>();
        a.add(assign);
        return ass.setValue(key_index, a);
    }

    @Override
    public Distrib getDistrib(Variable nvar) {
        if (this.getSize() == 1)
            return jdf_atomic.getDistrib(nvar);
        throw new SparseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }

    @Override
    public Distrib getDistrib(int index, Variable nvar) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid index");
        JDF myjdf = jdf.getValue(index);
        if (myjdf != null)
            return myjdf.getDistrib(nvar);
        return null;
    }

    @Override
    public int setDistrib(int key_index, Variable nonenum, Distrib d) {
        if (key_index >= getSize() || key_index < 0 || getSize() == 1)
            throw new SparseFactorRuntimeException("Invalid key index: outside map");
        JDF myjdf = jdf.getValue(key_index);
        if (myjdf == null) {
            myjdf = new JDF(nvars);
            jdf.setValue(key_index, myjdf);
        }
        myjdf.setDistrib(d, nonenum);
        return key_index; // 
    }

    @Override
    public int setDistrib(Variable nvar, Distrib d) {
        if (this.getSize() != 1)
            throw new SparseFactorRuntimeException("Table has variables that must be used to index access");
        jdf_atomic.setDistrib(d, nvar);
        return 0;
    }

    @Override
    public void setEmpty() {
        if (fac != null)
            fac.setEmpty();
        if (jdf != null)
            jdf.setEmpty();
        if (ass != null)
            ass.setEmpty();
    }

    /**
     * Get indices for all non-null entries.
     * @return the indices
     */
    public int[] getIndices() {
        return fac.getIndices();
    }
    

    /**
     * Identify each index that is linked to the specified key (which may
     * include "wildcards", indicated by null values).
     *
     * @param key
     * @return an array with all matching indices
     */
    @Override
    public int[] getIndices(Object[] key) {
        return fac.getIndices(key);
    }

    @Override
    public int getOccupied() {
        return fac.getSize();
    }
    
}

class SparseFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public SparseFactorRuntimeException(String string) {
        message = string;
    }
}
