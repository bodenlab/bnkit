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
import bn.Predef;
import bn.prob.GaussianDistrib;
import dat.EnumVariable;
import dat.Variable;

import java.util.*;

/**
 * Table for storing and retrieving doubles based on Enumerable keys.
 * 
 * The intended use is for "factors", representing the factorisation of (conditional) probabilities.
 * They form the backbone data structure in algorithms such as variable elimination and 
 * message-passing in belief propagation. The main operations that are supported include 
 * product, and variable marginalisation (summing-out or maxing-out).
 * These operations need to be efficient. DenseFactor is efficient when data are dense, 
 * i.e. when almost all possible key combinations are used. However, storage requirements 
 * grow exponentially with number of variables.
 * 
 * The class is sensitive to the order of variables, and exposes the user to 
 * details so caution must be exercised. Of specific importance is that the table stores the
 * natural logarithm of each factor, so that operations can be performed entirely in log space,
 * to avoid numerical issues, e.g. underflow.
 * 
 * TODO: Consider improving efficiency further to exploit the fact that now all variables are sorted in the constructor. (Currently, some code does not assume order.)
 *
 * @author mikael
 */
public class DenseFactor extends AbstractFactor {

    protected final double[] map; // the factors
    protected Set<Variable.Assignment>[] assigned = null; // the latent assignments associated with each factor, disabled by default
    protected JDF[] jdf = null; // the Joint Density Function associated with each factor, disabled by default
    protected int occupied = 0; // number of entries that are occupied
    
    /**
     * Construct a new table without any variables.
     * This type of factor is used when all variables are summed out. 
     * They can appear as part of products to "scale" its opposite.
     */
    public DenseFactor() {
        super();
        this.map = new double[1];
        // this.map[0] = 1.0; // log(0) = 1
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
    public DenseFactor(Variable... useVariables) {
        super(useVariables);
        // most of the time, there would be one or more enumerable variables so allocate a map
        // for all factor values
        if (this.nEVars > 0) {
            this.map = new double[getSize()];
            Arrays.fill(map, LOG0);
        } else { // sometimes there are no enumerable variables, so need only one entry
            this.map = new double[1];
            // this.map[0] = 1.0; // log(0) = 1
        }
        // if one or more non-enumerable variables, we allocate a matching map for JDFs
        if (this.nNVars > 0) {
            this.jdf = new JDF[this.map.length];
            for (int i = 0; i < this.jdf.length; i ++) 
                this.jdf[i] = new JDF(nvars);
        }
    }

    
    /**
     * Retrieve the log value of the factor, if table is without enumerable variables.
     * @return the only value of the factor
     */
    @Override
    public double getLogValue() {
        if (this.getSize() == 1)
            return this.map[0];
        throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }
    
    /**
     * Retrieve the JDF of the factor without enumerable variables.
     * @return probability distribution for variable
     */
    @Override
    public JDF getJDF() {
        if (this.getSize() == 1)
            return this.jdf[0];
        throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }
    
    /**
     * Retrieve the value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    @Override
    public double getLogValue(int index) {
        if (index >= getSize() || index < 0 || getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index: " + index + " in factor with " + getSize() + " entries");
        double value = map[index];
        return value;
    }

    /**
     * Retrieve the JDF
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @return probability distribution for variable
     */
    @Override
    public JDF getJDF(int index) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index");
        return this.jdf[index];
    }
    
    /**
     * Set the (only) log value associated with a table without enumerable variables.
     * @param value
     * @return 
     */
    @Override
    public int setLogValue(double value) {
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
        if (Double.isNaN(value) || value == LOG0) { // trying to set cell to NaN (should not happen) or to LOG0 (does happen)
            if (map[0] != LOG0) // there was a value previously
                occupied --;    // so table is now smaller
            map[0] = LOG0;      // revert to LOG0 even if NaN
            return 0;
        }
        if (value != LOG0 && map[0] == LOG0)
            occupied ++;
        map[0] = value;
        return 0;
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
    @Override
    public int setLogValue(int key_index, double value) {
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        if (Double.isNaN(value) || value == LOG0) { // trying to set cell to NaN (should not happen) or to LOG0 (does happen)
            if (map[key_index] != LOG0) // there was a value previously
                occupied --;    // so table is now smaller
            map[key_index] = LOG0;      // revert to LOG0 even if NaN
            return key_index;
        }
        if (value != LOG0 && map[key_index] == LOG0)
            occupied ++;
        map[key_index] = value;
        return key_index;
    }


    /**
     * Find out how many entries that occupied.
     * @return how many non-zero entries the factor table holds
     */
    @Override
    public int getOccupied() {
        return occupied;
    }
    
    /**
     * Set the JDF for a table without enumerable variables.
     * @param value JDF
     * @return the index (always 0)
     */
    @Override
    public int setJDF(JDF value) {
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
        jdf[0] = value;
        return 0;
    }
    
    /**
     * Associate the specified key-index with the given JDF. 
     *
     * @param key_index index for entry
     * @param value
     * @return the index at which the value was stored
     */
    @Override
    public int setJDF(int key_index, JDF value) {
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        jdf[key_index] = value;
        return key_index;
    }


    /**
     * Activate tracing of implied assignments.
     * @param status true if activated, false otherwise
     */
    @Override
    public void setTraced(boolean status) {
        if (status) {
            this.assigned = new Set[this.getSize()];
            for (int i = 0; i < this.assigned.length; i ++) 
                this.assigned[i] = new HashSet();
        } else {
            this.assigned = null;
        }
    }

    /**
     * Find out if tracing of implied assignments is active.
     * @return true if activated, false otherwise
     */
    @Override
    public boolean isTraced() {
        return this.assigned != null;
    }

     /**
     * Retrieve a collection of assignments from the specified entry.
     * @param key_index the index of the entry
     * @return set of assignments
     */
    @Override
    public Set<Variable.Assignment> getAssign(int key_index) {
        if (key_index >= assigned.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        return assigned[key_index];
    }

    /**
     * Retrieve a collection of assignments from the single entry.
     * @return set of assignments
     */
    @Override
    public Set<Variable.Assignment> getAssign() {
        if (this.getSize() == 1)
            return assigned[0];
        throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
    }

   /**
     * Associate the only entry with a collection of assignments.
     * @param assign assignment
     * @return index of entry
     */
    @Override
    public int addAssign(Collection<Variable.Assignment> assign) {
        if (!isTraced()) 
            throw new DenseFactorRuntimeException("Tracing is not enabled");
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
        if (assigned[0] == null)
            assigned[0] = new HashSet<>();
        assigned[0].addAll(assign);
        return 0;
    }

    /**
     * Associate the entry with a collection of assignments.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    @Override
    public int addAssign(int key_index, Collection<Variable.Assignment> assign) {
        if (!isTraced()) 
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        if (assigned[key_index] == null)
            assigned[key_index] = new HashSet<>();
        assigned[key_index].addAll(assign);
        return key_index;
    }
    
    /**
     * Tag the sole entry with an assignment.
     * @param assign assignment
     * @return index of entry
     */
    @Override
    public int addAssign(Variable.Assignment assign) {
        if (!isTraced()) 
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
        if (assigned[0] == null)
            assigned[0] = new HashSet<>();
        assigned[0].add(assign);
        return 0;
    }
    
    /**
     * Tag the entry with an assignment.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    @Override
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
    
    /**
     * Retrieve the distribution of the specified non-enumerable variable.
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    @Override
    public Distrib getDistrib(Variable nvar) {
        if (this.getSize() == 1)
            return this.jdf[0].getDistrib(nvar);
        throw new DenseFactorRuntimeException("This table must be accessed with a enumerable variable key");
    }
    
    /**
     * Retrieve the distribution of the specified non-enumerable variable, 
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    @Override
    public Distrib getDistrib(int index, Variable nvar) {
        if (index >= getSize() || index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid index");
        return this.jdf[index].getDistrib(nvar);
    }
    
   /**
     * Set the conditional distribution of a non-enumerable variable. 
     * Use with care as the method may not do all integrity checks.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    @Override
    public int setDistrib(int key_index, Variable nonenum, Distrib d) {
        if (key_index >= map.length || key_index < 0 || this.getSize() == 1)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        this.jdf[key_index].setDistrib(d, nonenum);
        return key_index; // 
    }

    /**
     * Set the distribution of a non-enumerable variable for a factor without enumerable variables. 
     * @param nvar the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry, always 0
     */
    @Override
    public int setDistrib(Variable nvar, Distrib d) {
        if (this.getSize() != 1)
            throw new DenseFactorRuntimeException("Table has variables that must be used to index access");
        this.jdf[0].setDistrib(d, nvar);
        return 0;
    }

    /**
     * Set all entries to log(0) to indicate that they are empty (represent a probability 0).
     */
    @Override
    public void setEmpty() {
        Arrays.fill(map, LOG0);
        if (isJDF()) {
            for (int i = 0; i < jdf.length; i ++) 
                this.jdf[i] = new JDF(nvars);
        }
        occupied = 0;
    }
    

    /**
     * Identify each "theoretical" index that is linked to the specified key (which may
     * include "wildcards", indicated by null values).
     * This function finds all indices that in theory match the key.
     *
     * @param key
     * @return an array with all matching indices
     */
    @Override
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


    @Override
    public Iterator<Integer> iterator() {
        return new FactorIterator();
    }

    private class FactorIterator implements Iterator<Integer> {
        private int index = 0;
        @Override
        public boolean hasNext() {
            int size = getSize();
            if (index < size) {
                do {
                    if (map[index] != LOG0)
                        return true;
                    index ++;
                } while (index < size);
            }
            return false;
        }

        @Override
        public Integer next() {
            if (hasNext())
                return index ++;
            throw new NoSuchElementException();
        }
    }

    public static void main(String[] args) {
        long seed = 1;
        Random random = new Random(seed);

        Factorize.VERBOSE = false;
        EnumVariable x1 = Predef.Boolean("X1");
        EnumVariable x2 = Predef.AminoAcid("AA1");
        EnumVariable y1 = Predef.Number(2);
        EnumVariable y2 = Predef.Nominal(new String[] {"a", "b", "c"});
        EnumVariable z1 = Predef.NucleicAcid("NA1");
        EnumVariable z2 = Predef.AminoAcid("AA2");
        Variable r1 = Predef.Real("R1");
        Variable r2 = Predef.Real("R2");
        AbstractFactor dt0 = new DenseFactor(x1,y1,r1,y2);
        for (int key_index = 0; key_index < dt0.getSize(); key_index ++) {
            if (random.nextInt(100) > 20)
                dt0.setValue(key_index, random.nextDouble());
            
            dt0.setDistrib(key_index, r1, new GaussianDistrib(random.nextGaussian(), random.nextDouble()));
        }
        dt0.display();
    
        AbstractFactor mt0 = Factorize.getMargin(dt0, y1, x1);
        mt0.display();
        
        AbstractFactor dt1 = new DenseFactor(y1,x1, r2);
        for (int key_index = 0; key_index < dt1.getSize(); key_index ++) {
            dt1.setValue(key_index, random.nextDouble());
            dt1.setDistrib(key_index, r2, new GaussianDistrib(random.nextGaussian(), random.nextDouble()));
        }
        dt1.display();
        
        AbstractFactor mt1 = Factorize.getMargin(dt1, x1);
        mt1.display();
        AbstractFactor xt1 = Factorize.getMaxMargin(dt1, x1);
        xt1.display();
        
        
        Factorize.getProduct(dt0, dt1);
        
        AbstractFactor dt2 = new DenseFactor(x2,y2,z1,x1);
        for (int key_index = 0; key_index < dt2.getSize(); key_index ++)
            if (random.nextInt(100) > 20)
                dt2.setValue(key_index, random.nextDouble());
        AbstractFactor mt2 = Factorize.getMargin(dt2, x2, y2);
        mt2.display();
        
        Factorize.getProduct(mt2, mt0);

        AbstractFactor dt3 = new DenseFactor(z1,x2,z2);
        for (int key_index = 0; key_index < dt3.getSize(); key_index ++)
            dt3.setValue(key_index, random.nextDouble());
        AbstractFactor dt4 = Factorize.getProduct(dt2, dt1);
        AbstractFactor dt5 = Factorize.getProduct(dt4, mt0);
        AbstractFactor dt6 = Factorize.getProduct(dt3, dt5);

        long startTime, endTime;

        startTime = System.nanoTime();
        Factorize.getProduct(dt0, dt1);
        endTime = System.nanoTime();
        System.out.println((endTime - startTime)/100000.0 + "ms");
        System.out.println(Factorize.getOverlap(dt0, dt1) + " : " + Factorize.getComplexity(dt0, dt1, true) + "\t" + Factorize.getComplexity(dt0, dt1, false) + " : " + dt0.getSize() + " v " + dt1.getSize());

        
        startTime = System.nanoTime();
        Factorize.getProduct(dt1, dt0);
        endTime = System.nanoTime();
        System.out.println((endTime - startTime)/100000.0 + "ms");
        System.out.println(Factorize.getOverlap(dt1, dt0) + " : " + Factorize.getComplexity(dt1, dt0, true) + "\t" + Factorize.getComplexity(dt1, dt0, false) + " : " + dt1.getSize() + " v " + dt0.getSize());

        startTime = System.nanoTime();
        Factorize.getProduct(mt1, dt2);
        endTime = System.nanoTime();
        System.out.println((endTime - startTime)/100000.0 + "ms");
        System.out.println(Factorize.getOverlap(mt1, dt2) + " : " + Factorize.getComplexity(mt1, dt2, true) + "\t" + Factorize.getComplexity(mt1, dt2, false) + " : " + mt1.getSize() + " v " + dt2.getSize());

        startTime = System.nanoTime();
        Factorize.getProduct(dt3, dt4);
        endTime = System.nanoTime();
        System.out.println((endTime - startTime)/100000.0 + "ms");
        System.out.println(Factorize.getOverlap(dt3, dt4) + " : " + Factorize.getComplexity(dt3, dt4, true) + "\t" + Factorize.getComplexity(dt3, dt4, false) + " : " + dt3.getSize() + " v " + dt4.getSize());

        startTime = System.nanoTime();
        Factorize.getProduct(dt5, dt6);
        endTime = System.nanoTime();
        System.out.println((endTime - startTime)/100000.0 + "ms");
        System.out.println(Factorize.getOverlap(dt5, dt6) + " : " + Factorize.getComplexity(dt5, dt6, true) + "\t" + Factorize.getComplexity(dt5, dt6, false) + " : " + dt5.getSize() + " v " + dt6.getSize());

        startTime = System.nanoTime();
        Factorize.getProduct(dt0, dt0);
        endTime = System.nanoTime();
        System.out.println((endTime - startTime)/100000.0 + "ms");
        System.out.println(Factorize.getOverlap(dt0, dt0) + " : " + Factorize.getComplexity(dt0, dt0, true) + "\t" + Factorize.getComplexity(dt0, dt0, false) + " : " + dt0.getSize() + " v " + dt0.getSize());
    }
}

class DenseFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public DenseFactorRuntimeException(String string) {
        message = string;
    }
}

