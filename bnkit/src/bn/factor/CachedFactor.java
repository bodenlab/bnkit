package bn.factor;

import bn.Distrib;
import bn.JDF;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Table for storing and retrieving doubles based ONLY on Enumerable keys, and where identities of variables are unknown.
 * The intended use is for "factors", representing the factorisation of (conditional) probabilities.
 * CachedFactors use the fact that the same factor/s can be used for many variables, which is especially useful for
 * phylogenetic analyses where similar relationships exists across sequence positions.
 *
 * They form the backbone data structure in algorithms such as variable elimination and
 * message-passing in belief propagation. The main operations that are supported include
 * product, and variable marginalisation (summing-out or maxing-out).
 * These operations need to be efficient.
 *
 * The table stores the natural logarithm of each factor, so that operations can be performed entirely in log space,
 * to avoid numerical issues, e.g. underflow.
 *
 * @author mikael
 */

public class CachedFactor extends AbstractFactor {

    public final static ConcurrentHashMap<Integer, CachedTable> cache = new ConcurrentHashMap();

    public CachedFactor() {
        super();
    }

    public CachedFactor(EnumVariable... useVariables) {
        super(useVariables);
        // factor values need space, but we should check if this map is available in cache first
        double[] map = new double[getSize()];
        Arrays.fill(map, LOG0);
    }

    /**
     * Retrieve the log value of the factor, if table is without enumerable variables.
     * @return the only value of the factor
     */
    public double getLogValue() {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Retrieve the JDF of the factor without enumerable variables.
     * @return probability distribution for variable
     */
    public JDF getJDF() {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Retrieve the log value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    public double getLogValue(int index) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Retrieve the JDF
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @return probability distribution for variable
     */
    public JDF getJDF(int index) {
        return null;
    }


    // Setters of cells in factor ---------------------------------------------------------

    /**
     * Set the only value associated with a table without enumerable variables.
     * @param value
     * @return
     */
    public int setLogValue(double value) {
        throw new CachedFactorRuntimeException("Not implemented");
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
    public int setLogValue(int key_index, double value) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Set the JDF for a table without enumerable variables.
     * @param value JDF
     * @return the index (always 0)
     */
    public int setJDF(JDF value) {
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
    }

    /**
     * Associate the specified key-index with the given JDF.
     *
     * @param key_index index for entry
     * @param value
     * @return the index at which the value was stored
     */
    public int setJDF(int key_index, JDF value) {
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
    }

    // Other setters...

    /**
     * Activate tracing of implied assignments.
     * @param status true if activated, false otherwise
     */
    public void setTraced(boolean status) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Find out if tracing of implied assignments is active.
     * @return true if activated, false otherwise
     */
    public boolean isTraced() {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Associate the only entry with a collection of assignments.
     * @param assign assignment
     * @return index of entry
     */
    public int addAssign(Collection<Variable.Assignment> assign) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Associate the entry with a collection of assignments.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    public int addAssign(int key_index, Collection<Variable.Assignment> assign) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Retrieve a collection of assignments from the single entry.
     * @return set of assignments
     */
    public Set<Variable.Assignment> getAssign() {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Retrieve a collection of assignments from entry.
     * @param key_index index of entry
     * @return set of assignments
     */
    public Set<Variable.Assignment> getAssign(int key_index) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Tag the sole entry with an assignment.
     * @param assign assignment
     * @return index of entry
     */
    public int addAssign(Variable.Assignment assign) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Tag the entry with an assignment.
     * @param key_index index of entry
     * @param assign assignment
     * @return index of entry
     */
    public int addAssign(int key_index, Variable.Assignment assign) {
        throw new CachedFactorRuntimeException("Not implemented");
    }

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
    public Distrib getDistrib(Variable nvar) {
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
    }

    /**
     * Retrieve the distribution of the specified non-enumerable variable,
     * conditioned on a particular instantiation of the enumerable variables in the table.
     * @param index the key index of the table, representing the instantiation
     * @param nvar non-enumerable variable
     * @return probability distribution for variable
     */
    public Distrib getDistrib(int index, Variable nvar) {
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
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
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
    }

    /**
     * Set the conditional distribution of a non-enumerable variable.
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(Object[] key, Variable nonenum, Distrib d) {
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
    }


    /**
     * Set the distribution of a non-enumerable variable for a factor without enumerable variables.
     * @param nvar the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry, always 0
     */
    public int setDistrib(Variable nvar, Distrib d) {
        throw new CachedFactorRuntimeException("Cached Factor table does not implement non-enumerable entries");
    }

    /**
     * Find out how many entries that are occupied.
     * @return  number of entries that are instantiated
     */
    public int getOccupied() {
        //        return occupied;
        throw new CachedFactorRuntimeException("Not implemented");
    }

    /**
     * Set all entries to 0.
     */
    public void setEmpty() {
        //         Arrays.fill(map, LOG0);
        throw new CachedFactorRuntimeException("Not implemented");
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
            throw new CachedFactorRuntimeException("Invalid key for EnumTable: key should be " + nEVars + " but is " + key.length + " values");
        int[] idx;          // size a function of domain sizes of null columns
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
        return new CachedFactor.FactorIterator();
    }

    private class FactorIterator implements Iterator<Integer> {
        private int index = 0;
        @Override
        public boolean hasNext() {
            int size = getSize();
            if (index < size) {
                do {
                    // TODO:
                    // if (map[index] != LOG0)
                    //    return true;
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
}

class CachedFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383971L;
    String message;

    public CachedFactorRuntimeException(String string) {
        message = string;
    }
}

class CachedTable {
    public double[] map;
    static int getHash(Enumerable[] enums, double[] map) {
        int hash = 7;
        for (int i = 0; i < enums.length; i++)
            hash = 31 * hash + enums[i].hashCode();
        for (int i = 0; i < map.length; i++)
            hash = 31 * hash + ((Double)map[i]).hashCode();
        return hash;
    }
}