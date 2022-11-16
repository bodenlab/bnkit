package bn.factor;

import bn.Distrib;
import bn.JDF;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * CachedFactors are factors that makes reference to a factor table with products that are quite possibly
 * shared by other factor tables. This implies that it stores and retrieves doubles based ONLY on Enumerable keys.
 * The sharing implies that identities of variables can be withheld.
 * The intended use is for "factors", representing the factorisation of (conditional) probabilities.
 * CachedFactors use the fact that the products are re-used, which is especially useful for
 * phylogenetic analyses where similar/same relationships exists across sequence positions.
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

    AbstractFactor.FactorMap mymap = null;
    double atomic = LOG0;
    protected Set<Variable.Assignment>[] assigned = null; // the latent assignments associated with each factor, disabled by default
    // trace-back pointers for each (parent) factor, disabled by default;
    // can work out latent assignments from this--so will eventually replace "assigned" above
    protected Map<AbstractFactor, int[]> traceback = null;

    public static int
            SOURCE_BNODE = 1,
            SOURCE_PRODUCT = 2,
            SOURCE_MARGIN = 4;

    private int key = 0; // hash key if generated, 0 if not yet initialised
    private final FactorCache cache;

    /**
     * Create a factor associated with a specified cache.
     * This does NOT imply that any of its data is IN the cache BUT it could be if deemed beneficial for efficient access and storage.
     * @param cache
     */
    private CachedFactor(FactorCache cache) {
        super();
        this.cache = cache;
        setFactorType(AbstractFactor.TYPE_CACHED);
    }
    /**
     * Create a factor associated with a specified cache, based on a list of variables.
     * This does NOT imply that any of its data is IN the cache BUT it could be if deemed beneficial for efficient access and storage.
     * @param cache
     */
    private CachedFactor(FactorCache cache, EnumVariable... useVariables) {
        super(useVariables);
        this.cache = cache;
        setFactorType(AbstractFactor.TYPE_CACHED);
    }

    /**
     * Create a new instance of CachedFactor based on its (enumerable) variables
     * Place potential cached data in the specified cache.
     * @param cache cache
     * @param useVariables
     * @return new CachedFactor
     */
    public static CachedFactor getFactor(FactorCache cache, EnumVariable... useVariables) {
        return new CachedFactor(cache, useVariables);
    }

    /**
     * Create a new instance of CachedFactor.
     * Place potential cached data in the specified cache (in this case, none?).
     * @param cache cache
     * @return new CachedFactor
     */
    public static CachedFactor getFactor(FactorCache cache) {
        return new CachedFactor(cache);
    }

    /**
     * Create a new instance of CachedFactor based on the product of two (existing) CachedFactors.
     * Place potential cached data in the same cache as that shared by the two existing.
     * @param f1
     * @param f2
     * @return
     */
    public static CachedFactor getProduct(CachedFactor f1, CachedFactor f2) {
        EnumVariable[] evars = Factorize.getConcat(f1.getEnumVars(), f2.getEnumVars());
        if (f1.cache == f2.cache) {
            CachedFactor cf = new CachedFactor(f1.cache, evars);
            cf.key = AbstractFactor.getFactorMapKey(new int[]{f1.key, f2.key});
            FactorMap fmap = f1.cache.get(cf.key);
            if (fmap != null)
                cf.mymap = fmap;
            // else, not currently in cache, but the "key" is set so that when values are assigned by setLogValues they are cached under it
            return cf;
        }
        return null;
    }

    public static CachedFactor getMargin(CachedFactor f, EnumVariable[] sumout) {
        EnumVariable[] evars_inX = f.getEnumVars();
        EnumVariable[] evars_inY = Factorize.getDifference(evars_inX, sumout);
        CachedFactor cf = new CachedFactor(f.cache, evars_inY);
        cf.key = AbstractFactor.getFactorMapKey(new int[] {f.key, AbstractFactor.getFactorMapKey(evars_inY)}); // generate a unique key for caching the marginalised factor
        FactorMap fmap = f.cache.get(cf.key);
        if (fmap != null)
            cf.mymap = fmap;
        // else, not currently in cache, but the "key" is set so that when values are assigned by setLogValues they are cached under it
        return cf;
    }

    /**
     * Factory method for CachedFactor when the variables are factorised from a parent-to-child link with a distance.
     * Since that distance will give the same factor values every time, with variables known, we can uniquely identify that factor from distance
     * @param cache
     * @param d distance, note that the effect it has will vary with evolutionary model
     * @param childvar main variable in node
     * @param parentvar parent variable for node
     * @return new factor table, which can be cached once values have been assigned
     */
    public static CachedFactor getByDistance(FactorCache cache, double d, EnumVariable childvar, EnumVariable parentvar) {
        EnumVariable[] useVariables;
        if (childvar == null) // parentvar must be instantiated
            useVariables = new EnumVariable[] {parentvar};
        else if (parentvar == null)
            useVariables = new EnumVariable[] {childvar};
        else
            useVariables = new EnumVariable[] {childvar, parentvar};
        CachedFactor cf = new CachedFactor(cache, useVariables);
        // next generate a unique key for caching the marginalised factor;
        // note if child variable is marginalised out, the distance is flipped to generate a distinct key from when the parent variable is child and present in another factor
        cf.key = AbstractFactor.getFactorMapKey(AbstractFactor.getFactorMapKey(cf.enums, (childvar == null ? -d : d)));
        FactorMap fmap = cache.get(cf.key);
        if (fmap != null)
            cf.mymap = fmap;
        // else, not currently in cache, but the "key" is set so that when values are assigned by setLogValues they are cached under it
        return cf;
    }

    /**
     * Factory method for CachedFactor when the variables are factorised from a parent-to-child link with a distance.
     * Since that distance will give the same factor values every time, with variables known, we can uniquely identify that factor from distance
     * @param cache
     * @param d distance, note that the effect it has will vary with evolutionary model
     * @param childvar main variable in node
     * @param parentvar parent variable for node
     * @return new factor table, which can be cached once values have been assigned
     */
    public static CachedFactor getByDistanceAndDesignation(FactorCache cache, double d, EnumVariable childvar, EnumVariable parentvar, Object value) {
        EnumVariable[] useVariables;
        if (childvar == null) { // parentvar must be instantiated
            useVariables = new EnumVariable[]{parentvar};
        } else if (parentvar == null) {
            useVariables = new EnumVariable[]{childvar};
        } else {
            useVariables = new EnumVariable[]{childvar, parentvar};
        }
        CachedFactor cf = new CachedFactor(cache, useVariables);
        // next generate a unique key for caching the marginalised factor;
        // note if child variable is marginalised out, the distance is flipped to generate a distinct key from when the parent variable is child and present in another factor
        cf.key = AbstractFactor.getFactorMapKey(AbstractFactor.getFactorMapKey(cf.enums, (childvar == null ? -d : d), value));
        FactorMap fmap = cache.get(cf.key);
        if (fmap != null)
            cf.mymap = fmap;
        // else, not currently in cache, but the "key" is set so that when values are assigned by setLogValues they are cached under it
        return cf;
    }

    /**
     * Flag if all factor values have been set.
     * @return true if no further calculation is required (e.g. when a cache has been used, or values have been assigned); false otherwise
     */
    public boolean isSet() {
        return (mymap != null);
    }

    /**
     * Retrieve the log value of the factor, if table is without enumerable variables.
     * @return the only value of the factor
     */
    public double getLogValue() {
        return atomic;
    }

    /**
     * Retrieve the log value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    public double getLogValue(int index) {
        if (mymap == null)
            throw new CachedFactorRuntimeException("Factor is unassigned");
        return mymap.get(index);
    }

    /**
     * Set and cache the factor.
     * Will generate a cache key based on either
     * the enumerables and the factor values, OR
     * pre-set key to represent an operation that
     * consistently reproduces it, e.g. a factor product between two factors (each with their own hash keys)
     * @param logmap
     */
    @Override
    public void setLogValues(double[] logmap) {
        FactorMap cached = new FactorMap(logmap);
        mymap = cached;
        if (key == 0) {
            key = AbstractFactor.getFactorMapKey(enums, logmap);
            cache.put(key, cached); // this first checks if map is already cached, then adds it if not
        } else {
            // key was set when the factor was created
            cache.put(key, cached); // this first checks if map is already cached, then adds it if not
        }
    }

    // Setters of cells in factor ---------------------------------------------------------

    /**
     * Set the only value associated with a table without enumerable variables.
     * @param value
     * @return
     */
    public void setLogValue(double value) {
        atomic = value;
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
     * Retrieve the JDF of the factor without enumerable variables.
     * @return probability distribution for variable
     */
    public JDF getJDF() {
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

    //
    // Note: MB copied from DenseFactor the trace of assignment functions below
    // Consider how these can be stored more efficiently
    //

    /**
     * Activate tracing of implied assignments.
     * @param status true if activated, false otherwise
     */
    @Override
    public void setTraced(boolean status) {
        if (status) {
            this.traceback = new HashMap<>();
//            this.assigned = new Set[this.getSize()];
//            for (int i = 0; i < this.assigned.length; i ++)
//                this.assigned[i] = new HashSet();
        } else {
            this.assigned = null;
            this.traceback = null;
        }
    }

    /**
     * Find out if tracing of implied assignments is active.
     * @return true if activated, false otherwise
     */
    @Override
    public boolean isTraced() {
        return this.assigned != null || this.traceback != null;
    }

    /**
     * Retrieve a collection of assignments from the specified entry.
     * @param key_index the index of the entry
     * @return set of assignments
     */
    @Override
    public Map<Variable, Object> getAssign(int key_index) {
        if (key_index < 0)
            throw new DenseFactorRuntimeException("Invalid key index: outside map");
        Map<Variable, Object> collect = new HashMap<>();
        if (nEVars > 0) {
            Object[] key = this.getKey(key_index);
            for (int i = 0; i < key.length; i++) {
                collect.put(evars[i], key[i]);
            }
        }
        if (traceback != null) {
            for (Map.Entry<AbstractFactor, int[]> entry : traceback.entrySet()) {
                AbstractFactor f = entry.getKey();
                int[] row = entry.getValue();
                Map<Variable, Object> back = f.getAssign(row[key_index]);
                collect.putAll(back);
            }
        } else if (assigned != null)
            return Variable.Assignment.toMap(assigned[key_index]);
        return collect;
    }

    /**
     * Retrieve a collection of assignments from the single entry.
     * @return set of assignments
     */
    @Override
    public Map<Variable, Object> getAssign() {
        if (this.getSize() == 1 && traceback != null)
            return getAssign(0);
        if (this.getSize() == 1 && assigned != null)
            return Variable.Assignment.toMap(assigned[0]);
        throw new CachedFactorRuntimeException("Table has variables that must be used to index access");
    }

    /**
     * Associate the only entry with a collection of assignments.
     * @param assign assignment
     * @return index of entry
     */
    @Override
    public int addAssign(Collection<Variable.Assignment> assign) {
        if (!isTraced())
            throw new CachedFactorRuntimeException("Tracing is not enabled");
        if (this.getSize() != 1)
            throw new CachedFactorRuntimeException("This table must be accessed with a enumerable variable key");
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
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        if (key_index >= getSize() || key_index < 0 || getSize() == 1)
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
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
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        if (this.getSize() != 1)
            throw new CachedFactorRuntimeException("This table must be accessed with a enumerable variable key");
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
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        if (key_index >= getSize() || key_index < 0 || getSize() == 1)
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        if (assigned[key_index] == null)
            assigned[key_index] = new HashSet<>();
        assigned[key_index].add(assign);
        return key_index;
    }

    /**
     * @param key_index   index to row in present factor
     * @param from_factor trace-back reference to factor
     * @param from_index  index to row for trace-back factor
     * @return
     */
    @Override
    public int addAssign(int key_index, AbstractFactor from_factor, int from_index) {
        if (!isTraced())
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        if (key_index >= getSize() || key_index < 0 || getSize() == 1)
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        int[] indices;
        if (traceback == null) {
            traceback = new HashMap<>();
            indices = new int[this.getSize()];
            traceback.put(from_factor, indices);
        } else {
            indices = traceback.get(from_factor);
            if (indices == null) {
                indices = new int[this.getSize()];
                traceback.put(from_factor, indices);
            }
        }
        indices[key_index] = from_index;
        return key_index;
    }

    /**
     * @param from_factor trace-back reference to factor
     * @param from_index  index to row for trace-back factor
     * @return
     */
    @Override
    public int addAssign(AbstractFactor from_factor, int from_index) {
        if (!isTraced())
            throw new CachedFactorRuntimeException("Invalid key index: outside map");
        if (getSize() != 1)
            throw new CachedFactorRuntimeException("This table must be accessed with a enumerable variable key");
        int[] indices;
        if (traceback == null) {
            traceback = new HashMap<>();
            indices = new int[this.getSize()];
            traceback.put(from_factor, indices);
        } else {
            indices = traceback.get(from_factor);
            if (indices == null) {
                indices = new int[this.getSize()];
                traceback.put(from_factor, indices);
            }
        }
        indices[0] = from_index;
        return 0;
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
                    if (mymap.get(index) != LOG0)
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
}

class CachedFactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383971L;
    String message;

    public CachedFactorRuntimeException(String string) {
        message = string;
    }
}
