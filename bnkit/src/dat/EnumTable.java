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
package dat;

import java.util.*;

/**
 * Table for storing and retrieving data of arbitrary types <E> based on
 * Enumerable keys
 *
 * @author mikael
 * @param <E>
 */
public class EnumTable<E> implements Iterable<Integer> {

    @Override
    public Iterator<Integer> iterator() {
        return map.keySet().iterator();
    }

    protected final Map<Integer, E> map;
    public final int nParents;
    protected final List<EnumVariable> parents;
    protected final EnumVariable[] pararr;
    protected final int[] period;
    protected final int[] step;
    protected final int[] domsize; // size of domain

    public EnumTable(EnumVariable... useParents) {
        this(EnumVariable.toList(useParents));
    }

    public EnumTable(Collection<EnumVariable> useParents) {
        this(useParents, new HashMap<Integer, E>());
    }

    private EnumTable(Collection<EnumVariable> useParents, Map<Integer, E> map) {
        this.map = map;
        this.parents = new ArrayList<>(useParents.size());
        this.pararr = new EnumVariable[useParents.size()];
        this.nParents = useParents.size();
        this.step = new int[this.nParents];
        this.period = new int[this.nParents];
        this.domsize = new int[this.nParents];
        int prod = 1;
        for (EnumVariable var : useParents) {
            this.parents.add(var);
        }
        for (int i = 0; i < nParents; i++) {
            int parent = nParents - i - 1;
            this.pararr[parent] = this.parents.get(parent);
            this.domsize[parent] = pararr[parent].size();
            this.step[parent] = prod;
            prod *= this.domsize[parent];
            this.period[parent] = prod;
        }
    }
    
    public EnumTable retrofit(List<EnumVariable> useParents) {
        if (this.nParents != useParents.size())
            throw new RuntimeException("Invalid retrofitting");
        for (int i = 0; i < this.nParents; i ++) {
            Variable p1 = this.getParents().get(i);
            Variable p2 = useParents.get(i);
            if (!p1.getDomain().equals(p2.getDomain()))
                throw new RuntimeException("Invalid retrofitting: " + p1.getName() + " does not share domain with " + p2.getName());
        }
        EnumTable dup = new EnumTable(useParents, this.getMap());
        return dup;
    }

    public boolean isEmpty() {
        return map.isEmpty();
    }
    
    public void setEmpty() {
        this.map.clear();
    }
    
    /**
     * Get the theoretical number of entries in this table. Note this number is
     * always greater or equal to the actual, populated number of entries.
     *
     * @return the size (number of entries)
     */
    public int getSize() {
        return period[0];
    }

    public List<EnumVariable> getParents() {
        return parents;
    }

    /**
     * Get the canonical names of the parent variables (names + "." + index)
     *
     * @return names of variables (in order)
     */
    public String[] getLabels() {
        String[] labels = new String[nParents];
        for (int i = 0; i < labels.length; i++) {
            labels[i] = this.pararr[i].toString();
        }
        return labels;
    }

    /**
     * Associate the specified key with the given value
     *
     * @param key
     * @param value
     * @return the index at which the value was stored
     */
    public int setValue(Object[] key, E value) {
        int index = this.getIndex(key);
        map.put(index, value);
        return index;
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
    public int setValue(int key_index, E value) {
        map.put(key_index, value);
        return key_index;
    }

    public int removeValue(int key_index) {
        map.remove(key_index);
        return key_index;
    }
    
    /**
     * Retrieve the index for the specified key
     *
     * @param key the values by which the index is identified (order the same as
     * when constructing the factor table)
     * @return the index for the instantiated key
     */
    public int getIndex(Object[] key) {
        int sum = 0;
        if (key.length != nParents) {
            throw new EnumTableRuntimeException("Invalid key: length is " + key.length + " not " + nParents);
        }
        for (int i = 0; i < nParents; i++) {
            if (key[i] == null) {
                throw new EnumTableRuntimeException("Null in key");
            }
            sum += (pararr[i].getIndex(key[i]) * step[i]);
        }
        return sum;
    }

    /**
     * Retrieve the index for the specified key
     *
     * @param key the values by which the index is identified (order the same as
     * when constructing the factor table)
     * @param vars enumerable variables on which to base the index (order is significant)
     * @return the index for the instantiated key
     */
    public static int getIndex(Object[] key, EnumVariable[] vars) {
        int sum = 0;
        if (key.length != vars.length) {
            throw new EnumTableRuntimeException("Invalid key: length is " + key.length + " not " + vars.length);
        }
        int prod = 1;
        for (int i = vars.length - 1; i >= 0; i --) {
            int step = prod;
            int domsize = vars[i].size();
            prod *= domsize;
            if (!vars[i].getDomain().isValid(key[i]))
                throw new EnumTableRuntimeException("Invalid key");
            sum += (vars[i].getIndex(key[i]) * step);
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
        int remain = index;
        if (index >= getSize() || index < 0) {
            throw new EnumTableRuntimeException("Invalid index");
        }
        Object[] key = new Object[nParents];
        for (int i = 0; i < nParents; i++) {
            int keyindex = remain / step[i];
            key[i] = pararr[i].getDomain().get(keyindex);
            remain -= keyindex * step[i];
        }
        return key;
    }

    /**
     * Check if there is a value assigned to the entry identified by the given index
     *
     * @param index the entry
     * @return true if assigned, false otherwise
     */
    public boolean hasValue(int index) {
        return map.containsKey(index);
    }

    /**
     * Retrieve the value of the entry identified by the given index
     *
     * @param index the entry
     * @return the value of the entry
     */
    public E getValue(int index) {
        E value = map.get(index);
        return value;
    }

    /**
     * Retrieve the value of the entry identified by the instantiated key
     *
     * @param key the entry
     * @return the value of the entry
     */
    public E getValue(Object[] key) {
        int index = getIndex(key);
        return getValue(index);
    }

    /**
     * Method to check if a key instance has the specified index. Must be fast
     * since it may be called frequently.
     *
     * @param key
     * @param index
     * @return true if the key instance maps to the index
     */
    public boolean isMatch(Object[] key, int index) {
        if (key.length != nParents || index < 0 || index >= getSize()) {
            throw new EnumTableRuntimeException("Invalid index or key");
        }
        int remain = index;
        for (int i = 0; i < nParents; i++) {
            if (key[i] != null) {
                int keyindex = pararr[i].getIndex(key[i]);
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
     * Retrieve an array of non-null values, with each element associated with
     * an entry matching the specified key. The key may contain values and
     * "wildcards" (with null matching any value for the corresponding parent)
     *
     * @param key the instantiated key (may contain null to indicate "any"
     * matches)
     * @return an array with values
     */
    public List<E> getValues(Object[] key) {
        boolean allWildcards = true;
        for (Object y : key) {
            if (y != null) {
                allWildcards = false;
                break;
            }
        }
        List<E> values = new ArrayList<>();
        if (allWildcards) {
            values.addAll(map.values());
        } else {
            for (Map.Entry<Integer, E> entry : map.entrySet()) {
                if (isMatch(key, entry.getKey())) {
                    values.add(entry.getValue());
                }
            }
        }
        return values;
    }

    /**
     * Get all values indiscriminately
     *
     * @return the collection of values associated with all keys
     */
    public Collection<E> getValues() {
        return map.values();
    }

    /**
     * Get all entries in the table
     *
     * @return the map entries, each an index and the associated value
     */
    public Set<Map.Entry<Integer, E>> getMapEntries() {
        return map.entrySet();
    }

    private Map<Integer, E> getMap() {
        return map;
    }
    
    /**
     * Identify each index that is linked to the specified key (which may
     * include "wildcards", indicated by null values).
     *
     * @param key
     * @return an array with all matching indices
     */
    public int[] getIndices(Object[] key) {
        boolean allWildcards = true;
        for (Object y : key) {
            if (y != null) {
                allWildcards = false;
                break;
            }
        }
        if (allWildcards) { // 
//            int[] ret = new int[this.getSize()];
//            for (int i = 0; i < ret.length; i++) {
//                ret[i] = i;
//            }
//            return ret;
            return getIndices();
        }
        List<Integer> indices = new ArrayList<>();
        for (Map.Entry<Integer, E> entry : map.entrySet()) {
            if (isMatch(key, entry.getKey())) {
                indices.add(entry.getKey());
            }
        }
        int[] ret = new int[indices.size()];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = indices.get(i);
        }
        return ret;
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
    public int[] getTheoreticalIndices(Object[] key) {
        if (key.length != nParents)
            throw new EnumTableRuntimeException("Invalid key for EnumTable: key should be " + nParents + " but is " + key.length + " values");
        int[] idx; // size a function of domain sizes of null columns
        int startentry = 0; // determined from non-null entries, where counting will start
        int tot = 1;
        for (int i = 0; i < key.length; i ++) {
            if (key[i] == null) {
                tot *= domsize[i];
            } else {
                int keyidx = pararr[i].getIndex(key[i]);
                startentry += keyidx * step[i];
            }
        }
        idx = new int[tot];
        if (tot == 1) // no null values
            idx[0] = startentry;
        else
            getTheoreticalIndicesRecursive(key, idx, 0, startentry, 0);
        return idx;
    }

    private synchronized int getTheoreticalIndicesRecursive(Object[] key, int[] idx, int my_idx, int my_tab, int parent) {
        if (key.length <= parent) 
            return -1;
        // parent is real
        if (key[parent] == null) {
            for (int i = 0; i < domsize[parent]; i ++) {
                int start = getTheoreticalIndicesRecursive(key, idx, my_idx, my_tab, parent + 1);
                if (start != -1) {
                    my_idx = start;
                    my_tab += step[parent];
                } else { // start == null meaning that we're at leaf
                    idx[my_idx ++] = my_tab;
                    my_tab += step[parent];
                }
            }
        } else {
            return getTheoreticalIndicesRecursive(key, idx, my_idx, my_tab, parent + 1);
        }
        return my_idx;
    }
    
    /**
     * Get indices for all non-null entries.
     * @return the indices
     */
    public int[] getIndices() {
        Set<Integer> all = map.keySet();
        int[] arr = new int[all.size()];
        int i = 0;
        for (Integer y : all) 
            arr[i ++] = y;
        return arr;
    }
    

    /**
     * Takes an entry index of the current table and "masks" out a subset of
     * parents, to determine the index in a table with only the parents that are
     * not masked. Time complexity is O(3n) where n is the number of parents in
     * the current table. (This computation could be done marginally more
     * efficiently.)
     *
     * @param origindex
     * @param maskMe
     * @return index in new more compact table
     */
    public int maskIndex(int origindex, Collection<EnumVariable> maskMe) {
        int origremain = origindex;
        int sum = 0;
        int jn = 0;
        int[] newstep = new int[nParents - maskMe.size()];
        int[] newvale = new int[nParents - maskMe.size()];
        for (int i = 0; i < nParents; i++) {
            if (!maskMe.contains(parents.get(i))) // if NOT masked-out
            {
                newvale[jn++] = parents.get(i).size();
            }
        }
        jn = newstep.length - 1;
        int prod = 1;
        for (int i = nParents - 1; i >= 0; i--) {
            if (!maskMe.contains(parents.get(i))) { // if NOT masked-out
                newstep[jn] = prod;
                prod *= newvale[jn--];
            }
        }
        jn = 0;
        for (int i = 0; i < nParents; i++) {
            int key = origremain / step[i];
            origremain -= key * step[i];
            if (!maskMe.contains(parents.get(i))) // if NOT masked-out
            {
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
        return EnumTable.getKey(this.parents, evid);
    }

    /**
     * Create indices that allow quick cross-referencing between this table's variables and a specified list.
     * The parameters need to be allocated as they are filled "in place."
     * @param xcross indices to go from X to Y, i.e. at [x] you find y (or -1 if no mapping)
     * @param yvars variables in a table Y
     * @param ycross indices to go from Y to X, i.e. at [y] you find x (or -1 if no mapping)
     * @return number of positions that overlap
     */
    public  int crossReference(int[] xcross, Variable[] yvars, int[] ycross) {
        Variable[] xvars = pararr;
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

    public void display() {
        System.out.print("Idx ");
        for (int j = 0; j < this.nParents; j++) {
            System.out.print(String.format("[%10s]", this.parents.get(j).toString()));
        }
        System.out.println(" P");
        for (int i = 0; i < this.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = this.getKey(i);
            Object val = this.getValue(i);
            for (Object key1 : key) {
                System.out.print(String.format(" %-10s ", key1.toString()));
            }
            if (val != null) {
                System.out.println(String.format(" %5s", this.getValue(i).toString()));
            } else {
                System.out.println(" null ");
            }

        }
    }

    /**
     * Get a key from a partial assignment of variables defined for the table.
     * @param <T>
     * @param vars variables in order
     * @param evid the evidence
     * @return the key that encodes the values for the provided variables
     */
    public static <T extends Variable> Object[] getKey(List<T> vars, Variable.Assignment[] evid) {
        Object[] key = new Object[vars.size()]; // key is initialised to nulls by default
        if (key.length <= evid.length) { // for efficiency, we check what to iterate over
            for (Variable.Assignment e : evid) {
                try {
                    EnumVariable evar = (EnumVariable) e.var;
                    int var_index = vars.indexOf(evar);
                    if (var_index >= 0)
                        key[var_index] = e.val;
                } catch (ClassCastException exception) {
                    // ignore non-enumerable variables
                }
            }
        } else { // evidence is longer than key, so let's iterate over key
            for (int i = 0; i < key.length; i ++) {
                Variable var = vars.get(i);
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
            throw new EnumTableRuntimeException("Invalid operation since keys are of difference lengths (" + target.length + " vs " + source.length + ")");
        for (int i = 0; i < target.length; i ++)
            if (source[i] != null)
                target[i] = source[i];
        return target;
    }
    
}

class EnumTableRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public EnumTableRuntimeException(String string) {
        message = string;
    }
}
