package dat.phylo;

import asr.ASRRuntimeException;
import dat.file.Newick;

import java.io.IOException;
import java.util.*;

/**
 * Class for representing the bare-bones of a tree (as defined by Tree via the TreeTopology interface) and
 * defining data structures and methods to work with enumerable values that can be assigned to tips and branch points.
 * This class is the vehicle for a unit of information to be transferred to inference algorithms.
 */
public class TreeInstance {

    /**
     * Turn on to print instantiations of branch points with an additional index (making it unique).
     */
    public static boolean LABEL_INCLUDES_INDEX = false;

    private final IdxTree  tree;        // index representation of a tree, reference to branch points is via an index
    private final Object[] instance;    // for each element, contains the value of a leaf or branch point, or null (un-instantiated)
    private Object[] possible;    // all possible values that a leaf or branch point can take

    private static Object[] map2arr(int n, Map<Integer, Object> assignmap) {
        Object[] assign = new Object[n];
        for (Map.Entry<Integer, Object> entry : assignmap.entrySet())
            assign[entry.getKey()] = entry.getValue();
        return assign;
    }

    /**
     * Create a barebone structure to infer labels on a tree.
     * @param tree the definition of the tree topology
     * @param assignmap the values that should be assigned (input to inference)
     */
    public TreeInstance(IdxTree tree, Map<Integer, Object> assignmap) {
        this(tree, map2arr(tree.getSize(), assignmap));
    }

    /**
     * Create a barebone structure to infer labels on a tree.
     * @param tree the definition of the tree topology
     * @param assign the values that should be assigned (input to inference) in order of the
     *               branchpoints of the tree (index as defined in TreeTopology class)
     */
    public TreeInstance(IdxTree tree, Object[] assign) {
        this.tree = tree;
        this.instance = assign;
        Set s = new HashSet();
        for (int i = 0; i < instance.length; i ++) {
            if (instance[i] != null) //
                s.add(instance[i]);
        }
        this.possible = new Object[s.size()];
        s.toArray(this.possible);
    }

    /**
     * Create a new instance for a tree based on the values in another instance
     * (by reference to branch point IDs). No value is copied if an ID is not matched.
     * @param tree the tree from which an instance is created
     * @param from the instance from which values are taken
     * @return a new tree instance
     */
    public static TreeInstance copyByTree(IdxTree tree, TreeInstance from) {
        BranchPoint[] bp_new = tree.bpoints;
        Object[] vals_new = new Object[tree.bpoints.length];
        IdxTree src = from.getTree();
        BranchPoint[] bp_old = src.bpoints;
        for (int j = 0; j < bp_old.length; j ++) {
            Object val = from.getInstance(j);
            if (val != null) { // instantiated, so copy if included in new tree
                for (int i = 0; i < bp_new.length; i++) {
                    if (bp_new[i].getID().equals(bp_old[j].getID())) {
                        vals_new[i] = val;
                        break; // no point in continuing to loop
                    }
                }
            }
        }
        return new TreeInstance(tree, vals_new);
    }

    /**
     * Determine a standard code for each symbol used by TreeInstance, use it to recode the existing values, and
     * return the key with which the TreeInstance can be converted back to its old symbol set.
     * Modifies the TreeInstance "in-place".
     * @param recodeNullAtLeaves if true, null values are replaced with a distinct symbol at leaves, else they are replaced everywhere
     * @return the key {@see TreeInstance#decode()}
     */
    public synchronized Object[] encode(boolean recodeNullAtLeaves) {
        Object[] recoverkey = new Object[possible.length + (recodeNullAtLeaves ? 1 : 0)];
        Object[] newpossible = new Object[possible.length + (recodeNullAtLeaves ? 1 : 0)];
        for (int i = 0; i < possible.length; i ++) {
            recoverkey[i] = possible[i];
            newpossible[i] = i;
        }
        if (recodeNullAtLeaves) {
            recoverkey[possible.length] = null;
            newpossible[possible.length] = possible.length;
        }
        for (int i = 0; i < instance.length; i ++) {
            if (instance[i] != null) {
                instance[i] = getIndexByValue(instance[i]);
            } else { // instance[i] == null
                if (recodeNullAtLeaves) {
                    if (tree.isLeaf(i))
                        instance[i] = possible.length; // null converted to symbol
                }
            }
        }
        this.possible = newpossible;
        return recoverkey;
    }

    /**
     * Determine a standard code for each symbol used by TreeInstance, use it to recode the existing values, and
     * return the key with which the TreeInstance can be converted back to its old symbol set.
     * Modifies the TreeInstance "in-place".
     * @return the key {@see TreeInstance#decode()}
     */
    public Object[] encode() {
        return encode(true);
    }


    /**
     * Use the supplied key to convert the existing values to those provided through the key.
     * The key is assumed to have been generated by {@see TreeInstance#encode()}
     * @param key
     */
    public void decode(Object[] key) {
        if (key.length < 1)
            return;
        boolean encodeNullAtLeaves = key[key.length - 1] == null;
        Object[] newpossible = new Object[possible.length - (encodeNullAtLeaves ? 1 : 0)];
        for (int i = 0; i < newpossible.length; i ++) {
            newpossible[i] = key[i];
        }
        for (int i = 0; i < instance.length; i ++) {
            if (instance[i] != null) {
                Integer nom = (Integer) instance[i];
                instance[i] = key[nom];
            }
        }
        this.possible = newpossible;
    }

    /**
     * Retrieve the tree associated with the instance
     * @return tree
     */
    public IdxTree getTree() {
        return tree;
    }

    /**
     * Retrieve the value at a branch point
     * @param index the index of the branch point in tree
     * @return the value at the index
     */
    public Object getInstance(int index) {
        return instance[index];
    }

    /**
     * Retrieve the values at all branch points, indexed by the associated IdxTree
     * @return instance values (assigned or inferred)
     */
    public Object[] getInstance() {
        return instance;
    }

    /**
     * The branch points can be assigned n defined values. This function determines the value based on the index i in 0 to n-1.
     * @param idx4val index of value
     * @return the value
     */
    public Object getValueByIndex(int idx4val) {
        return possible[idx4val];
    }

    /**
     * The branch points can be assigned n defined values. This function determines the the index i in 0 to n-1 based on the value.
     * @param val the value
     * @return the index, or -1 if the value is invalid
     */
    public int getIndexByValue(Object val) {
        for (int i = 0; i < possible.length; i++) {
            if (val.equals(possible[i]))
                return i;
        }
        return -1;
    }

    /**
     * Get how many valid values there are
     * @return the number of distinct values
     */
    public int getNPossibleValues() {
        return possible.length;
    }

    /**
     * Get an array containing all possible values that can instantiate branchpoints, excluding null
     * @return possible values as an array
     */
    public Object[] getPossible() {
        return possible;
    }

    /**
     * Get the total number of branch points in the tree, including leaves
     * @return the number of branch points
     */
    public int getSize() {
        return tree.getSize();
    }

    /**
     * Get the indices for all children of the branch point
     * @param index the branch point index
     * @return a list of all children, which is empty if the branch point is a leaf
     */
    public int[] getChildren(int index) {
        return tree.getChildren(index);
    }

    /**
     * String representation of the instantiated tree on the Newick format.
     * @return string representation
     */
    public String toString() {
        return toString(0) + ";";
    }

    /**
     * String representation of the node and its children (recursively) on the Newick format.
     * @return string representation
     */
    public String toString(int bpidx) {
        StringBuilder sb = new StringBuilder();
        int[] children = tree.getChildren(bpidx);
        int cnt = 0;
        for (int child : children) {
            try {
                sb.append(toString(child) + ":" + tree.getDistance(child));
            } catch (TreeRuntimeException e) { // does not have distances
                sb.append(toString(child));
            }
            if (++ cnt < children.length)
                sb.append(",");
        }
        Object y = getInstance(bpidx);
        String istr = y == null ? "?" : (LABEL_INCLUDES_INDEX ? "#" + bpidx + "_" + y.toString() : y.toString());
        if (children.length < 1)
            return istr;
        else
            return "(" + sb.toString() + ")" + istr;
    }

    /**
     * Save instantiation of tree in Newick format.
     * If multiple trees, each will be given an index in the file extension.
     * @param filename
     * @throws IOException
     */
    public void save(String filename) throws IOException {
        Newick.save(this.tree, filename, this.instance);
    }

    public static class BlanketTreeDecor<E> implements TreeDecor<E> {
        E val = null;
        public BlanketTreeDecor(int nNodes, E val) {
            this.val = val;
        }

        @Override
        public E getDecoration(int idx) {
            return val;
        }

        @Override
        public void decorate(TreeInstance ti) {
            // nothing needs to be done
        }
    }
}
