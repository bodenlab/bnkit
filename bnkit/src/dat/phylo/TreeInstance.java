package dat.phylo;

import java.util.*;

/**
 * Class for representing the bare-bones of a tree (as defined by Tree via the TreeTopology interface) and
 * defining data structures and methods to work with enumerable values that can be assigned to tips and branch points.
 * This class is the vehicle for a unit of information to be transferred to inference algorithms.
 * Currently, maximum parsimony is an inner class as one inference algorithm that is designed to work on
 * TreeInstance instantiations.
 */
public class TreeInstance {

    /**
     * Turn on to print instantiations of branch points with an additional index (making it unique).
     */
    public static boolean LABEL_INCLUDES_INDEX = false;

    private final IdxTree  tree;        // index representation of a tree, reference to branch points is via an index
    private final Object[] instance;    // for each element, contains the value of a leaf or branch point, or null (un-instantiated)
    private final Object[] possible;    // all possible values that a leaf or branch point can take

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
            s.add(instance[i]);
        }
        this.possible = new Object[s.size()];
        s.toArray(this.possible);
    }

    public IdxTree getTree() {
        return tree;
    }

    public Object getInstance(int index) {
        return instance[index];
    }

    public Object getValueByIndex(int idx4val) {
        return possible[idx4val];
    }

    public int getIndexByValue(Object val) {
        for (int i = 0; i < possible.length; i++) {
            if (val.equals(possible[i]))
                return i;
        }
        return -1;
    }

    public int getNPossibleValues() {
        return possible.length;
    }

    public int getSize() {
        return tree.getSize();
    }

    public int[] getChildren(int index) {
        return tree.getChildren(index);
    }

    /**
     * String representation of the instantiated tree on the Newick format.
     * @return string representation
     */
    public String toString() {
        return toString(0);
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
            sb.append(toString(child));
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


}
