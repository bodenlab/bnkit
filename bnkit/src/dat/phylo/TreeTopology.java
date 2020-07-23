package dat.phylo;

/**
 * Interface to define a tree that can be queried, but not modified.
 * All branchpoints are numbered from 0 to n - 1, where 0 is root (by definition)
 * and n is the total number of branchpoints including leaves.
 * Each branchpoint except root will have exactly one parent branchpoint.
 * Each branchpoint will have any number of children, with leaves having 0 children.
 */
public interface TreeTopology {
    /**
     * Determine the number of branchpoints (including leaves) in the tree.
     * @return
     */
     int getSize();

    /**
     * Determine the parent index of branchpoint.
     * @param idx branchpoint index: 0 for root
     * @return the index of the parent or -1 if no parent (i.e. for root)
     */
     int getParent(int idx);

    /**
     * Determine the indices of the children of given branchpoint
     * @param idx branchpoint index: 0 for root
     * @return the indices of all children as an array if int, empty array if root
     */
     int[] getChildren(int idx);

    /**
     * Retrieve the identifier/label of the branch point, if available
     * @param idx the branch point index
     * @return label if available, else null; note that branch points are not always labelled
     */
     Object getLabel(int idx);

}
