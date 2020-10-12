package dat.phylo;

import java.util.*;

/**
 * Indexed (light-weight) view of tree, or set of trees.
 * Not meant to be edited, instead create new instances.
 *
 * @author mikael
 */
public class IdxTree implements Iterable<Integer> {
    protected final Map<Object, Integer> index; // map label to index; FIXME: need to index ancestor@position where position is relevant only to duplications
    protected final BranchPoint[] bpoints;  // store branch points in order of their index
    protected final int[] parent;           // [child idx] = parent idx, parent idx is -1 if no parent (i.e. root)
    protected final int[][] children;       // [parent idx] = [child 1][child 2]...[child n] where n is #children for parent
    protected final double[] distance;      // [child idx] = distance to parent

    /**
     * Constructor for internal (factory method) use only.
     * @param n number of branch points
     * @param distances if true, create storage for distances
     */
    private IdxTree(int n, boolean distances) {
        bpoints = new BranchPoint[n];
        parent = new int[bpoints.length];
        children = new int[bpoints.length][];
        index = new HashMap<>();
        distance = (distances ? new double[bpoints.length] : null);
    }

    /**
     * Create an "index tree", based on branch points which each have pointers to children
     * @param bpointarr
     */
    public IdxTree(BranchPoint[] bpointarr) {
        bpoints = bpointarr;
        parent = new int[bpoints.length];
        children = new int[bpoints.length][];
        index = new HashMap<>();
        double[] distances = new double[bpoints.length]; // temp array, may be used permanently
        // Check branch points, and index them
        boolean distances_found = true; // assume true until a non-ancestor does not have a distance
        Integer count = 0;
        for (int i = 0; i < bpoints.length; i ++) {
            BranchPoint bp = bpoints[i];
            children[i] = new int[bp.getChildren().size()];
            if (bp.getID() == null) // if leaf, and label is not set, OR if ancestor and ancestor ID is not set
                bp.setAncestor(count ++);
            if (index.containsKey(bp.getID()))
                throw new TreeRuntimeException("Branch points form cycle; duplicate found: label is " + bp.getLabel() + ", identifier is " + bp.getID());
            else {
                index.put(bp.getID(), i);
                BranchPoint ancestor = bp.getParent();
                if (ancestor != null) { // root node, possibly internal
                    Double dist = null;
                    try {
                        if (ancestor.getParent() != null)
                            distances[i] = ancestor.getDistance();
                    } catch (RuntimeException e) {
                        distances_found = false;
                    }
                    Integer paridx = index.get(ancestor.getID());
                    if (paridx == null)
                        throw new TreeRuntimeException("Branch points are out of order, so parent indexing failed at: " + bp.getID());
                    parent[i] = paridx;
                } else
                    parent[0] = -1;
            }
        }
        if (distances_found)
            this.distance = distances;
        else
            this.distance = null;
        for (int i = 0; i < bpoints.length; i ++) {
            for (int j = 0; j < children[i].length; j ++) {
                Integer idx = index.get(bpoints[i].getChildren().get(j).getID());
                if (idx != null)
                    children[i][j] = idx;
                else
                    throw new TreeRuntimeException("Invalid branch points, so child indexing failed at: " + bpoints[i].getID());
            }
        }
    }

    /**
     * Create a new IdxTree which is a reduced version of the source tree by excluding indices for nominated branch points, and
     * adjusting parent-child relations, as well as distances, which are assumed to be additive.
     * Note that the source branch points are not cloned, but are still referenced by the new tree
     * @return new tree with nominated branch points excluded
     */
    public static IdxTree createPrunedTree(IdxTree source, Set<Integer> pruneMe) {
        // Create the bare-bones of the new pruned tree
        IdxTree pruned = new IdxTree(source.getSize() - pruneMe.size(), source.distance != null);
        // dstidx counts branch points (kept from the source tree and) added to the new tree
        int dstidx = 0; // also used as the index for destination tree
        int[] match = new int[source.bpoints.length]; // a map for matching the source index to the destination index, -1 means no match
        // the source tree is traversed depth-first, srcidx is the index for each branch point
        for (int srcidx = 0; srcidx < source.bpoints.length; srcidx ++) {
            if (!pruneMe.contains(srcidx)) {    // branch point definitely to be transferred to new tree
                pruned.bpoints[dstidx] = source.bpoints[srcidx]; // transfer branch point data "by reference"
                match[srcidx] = dstidx;         // update source-to-dest mapping, to use below
                dstidx += 1;
            } else
                match[srcidx] = -1;             // update source-to-dest mapping, to indicate no match
        }
        for (int srcidx = 0; srcidx < source.bpoints.length; srcidx ++) { // loop again, this time to transfer parent and child pointers
            dstidx = match[srcidx];
            if (dstidx >= 0) { // this branch point is retained in new, destination tree
                if (pruned.distance != null) pruned.distance[dstidx] = source.distance[srcidx]; // transfer distance
                if (source.children[srcidx] != null) {                          // transfer children?
                    List<Integer> keep = new ArrayList<>();                     // collect what children that should remain
                    for (int j = 0; j < source.children[srcidx].length; j++) {  // by iterating through existing
                        if (match[source.children[srcidx][j]] >= 0)             // this child maps to a branch point that remains
                            keep.add(match[source.children[srcidx][j]]);        // so cache it
                    }
                    int[] replace = new int[keep.size()];                       // allocate array for storing those that remin
                    for (int j = 0; j < keep.size(); j++)                       // copy the pointers
                        replace[j] = keep.get(j);
                    pruned.children[dstidx] = replace;                   // place in index
                }
                if (source.parent[srcidx] >= 0) { // the original has a parent
                    // transfer parent pointer if it remains
                    pruned.parent[dstidx] = match[source.parent[srcidx]];
                } else
                    pruned.parent[dstidx] = -1;
            }
        }
        return pruned;
    }

    /**
     * Isolate a subtree rooted by a given node, and create a separate instance for it.
     * Note that the source branch points are not cloned, but are still referenced by the new tree
     * @param source source tree
     * @param rootidx the root branch point for the new, isolated tree by reference to the source tree indices
     * @return new tree instance
     */
    public static IdxTree createSubtree(IdxTree source, int rootidx) {
        // Create the bare-bones of the new isolated tree
        Set<Integer> subidxs = source.getSubtreeIndices(rootidx);
        IdxTree pruned = new IdxTree(subidxs.size(), source.distance != null);
        // dstidx counts branch points (kept from the source tree and) added to the new tree
        int[] match = new int[source.bpoints.length]; // a map for matching the source index to the destination index, -1 means no match
        for (int srcidx = 0; srcidx < source.bpoints.length; srcidx ++) match[srcidx] = -1; // default before mapping
        // the source tree is traversed depth-first, srcidx is the index for each branch point
        int dstidx = 0; // also used as the index for destination tree
        List<Integer> sorted_subidxs = new ArrayList<>(subidxs);
        Collections.sort(sorted_subidxs);
        for (int srcidx : sorted_subidxs) { // branch point definitely to be transferred to new tree
            pruned.bpoints[dstidx] = source.bpoints[srcidx]; // transfer branch point data "by reference"
            match[srcidx] = dstidx;         // update source-to-dest mapping, to use below
            dstidx += 1;
        }
        for (int srcidx : sorted_subidxs) { // loop again, this time to transfer parent and child pointers
            dstidx = match[srcidx];
            if (pruned.distance != null) pruned.distance[dstidx] = source.distance[srcidx]; // transfer distance
            if (source.children[srcidx] != null) {                          // transfer children?
                pruned.children[dstidx] = new int[source.children[srcidx].length];  // allocate array for storing those that remain
                for (int j = 0; j < source.children[srcidx].length; j++)    // by iterating through existing
                    pruned.children[dstidx][j] = match[source.children[srcidx][j]];
            }
            if (source.parent[srcidx] >= 0 && dstidx > 0) { // the original has a parent
                // transfer parent pointer if it remains
                pruned.parent[dstidx] = match[source.parent[srcidx]];
            } else
                pruned.parent[dstidx] = -1;
        }
        return pruned;
    }

    /**
     * Create an extended tree by duplicating a subtree from a nominated ancestor branch point.
     * @return tree with additional indices created for duplicated branch points
     */
    public static IdxTree createDuplicatedSubtree() {
        // FIXME: implement; note that "index" map needs to map from ancestor@position to recover duplicated branch points
        throw new RuntimeException("Not implemented");
    }

    @Override
    public Iterator<Integer> iterator() {
        return new IndexIterator();
    }

    public Iterator<Integer> getBreadthFirstIterator() {
        return new IndexIterator(IndexIterator.MODE_BREADTH_FIRST);
    }

    public Iterator<Integer> getDepthFirstIterator() {
        return new IndexIterator(IndexIterator.MODE_DEPTH_FIRST);
    }

    public TreeInstance getInstance(Object[] labels, Object[] values) {
        assert labels.length == values.length;
        Map<Object, Object> assign = new HashMap<>();
        for (int n = 0; n < labels.length; n ++)
            assign.put(labels[n], values[n]);
        return getInstance(assign);
    }

    /**
     * Instantiate the tree with values, so that other values can be inferred.
     * @param values selected identifiers with corresponding value
     * @return an instance of the tree that can be processed
     */
    public TreeInstance getInstance(Map<Object, Object> values) {
        Map<Integer, Object> idxvals = new HashMap<>();
        for (Map.Entry<Object, Object> entry : values.entrySet()) {
            Object val = entry.getValue();
            if (val != null) {
                Integer idx = index.get(entry.getKey());
                if (idx != null)
                    idxvals.put(idx, val);
            }
        }
        return new TreeInstance(this, idxvals);
    }

    /**
     * Determine the number of branchpoints (including leaves) in the tree.
     * @return number of branchpoints including root and leaves
     */
    public int getSize() {
        return bpoints.length;
    }

    /**
     * Determine the number of leaves in the tree.
     * @return the number of leaves/tips in the tree
     */
    public int getNLeaves() {
        int count = 0;
        for (BranchPoint bp : bpoints)
            count += bp.isLeaf() ? 1 : 0;
        return count;
    }

    /**
     * Determine if this branch point is a leaf (of the tree), i.e. does not have children.
     * @param idx
     * @return
     */
    public boolean isLeaf(int idx) {
        if (idx >= 0 && idx < getSize()) {
            if (children[idx] == null)
                return true;
            if (children[idx].length == 0)
                return true;
            return false;
        }
        throw new TreeRuntimeException("Invalid branch point index: " + idx);
    }

    /**
     * Determine the parent index of branchpoint.
     * @param idx branchpoint index: 0 for root
     * @return the index of the parent or -1 if no parent (i.e. for root)
     */
    public int getParent(int idx) {
        if (idx >= 0 && idx < getSize())
            return parent[idx];
        throw new TreeRuntimeException("Invalid branch point index: " + idx);
    }

    /**
     * Determine the indices of the children of given branchpoint
     * @param idx branchpoint index: 0 for root
     * @return the indices of all children as an array if int, empty array if root
     */
    public int[] getChildren(int idx) {
        if (idx >= 0 && idx < getSize())
            return children[idx];
        throw new TreeRuntimeException("Invalid branch point index: " + idx);
    }

    /**
     * Get the label of the node at the specified index. The label is either user specified (if leaf node) or
     * an automatically generated ancestor tag (an Integer instance, incrementing from 0).
     * @param idx
     * @return
     */
    public Object getLabel(int idx) {
        if (idx >= 0 && idx < getSize())
            return bpoints[idx].getID();
        throw new TreeRuntimeException("Invalid branch point index: " + idx);
    }

    /**
     * Retrieve distance to parent from the branch point.
     * @param idx index of branch point
     * @return distance to parent
     * @throws TreeRuntimeException if the tree does not have distances
     */
    public double getDistance(int idx) {
        if (this.distance != null)
            return this.distance[idx];
        else
            throw new TreeRuntimeException("No distance associated with this tree");
    }

    /**
     * Retrieve the branch point for a given index
     * @param idx
     * @return
     */
    public BranchPoint getBranchPoint(int idx) {
        if (idx >= 0 && idx < getSize())
            return bpoints[idx];
        throw new TreeRuntimeException("Invalid branch point index: " + idx);
    }

    /**
     * Get main root node of tree. Note that an IdxTree may have several subtrees, each with separate roots.
     * For that use IdxTree.getRoots()
     *
     * @return the root of the tree
     */
    public BranchPoint getRoot() {
        return bpoints[0];
    }

    /**
     * Determine the depth of the branch point (root is at 0)
     * @param idx
     * @return
     */
    public int getDepth(int idx) {
        if (idx == 0)
            return 0;
        else
            return getDepth(getParent(idx)) + 1;
    }

    /**
     * Determine size of a subtree, i.e. the number of branch points including leaves and internals.
     * @param idx
     * @return
     */
    public int getSubtreeSize(int idx) {
        int cnt = 1;
        for (int i : getChildren(idx)) {
            cnt += getSubtreeSize(i);
        }
        return cnt;
    }

    /**
     * Retrieve all indices that form part of the subtree under specified (root) index.
     * @param idx
     * @return
     */
    public Set<Integer> getSubtreeIndices(int idx) {
        Set<Integer> childset = new HashSet<>();
        for (int i : getChildren(idx)) {
            childset.addAll(getSubtreeIndices(i));
        }
        childset.add(idx);
        return childset;
    }

    public String getSubtreeString(int idx) {
        int[] ch = getChildren(idx);
        if (ch != null) {
            if (ch.length > 0) {
                StringBuilder sb = new StringBuilder();
                sb.append("(");
                for (int i = 0; i < ch.length; i ++)
                    sb.append(getSubtreeString(ch[i]) + (i < ch.length - 1 ? "," : ""));
                sb.append(")");
                return sb.toString();
            }
        }
        return getLabel(idx).toString();
    }

    @Override
    public String toString() {
        int[] r = getRoots();
        StringBuilder sb = new StringBuilder();
        for (int root : r) {
            sb.append(getSubtreeString(root) + "; ");
        }
        return sb.toString();
    }

    /**
     * Retrieve indices of all roots, i.e. branch points with no parents
     * @return indices to roots
     */
    public int[] getRoots() {
        int cnt = 0;
        for (int idx = 0; idx < bpoints.length; idx ++)
            if (parent[idx] < 0)
                cnt += 1;
        int[] ret = new int[cnt];
        cnt = 0;
        for (int idx = 0; idx < bpoints.length; idx ++)
            if (parent[idx] < 0)
                ret[cnt ++] = idx;
        return ret;
    }

    public int[] getLeaves() {
        int[] ret = new int[getNLeaves()];
        int cnt = 0;
        for (int idx = 0; idx < bpoints.length; idx ++)
            if (isLeaf(idx))
                ret[cnt ++] = idx;
        return ret;
    }

    public int[] getAncestors() {
        int[] ret = new int[getSize() - getNLeaves()];
        int cnt = 0;
        for (int idx = 0; idx < bpoints.length; idx ++)
            if (!isLeaf(idx))
                ret[cnt ++] = idx;
        return ret;
    }

    /**
     * Retrieve the index for a specified label (either user specified label or automatically assigned
     * @param label name or ancestor ID
     * @return index of node, or -1 if not found
     */
    public int getIndex(Object label) {
        for (int i : this) {
            if (this.getLabel(i).equals(label))
                return i;
        }
        return -1;
    }

    /**
     * Class to iterate through branch point indices.
     * By default a depth-first like traversal is used, but the implementation may change.
     * Instead, to be sure, use explicit flags to set traversal order.
     */
    public class IndexIterator implements java.util.Iterator<Integer> {
        public final static int MODE_BREADTH_FIRST = 0x01;
        public final static int MODE_DEPTH_FIRST = 0x02;
        public final static int MODE_RIGHT_FIRST = 0x04; // default is LEFT FIRST
        private int current = 0;
        private final int MODE;
        private LinkedList<Integer> queue = null;
        public IndexIterator() { /* initialize the Iterator */
            this(0);
        }
        private int updateQueue() {
            if (queue == null) {
                queue = new LinkedList<>();
            } else {
                if (queue.size() > 0)
                    current = queue.removeFirst();
                else
                    current = -1;
            }
            if (current == -1)
                return -1;
            int[] to_append = getChildren(current);
            if ((MODE & MODE_BREADTH_FIRST) > 0) {
                for (int i = 0; i < to_append.length; i++)
                    queue.addLast((MODE & MODE_RIGHT_FIRST) == 0 ? to_append[i] : to_append[to_append.length - i - 1]);
            } else if ((MODE & MODE_DEPTH_FIRST) > 0) {
                for (int i = 0; i < to_append.length; i++)
                    queue.addFirst((MODE & MODE_RIGHT_FIRST) > 0 ? to_append[i] : to_append[to_append.length - i - 1]);
            }
            return current;
        }
        public IndexIterator(int mode) { /* initialize the Iterator */
            this.MODE = mode;
            if (MODE == MODE_BREADTH_FIRST) {

            }
        }
        public boolean hasNext() { /* see if any elements remain */
            if (current < getSize() && current >= 0) {
                if (queue != null)
                    if (queue.size() == 0)
                        return false;
                return true;
            }
            return false;
        }
        public Integer next() { /* return next element, or throw exception if none */
            if (hasNext()) {
                if ((MODE & MODE_BREADTH_FIRST) > 0) {
                    return updateQueue();
                } else if ((MODE & MODE_DEPTH_FIRST) > 0) {
                    return updateQueue();
                } else
                    return current++;
            }
            throw new NoSuchElementException();
        }
    }

    /**
     * Helper function to extract the indices with specified value
     * @param findMe value to search for
     * @param instances an array representing instantiations across a tree
     * @return the indices which have the specified value
     */
    public static Set<Integer> getIndices(Object findMe, Object[] instances) {
        Set<Integer> collect = new HashSet<>();
        for (int i = 0; i < instances.length; i ++) {
            if (instances[i] == findMe)
                collect.add(i);
        }
        return collect;
    }

}
