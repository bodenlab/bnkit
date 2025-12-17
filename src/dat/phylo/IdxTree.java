package dat.phylo;

import bn.Distrib;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;
import dat.Enumerable;
import dat.file.Newick;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;
import smile.stat.distribution.ExponentialFamilyMixture;
import smile.stat.distribution.GammaDistribution;
import smile.stat.distribution.Mixture;
import stats.RateModel;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Indexed (light-weight) view of tree, or set of trees.
 * Not meant to be edited, instead create new instances.
 *
 * @author mikael
 */
public class IdxTree implements Iterable<Integer> {
    protected final Map<Object, Integer> index; // map label to index; FIXME: need to index ancestor@position where position is relevant only to duplications
    protected final BranchPoint[] bpoints;      // store branch points in order of their index
    protected final int[] parent;               // [child idx] = parent idx, parent idx is -1 if no parent (i.e. root)
    protected final int[][] children;           // [parent idx] = [child 1][child 2]...[child n] where n is #children for parent
    protected final double[] distance;          // [child idx] = distance to parent

    /**
     * Constructor for internal (factory method) use only.
     * @param n number of branch points
     * @param usedistances if true, create storage for distances
     */
    private IdxTree(int n, boolean usedistances) {
        bpoints = new BranchPoint[n];
        parent = new int[bpoints.length];
        children = new int[bpoints.length][];
        index = new HashMap<>();
        distance = (usedistances ? new double[bpoints.length] : null);
    }

    /**
     * Calculates the depth of each branch point in the tree.
     * The depth is defined as the number of edges from the root to the branch point.
     *
     * @return an array where each element represents the depth of the corresponding branch point
     */
    public int[] getDepth() {
        int[] depth = new int[parent.length];
        for (int i = 0; i < parent.length; i ++) { // go through all indices
            int myparent = parent[i];   // parent of the current index
            if (myparent == -1)         // if root, depth is 0
                depth[i] = 0;
            else
                depth[i] = depth[myparent] + 1;
        }
        return depth;
    }

    /**
     * Calculates the width at each branch point in the tree.
     * The width is defined as the number of leaves under the branch point.
     *
     * @return an array where each element represents the width of the corresponding branch point
     */
    public int[] getWidth() {
        int[] width = new int[parent.length];
        for (int i = 0; i < parent.length; i++) {
            width[i] = countLeaves(i);
        }
        return width;
    }

    /**
     * Calculates the whole-number log2 width at each branch point in the tree.
     * The width is defined as the number of leaves under the branch point.
     *
     * @return an array where each element is the log2 of the width of the corresponding branch point
     */
    public double[] getLog2Width() {
        double[] width = new double[parent.length];
        for (int i = 0; i < parent.length; i++) {
            int nleaves = countLeaves(i);
            width[i] = Math.log(nleaves)/Math.log(2);
        }
        return width;
    }

    public double[] getLog2WidthThresholds(int nbins, double[] log2width) {
        double[] log2sorted = Arrays.copyOf(log2width, log2width.length);
        Arrays.sort(log2sorted);
        double[] thresholds = new double[Math.max(nbins, 2)];
        thresholds[0] = 0; // all leaf nodes are in the first bin
        thresholds[1] = 1; // all ancestors with two leaves are in the second bin
        // if only two bins, we are done because we include all ancestors in the second bin
        if (nbins > 2) { // three bins or more, we'll have to figure out what other thresholds to us
            int start = 0;
            while (log2sorted[start] == 0 || log2sorted[start] == 1)
                start++;
            int binsize = (log2width.length - start) / (nbins - 2);
            for (int t = 2; t < thresholds.length; t++)
                thresholds[t] = log2sorted[start + (t - 2) * binsize];
        }
        return thresholds;
    }

    public int[] getBinned(double[] thresholds, double[] values2bin) {
        int[] binned = new int[parent.length];
        for (int i = 0; i < parent.length; i++) {
            for (int j = 0; j < thresholds.length; j++) {
                if (values2bin[i] <= thresholds[j]) {
                    binned[i] = j;
                    break;
                }
            }
        }
        return binned;
    }

    public int[] getBinned(int nbins, double[] values2bin) {
        double[] thresholds = getLog2WidthThresholds(nbins, values2bin);
        return getBinned(thresholds, values2bin);
    }

    /**
     * Recursively counts the number of leaves under a given branch point.
     *
     * @param idx the index of the branch point
     * @return the number of leaves under the branch point
     */
    private int countLeaves(int idx) {
        if (isLeaf(idx)) {
            return 1;
        }
        int count = 0;
        for (int child : getChildren(idx)) {
            count += countLeaves(child);
        }
        return count;
    }

    /**
     * Get the indices of the leaves in the tree that are connected to the given leaf;
     * limit the search to a given level (0 for root, 1 for the next level of sub-trees, etc.)
     *
     * @param leafcnt the left-to-right count of the leaf
     * @param level the level of the search
     * @param dmat  the distance matrix
     * @return an array of count indices of the connected leaves
     */
    public static int[] getConnected(int leafcnt, int level, Double[][] dmat) {
        List<Integer> idxs = new ArrayList<>();
        int skip = 0; // number of skipped columns (corresponds to level)
        int select = -1;
        for (int b = 0; b < dmat[leafcnt].length; b ++) {
            if (skip == level && dmat[leafcnt][b] != null) {
                select = b;
                break;
            }
            if (dmat[leafcnt][b] != null)
                skip += 1;
        }
        if (select > 0) { // valid
            for (int j = 0; j < dmat.length; j ++) {
                if (dmat[j][select] != null && j != leafcnt)
                    idxs.add(j);
            }
        }
        int[] connected = new int[idxs.size()];
        int cnt = 0;
        for (Integer i : idxs)
            connected[cnt ++] = i;
        return connected;
    }

    /**
     * Calculate the root-to-leaf distances for each leaf in the tree.
     * @param dmat the distance matrix
     * @return  an array of distances from the root to each leaf
     */
    protected static double[] getLeafDistances(Double[][] dmat) {
        double[] leafdist = new double[dmat.length];
        for (int i = 0; i < dmat.length; i ++) { // each row is a leaf
            for (int j = 0; j < dmat[i].length; j ++) // each col is a branch (with a node idx)
                leafdist[i] += dmat[i][j] == null ? 0 : dmat[i][j];
        }
        return leafdist;
    }

    /**
     * Determine a Gaussian distribution that models the leaf-to-root distances of the tree
     * @param dmat
     * @return
     */
    protected static GaussianDistrib getLeafDistanceDistribution(Double[][] dmat) {
        double[] leafdist = getLeafDistances(dmat);
        return GaussianDistrib.estimate(leafdist);   // estimate distribution
    }

    /**
     * Determine a Gaussian distribution that models the leaf-to-root distances of the tree
     * @return
     */
    public GaussianDistrib getLeaf2RootDistrib() {
        Double[][] dmat = getDistance2RootMatrix();
        double[] leafdist = getLeafDistances(dmat);
        return GaussianDistrib.estimate(leafdist);   // estimate distribution
    }

    /**
     * Constructs a matrix of distances from each leaf node to the root.
     * The matrix is represented as a 2D array where each row corresponds to a leaf node,
     * and each column contains the distance from that leaf node to the root.
     *
     * @return a 2D array where each element represents the distance from a leaf node to the root
     */
    public Double[][] getDistance2RootMatrix() {
        int[] leaves = getLeaves();
        Double[][] distmat = new Double[leaves.length][distance.length];
        for (int cnt = 0; cnt < leaves.length; cnt ++) {
            int idx = leaves[cnt];
            int parent = getParent(idx);
            distmat[cnt][idx] = getDistance(idx);
            while (parent != -1) {
                distmat[cnt][parent] = getDistance(parent);
                parent = getParent(parent);
            }
        }
        return distmat;
    }

    /**
     * Sets the distances from each leaf node to the root based on the given matrix.
     * The matrix is represented as a 2D array where each row corresponds to a leaf node,
     * and each column contains the distance from that leaf node to the root.
     *
     * @param distmat a 2D array where each element represents the distance from a leaf node to the root
     */
    public void setDistance2RootMatrix(Double[][] distmat) {
        int[] leaves = getLeaves();
        for (int cnt = 0; cnt < leaves.length; cnt ++) {
            int idx = leaves[cnt];
            int parent = getParent(idx);
            setDistance(idx, distmat[cnt][idx]);
            while (parent != -1 && distmat[cnt][parent] != null) {
                setDistance(parent, distmat[cnt][parent]);
                parent = getParent(parent);
            }
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        IdxTree integers = (IdxTree) o;
        return Objects.equals(index, integers.index) && Arrays.equals(parent, integers.parent) && Arrays.equals(distance, integers.distance);
    }

    @Override
    public int hashCode() {
        String[] labels = new String[this.getNLeaves()];
        int cnt = 0;
        for (int i : this.getLeaves())
            labels[cnt ++] = getLabel(i).toString();
        int result = Objects.hash(labels);
        result = 31 * result + Arrays.hashCode(parent);
        result = 31 * result + Arrays.hashCode(distance);
        return result;
    }

    public static IdxTree fromJSON(JSONObject json) {
        int n = json.getInt("Branchpoints");
        JSONArray jdists = json.getJSONArray("Distances");
        IdxTree tree = new IdxTree(n, jdists != null);
        JSONArray jlabs = json.getJSONArray("Labels");
        JSONArray jpars = json.getJSONArray("Parents");
        for (int p = 0; p < n; p ++) {
            try {
                tree.parent[p] = jpars.getInt(p);
                if (jdists != null)
                    tree.distance[p] = jdists.getDouble(p);
                tree.index.put(jlabs.getString(p), p);
            } catch (JSONException e) {
                throw new TreeRuntimeException("Invalid JSON format: " + e.getMessage());
            }
        }
        for (int p = 0; p < n; p ++) {
            // check through all children if they have p as parent
            List<Integer> ch = new ArrayList<>();
            for (int c = 0; c < n; c ++) {
                if (tree.parent[c] == p) // if c has p as parent add it to list
                    ch.add(c);
            }
            // convert list to array
            tree.children[p] = new int[ch.size()];
            for (int i = 0; i < ch.size(); i ++)
                tree.children[p][i] = ch.get(i);
        }
        // finally fix the branchpoints
        for (int p = 0; p < n; p ++) {
            BranchPoint bp = new BranchPoint(jlabs.getString(p));
            if (tree.distance != null)
                bp.setDistance(tree.distance[p]);
            tree.bpoints[p] = bp;
        }
        Integer ancid = 0;
        for (int p = 0; p < n; p ++) {
            BranchPoint bp = tree.bpoints[p];
            if (tree.parent[p] >= 0)
                bp.setParent(tree.bpoints[tree.parent[p]]);
            for (int c = 0; c < tree.children[p].length; c ++)
                bp.addChild(tree.bpoints[tree.children[p][c]]);
            if (bp.isParent())
                bp.setAncestor(ancid ++);
        }
        return tree;
    }

    /**
     * Convert instance to JSON
     * @return
     */
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        json.put("Branchpoints", bpoints.length);
        String[] labels = new String[bpoints.length];
        for (Map.Entry<Object, Integer> entry : index.entrySet()) {
            //labels[entry.getValue()] = entry.getKey().toString();
            labels[entry.getValue()] = (String) this.getBranchPoint(entry.getValue()).getLabel();
        }
        JSONArray labarr = new JSONArray();
        JSONArray pararr = new JSONArray();
        JSONArray distarr = new JSONArray();
        for (int i = 0; i < bpoints.length; i ++) {
            labarr.put(labels[i]);
            pararr.put(parent[i]);
            if (distance != null)
                distarr.put(distance[i]);
        }
        json.put("Labels", labarr);
        json.put("Parents", pararr);
        if (distance != null)
            json.put("Distances", distarr);
        return json;
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
                    try {
                        if (bp.getParent() != null)
                            distances[i] = bp.getDistance();
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
     * Create an index from the current tree to create a new pruned version, by excluding branch points nominated in the parameter pruneMe
     * @param pruneMe indices of branch points to be removed
     * @return index with nominated branch points excluded
     */
    public int[] getPrunedIndex(Set<Integer> pruneMe, boolean removeOrphans) {
        // dstidx counts branch points (kept from the source tree and) added to the new tree
        int dstidx = 0; // also used as the index for destination tree
        int[] match = new int[bpoints.length]; // a map for matching the source index to the destination index, -1 means no match

        if (removeOrphans) {
            boolean[] orphan = new boolean[bpoints.length]; // map for orphans (nodes with no extant descendant)
            Arrays.fill(orphan, true);                  // assume all positions are orphans
            for (int srcidx = 0; srcidx < bpoints.length; srcidx++) {
                if (this.isLeaf(srcidx)) { // at leaf so traverse upwards to mark non-orphan status
                    int top = -1, idx = srcidx;
                    while (!pruneMe.contains(idx) && orphan[idx]) { // stop once hitting a pruned node or one already marked non-orphan
                        orphan[idx] = false;        // mark it non-orphan
                        top = idx;
                        idx = this.getParent(idx);  // go to parent
                        if (idx == -1)              // hitting a node without parent
                            break;
                    }
                    // now flood all descendants (as non-orphans) because we've hit the top of a non-orphaned sub-tree
                    if (top != -1) {
                        markMyChildrenAsNonOrphans(top, pruneMe, orphan);
                    }
                }
            }
            for (int srcidx = 0; srcidx < bpoints.length; srcidx++) {
                if (orphan[srcidx])
                    pruneMe.add(srcidx);
            }
        } // any orphans are now included in the pruneMe set

        // build the match index (to go from source tree idx to pruned position-specific tree
        // the source tree is traversed depth-first, srcidx is the index for each branch point
        for (int srcidx = 0; srcidx < bpoints.length; srcidx ++) {
            if (!pruneMe.contains(srcidx)) {    // branch point is not in the nominated kill-list, so WILL BE transferred to new tree
                match[srcidx] = dstidx;         // update source-to-dest mapping, to use below
                dstidx += 1;
            } else {
                match[srcidx] = -1;             // update source-to-dest mapping, to indicate no match
            }
        }
        return match;
    }

    private void markMyChildrenAsNonOrphans(int idx, Set<Integer> pruned, boolean[] orphan) {
        if (isLeaf(idx) || pruned.contains(idx))
            return;
        int[] children = getChildren(idx);
        for (int c : children) {
            if (orphan[c] == false || pruned.contains(c))
                continue;
            orphan[c] = false;
            markMyChildrenAsNonOrphans(c, pruned, orphan);
        }
    }
    /**
     *
     * @param source
     * @param prunedIndex
     * @return
     */
    public static IdxTree createPrunedTree(IdxTree source, int[] prunedIndex) {
        int cnt = 0;
        for (int i : prunedIndex)
            if (i != -1)
                cnt += 1;
        // Create the bare-bones of the new pruned tree
        IdxTree pruned = new IdxTree(cnt, source.distance != null);
        for (int srcidx = 0; srcidx < source.bpoints.length; srcidx ++) { // loop again, this time to transfer parent and child pointers
            int dstidx = prunedIndex[srcidx];
            if (dstidx >= 0) { // this branch point is retained in new, destination tree
                pruned.bpoints[dstidx] = source.bpoints[srcidx]; // transfer branch point data "by reference"
                if (pruned.distance != null)
                    pruned.distance[dstidx] = source.distance[srcidx]; // transfer distance
                if (source.children[srcidx] != null) {                          // transfer children?
                    List<Integer> keep = new ArrayList<>();                     // collect what children that should remain
                    for (int j = 0; j < source.children[srcidx].length; j++) {  // by iterating through existing
                        if (prunedIndex[source.children[srcidx][j]] >= 0)             // this child maps to a branch point that remains
                            keep.add(prunedIndex[source.children[srcidx][j]]);        // so cache it
                    }
                    int[] replace = new int[keep.size()];                       // allocate array for storing those that remin
                    for (int j = 0; j < keep.size(); j++)                       // copy the pointers
                        replace[j] = keep.get(j);
                    pruned.children[dstidx] = replace;                   // place in index
                }
                if (source.parent[srcidx] >= 0) { // the original has a parent
                    // transfer parent pointer if it remains
                    pruned.parent[dstidx] = prunedIndex[source.parent[srcidx]];
                } else
                    pruned.parent[dstidx] = -1;
            }
        }
        return pruned;
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
                    int[] replace = new int[keep.size()];                       // allocate array for storing those that remain
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

    public Set<Integer> getIndicesOfOrphanedTrees() {
        int[] roots = this.getRoots();
        Set<Integer> prune = new HashSet<>();
        Set<Integer> members = new HashSet<>();
        for (int r : roots) {
            for (int i : getSubtreeIndices(r)) {
                if (getBranchPoint(i).isLeaf()) {
                    members.clear();
                    break;
                } else
                    members.add(i);
            }
            prune.addAll(members); // will only be non-empty if no members is a proper leaf
            members.clear();
        }
        return prune;
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
                Integer idx = getIndex(entry.getKey());
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
     * Determine the number of ancestors in the tree.
     * @return the number of ancestors in the tree
     */
    public int getNParents() {
        int count = 0;
        for (int idx = 0; idx < children.length; idx ++)
            if (children[idx] != null)
                count += (children[idx].length > 0) ? 1 : 0;
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
     * Determine if this branch point is a parent, i.e. has children.
     * @param idx
     * @return
     */
    public boolean isParent(int idx) {
        if (idx >= 0 && idx < getSize()) {
            if (children[idx] == null)
                return false;
            if (children[idx].length > 0)
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
     * Determine if branchpoint is connected to others, i.e. have either children or parents.
     * @param idx branchpoint index: 0 for root
     * @return true if connected, false otherwise
     */
    public boolean isConnected(int idx) {
        if (this.children[idx].length > 0 || this.parent[idx] != -1)
            return true;
        else
            return false;
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
     * Set distance to parent from the branch point.
     * @param idx index of branch point
     * @param dist distance
     */
    public void setDistance(int idx, double dist) {
        distance[idx] = dist;
        bpoints[idx].setDistance(dist);
    }

    /**
     * Determine the distance to root from the branch point
     * @param idx branchpoint index
     * @return additive distance to root (of a subtree in which the branchpoint is located)
     */
    public double getDistanceToRoot(int idx) {
        if (idx == 0)
            return 0;
        else
            return getDistanceToRoot(getParent(idx)) + distance[idx];
    }

    /**
     * Retrieve the branch point for a given index
     * @param idx branchpoint index
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

    /**
     * Text representation (Newick) of subtree with the specified node index as root
     * @param idx index of the root of the subtree
     * @return text string representation on Newick format
     */
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
     * Retrieve indices of all roots, i.e. branch points with no parents (but children, TODO)
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

    /**
     * Get the names of all leaf nodes in this tree.
     * @return array with names (in order of increasing branch point index)
     */
    public String[] getNames() {
        int[] idxs = getLeaves();
        String[] ret = new String[idxs.length];
        for (int i = 0; i < ret.length; i ++) {
            ret[i] = (String)this.getLabel(idxs[i]);
        }
        return ret;
    }

    /**
     * Get the labels of all nodes in this tree, as indexed by the tree.
     * @return array with labels (in order of increasing branch point index, from 0 to N-1, where N is the number of nodes in the tree)
     */
    public Object[] getLabels() {
        Object[] ret = new Object[getSize()];
        for (int i : this) {
            ret[i] = this.getLabel(i);
        }
        return ret;
    }


    /**
     * Get the branch point indices for all leaf nodes in this tree.
     * @return the subset of indices that identify leaves/extant species
     */
    public int[] getLeaves() {
        int[] ret = new int[getNLeaves()];
        int cnt = 0;
        for (int idx = 0; idx < bpoints.length; idx ++)
            if (isLeaf(idx))
                ret[cnt ++] = idx;
        return ret;
    }

    /**
     * Determines the leaf states of the tree.
     * A leaf state is true if the branch point is a leaf (i.e., has no children), and false otherwise.
     *
     * @return a boolean array where each element represents the leaf state of the corresponding branch point
     */
    public boolean[] getLeafStates() {
        boolean[] state = new boolean[children.length];
        for (int i = 0; i < children.length; i ++)  // go through all indices
            state[i] = isLeaf(i);
        return state;
    }


    /**
     * Get the branch point indices for all nodes in this tree that are NOT leaves,
     * i.e. ancestors.
     * @return the indices of all ancestor nodes
     */
    public int[] getAncestors() {
        int[] ret = new int[getSize() - getNLeaves()];
        int cnt = 0;
        for (int idx = 0; idx < bpoints.length; idx ++)
            if (!isLeaf(idx))
                ret[cnt ++] = idx;
        return ret;
    }

    /**
     * Retrieve all indices of leaves that are part of the subtree under specified (root) index.
     * @param bpidx specified (root) index
     * @return an array with indices of all leaves
     */
    public int[] getLeaves(int bpidx) {
        Set<Integer> subtree = getSubtreeIndices(bpidx);
        Set<Integer> leaves = new HashSet<>();
        for (int idx : subtree)
            if (isLeaf(idx))
                leaves.add(idx);
        int[] leafarr = new int[leaves.size()];
        int i = 0;
        for (int idx : leaves)
            leafarr[i ++] = idx;
        return leafarr;
    }

    /**
     * Retrieve the index for a specified label (either user specified label or automatically assigned
     * @param label name or ancestor ID
     * @return index of tree node, or -1 if not found
     */
    public int getIndex(Object label) {
        if (label instanceof String) {
            for (int i : this) {
                if (this.getBranchPoint(i).getLabel().toString().equals(label))
                    return i;
                if (this.getLabel(i).toString().equals(label))
                    return i;
            }
            if (((String) label).startsWith("N")) {
                String shortened = ((String) label).substring(1);
                try {
                    Object id = Integer.parseInt(shortened);
                    return getIndex(id);
                } catch (NumberFormatException e) {
                    // nope, not an ancestor ID
                }
            }
        } else {
            for (int i : this) {
                if (this.getLabel(i).equals(label))
                    return i;
            }
        }
        return -1;
    }

    /**
     * Calculate the mean distance between all direct ancestor-child pairs.
     * @return the mean distance between all direct ancestor-child pairs
     */
    public double getMeanDistance() {
        double d = 0;
        try {
            for (int i = 0; i < bpoints.length; i ++)
                d += bpoints[i].getDistance();
        } catch (TreeRuntimeException e) {
            throw new TreeRuntimeException("One or more branchpoints do not have distances assigned");
        }
        return d;
    }

    /**
     * Calculate the mean distance for any leaf to their root
     * @return mean distance
     */
    public double getMeanDistanceToRoot() {
        double sum = 0;
        int[] leaves = getLeaves();
        for (int idx : leaves)
            sum += getDistanceToRoot(idx);
        return(sum / leaves.length);
    }
    /**
     * Adjust distances from extants to root
     * @param target2root
     */
    public void adjustDistances(double target2root) {
        double mu = getMeanDistanceToRoot();
        double multiplier = target2root / mu;
        for (int idx = 1; idx < bpoints.length; idx ++)
            setDistance(idx, getDistance(idx)*multiplier);
    }

    /**
     * Collect distances from all branches in the tree (except for the root, and those which are 0)
     * @return the distances in an array
     */
    public double[] getValidDistances() {
        List<Double> dists = new ArrayList<>();
        for (int idx : this) {
            if (idx == 0)
                continue;
            Double dist = getDistance(idx);
            if (dist > 0)
                dists.add(dist);
        }
        double[] dists_arr = new double[dists.size()];
        int k = 0;
        for (double d : dists)
            dists_arr[k ++] = d;
        return dists_arr;
    }

    /**
     * Fits the root-to-leaf distances to match a target distribution by swapping their order,
     * not changing the composition.
     *
     * This method iteratively adjusts the distances in the tree to better fit the specified target distribution.
     * It performs a specified number of iterations, where in each iteration it samples pairs of leaves and swaps
     * their distances to improve the fit to the target distribution.
     *
     * @param niter the number of iterations to perform
     * @param target the target distribution to fit the distances to
     */
    public void fitDistances(int niter, Distrib target, long seed) {
        Random rand = new Random(seed);
        Double[][] dmat = getDistance2RootMatrix();
        int[][] connected = new int[this.getNLeaves()][];
        for (int i = 0; i < connected.length; i++)
            connected[i] = IdxTree.getConnected(i, 1, dmat);
        for (int i = 0; i < niter; i++) {
            // Step 1: sample a pair of leaves j1 v j2, (TODO: currently "greedy" not stochastic)
            // where j1 is "good" (root distance does not need to be changed) and
            // j2 is "bad" (root distance should be changed, as a consequence of swapping distances in j1)
            double[] leafdists = IdxTree.getLeafDistances(dmat);
            double[][] logodds = new double[leafdists.length][];
            double sumodds = 0;
            int nomj1 = 0, nomj2 = 0;
            for (int j1 = 0; j1 < leafdists.length; j1++) {
                double p1 = target.get(leafdists[j1]);
                int[] conn = connected[j1];
                logodds[j1] = new double[conn.length];
                int cnt = 0;
                for (int j2 : conn) {
                    double p2 = target.get(leafdists[j2]);
                    logodds[j1][cnt] = Math.max(Math.log(p1 / p2), 0);
                    sumodds += logodds[j1][cnt];
                    cnt ++;
                }
            }
            double toss = rand.nextDouble();
            double accum = 0;
            boolean found = false;
            for (int j1 = 0; j1 < logodds.length; j1 ++) {
                for (int cnt = 0; cnt < logodds[j1].length; cnt ++) {
                    int j2 = connected[j1][cnt];
                    double normalised = logodds[j1][cnt] / sumodds;
                    if (toss < (normalised + accum)) {
                        nomj1 = j1;
                        nomj2 = j2;
                        found = true;
                        break;
                    }
                    accum += normalised;
                }
                if (found)
                    break;
            }

            // Step 2: identify what branches b1 v b2 that should be swapped in j1,
            // so that ...
            Set<Integer> shared = new HashSet<>(); // branches shared between j1 and j2
            Set<Integer> unique = new HashSet<>(); // branches unique to j1
            for (int b = 1; b < dmat[nomj1].length; b ++) {
                if (dmat[nomj1][b] != null) {
                    if (dmat[nomj2][b] != null)
                        shared.add(b);
                    else
                        unique.add(b);
                }
            }
            // Step 3: make a choice,
            Integer favb1 = null, favb2 = null;
            // two possibilities...
            if (leafdists[nomj2] < leafdists[nomj1]) { // increase j2 relative j1:
                // pick largest unique to j1
                for (Integer b1 : unique)
                    if (favb1 == null)
                        favb1 = b1;
                    else if (dmat[nomj1][b1] > dmat[nomj1][favb1])
                        favb1 = b1;
                // pick smallest shared between j1 and j2
                for (Integer b2 : shared)
                    if (favb2 == null)
                        favb2 = b2;
                    else if (dmat[nomj1][b2] < dmat[nomj1][favb2])
                        favb2 = b2;
            } else { // decrease j2 relative j1
                // pick smallest unique to j1
                for (Integer b1 : unique)
                    if (favb1 == null)
                        favb1 = b1;
                    else if (dmat[nomj1][b1] < dmat[nomj1][favb1])
                        favb1 = b1;
                // pick largest shared between j1 and j2
                for (Integer b2 : shared)
                    if (favb2 == null)
                        favb2 = b2;
                    else if (dmat[nomj1][b2] > dmat[nomj1][favb2])
                        favb2 = b2;
            }
            // Step 4: make the swap in the distance matrix
            ;
            double x1 = dmat[nomj1][favb1];
            double x2 = dmat[nomj1][favb2];
            for (int leaf = 0; leaf < dmat.length; leaf ++) {
                if (dmat[leaf][favb1] != null)
                    dmat[leaf][favb1] = x2;
                if (dmat[leaf][favb2] != null)
                    dmat[leaf][favb2] = x1;
            }
        }
        setDistance2RootMatrix(dmat);
    }

    /**
     * Calculates the parameters of a Gamma distribution based on the distances between branch points in the tree.
     *
     * This method collects all positive distances between branch points, then calculates the shape and scale parameters
     * of the Gamma distribution that best fits these distances.
     *
     * @return an array containing the alpha and beta parameters of the Gamma distribution
     */
    public double[] getGammaParams() {
        double[] dists_arr = getValidDistances();
        double alpha1 = GammaDistrib.getAlpha(dists_arr);
        double beta1 = GammaDistrib.getBeta(dists_arr, alpha1);
        return new double[] {alpha1, beta1};
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

    /**
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture
     * @return a sample from the mixture
     */
    public static double sampleMixture(GammaDistrib[] mixture, double[] priors) {
        GammaDistrib.Mixture mix = new GammaDistrib.Mixture(mixture, priors);
        //mix.s
        double r = Math.random();
        double sum = 0.0;
        for (int i = 0; i < mixture.length; i++) {
            sum += priors[i];
            if (r < sum)
                return mixture[i].sample();
        }
        return mixture[mixture.length - 1].sample();
    }

    /**
     * @param mixture a mixture of Gamma distributions, a uniform prior is assumed
     * @param r a value to evaluate the mixture at
     * @return distribution of the components that generated the value r
     * */
    public static EnumDistrib getComponentDistribution(GammaDistrib[] mixture, double r) {
        double[] lhoods = new double[mixture.length];
        EnumDistrib ed = new EnumDistrib(new Enumerable(mixture.length));
        for (int i = 0; i < mixture.length; i++)
            lhoods[i] = mixture[i].get(r) + 1e-10;
        ed.set(lhoods);
        return ed;
    }
    /**
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param r a value to evaluate the mixture at
     * @return distribution of the components that generated the value r
     * */
    public static EnumDistrib getComponentDistribution(GammaDistrib[] mixture, double[] priors, double r) {
        if (priors == null)
            return getComponentDistribution(mixture, r);
        double[] lhoods = new double[mixture.length];
        EnumDistrib ed = new EnumDistrib(new Enumerable(mixture.length));
        for (int i = 0; i < mixture.length; i++)
            lhoods[i] = mixture[i].get(r) * priors[i] + 1e-10;
        ed.set(lhoods);
        return ed;
    }

    /**
     * Samples the component that has generated the given value from a mixture of Gamma distributions.
     *
     * This method evaluates the mixture of Gamma distributions at the given value and samples
     * the component that generated the value and returns its index.
     *
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param r a value to evaluate the mixture at
     * @return sample (the index of) the component that has generated the value r
     * */
    public static int sampleComponent(GammaDistrib[] mixture, double[] priors, double r) {
        return (Integer)getComponentDistribution(mixture, priors, r).sample();
    }

    /**
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param r a value to evaluate the mixture at
     * @return sample the component that has generated the value r
     * */
    public static int getMaxComponent(GammaDistrib[] mixture, double[] priors, double r) {
        return (Integer)getComponentDistribution(mixture, priors, r).getMaxIndex();
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     *
     * This method samples the components of a mixture of Gamma distributions based on the given rates, assuming uniform priors.
     * It performs a specified number of samples and returns an estimated distribution of the components.
     *
     * @param mixture a mixture of Gamma distributions
     * @param rates the values to evaluate the mixture for
     * @param nsamples number of samples from which the estimate is based
     * @return an estimate of the distribution of the components that generated the rates
     * */
    public static EnumDistrib estComponentDistribution(GammaDistrib[] mixture, double[] rates, int nsamples) {
        return estComponentDistribution(mixture, null, rates, nsamples);
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     *
     * This method samples the components of a mixture of Gamma distributions based on the given rates and priors.
     * It performs a specified number of samples and returns an estimated distribution of the components.
     *
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param rates the values to evaluate the mixture for
     * @param nsamples number of samples from which the estimate is based
     * @return an estimate of the distribution of the components that generated the rates
     * */
    public static EnumDistrib estComponentDistribution(GammaDistrib[] mixture, double[] priors, double[] rates, int nsamples) {
        if (rates == null || rates.length == 0)
            throw new RuntimeException("Invalid rates: cannot be empty");
        int draws = 0;
        int[] samples = new int[mixture.length];
        while (draws < nsamples) {
            for (int i = 0; i < rates.length; i++) {
                int comp = sampleComponent(mixture, priors, rates[i]);
                samples[comp] += 1;
                draws += 1;
            }
        }
        double[] p = new double[mixture.length];
        for (int k = 0; k < samples.length; k++)
            p[k] = (double)samples[k] / (double)draws;
        EnumDistrib ed = new EnumDistrib(new Enumerable(mixture.length));
        ed.set(p);
        return ed;
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     * This method uses expectation maximization to estimate the distribution of the components that generated the rates.
     *
     * @param ncomponents number of components in the mixture
     * @return a mixture of Gamma distributions fitted to the distances in the tree; null if EM fails
     */
    public static GammaDistrib.Mixture getGammaMixture(IdxTree tree, int ncomponents) {
        return getGammaMixture(tree, ncomponents, System.currentTimeMillis());
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     * This method uses expectation maximization to estimate the distribution of the components that generated the rates.
     *
     * @param ncomponents number of components in the mixture
     * @return a mixture of Gamma distributions fitted to the distances in the tree; null if EM fails
     */
    public static GammaDistrib.Mixture getGammaMixture(IdxTree tree, int ncomponents, long seed) {
        double[] data = tree.getValidDistances();
        GammaDistrib.Mixture gdm = GammaDistrib.Mixture.fitMLE(data, ncomponents, seed);
        return gdm;
    }

    public static IdxTree generateTreeFromMixture(IdxTree tree1,  int NCOMP, long SEED, int NITER) {
        // First we process an already loaded tree or synthesise a new tree
        GaussianDistrib gds1 = IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix());

        // Now we are estimating a mixture of Gamma distributions from the first/source tree (loaded or synthesised)
        GammaDistrib.Mixture mixture = getGammaMixture(tree1, NCOMP, SEED);

        // 2. Generate a new tree based on the above mixture of Gamma distributions
        //   a) assume that branch lengths are uniformly distributed across topology
        Tree tree2 = Tree.Random(tree1.getNLeaves(), mixture, 2, 2, SEED);

        //   b) adjust the placement of distances to better fit the original distribution
        tree2.fitDistances(NITER, gds1, SEED + 202);

        return tree2;
    }

    public static IdxTree generateTreeFromDistrib(RateModel distmodel, int NLEAVES, long SEED) {
        return generateTreeFromDistrib(distmodel, null, NLEAVES, SEED, 0);
    }

    public static IdxTree shuffleWithLeaf2RootDistrib(IdxTree tree1, Distrib dist, long SEED, int NITER)  {
        try {
            GaussianDistrib gds1 = (GaussianDistrib) dist;
            //   b) adjust the placement of distances to better fit the original distribution
            tree1.fitDistances(NITER, gds1, SEED + 202);
            return tree1;
        } catch (ClassCastException e) {
            throw new RuntimeException("Leaf-2-root distance distribution is not a Gaussian");
        }
    }

    public static IdxTree generateTreeFromDistrib(RateModel distmodel, Distrib leaf2root, int NLEAVES, long SEED, int NITER) {
        // 2. Generate a new tree based on the above mixture of Gamma distributions
        //   a) assume that branch lengths are uniformly distributed across topology
        Tree tree = null;
        tree = Tree.Random(NLEAVES, distmodel, 2,2, SEED);
        if (leaf2root == null)
            return tree;
        try {
            return shuffleWithLeaf2RootDistrib(tree, leaf2root, SEED, NITER);
        } catch (ClassCastException e) {
            throw new RuntimeException("Leaf-2-root distance distribution is not a Gaussian");
        }
    }

    public static void main(String[] args) {
        long SEED = System.currentTimeMillis(); // random seed
        int NCOMP = 3; // number of components in Gamma mixture to model branch distances
        int NITER = 100; // max number of iterations to estimate mixture distribution (EM)
        IdxTree tree1 = null;
        for (int argc = 0; argc < args.length; argc++) {
            switch (args[argc]) {
                case "--ncomp":
                    NCOMP = Integer.parseInt(args[++argc]);
                    break;
                case "--seed":
                    SEED = Long.parseLong(args[++argc]);
                    break;
                case "--niter":
                    NITER = Integer.parseInt(args[++argc]);
                    break;
                case "--nwk":
                    try {
                        tree1 = Newick.load(args[++argc]);
                    } catch (IOException e) {
                        System.err.println("Error loading " + args[argc] + ": " + e.getMessage());
                        System.exit(1);
                    }
                default:
                    break;
            }
        }
        if (tree1 != null) {
            // What this function does: IdxTree ntree = generateTreeFromMixture(tree, NCOMP, SEED, NITER);
            // Namely... process an already loaded tree or synthesise a new tree
            GaussianDistrib gds1 = IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix());
            // Then estimate a mixture of Gamma distributions from the first/source tree (loaded or synthesised)
            RateModel distmodel = getGammaMixture(tree1, NCOMP, SEED);
            // Then generate a new tree based on the above mixture of Gamma distributions
            Tree tree2 = Tree.Random(tree1.getNLeaves(), distmodel, 2, 2, SEED);
            // Finally, adjust the placement of distances to better fit the original distribution
            tree2.fitDistances(NITER, gds1, SEED + 202);
            System.out.println("--dist-distrib " + distmodel.getTrAVIS() + "\n--leaf2root-distrib " + gds1.getTrAVIS());
        }
    }
}
