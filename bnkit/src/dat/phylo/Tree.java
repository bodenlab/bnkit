/*
 * binfkit -- software for bioinformatics
 * Copyright (C) 2014  M. Boden et al.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package dat.phylo;

import dat.EnumSeq;
import dat.EnumSeq.Alignment;
import dat.Enumerable;
import dat.file.Newick;

import java.io.IOException;
import java.util.*;

/**
 * Class to represent a phylogenetic tree, refactored from old PhyloTree (now deprecated).
 * Rooted, multifurcating tree for representing phylogenetic relationships.
 * Functionality includes labeling and traversing branchpoints; reading and writing to Newick format;
 * Programmers should note that almost all functionality is implemented through recursion.
 *
 * The current design separates the tree topology (with branch points and their labels, represented by this class)
 * from instantiations (values assigned to tips and internal branch points, represented by TreeInstance,
 * of which several can be based on the same topology).
 *
 * @author mikael
 */
public class Tree implements TreeTopology {

    private final BranchPoint root; // the root of the tree, all descendants linked by pointer in BranchPoint class
    private final Map<Object, Integer> index; // map label to index
    private final BranchPoint[] bpoints; // store branch points in order of their index
    private final int[] parent; // [child idx] = parent idx, parent idx is -1 if no parent (i.e. root)
    private final int[][] children; // [parent idx] = [child 1][child 2]...[child n] where n is #children for parent

    /**
     * Constructor for tree from a root with all nodes connected off that.
     * Use factory methods to construct trees by assembling BranchPoints relative a root.
     */
    public Tree(BranchPoint root) {
        this.root = root;
        // assume all branch points are OK and store them
        List<Tree.BranchPoint> branchPoints = root.getSubtree();
        bpoints = new BranchPoint[branchPoints.size()];
        branchPoints.toArray(bpoints);
        // Check branch points, and index them
        index = new HashMap<>();
        parent = new int[bpoints.length];
        children = new int[bpoints.length][];
        int maxchild = 0;
        Integer count = 0;
        for (int i = 0; i < bpoints.length; i ++) {
            BranchPoint bp = branchPoints.get(i);
            children[i] = new int[bp.getChildren().size()];
            if (bp.getID() == null) // if leaf, and label is not set, OR if ancestor and ancestor ID is not set
                bp.setAncestor(count ++);
            if (index.containsKey(bp.getID()))
                throw new TreeRuntimeException("Branch points form cycle; duplicate found: label is " + bp.getLabel() + ", identifier is " + bp.getID());
            else {
                index.put(bp.getID(), i);
                if (i > 0) { // not root
                    BranchPoint ancestor = bp.getParent();
                    Integer paridx = index.get(ancestor.getID());
                    if (paridx == null)
                        throw new TreeRuntimeException("Branch points are out of order, so parent indexing failed at: " + bp.getID());
                    parent[i] = paridx;
                } else
                    parent[0] = -1;
            }
        }
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
     * Retrieve the branch point for a given index
     * @param idx
     * @return
     */
    public BranchPoint getBranchPoint(int idx) {
        if (idx >= 0 && idx < getSize())
            return bpoints[idx];
        throw new TreeRuntimeException("Invalid branch point index: " + idx);
    }

    public static Tree load(String filename, String format) throws IOException {
        if (format.equalsIgnoreCase("newick") || format.equalsIgnoreCase("nwk"))
            return Newick.load(filename);
        else
            throw new IOException("Unknown format: " + format);
    }

    public void save(String filename, String format) throws IOException {
        if (format.equalsIgnoreCase("newick") || format.equalsIgnoreCase("nwk"))
            Newick.save(this, filename, Newick.MODE_DEFAULT);
        else if (format.equalsIgnoreCase("ancestor") || format.equalsIgnoreCase("anwk"))
            Newick.save(this, filename, Newick.MODE_ANCESTOR);
        else
            throw new IOException("Unknown format: " + format);
    }
    /**
     * String representation in the Newick format.
     *
     * @return string representation of tree
     */
    @Override
    public String toString() {
        return root.toString();
    }

    /**
     * Get root node of tree.
     *
     * @return the root of the tree
     */
    public BranchPoint getRoot() {
        return root;
    }

    /**
     * Find the node with the specified label.
     *
     * @param content label or label
     * @return matching node, or null if not found
     */
    public BranchPoint find(Object content) {
        return root.find(content);
    }

    public void setInternalLabels() {
        root.setInternalLabels(0);
    }

    public static void main(String[] args) {
        Tree phyloTree = Newick.parse("((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2");
        System.out.println(phyloTree.root);
        phyloTree.setInternalLabels();
        System.out.println(phyloTree.root);
        try {
            Tree edge1 = Newick.load("/Users/mikael/simhome/ASR/edge1.nwk");
            System.out.println(edge1);
            Alignment aln = new Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/gap1.aln", Enumerable.aacid));
            TreeInstance.Parsimony tip = edge1.getInstance(aln.getNames(), aln.getGapColumn(1)).new Parsimony();
            System.out.println(tip);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Class for branch points that make up tree.
     * Note recursive definition. Supports any branching factor.
     */
    public static class BranchPoint {
        private List<BranchPoint> children = new ArrayList<>(); // the children of this branch point
        private Object label = null;               // arbitrary label of branch point
        private Integer ancestor = null;            // ancestor ID (a count which conventionally starts with 0 at root)
        private Double dist = null;         // optional distance (from this node to its parent)
        private BranchPoint parent = null;         // link to parent

        /**
         * Construct a branch point
         */
        public BranchPoint() {
        }

        /**
         * Construct node from label/label.
         *
         * @param label label/label
         */
        public BranchPoint(String label) {
            this.label = label;
        }

        /**
         * Construct node from label/label, and parent at distance.
         *
         * @param label  label/label
         * @param parent parent node
         * @param dist   distance to parent from this node
         */
        public BranchPoint(String label, BranchPoint parent, Double dist) {
            this.label = label;
            this.parent = parent;
            this.dist = dist;
        }

        /**
         * Construct a branch point for a parent based on existing children
         *
         * @param label
         * @param children
         */
        public BranchPoint(String label, BranchPoint... children) {
            this.label = label;
            for (BranchPoint child : children)
                this.children.add(child);
        }

        public boolean isLeaf() {
            if (children == null)
                return true;
            return children.isEmpty();
        }

        /**
         * Retrieve identifier for branch point.
         * This is the label (leaf) or ancestor counter (internal).
         * Should be unique to be used in a tree.
         *
         * @return identifier
         */
        public Object getID() {
            return (isLeaf() ? label : ancestor);
        }

        public Object getLabel() {
            return label;
        }

        public void setLabel(String label) {
            this.label = label;
        }

        public void setAncestor(Integer id) {
            this.ancestor = id;
        }

        public Integer getAncestor() {
            return ancestor;
        }

        public void setParent(BranchPoint parent) {
            this.parent = parent;
        }

        public BranchPoint getParent() {
            return parent;
        }

        /**
         * String representation of the node and its children (recursively) that uses the Newick format.
         *
         * @return string representation
         */
        public String toString() {
            StringBuilder sb = new StringBuilder();
            String dstr = null;
            int nchildren = children.size();
            int cnt = 0;
            for (BranchPoint child : children) {
                sb.append(child.toString());
                if (++cnt < nchildren)
                    sb.append(",");
            }
            if (dist != null)
                dstr = ":" + dist.toString();
            if (nchildren < 1)
                return label.toString() + ((dist != null) ? (dstr) : (""));
            else
                return "(" + sb.toString() + ")" + label.toString() + ((dist != null) ? (dstr) : (""));
        }

        /**
         * Add child to node.
         *
         * @param child
         */
        public void addChild(BranchPoint child) {
            children.add(child);
        }


        /**
         * Remove child from node.
         *
         * @param child
         */
        public void removeChild(BranchPoint child) {
            children.remove(child);
        }


        /**
         * Retrieve all the children of the node.
         *
         * @return
         */
        public List<BranchPoint> getChildren() {
            return children;
        }

        public List<BranchPoint> getSubtree() {
            if (this.isLeaf())
                return Collections.singletonList(this);
            List<BranchPoint> below = new ArrayList<>();
            below.add(this);
            for (BranchPoint child : getChildren()) {
                below.addAll(child.getSubtree());
            }
            return below;
        }

        /**
         * Retrieve and return all direct and indirect descendants of this branch point (excluding itself)
         *
         * @return
         */
        public List<BranchPoint> getDescendants() {
            if (this.isLeaf())
                return Collections.EMPTY_LIST;
            List<BranchPoint> below = new ArrayList<>();
            for (BranchPoint child : getChildren()) {
                below.add(child);
                below.addAll(child.getDescendants());
            }
            return below;
        }

        /**
         * Find node by label/label.
         * Searches the tree recursively using the current node as root.
         *
         * @param label
         * @return the node that contains the specified label, or null if not found
         */
        public BranchPoint find(Object label) {
            if (this.label.equals(label))
                return this;
            else {
                for (BranchPoint child : children) {
                    BranchPoint branchPoint = child.find(label);
                    if (branchPoint != null)
                        return branchPoint;
                }
                return null;
            }
        }

        private int setInternalLabels(int count) {
            if (isLeaf()) {     // no need to label anything, or to increase count
                return count;
            } else {            // this is an ancestor/internal node, so label it and investigate children
                this.setLabel("N" + count++);
                for (BranchPoint child : children) {
                    count = child.setInternalLabels(count);
                }
                return count;
            }
        }

        /**
         * Set the distance for this node (from this node to its parent)
         *
         * @param dist the distance
         */
        public void setDistance(double dist) {
            this.dist = dist;
        }

        /**
         * Retrieve the distance of this node (from this node to its parent)
         *
         * @return the distance
         */
        public double getDistance() {
            if (this.dist == null)
                throw new RuntimeException("Node " + this + " with content " + label + " does not have a distance");
            return this.dist;
        }
    }

    public class TreeRuntimeException extends RuntimeException {

        public TreeRuntimeException(String errmsg) {
            super(errmsg);
        }
    }

}

