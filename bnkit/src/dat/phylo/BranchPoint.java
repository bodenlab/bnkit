package dat.phylo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Class for branch points that make up tree.
 * Note recursive definition. Supports any branching factor.
 */
public class BranchPoint {
    private List<BranchPoint> children = new ArrayList<>(); // the children of this branch point
    private Object label = null;                // arbitrary label of branch point
    private Integer ancestor = null;            // ancestor ID (a count that conventionally starts with 0 at root)
    private Double dist = null;                 // optional distance (from this node to its parent)
    private BranchPoint parent = null;          // link to parent

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

    public int setInternalLabels(int count) {
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