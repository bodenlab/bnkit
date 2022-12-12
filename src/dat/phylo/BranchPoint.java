package dat.phylo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Class for branch points that make up tree.
 * Note recursive definition. Supports any branching factor.
 * Should be used with great care, since pointers may be used in parallel with parent-child indices in index trees.
 */
public class BranchPoint {
    private List<BranchPoint> children = new ArrayList<>(); // the children of this branch point
    private Object label = null;                // arbitrary label of branch point
    private Integer ancestor = null;            // ancestor ID (a count that conventionally starts with 0 at root)
    private Double dist = null;                 // optional distance (from this node to its parent)
    private BranchPoint parent = null;          // link to parent

    /**
     * Construct an empty or disconnected branch point
     */
    private BranchPoint() {
    }

    /**
     * Construct a disconnected branch point but with label.
     *
     * @param label label
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
     * @// TODO: 21/9/20 if parent is specified, link it to this branch point as child or will this create issues?
     */
    public BranchPoint(String label, BranchPoint parent, Double dist) {
        this.label = label;
        this.parent = parent;
        this.dist = dist;
    }

    /**
     * Construct a branch point for a parent based on existing children.
     * Note that this will overwrite the parent of those existing children.
     *
     * @param label
     * @param children
     */
    public BranchPoint(String label, BranchPoint... children) {
        this.label = label;
        for (BranchPoint child : children) {
            this.children.add(child);
            child.setParent(this);
        }
    }

    /**
     * Check if this branch point is a leaf, i.e. has no children of its own
     * @return
     */
    public boolean isLeaf() {
        if (children == null)
            return true;
        return children.isEmpty();
    }

    /**
     * Check if this branch point is a parent, i.e. has children of its own
     * @return
     */
    public boolean isParent() {
        return !children.isEmpty();
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

    /**
     * Retrieve the label of the branch point
     * @return the label
     */
    public Object getLabel() {
        return label;
    }

    /**
     * Set the label of the branch point
     * @param label the label
     */
    public void setLabel(String label) {
        this.label = label;
    }

    /**
     * Set the ancestor identifier of the branch point.
     * Would conventionally be an integer increment (starting with 0 at root)
     * By default the branch point has none (i.e. null) which should/could indicate a leaf.
     * @param id
     */
    public void setAncestor(Integer id) {
        this.ancestor = id;
    }

    /**
     * Get the ancestor identifier.
     * Would conventionally be an integer increment (starting with 0 at root)
     * By default the branch point has none (i.e. null) which should/could indicate a leaf.
     * @return ancestor identifier
     */
    public Integer getAncestor() {
        return ancestor;
    }

    /**
     * Set the parent branch point pointer of the present branch point
     * Note that this does not automatically assign the present node as a child of that specified parent
     * @param parent
     */
    public void setParent(BranchPoint parent) {
        this.parent = parent;
    }

    /**
     * Get the parent branch point of the present branch point
     * @return
     */
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
     * Note that this does not automatically assign a parent to that child
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

    /**
     * Retrieve all the branch points at and below the current branch point
     * @return a list of branch points (including the present)
     */
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
     * Find node by label.
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

    /**
     * Set the labels for all internal, i.e. non-leaf branch points,
     * by the convention "N" then the number incrementing from given count for root,
     * by depth-first search.
     * @param count the start count, e.g. 0 for the ultimate root
     * @return
     */
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
            throw new TreeRuntimeException("Node " + this + " with content " + label + " does not have a distance");
        return this.dist;
    }
}