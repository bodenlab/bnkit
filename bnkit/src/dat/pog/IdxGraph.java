package dat.pog;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Indexed graph data structure.
 * Based on "adjacency matrix" but uses efficient representation using BitSet.
 * The idea is to keep it thread-safe.
 */
public class IdxGraph {
    private Node[] nodes = null; // indexed by internal node ID i in {0..N}, node does not exist if nodes[i] = null
    final private boolean undirected;
    private BitSet[] edgesForward = null;
    private BitSet[] edgesBackward = null;
    private Map<Integer, Edge> edges = null;
    private BitSet startNodes = null;
    private BitSet endNodes = null;

    public static IdxGraph makeTerminatedDAG(int nNodes) {
        return new IdxGraph(nNodes, false, true);
    }

    public static IdxGraph makeUndirectedGraph(int nNodes) {
        return new IdxGraph(nNodes, true, false);
    }

    /**
     * Create graph.
     * @param nNodes number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
     * @param undirected if true, edges are not directed, else
     * @param terminated if true, a virtual start and a virtual end nodes are created to enable edges to mark multiple entry and exit nodes
     */
    public IdxGraph(int nNodes, boolean undirected, boolean terminated) {
        this.undirected = undirected;
        this.nodes = new Node[nNodes];
        this.edgesForward = new BitSet[nNodes];
        if (!undirected)
            this.edgesBackward = new BitSet[nNodes];
        this.edges = new HashMap<>();
        if (terminated) {
            this.startNodes = new BitSet(nNodes);
            this.endNodes = new BitSet(nNodes);
        }
    }

    /**
     * Maximum number of nodes.
     * @return the maximum number of nodes that this graph can include.
     */
    public int size() {
        return this.nodes.length;
    }

    /**
     * Return whether graph is terminated
     * @return true if the graph has virtual termination nodes, allowing edges to indicate multiple start and end nodes
     */
    public boolean isTerminated() {
        return (startNodes != null && endNodes != null);
    }

    /**
     * Check if the graph is directed
     * @return true if directed, false if undirected
     */
    public boolean isDirected() {
        return !undirected;
    }

    /**
     * Check if index is a valid index for a node; does not include virtual termination nodes and
     * does not indicate if the index is indeed occupied by a node
     * @param idx index
     * @return true if index a valid node index, else false
     */
    public boolean isIndex(int idx) {
        return (idx >= 0 && idx < this.size());
    }

    /**
     * Check if index is occupied by a node; does not include virtual termination nodes
     * @param idx index
     * @return true if index occupied by a node, else false
     */
    public boolean isNode(int idx) {
        if (isIndex(idx))
            return (nodes[idx] != null);
        return false;
    }

    /**
     * Check if the indices map to an existing edge.
     * @param from the source node index (note: -1 for terminal edge)
     * @param to the destination node index (note: N for terminal edge, where N is the size of the possible/valid node-set)
     * @return true if edge exists, else false
     */
    public boolean isEdge(int from, int to) {
        if (isNode(from) && isNode(to)) {
            if (undirected) {
                if (to < from) {
                    int bigger = from;
                    from = to;
                    to = bigger;
                }
            }
            return edgesForward[from].get(to); // note edgesBackward is only non-redundant if directed
        } else if (from == -1 && isNode(to) && isTerminated()) {
            return startNodes.get(to);
        } else if (isNode(from) && to == this.size() && isTerminated()) {
            return endNodes.get(from);
        }
        return false;
    }

    /**
     * Determine the canonical index for an edge, including terminal edges.
     * @param from the source node index (note: -1 for terminal edge)
     * @param to the destination node index (note: N for terminal edge, where N is the size of the possible/valid node-set)
     * @return a unique index for the edge
     */
    public Integer getEdgeIndex(int from, int to) {
        if (isNode(from) && isNode(to)) {
            if (undirected && to < from) {
                int bigger = from;
                from = to;
                to = bigger;
            }
            return Integer.valueOf(from * size() + to);
        } else if (from == -1 && isNode(to) && isTerminated()) {
            return -Integer.valueOf(to);
        } else if (isNode(from) && to == this.size() && isTerminated()) {
            return Integer.valueOf(from * size() + to);
        }
        throw new RuntimeException("Cannot be an edge: " + from + " -- " + to);
    }

    /**
     * Retrieve the node instance for a specified index
     * @param idx index of node
     * @return the node instance
     */
    public Node getNode(int idx) {
        if (isIndex(idx))
            return nodes[idx];
        else
            throw new InvalidIndexRuntimeException("Index outside bounds: " + idx);
    }

    /**
     * Modify the graph by adding a node at a specified index.
     * Note: does not connect the node.
     * @param nid node index, which must be valid, i.e. between 0 and N - 1, where N is the size of the possible/valid node-set
     */
    public synchronized void addNode(int nid) {
        addNode(nid, new DefaultNode());
    }

    /**
     * Modify the graph by adding an instance of a node at a specified index.
     * Note: does not connect the node.
     * @param nid node index, which must be valid, i.e. between 0 and N - 1, where N is the size of the possible/valid node-set
     * @param node the instance of the node, which belongs to class Node
     */
    public synchronized void addNode(int nid, Node node) {
        if (isIndex(nid)) {
            this.nodes[nid] = node;
            this.edgesForward[nid] = new BitSet(size());
            if (!undirected)
                this.edgesBackward[nid] = new BitSet(size());
        } else {
            throw new InvalidIndexRuntimeException("Index outside bounds: " + nid);
        }
    }

    /**
     * Modify the graph by removing a node at a specified index.
     * Note: will disconnect the node (i.e. remove any exiting edge involving the node)
     * @param nid node index
     */
    public synchronized void removeNode(int nid) {
        if (isIndex(nid)) {
            this.nodes[nid] = null;
            this.edgesForward[nid] = null;
            if (!undirected)
                this.edgesBackward[nid] = null;
        } else {
            throw new InvalidIndexRuntimeException("Index outside bounds: " + nid);
        }
    }

    /**
     * Retrieve the edge instance for a pair of node indices
     * @param from source node index
     * @param to target node index
     * @return edge instance if exists, else null
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public Edge getEdge(int from, int to) {
        if (isNode(from) && isNode(to)) {
            return this.edges.get(getEdgeIndex(from, to));
        } else if (from == -1 && isNode(to) && isTerminated()) {
            return this.edges.get(getEdgeIndex(from, to));
        } else if (isNode(from) && to == this.size() && isTerminated()) {
            return this.edges.get(getEdgeIndex(from, to));
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve edge between non-existent node/s: " + from + " or " + to);
        }
    }

    /**
     * Modify the graph by adding an edge between two existing nodes.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     */
    public synchronized void addEdge(int from, int to) {
        addEdge(from, to, new DefaultEdge());
    }

    /**
     * Modify the graph by adding an instance of an edge between two existing nodes.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @param edge the instance of the edge (optional)
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized void addEdge(int from, int to, Edge edge) {
        if (isNode(from) && isNode(to)) {
            if (!undirected) { // directed, we set both forward and backward pointers accordingly
                this.edgesBackward[to].set(from);
                this.edgesForward[from].set(to);
            } else { // un-directed, which means we record the same edge in BOTH directions
                this.edgesForward[to].set(from);
                this.edgesForward[from].set(to);
            }
            this.edges.put(getEdgeIndex(from, to), edge);
        } else if (from == -1 && isNode(to) && isTerminated()) {
            this.startNodes.set(to);
            this.edges.put(getEdgeIndex(from, to), edge);
        } else if (isNode(from) && to == this.size() && isTerminated()) {
            this.endNodes.set(from);
            this.edges.put(getEdgeIndex(from, to), edge);
        } else {
            throw new InvalidIndexRuntimeException("Cannot add edge \"" + edge.toString() + "\" between non-existent node/s: " + from + " or " + to);
        }
    }

    /**
     * Modify the graph by adding an instance of an edge between two existing nodes.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized void removeEdge(int from, int to) {
        if (isNode(from) && isNode(to)) {
            if (!undirected)
                this.edgesBackward[to].set(from, false);
            else
                this.edgesForward[to].set(from, false);
            this.edgesForward[from].set(to, false);
            this.edges.remove(getEdgeIndex(from, to));
        } else if (from == -1 && isNode(to) && isTerminated()) {
            this.startNodes.set(to, false);
            this.edges.remove(getEdgeIndex(from, to));
        } else if (isNode(from) && to == this.size() && isTerminated()) {
            this.endNodes.set(from, false);
            this.edges.remove(getEdgeIndex(from, to));
        } else {
            throw new InvalidIndexRuntimeException("Cannot remove edge between non-existent node/s: " + from + " or " + to);
        }
    }

    /**
     * Get indices of all nodes that can be reached from the given node index.
     * If graph is undirected, then it does not matter what the source node was when added.
     * If graph is directed, then the default is to only return the edges looking forward (from the source node when added).
     * @param idx
     * @return
     */
    public int[] getNodeIndices(int idx) {
        return getNodeIndices(idx, true);
    }

    /**
     * Get indices of all nodes that can be reached relative the given node index
     * @param idx
     * @param direction if true, forward, else backward
     * @return
     */
    public int[] getNodeIndices(int idx, boolean direction) {
        int count = 0;
        if (isNode(idx)) {
            BitSet bs = (direction || undirected) ? edgesForward[idx] : edgesBackward[idx];
            int[] indices = new int[bs.cardinality()];
            for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else if (idx == -1 && (direction || undirected)) {
            int[] indices = new int[startNodes.cardinality()];
            for (int i = startNodes.nextSetBit(0); i >= 0; i = startNodes.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else if (idx == size() && (!direction || undirected)) {
            int[] indices = new int[endNodes.cardinality()];
            for (int i = endNodes.nextSetBit(0); i >= 0; i = endNodes.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve edges from non-existent/invalid node: " + idx);
        }
    }

}

/**
 *
 */
class InvalidIndexRuntimeException extends RuntimeException {
    public InvalidIndexRuntimeException(String errmsg) {
        super(errmsg);
    }
}
/**
 * Default, place-holder implementation of a node
 */
class DefaultNode implements Node {
}

/**
 * Default, place-holder implementation of an edge
 */
class DefaultEdge implements Edge {
}