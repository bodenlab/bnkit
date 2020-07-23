package dat.pog;

import java.util.*;

/**
 * Indexed graph data structure.
 * Based on "adjacency matrix" but uses efficient representation using BitSet.
 * Beware: The plan is to keep this class thread-safe.
 */
public class IdxEdgeGraph<E extends Edge> extends IdxGraph {
    private Map<Integer, E> edges = null;

    /**
     * Create graph based on indexed nodes, edges represented by pairs of indices. Node and edge data can be attached.
     * The data structure allows for multiple "entry" and "exit" nodes, as would be required in some graphs with node order (partial order graphs for instance).
     * @param nNodes maximum number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
     * @param undirected if true, edges are "undirected", else "directed"
     * @param terminated if true, a virtual start and a virtual end node are created to enable edges to mark multiple entry and exit nodes
     */
    public IdxEdgeGraph(int nNodes, boolean undirected, boolean terminated) {
        super(nNodes, undirected, terminated);
        this.edges = new HashMap<>();
    }

    /**
     * Retrieve the edge instance for a pair of node indices
     * @param from source node index (use -1 for terminal start edge)
     * @param to target node index (use N for a terminal end edge, where N is the number of possible/valid nodes)
     * @return edge instance if exists, else null
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public E getEdge(int from, int to) {
        if (isNode(from) && isNode(to)) {
            return this.edges.get(getEdgeIndex(from, to));
        } else if (from == -1 && isNode(to) && isTerminated()) {
            return this.edges.get(getEdgeIndex(from, to));
        } else if (isNode(from) && to == this.maxsize() && isTerminated()) {
            return this.edges.get(getEdgeIndex(from, to));
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve edge between non-existent node/s: " + from + " or " + to);
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
        super.removeEdge(from , to);
        this.edges.remove(getEdgeIndex(from, to));
    }

    /**
     * Modify the graph by adding an instance of an edge between two existing nodes.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @param edge the instance of the edge (optional)
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized boolean addEdge(int from, int to, E edge) {
        if (addEdge(from, to)) {
            this.edges.put(getEdgeIndex(from, to), edge);
            return true;
        } else
            return false;
    }

    public synchronized boolean addTerminalEdge(int from, E edge) {
        return addEdge(from, maxsize(), edge);
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
