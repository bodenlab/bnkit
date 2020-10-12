package dat.pog;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class IdxGraph {
    protected final boolean directed;
    public String nodeDOT = "style=\"rounded,filled\", shape=box, fixedsize=true";
    public String edgeDOT = ""; //""style=\"bold\"";
    protected Node[] nodes;    // indexed by internal node ID i in {0..N}, node i does not exist if nodes[i] = null
    protected int nNodes;
    protected BitSet[] edgesForward;
    protected BitSet[] edgesBackward = null;
    protected BitSet startNodes = null;
    protected BitSet endNodes = null;
    protected String name = null;

    public IdxGraph(int nNodes, boolean undirected, boolean terminated) {
        this.nodes = new Node[nNodes];
        //for (int i = 0; i < nNodes; i ++) this.nodes.add(null); // set all node indices to indicate "no node"
        this.nNodes = nNodes;
        this.directed = !undirected;
        this.edgesForward = new BitSet[nNodes];
        if (directed)
            this.edgesBackward = new BitSet[nNodes];
        if (terminated) {
            this.startNodes = new BitSet(nNodes);
            this.endNodes = new BitSet(nNodes);
        }
    }

    public static IdxGraph makeTerminatedDAG(int nNodes) {
        return new IdxEdgeGraph(nNodes, false, true);
    }

    public static IdxGraph makeUndirectedGraph(int nNodes) {
        return new IdxEdgeGraph(nNodes, true, false);
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        if (name == null)
            return "";
        return name;
    }

    /**
     * Maximum number of nodes.
     * @return the maximum number of nodes that this graph can include.
     */
    public int maxsize() {
        return this.nNodes;
    }

    /**
     * Determine actual number of nodes that are defined for graph
     * @return number of nodes
     */
    public int size() {
        int cnt = 0;
        for (int i = 0; i < this.nNodes; i ++)
            cnt += this.nodes[i] != null ? 1 : 0;
        return cnt;
    }

    public synchronized int getFreeIndex() {
        for (int i = 0; i < this.nNodes; i ++)
            if (this.nodes[i] == null)
                return i;
        throw new RuntimeException("There are no free indices in graph");
    }

    /**
     * Determine whether graph is terminated
     * @return true if the graph has virtual termination nodes, allowing edges to indicate multiple start and end nodes
     */
    public boolean isTerminated() {
        return (startNodes != null && endNodes != null);
    }

    /**
     * Determine whether the graph is directed
     * @return true if directed, false if undirected
     */
    public boolean isDirected() {
        return directed;
    }

    /**
     * Check if index is a valid index for a node; does not include virtual termination nodes and
     * does not indicate if the index is indeed occupied by a node
     * @param idx index
     * @return true if specified index is a valid node index, else false
     */
    public boolean isIndex(int idx) {
        return (idx >= 0 && idx < this.maxsize());
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
            if (!isDirected()) {
                if (to < from) {
                    int bigger = from;
                    from = to;
                    to = bigger;
                }
            }
            return edgesForward[from].get(to); // note edgesBackward is only non-redundant if directed
        } else if (from == -1 && isNode(to) && isTerminated()) {
            return startNodes.get(to);
        } else if (isNode(from) && to == this.maxsize() && isTerminated()) {
            return endNodes.get(from);
        }
        return false;
    }

    public int getEdgeCount() {
        int n = 0;
        for (int i = 0; i < maxsize(); i ++) {
            if (isNode(i))
                n += edgesForward[i].cardinality();
        }
        if (!isDirected())
            n /= 2;
        if (isTerminated()) {
            n += startNodes.cardinality();
            n += endNodes.cardinality();
        }
        return n;
    }

    /**
     * Determine if node is a valid entry point in a terminated graph.
     * @param idx node index
     * @return true if index is a valid start node, else false
     */
    public boolean isStartNode(int idx) {
        if (isTerminated())
            return startNodes.get(idx);
        return false;
    }

    /**
     * Determine if node is a valid exit point in a terminated graph.
     * @param idx node index
     * @return true if index is a valid exit node, else false
     */
    public boolean isEndNode(int idx) {
        if (isTerminated())
            if (isNode(idx))
                return endNodes.get(idx);
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
            if (!isDirected() && to < from) {
                int bigger = from;
                from = to;
                to = bigger;
            }
            return from * maxsize() + to;
        } else if (from == -1 && isNode(to) && isTerminated())
            return -to;
        else if (isNode(from) && to == this.maxsize() && isTerminated())
            return from * maxsize() + to;
        throw new InvalidIndexRuntimeException("Cannot be an edge: " + from + " -- " + to);
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
     * Modify the graph by adding an instance of a node at a specified index.
     * Note: does not connect the node.
     * @param nid node index, which must be valid, i.e. between 0 and N - 1, where N is the size of the possible/valid node-set
     * @param node the instance of the node, which belongs to class Node
     */
    public synchronized int addNode(int nid, Node node) {
        if (isIndex(nid)) {
            this.nodes[nid] = node;
            this.edgesForward[nid] = new BitSet(maxsize());
            if (isDirected())
                this.edgesBackward[nid] = new BitSet(maxsize());
            return nid;
        } else {
            throw new InvalidIndexRuntimeException("Index outside bounds: " + nid);
        }
    }

    public synchronized int addNode(Node node) {
        try {
            int idx = this.getFreeIndex();
            return addNode(idx, node);
        } catch (RuntimeException e) {
            throw new RuntimeException("Cannot add new node: no index available");
        }
    }

    /**
     * Modify the graph by removing a node at a specified index.
     * Note: will disconnect the node (i.e. remove any existing edge involving the node)
     * @param nid node index
     */
    public synchronized void removeNode(int nid) {
        if (isIndex(nid)) {
            for (int next : getNodeIndices(nid, true))
                this.removeEdge(nid, next);
            for (int prev : getNodeIndices(nid, false))
                this.removeEdge(prev, nid);
            this.nodes[nid] = null;
            this.edgesForward[nid] = null;
            if (isDirected())
                this.edgesBackward[nid] = null;
            if (isTerminated()) {
                this.startNodes.set(nid, false);
                this.endNodes.set(nid, false);
            }
        } else {
            throw new InvalidIndexRuntimeException("Index outside bounds: " + nid);
        }
    }

    /**
     * Modify the graph by adding an instance of an edge between two existing nodes.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized boolean addEdge(int from, int to) {
        if (isNode(from) && isNode(to)) {
            if (isDirected()) { // directed, we set both forward and backward pointers accordingly
                this.edgesBackward[to].set(from);
                this.edgesForward[from].set(to);
            } else { // un-directed, which means we record the same edge in BOTH directions
                this.edgesForward[to].set(from);
                this.edgesForward[from].set(to);
            }
        } else if (from == -1 && isNode(to) && isTerminated()) {
            this.startNodes.set(to);
        } else if (isNode(from) && to == this.maxsize() && isTerminated()) {
            this.endNodes.set(from);
        } else
            return false;
        return true;
    }

    public synchronized boolean addTerminalEdge(int from) {
        return addEdge(from, maxsize());
    }

    /**
     * Modify the graph by adding an instance of an edge between two existing nodes.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized void removeEdge(int from, int to) {
        if (from == -1 && to == 0)
            from = -1;
        if (isNode(from) && isNode(to)) {
            if (isDirected())
                this.edgesBackward[to].set(from, false);
            else
                this.edgesForward[to].set(from, false);
            this.edgesForward[from].set(to, false);
        } else if (from == -1 && isNode(to) && isTerminated()) {
            this.startNodes.set(to, false);
        } else if (isNode(from) && to == this.maxsize() && isTerminated()) {
            this.endNodes.set(from, false);
        } else {
            throw new InvalidIndexRuntimeException("Cannot remove edge between non-existent node/s: " + from + " or " + to);
        }
    }

    /**
     * Get indices of all nodes that can be reached in ONE step from the given node index.
     * If graph is undirected, then it does not matter what the source node was when added.
     * If graph is directed, then the default is to only return the edges looking forward (from the source node when added).
     * @param idx
     * @return
     * @throws InvalidIndexRuntimeException if the index is invalid
     */
    public int[] getNodeIndices(int idx) {
        return getNodeIndices(idx, true);
    }

    /**
     * Get indices of all nodes that can be reached in ONE step relative the given node index
     * @param idx index to node
     * @param look_forward if true, look forward, else backward
     * @return
     * @throws InvalidIndexRuntimeException if the index is invalid
     */
    public int[] getNodeIndices(int idx, boolean look_forward) {
        int count = 0;
        if (isNode(idx)) {
            BitSet bs = (look_forward || !isDirected()) ? edgesForward[idx] : edgesBackward[idx];
            int[] indices = new int[bs.cardinality()];
            if (indices.length > 1)
                count = 0;
            for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else if (idx == -1 && (look_forward || !isDirected())) {
            int[] indices = new int[startNodes.cardinality()];
            for (int i = startNodes.nextSetBit(0); i >= 0; i = startNodes.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else if (idx == maxsize() && (!look_forward || !isDirected())) {
            int[] indices = new int[endNodes.cardinality()];
            for (int i = endNodes.nextSetBit(0); i >= 0; i = endNodes.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve edges from non-existent/invalid node: " + idx);
        }
    }


    /**
     * Generate a text string that describes the graph.
     * @return text description
     */
    public String toString() {
        StringBuffer buf = new StringBuffer();
        char direc = isDirected() ? '>' : '-';
        buf.append("{");
        if (isTerminated()) {
            int[] idxs = this.getNodeIndices(-1, true);
            for (int i = 0; i < idxs.length; i ++)
                buf.append(idxs[i] + (i < idxs.length - 1 ? "," : ""));
            buf.append("-" + direc + "[N=" + this.size() + "|E=" + this.getEdgeCount() + "]-" + direc);
            idxs = this.getNodeIndices(this.maxsize(), false);
            for (int i = 0; i < idxs.length; i ++)
                buf.append(idxs[i] + (i < idxs.length - 1 ? "," : ""));
        } else { // not terminated
            buf.append("[N=" + this.size() + "|E=" + this.getEdgeCount() + "]-" + direc);
        }
        buf.append("}");
        return buf.toString();
    }


    /**
     * Generate a text string that describes the graph, following the DOT format.
     * https://www.graphviz.org/pdf/dotguide.pdf
     *
     *   node  [style="rounded,filled,bold", shape=box, fixedsize=true, width=1.3, fontname="Arial"];
     *   edge  [style=bold, fontname="Arial", weight=100]
     * @return
     */
    public String toDOT() {
        StringBuffer buf = new StringBuffer();
        if (!isDirected())
            buf.append("graph " + getName() + " {\nnode [" + nodeDOT + "];\n");
        else
            buf.append("digraph " + getName() + " {\nrankdir=\"LR\";\nnode [" + nodeDOT + "];\n");
        for (int i = 0; i < nodes.length; i ++) {
            Node n = nodes[i];
            if (n != null) {
                buf.append(Integer.toString(i) + " [" + n.toDOT() + "];\n");
            }
        }
        if (isTerminated()) {
            buf.append("_start [label=\"S\",style=bold,fontcolor=red,width=0.25,fillcolor=gray,penwidth=0];\n");
            buf.append("_end [label=\"E\",style=bold,fontcolor=red,width=0.25,fillcolor=gray,penwidth=0];\n");
            buf.append("{rank=source;_start;}\n{rank=sink;_end;}\n");
        }
        buf.append("edge [" + edgeDOT + "];\n");
        if (isTerminated() && isDirected()) {
            for (int i = 0; i < startNodes.length(); i++) {
                if (startNodes.get(i)) {
                    buf.append("_start -> " + i + "\n");
                }
            }
        }
        for (int from = 0; from < edgesForward.length; from ++) {
            if (isNode(from)) {
                for (int to = isDirected() ? 0 : from; to < edgesForward[from].length(); to++) {
                    if (edgesForward[from].get(to)) {
                        if (!isDirected())
                            buf.append(from + " -- " + to + "\n");
                        else
                            buf.append(from + " -> " + to + "\n");
                    }
                }
            }
        }
        if (isTerminated() && isDirected()) {
            for (int i = 0; i < endNodes.length(); i++) {
                if (endNodes.get(i)) {
                    buf.append(i + " -> _end\n");
                }
            }
        }
        buf.append("}\n");
        return buf.toString();
    }

    public static class DefaultGraph extends IdxGraph {

        /**
         * Create graph.
         *
         * @param nNodes     number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
         * @param undirected if true, edges are not directed, else
         * @param terminated if true, a virtual start and a virtual end nodes are created to enable edges to mark multiple entry and exit nodes
         */
        public DefaultGraph(int nNodes, boolean undirected, boolean terminated) {
            super(nNodes, undirected, terminated);
        }

        /**
         * Modify the graph by adding a node at a specified index.
         * Note: does not connect the node.
         * @param nid node index, which must be valid, i.e. between 0 and N - 1, where N is the size of the possible/valid node-set
         */
        public synchronized int addNode(int nid) {
            Node node = new Node();
            return super.addNode(nid, node);
        }

        /**
         * Modify the graph by adding a node at an unspecified specified index (first free).
         * Note: does not connect the node.
         */
        public synchronized int addNode() {
            Node node = new Node();
            return super.addNode(node);
        }

    }

    public void saveToDOT(String filename) throws IOException {
        FileWriter fwriter=new FileWriter(filename);
        BufferedWriter writer=new BufferedWriter(fwriter);
        writer.write(toDOT());
        writer.close();
        fwriter.close();
    }

    public static void saveToDOT(String filename, IdxGraph... graphs) throws IOException {
        FileWriter fwriter=new FileWriter(filename);
        BufferedWriter writer=new BufferedWriter(fwriter);
        int cnt = 1;
        for (IdxGraph g : graphs) {
            writer.write(g.toDOT());
            writer.newLine();
        }
        writer.close();
        fwriter.close();
    }
}
