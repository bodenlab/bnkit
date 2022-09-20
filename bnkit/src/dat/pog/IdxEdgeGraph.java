package dat.pog;

import asr.ASRRuntimeException;
import asr.GRASP;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

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

    @Override
    public int hashCode() {
        int result = super.hashCode();
        Integer[] edgeidxs = new Integer[edges.keySet().size()];
        edges.keySet().toArray(edgeidxs);
        Arrays.sort(edgeidxs);
        result = 31 * result + Arrays.hashCode(edgeidxs);
        return result;
    }

    /**
     * Create a JSON representation of the instance
     * See {@link IdxGraph#toJSON}.
     * @return JSON object of this instance
     */
    public JSONObject toJSON() {
        JSONObject json = super.toJSON();
        List<JSONObject> edgelist = new ArrayList<>();
        JSONArray edgeidxs = new JSONArray();
        JSONArray edgeinst = new JSONArray();
        Map<Integer, E> edgemap = getEdges();
        Class edgeclass = null;
        for (Map.Entry<Integer, E> entry : edgemap.entrySet()) { //
            int edgeidx = entry.getKey();
            int from = getFrom(edgeidx);
            int to = getTo(edgeidx);
            edgeidxs.put(new JSONArray(new int[] {from, to}));
            E e = entry.getValue();
            if (edgeclass == null)
                edgeclass = e.getClass();
            else if (e.getClass() != edgeclass)
                throw new ASRRuntimeException("Mixing edges");
            edgeinst.put(e.toJSON());
        }
        if (edgeclass != null) {
            json.put("Edgeindices", edgeidxs);
            json.put("Edges", edgeinst);
            json.put("Edgetype", edgeclass);
        }
        return json;
    }


    /**
     * Retrieve all the edges by ref to index; these include disabled edges
     * @return all edges
     */
    public Map<Integer, E> getEdges() {
        return edges;
    }

    /**
     * Retrieve the edge instance for a pair of node indices
     * @param from source node index (use -1 for terminal start edge)
     * @param to target node index (use N for a terminal end edge, where N is the number of possible/valid nodes)
     * @return edge instance if exists, else null
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public E getEdge(int from, int to) {
        boolean valid_from = isNode(from);
        boolean valid_to = isNode(to);
        if (isNode(from) && isNode(to)) {
            int eidx = getEdgeIndex(from, to);
            return this.edges.get(eidx);
        } else if (from == -1 && isNode(to) && isTerminated()) {
            int eidx = getEdgeIndex(from, to);
            return this.edges.get(eidx);
        } else if (isNode(from) && to == this.maxsize() && isTerminated()) {
            int eidx = getEdgeIndex(from, to);
            return this.edges.get(eidx);
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve edge between non-existent node/s: " + from + (valid_from ? " (valid) -- " : " (invalid) -- ") + to + (valid_to ? " (valid)" : " (invalid)"));
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
        if (isEdge(from, to)) {
            super.removeEdge(from, to);
            this.edges.remove(getEdgeIndex(from, to));
        }
    }

    /**
     * Modify the graph by disabling an edge between two existing nodes.
     * If there is an edge, it will be stashed to be enabled at a later stage.
     * Exercise some caution with this functionality since edges are enabled at the level of IdxGraph, but their instances are stored in IdxEdgeGraph.
     *
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @return true if an edge exists, and then is disabled, false if no edge could be found
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized boolean disableEdge(int from, int to) {
        if (isEdge(from, to)) {
            super.removeEdge(from, to);
            return true;
        }
        return false;
    }

    /**
     * Modify the graph by adding an instance of an edge between two existing nodes.
     * If the edge already exists between the same node indices, it will be replaced.
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

    /**
     * Modify the graph by enabling an instance of an edge between two existing nodes.
     * Requires that the edge already exists between the same node indices.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     * @param from the source node index of the edge, -1 for a terminal start edge
     * @param to the destination node index of the edge, N for a terminal end edge, where N is the number of possible/valid nodes
     * @return true if the edge instance was already there, so enabling could be completed, false if no edge instance was available
     * @throws InvalidIndexRuntimeException if either node index is invalid
     */
    public synchronized boolean enableEdge(int from, int to) {
        if (edges.get(getEdgeIndex(from, to)) != null) {
            super.addEdge(from, to);
            return true;
        } else
            return false;
    }

    public synchronized boolean addTerminalEdge(int from, E edge) {
        return addEdge(from, maxsize(), edge);
    }

    /**
     * Generate a text string that describes the graph, following the DOT format.
     * https://www.graphviz.org/pdf/dotguide.pdf
     *
     *   node  [style="rounded,filled,bold", shape=box, fixedsize=true, width=1.3, fontname="Arial"];
     *   edge  [style=bold, fontname="Arial", weight=100]
     * @return
     */
    @Override
    public String toDOT() {
        StringBuffer buf = new StringBuffer();
        if (!isDirected())
            buf.append("graph " + getName() + " {\nnode [" + nodeDOT + "];\n");
        else
            buf.append("digraph " + getName() + " {\nrankdir=\"LR\";\nnode [" + nodeDOT + "];\n");
        for (int i = 0; i < nodes.length; i ++) {
            Node n = nodes[i];
            if (n != null) {
//                if (n.getLabel() == null) {
//                    n.setLabel(Integer.toString(i));
//                    buf.append(Integer.toString(i) + " [" + n.toDOT() + "];\n");
//                    n.setLabel(null);
//                } else
                    buf.append(Integer.toString(i) + " [" + n.toDOT() + "];\n");
            }
        }
        if (isTerminated()) {
            buf.append("_start [label=\"S(" + getName() + ")\",style=bold,fontcolor=red,fillcolor=gray,penwidth=0];\n");
            buf.append("_end [label=\"E(" + getName() +")\",style=bold,fontcolor=red,fillcolor=gray,penwidth=0];\n");
            buf.append("{rank=source;_start;}\n{rank=sink;_end;}\n");
        }

        buf.append("edge [" + edgeDOT + "];\n");
        // _start -> 0 -> 1 -> 2 -> 3 -> 4 -> 5 -> 6 -> _end [style=invis]
        if (isIndexTopologicalOrder()) {
            if (isTerminated() && isDirected())
                buf.append("_start -> ");
            boolean cropped = false;
            for (int from = 0; from < edgesForward.length; from ++) {
                if (isNode(from)  && isDirected()) {
                    buf.append(from + (from < edgesForward.length - 1 ? " -> " : " "));
                    cropped = false;
                } else
                    cropped = true;
            }
            if (isTerminated()) {
                if (!cropped)
                    buf.append("-> ");
                buf.append("_end [style=invis]\n");
            } else
                buf.append("[style=invis]\n");
        }
        if (isTerminated() && isDirected()) {
            for (int i = 0; i < startNodes.length(); i++) {
                if (startNodes.get(i)) {
                    E edge = getEdge(-1, i);
                    buf.append("_start -> " + i + (edge == null ? "\n" : "[" + edge.toDOT() + "]\n"));
                }
            }
        }
        for (int from = 0; from < edgesForward.length; from ++) {
            if (isNode(from)) {
                for (int to = isDirected() ? 0 : from; to < edgesForward[from].length(); to++) {
                    if (edgesForward[from].get(to)) {
                        E edge = getEdge(from, to);
                        if (!isDirected())
                            buf.append(from + " -- " + to + (edge == null ? "\n" : "[" + edge.toDOT() + "]\n"));
                        else
                            buf.append(from + " -> " + to + (edge == null ? "\n" : "[" + edge.toDOT() + "]\n"));
                    }
                }
            }
        }
        if (isTerminated() && isDirected()) {
            for (int i = 0; i < endNodes.length(); i++) {
                if (endNodes.get(i)) {
                    E edge = getEdge(i, nodes.length);
                    buf.append(i + " -> _end" + (edge == null ? "\n" : "[" + edge.toDOT() + "]\n"));
                }
            }
        }
        buf.append("}\n");
        return buf.toString();
    }

    /**
     * Generate a matrix that describes the graph, where
     * row means "from" index
     * col means "to" index
     * @return matrix representation of graph
     */
    @Override
    public int[][] toMatrix() {
        int[][] m = new int[this.nNodes+(isTerminated()?2:0)][this.nNodes + (isTerminated()?2:0)];
        if (isTerminated()) {
            for (int i = 0; i < startNodes.length(); i++) {
                if (startNodes.get(i)) {
                    Edge e = getEdge(-1, i);
                    if (e==null && isEdge(-1, i))
                        m[0][i + 1] = 1;
                    else {
                        try {
                            POGraph.BidirEdge bidir = (POGraph.BidirEdge) e;
                            if (bidir.isForward())
                                m[0][i + 1] = 1;
                            if (bidir.isBackward())
                                m[i + 1][0] = 1;
                        } catch (ClassCastException e1) {
                            try {
                                POGraph.StatusEdge edge = (POGraph.StatusEdge) e;
                                if (edge.getReciprocated())
                                    m[i + 1][0] = 1;
                                else
                                    m[0][i + 1] = 1;
                            } catch (ClassCastException e2) {
                                try {
                                    SeqEdge edge = (SeqEdge) e;
                                    if (edge != null) {
                                        if (edge != null)
                                            m[0][i + 1] = edge.getSeqs().size();
                                    }
                                } catch (ClassCastException e3) {
                                    ;
                                }
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < endNodes.length(); i++) {
                if (endNodes.get(i)) {
                    Edge e = getEdge(i, nNodes);
                    if (e==null && isEdge(i, nNodes))
                        m[i + 1][nNodes + 1] = 1;
                    else {
                        try {
                            POGraph.BidirEdge bidir = (POGraph.BidirEdge) e;
                            if (bidir != null) {
                                if (bidir.isForward())
                                    m[i + 1][nNodes + 1] = 1;
                                if (bidir.isBackward())
                                    m[nNodes + 1][i + 1] = 1;
                            }
                        } catch (ClassCastException e1) {
                            try {
                                POGraph.StatusEdge edge = (POGraph.StatusEdge) e;
                                if (edge != null) {
                                    if (edge.getReciprocated())
                                        m[nNodes + 1][i + 1] = 1;
                                    else
                                        m[i + 1][nNodes + 1] = 1;
                                }
                            } catch (ClassCastException e2) {
                                try {
                                    SeqEdge edge = (SeqEdge) e;
                                    if (edge != null) {
                                        if (edge != null)
                                            m[i + 1][nNodes + 1] = edge.getSeqs().size();
                                    }
                                } catch (ClassCastException e3) {
                                    ;
                                }

                            }
                        }
                    }
                }
            }
            for (int from = 0; from < edgesForward.length; from++) {
                if (isNode(from)) {
                    for (int to = isDirected() ? 0 : from; to < edgesForward[from].length(); to ++) {
                        if (edgesForward[from].get(to)) {
                            Edge e = getEdge(from, to);
                            if (e==null && isEdge(from, to))
                                m[from + 1][to + 1] = 1;
                            else {
                                try {
                                    POGraph.BidirEdge bidir = (POGraph.BidirEdge) e;
                                    if (bidir != null) {
                                        if (bidir.isForward())
                                            m[from + 1][to + 1] = 1;
                                        if (bidir.isBackward())
                                            m[to + 1][from + 1] = 1;
                                    }
                                } catch (ClassCastException e1) {
                                    try {
                                        POGraph.StatusEdge edge = (POGraph.StatusEdge) e;
                                        if (edge != null) {
                                            if (edge.getReciprocated())
                                                m[to + 1][from + 1] = 1;
                                            else
                                                m[from + 1][to + 1] = 1;
                                        }
                                    } catch (ClassCastException e2) {
                                        try {
                                            SeqEdge edge = (SeqEdge) e;
                                            if (edge != null) {
                                                if (edge != null)
                                                    m[from + 1][to + 1] = edge.getSeqs().size();
                                            }
                                        } catch (ClassCastException e3) {
                                            ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else { // not terminated
            for (int from = 0; from < edgesForward.length; from++) {
                if (isNode(from)) {
                    for (int to = isDirected() ? 0 : from; to < edgesForward[from].length(); to ++) {
                        if (edgesForward[from].get(to)) {
                            Edge e = getEdge(from, to);
                            if (e==null && isEdge(from, to))
                                m[from + 1][to + 1] = 1;
                            else {
                                try {
                                    POGraph.BidirEdge bidir = (POGraph.BidirEdge) e;
                                    if (bidir != null) {
                                        if (bidir.isForward())
                                            m[from + 1][to + 1] = 1;
                                        if (bidir.isBackward())
                                            m[to + 1][from + 1] = 1;
                                    }
                                } catch (ClassCastException e1) {
                                    try {
                                        POGraph.StatusEdge edge = (POGraph.StatusEdge) e;
                                        if (edge != null) {
                                            if (edge.getReciprocated())
                                                m[to + 1][from + 1] = 1;
                                            else
                                                m[from + 1][to + 1] = 1;
                                        }
                                    } catch (ClassCastException e2) {
                                        try {
                                            SeqEdge edge = (SeqEdge) e;
                                            if (edge != null) {
                                                if (edge != null)
                                                    m[from + 1][to + 1] = edge.getSeqs().size();
                                            }
                                        } catch (ClassCastException e3) {
                                            ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return m;
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
