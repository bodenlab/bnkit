package dat.pog;

import asr.ASRException;
import asr.ASRRuntimeException;
import asr.GRASP;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class IdxGraph {
    protected final boolean directed;
    public String nodeDOT = "style=\"rounded,filled\", shape=box, fixedsize=true";
    public String edgeDOT = ""; //""style=\"bold\"";
    protected Node[] nodes;    // indexed by internal node ID i in {0..N}; FIXME: node i does not exist if nodes[i] = null
    protected boolean[] allnodes;
    protected int nNodes;
    protected BitSet[] edgesForward;
    protected BitSet[] edgesBackward = null;
    protected BitSet startNodes = null;
    protected BitSet endNodes = null;
    protected String name = null;

    public IdxGraph(int nNodes, boolean undirected, boolean terminated) {
        this.nodes = new Node[nNodes];
        this.allnodes = new boolean[nNodes]; // initially NO nodes
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

    /**
     * Compare two graphs logically
     * Implemented so that two graphs can be compared "logically"
     * @param obj
     * @return
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null)
            return false;
        if (!(obj instanceof IdxGraph))
            return false;
        IdxGraph other = (IdxGraph) obj;
        if ((name == null && other.name != null) || (name != null && other.name == null))
            return false;
        if (name != null && other.name != null) {
            if (!name.equals(other.name))
                return false;
        }
        if (directed != other.isDirected())
            return false;
        if (nNodes != other.nNodes)
            return false;
        for (int i = 0; i < nNodes; i ++) {
            if (nodes[i] != null)
                if (!nodes[i].equals(other.nodes[i]))
                    return false;
        }
        if (isTerminated()) {
            if (!startNodes.equals(other.startNodes))
                return false;
            if (!endNodes.equals(other.endNodes))
                return false;
        }
        if (edgesForward.length != other.edgesForward.length)
            return false;
        for (int i = 0; i < edgesForward.length; i ++) {
            if (edgesForward[i] == null && other.edgesForward[i] == null)
                continue;
            if (edgesForward[i] != null)
                if (!edgesForward[i].equals(other.edgesForward[i]))
                    return false;
        }
        if (isDirected()) {
            if (edgesBackward.length != other.edgesBackward.length)
                return false;
            for (int i = 0; i < edgesBackward.length; i ++) {
                if (edgesBackward[i] == null && other.edgesBackward[i] == null)
                    continue;
                if (edgesBackward[i] != null)
                    if (!edgesBackward[i].equals(other.edgesBackward[i]))
                        return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(directed, nNodes, startNodes, endNodes, name);
        result = 31 * result + Arrays.hashCode(edgesForward);
        result = 31 * result + Arrays.hashCode(edgesBackward);
        return result;
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
            cnt += this.allnodes[i] ? 1 : 0;
        return cnt;
    }

    /**
     * Check if the index 0, 1, 2, ..., N is also a valid topological order
     * @return true if index is a valid topological ordering, false if it isn't
     */
    public boolean isIndexTopologicalOrder() {
        for (int i = 0; i < edgesForward.length; i ++) {
            BitSet next = edgesForward[i];
            if (next != null) {
                for (int j = next.nextSetBit(0); j >= 0; j = next.nextSetBit(j + 1))
                    if (j < i)
                        return false;
            }
        }
        return true;
    }

    public synchronized int getFreeIndex() {
        for (int i = 0; i < this.nNodes; i ++)
            if (this.allnodes[i] == false)
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
            return allnodes[idx];
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
        int supersize = maxsize() + 1;
        if (isNode(from) && isNode(to)) {
            if (!isDirected() && to < from) {
                int bigger = from;
                from = to;
                to = bigger;
            }
            return from * supersize + to;
        } else if (from == -1 && isNode(to) && isTerminated())
            return -to;
        else if (isNode(from) && to == this.maxsize() && isTerminated())
            return from * supersize + to;
        throw new InvalidIndexRuntimeException("Cannot be an edge: " + from + " -- " + to);
    }

    /**
     * Recover the "from" index from an edge index generated by {@link #getEdgeIndex(int from, int to)}
     * @param edgeidx
     * @return
     */
    public int getFrom(Integer edgeidx) {
        int supersize = maxsize() + 1;
        int from;
        if (isTerminated()) {
            if (edgeidx <= 0)
                from = -1;
            else
                from = edgeidx / supersize;
            return from;
        }
        from = edgeidx / supersize;
        return isDirected() ? from : Math.min(from, edgeidx % supersize);
    }

    /**
     * Recover the "to" index from an edge index generated by {@link #getEdgeIndex(int from, int to)}
     * @param edgeidx
     * @return
     */
    public int getTo(Integer edgeidx) {
        int supersize = maxsize() + 1;
        int to;
        if (isTerminated()) {
            if (edgeidx <= 0)
                to = -edgeidx;
            else
                to = edgeidx % supersize;
            return to;
        }
        to = edgeidx % supersize;
        return isDirected() ? to : Math.max(edgeidx / supersize ,to);
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
            this.allnodes[nid] = true;
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
            this.allnodes[nid] = false;
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
     * If the edge already exists it will be replaced.
     * If the graph is terminated, an edge can be added from a virtual start, or to a virtual end node.
     *
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
     * Modify the graph by removing an edge between two existing nodes.
     * If the graph is terminated, an edge can be removed from a virtual start, or to a virtual end node.
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
     * @return an array with all indices; if no indices, the array will be length 0
     * @throws InvalidIndexRuntimeException if the index is invalid
     */
    public int[] getNodeIndices(int idx) {
        return getNodeIndices(idx, true);
    }

    /**
     * Get indices of all nodes that can be reached in ONE step relative the given node index
     * @param idx index to node
     * @param look_forward if true, look forward, else backward
     * @return an array with all indices; if no indices, the array will be length 0
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
        } else if (idx == -1 && (look_forward && isDirected())) {
            int[] indices = new int[startNodes.cardinality()];
            for (int i = startNodes.nextSetBit(0); i >= 0; i = startNodes.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else if (idx == maxsize() && (!look_forward && isDirected())) {
            int[] indices = new int[endNodes.cardinality()];
            for (int i = endNodes.nextSetBit(0); i >= 0; i = endNodes.nextSetBit(i + 1))
                indices[count ++] = i;
            return indices;
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve edges from non-existent/invalid node: " + idx);
        }
    }

    /**
     * Get number of nodes that can be reached in ONE step relative the given node index
     * @param idx index of node
     * @param look_forward if true, look forward, else backward
     * @return cardinality of node
     * @throws InvalidIndexRuntimeException if the node index is invalid
     */
    public int getCardinality(int idx, boolean look_forward) {
        if (isNode(idx)) {
            BitSet bs = (look_forward || !isDirected()) ? edgesForward[idx] : edgesBackward[idx];
            return bs.cardinality();
        } else if (idx == -1 && (look_forward || !isDirected())) {
            return startNodes.cardinality();
        } else if (idx == maxsize() && (!look_forward || !isDirected())) {
            return endNodes.cardinality();
        } else {
            throw new InvalidIndexRuntimeException("Cannot retrieve cardinality of non-existent/invalid node: " + idx);
        }
    }

    private static int[] concat(int[] array1, int[] array2) {
        int[] result = Arrays.copyOf(array1, array1.length + array2.length);
        System.arraycopy(array2, 0, result, array1.length, array2.length);
        return result;
    }

    /**
     * Determine a topological ordering of the nodes. This function can only be used on directed and terminated graphs.
     * Note that more than one topological order can exist for bi- and multi-furcating POGs
     * @return an array with the indices of the POG in a topological order
     * TODO: deprecate this in favour of getTopologicalDepthFirst. Currently they differ in that this function includes -1, which seems unreasonable at first glance...
     */
    public int[] getTopologicalOrder() {
        if (!(isDirected() && isTerminated()))
            throw new RuntimeException("Topological order cannot be determined when graph is undirected");
        int[] ret = new int[this.size() + 1];
        boolean[] added = new boolean[this.nNodes];
        int[] inedges = new int[this.nNodes];
        int add_ptr = 0; // pointer to position for adding index
        ret[add_ptr ++] = -1; // first node
        // add all nodes with NO edges leading IN, i.e. zero backward edges
        for (int i = 0; i < this.nNodes; i++) {
            try {
                if (getCardinality(i, false) == 0) {
                    ret[add_ptr ++] = i;
                    added[i] = true;
                }
            } catch (InvalidIndexRuntimeException e) {
                ;
            }
        }
        // go through all nodes, first those added above, incrementally adding those that are beyond/forward
        int exp_ptr = 0; // pointer to position to expand
        while (exp_ptr < add_ptr) {
            int current = ret[exp_ptr ++];
            int[] next = getNodeIndices(current, true);
            // any node in the just expanded list could be the start of a linear stretch of nodes
            for (int idx : next) {
                if (!added[idx]) {
                    inedges[idx] = getCardinality(idx, false) + (isStartNode(idx) ? 1 : 0); // check how many paths lead to it
                    added[idx] = true;
                }
                if (inedges[idx] == 1) { // only add it, when exactly one path leads to it, or when we have expanded the last path to it
                    ret[add_ptr++] = idx;
                }
                if (inedges[idx] > 1) // more paths will lead to this node index, so count down to catch it later
                    inedges[idx] --;
            }
        }
        return ret;
    }

    /**
     * Determine a topological ordering of the nodes by depth-first search, which maximises continuity of paths.
     * This function can only be used on directed and terminated graphs.
     * Note that more than one topological order can exist for bi- and multi-furcating POGs
     * @return an array with the indices of the POG in a topological order
     */
    public int[] getTopoSortDepthFirst() {
        if (!(isDirected() && isTerminated()))
            throw new RuntimeException("Topological order cannot be determined when graph is undirected");
        TopologicalSort tsort = new TopologicalSort(this);
        return tsort.getSorted();
    }

    /**
     * Topological sort by depth-first search
     * https://en.wikipedia.org/wiki/Topological_sorting
     * Cormen, Thomas H.; Leiserson, Charles E.; Rivest, Ronald L.; Stein, Clifford (2001),
     * "Section 22.4: Topological sort", Introduction to Algorithms (2nd ed.), MIT Press and McGraw-Hill, pp. 549–552, ISBN 0-262-03293-7
     */
    private static class TopologicalSort {
        final boolean[] permanent;
        final boolean[] temporary;
        final LinkedList<Integer> sorted;
        final IdxGraph graph;

        TopologicalSort(IdxGraph graph) {
            this.graph = graph;
            permanent = new boolean[graph.nNodes];
            temporary = new boolean[graph.nNodes];
            sorted = new LinkedList<>();
        }

        private boolean visit(int n) {
            if (permanent[n]) return false;
            if (temporary[n]) throw new RuntimeException("Graph is not a valid DAG, hence topological sort not defined");
            temporary[n] = true;
            int[] next = graph.getNodeIndices(n, true);
            for (int m : next) {
                visit(m);
            }
            temporary[n] = false;
            permanent[n] = true;
            sorted.addFirst(n);
            return true;
        }

        /**
         * Determine the topological sort, starting at the start nodes, using depth-first search.
         * Any nodes that are not reachable from the start nodes will not be contained in the sort.
         * As an aside, unreachable/disconnected nodes can be included if the loop is extended to all those recorded in the graph.
         * @return an array with the indices of nodes ordered by topological sort
         */
        int[] getSorted() {
            int[] starts = graph.getStarts();
            for (int j = 0; j < starts.length; j++) {
                if (!permanent[starts[j]]) {
                    visit(starts[j]);
                }
            }
            int[] ret = new int[sorted.size()];
            int i = 0;
            for (Integer idx : sorted)
                ret[i ++] = idx;
            return ret;
        }

    }
    /**
     * Determine the linear order of nodes that start at given index, and end before cardinality exceeds one
     * @param idx start node
     * @param look_forward increment order going forward (true) or backward (false)
     * @return the node indices in linear order
     */
    public int[] getOrdered(int idx, boolean look_forward) {
        if (!isDirected())
            throw new RuntimeException("Topological order undefined when graph is undirected");
//        if (getCardinality(idx, look_forward) != 1 && getCardinality(idx, !look_forward) != 1)
//            return new int[0];
        List<Integer> path = new ArrayList<>();
        while (isPath(idx)) {
            path.add(idx);
            idx = look_forward ? this.edgesForward[idx].nextSetBit(0) : this.edgesBackward[idx].nextSetBit(0);
        }
        int[] ret = new int[path.size()];
        for (int i = 0; i < path.size(); i ++)
            ret[i] = path.get(i);
        return ret;
    }

    /**
     * Get indices that can be found and linearly ordered from each of the given starting indices.
     * @param expand_us array of starting indices
     * @param look_forward direction of search
     * @return one array of ordered indices for each starting index;
     * note the last index in each array will potentially be connected to other parts in the graph
     */
    public int[][] getOrdered(int[] expand_us, boolean look_forward) {
        if (!isDirected())
            throw new RuntimeException("Topological order undefined when graph is undirected");
        int[][] arr = new int[expand_us.length][];
        int ntot = 0;
        for (int i = 0; i < expand_us.length; i ++) {
            arr[i] = getOrdered(expand_us[i], look_forward);
            // last node here, will potentially be expanded next, save ...
        }
        return arr;
    }


    /**
     * Determine if node is on a linear path, i.e. has only one entry edge and one exit edge.
     * @param idx node index
     * @return true if node is part of a linear sequence of nodes
     */
    public boolean isPath(int idx) {
        if (idx == -1 || idx == nNodes) // termination does not form part of any path
            return false;
        int cntForward = edgesForward[idx].cardinality() + (isEndNode(idx) ? 1 : 0);
        int cntBackward = edgesBackward[idx].cardinality() + (isStartNode(idx) ? 1 : 0);
        return (cntForward == 1 && cntBackward == 1);
    }

    /**
     * Retrieve all indices for nodes that start this graph
     * @return
     */
    public int[] getStarts() {
        if (!isTerminated())
            throw new RuntimeException("Graph must be terminated to have start nodes");
        List<Integer> list = new ArrayList<>();
        for (int i = startNodes.nextSetBit(0); i >= 0; i = startNodes.nextSetBit(i+1)) {
            list.add(i);
            if (i == Integer.MAX_VALUE)
                break; // or (i+1) would overflow
        }
        int[] ret = new int[list.size()];
        for (int i = 0; i < list.size(); i ++)
            ret[i] = list.get(i);
        return ret;
    }

    /**
     * Retrieve all indices for nodes that end this graph
     * @return
     */
    public int[] getEnds() {
        if (!isTerminated())
            throw new RuntimeException("Graph must be terminated to have end nodes");
        List<Integer> list = new ArrayList<>();
        for (int i = endNodes.nextSetBit(0); i >= 0; i = endNodes.nextSetBit(i+1)) {
            list.add(i);
            if (i == Integer.MAX_VALUE)
                break; // or (i+1) would overflow
        }
        int[] ret = new int[list.size()];
        for (int i = 0; i < list.size(); i ++)
            ret[i] = list.get(i);
        return ret;
    }

    /**
     * Get all indices that can be found and linearly ordered before expanding again.
     * @param idx
     * @param look_forward
     * @return
     */
    public int[] getOrderedWrap(int idx, boolean look_forward) {
        if (!isDirected())
            throw new RuntimeException("Topological order undefined when graph is undirected");
        int[] expand_us = null;
        if (isTerminated() && idx == -1 && look_forward){
            expand_us = getStarts();
        } else if (isTerminated() && idx == nNodes && !look_forward) {
            expand_us = getEnds();
        } else {
            expand_us = getNodeIndices(idx, look_forward);
        }
        int[][] arr = getOrdered(expand_us, look_forward);
        int ntot = 0;
        for (int i = 0; i < arr.length; i ++)
            ntot += arr[i].length;
        int[] ret = new int[ntot];
        int cnt = 0;
        for (int i = 0; i < arr.length; i ++)
            for (int j = 0; j < arr[i].length; j ++)
                ret[cnt ++] = arr[i][j];
        return ret;
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
            if (n != null)
                buf.append(i + " [" + n.toDOT() + "];\n");
            else
                buf.append(i + ";\n");
        }
        if (isTerminated()) {
            buf.append("_start [label=\"S(" + getName() + ")\",style=bold,fontcolor=red,fillcolor=gray,penwidth=0];\n");
            buf.append("_end [label=\"E(" + getName() +")\",style=bold,fontcolor=red,fillcolor=gray,penwidth=0];\n");
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

    /**
     * Generate a matrix that describes the graph, where
     * row means "from" index
     * col means "to" index
     * @return
     */
    public int[][] toMatrix() {
        int[][] m = new int[this.nNodes+(isTerminated()?2:0)][this.nNodes + (isTerminated()?2:0)];
        if (isTerminated()) {
            for (int i = 0; i < startNodes.length(); i++)
                m[0][i + 1] = startNodes.get(i) ? 1 : 0;
            for (int i = 0; i < endNodes.length(); i++)
                m[i + 1][nNodes + 1] = endNodes.get(i) ? 1 : 0;
            for (int from = 0; from < edgesForward.length; from++) {
                if (isNode(from)) {
                    for (int to = isDirected() ? 0 : from; to < edgesForward[from].length(); to++)
                        m[from + 1][to + 1] = edgesForward[from].get(to) ? 1 : 0;
                }
            }
        } else { // not terminated
            for (int from = 0; from < edgesForward.length; from++) {
                if (isNode(from)) {
                    for (int to = isDirected() ? 0 : from; to < edgesForward[from].length(); to++)
                        m[from][to] = edgesForward[from].get(to) ? 1 : 0;
                }
            }
        }
        return m;
    }

    public String toMatrixString() {
        StringBuilder sb = new StringBuilder();
        int[][] m = toMatrix();
        for (int r = 0; r < m.length; r++) {
            for (int c = 0; c < m[r].length; c++)
                sb.append(String.format("%4d ", m[r][c]));
            sb.append("\n");
        }
        return sb.toString();
    }

    /*
    ${\tt\left(\begin{array}{ccccccccc}
      0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
      1 & 0 & 1 & 1 & 1 & 0 & 0 & 0 & 0\\
      0 & 1 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
      0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
      0 & 1 & 1 & 1 & 0 & 1 & 0 & 0 & 0\\
      0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 0\\
      0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0\\
      0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1\\
      0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0
    \end{array}\right)}$
     */
    public String toLaTeXString(String label) {
        StringBuilder sb = new StringBuilder();
        int[][] m = toMatrix();
        sb.append("$"+ label + "{\\tt\\left(\\begin{array}{");
        if (m.length>0)
            for (int c = 0; c < m[0].length; c++)
                sb.append("c");
        sb.append("}\n");
        for (int r = 0; r < m.length; r++) {
            for (int c = 0; c < m[r].length; c++) {
                sb.append(String.format("%4d ", m[r][c]));
                if (c < m[r].length - 1)
                    sb.append("& ");
            }
            sb.append("\\\\ \n");
        }
        sb.append("\\end{array}\\right)}$\n");
        return sb.toString();
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

    /**
     * Put the class specific info from JSON into instance
     * @param json
     */
    protected void useJSON(JSONObject json) {
        try {
            String name = json.optString("Name", null);
            if (name != null) this.setName(name);
            JSONArray narr = json.getJSONArray("Indices");
            JSONArray earr = json.getJSONArray("Adjacent");
            if (narr.length() != earr.length())
                throw new RuntimeException("IdxGraph format is wrong in JSON object: mismatch in number of Indices and Adjacent");
            // first, retrieve indices for nodes
            int[] indices = new int[narr.length()];
            for (int i = 0; i < narr.length(); i ++)
                indices[i] = narr.getInt(i);
            JSONArray nodelist = json.getJSONArray("Nodes");
            try {
                Class nodetype = (Class) json.get("Nodetype");
                // second, retrieve actual nodes
                for (int i = 0; i < nodelist.length(); i++) {
                    JSONObject jnode = nodelist.getJSONObject(i);
                    if (nodetype.isAssignableFrom(SymNode.class))
                        this.addNode(indices[i], SymNode.fromJSON(jnode));
                    else if (nodetype.isAssignableFrom(EnumNode.class))
                        this.addNode(indices[i], EnumNode.fromJSON(jnode));
                    else
                        this.addNode(indices[i], Node.fromJSON(jnode));
                }
            } catch (ClassCastException e) {
                String nodetype = json.getString("Nodetype");
                // second, retrieve actual nodes
                for (int i = 0; i < nodelist.length(); i++) {
                    JSONObject jnode = nodelist.getJSONObject(i);
                    if (nodetype.equals(SymNode.class.toString()))
                        this.addNode(indices[i], SymNode.fromJSON(jnode));
                    else if (nodetype.equals(EnumNode.class.toString()))
                        this.addNode(indices[i], EnumNode.fromJSON(jnode));
                    else
                        this.addNode(indices[i], Node.fromJSON(jnode));
                }
            }
            // add edges next; connected nodes need to have been added previously
            for (int i = 0; i < narr.length(); i ++) {
                int idx = indices[i];
                JSONArray edgelist = earr.getJSONArray(i);
                for (int j = 0; j < edgelist.length(); j ++)
                    this.addEdge(idx, edgelist.getInt(j));
            }
            if (isTerminated()) {
                // add edges that start [from -1] and end [to n] graph
                JSONArray starts = json.getJSONArray("Starts");
                for (int i = 0; i < starts.length(); i ++)
                    addEdge(-1, starts.getInt(i));
                JSONArray ends = json.getJSONArray("Ends");
                for (int i = 0; i < ends.length(); i ++)
                    addTerminalEdge(ends.getInt(i));
            }
        } catch (JSONException e) {
            throw new ASRRuntimeException("IdxGraph format is wrong in JSON object: " + e.getMessage());
        }
    }

    /**
     * Decode JSON object into an instance of IdxGraph
     * @param json
     * @return
     * @throws RuntimeException if there is a formatting problem
     */
    public static IdxGraph fromJSON(JSONObject json) {
        try {
            String version = json.optString("GRASP_version", null);
            Class datatype = (Class)json.get("Datatype");
            if (version != null)
                if (!(version.equals(GRASP.VERSION) || version.equals("30-July-2022")))
                    throw new ASRRuntimeException("Invalid version: " + version);
            if (datatype != null)
                if (!(IdxGraph.class.isAssignableFrom(datatype)))
                    throw new ASRRuntimeException("Invalid datatype: " + datatype);
            int n = json.getInt("Size");
            boolean terminated = json.getBoolean("Terminated");
            boolean directed = json.getBoolean("Directed");
            IdxGraph g = new IdxGraph(n, !directed, terminated);
            g.useJSON(json);
            return g;
        } catch (JSONException e) {
            throw new ASRRuntimeException("IdxGraph format is wrong in JSON object: " + e.getMessage());
        }
    }

    /**
     * Create a JSON representation of the instance
     * @return JSON object of this instance
     */
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        json.put("Datatype", this.getClass());
        json.put("GRASP_version", GRASP.VERSION);
        json.put("Name", getName().length() == 0 ? null : getName());
        json.put("Size", this.nNodes);
        json.put("Terminated", isTerminated());
        json.put("Directed", isDirected());
        if (this.isTerminated()) {
            json.put("Starts", new JSONArray(this.getStarts()));
            json.put("Ends", new JSONArray(this.getEnds()));
        }
        JSONArray narr = new JSONArray();
        JSONArray earr = new JSONArray();
        Class nodetype = null;
        List<JSONObject> nodelist = new ArrayList<>();
        for (int idx = 0; idx < nNodes; idx ++) {
            if (nodes[idx] != null) {
                if (nodetype == null)
                    nodetype = nodes[idx].getClass();
                else if (nodetype != nodes[idx].getClass()) {
                    throw new ASRRuntimeException("Invalid mix of nodetypes " + nodetype + " and " + nodes[idx].getClass().toString());
                }
                JSONObject node = nodes[idx].toJSON();
                narr.put(idx);
                int[] edges = getNodeIndices(idx, true);
                earr.put(new JSONArray(edges));
                nodelist.add(node);
            }
        }
        json.put("Nodetype", nodetype);
        json.put("Indices", narr);
        json.put("Adjacent", earr);
        json.put("Nodes", new JSONArray(nodelist.toArray()));
        return json;
    }

    /**
     * Create a JSON array from a collection of graphs, each a IdxGraph or a sub-class thereof.
     * @param graphs
     * @return a series of JSON objects, forming a JSONArray (which is also a JSON object)
     */
    public static <E extends IdxGraph> JSONArray toJSONArray(Collection<E> graphs) {
        JSONArray narr = new JSONArray();
        int cnt = 0;
        for (IdxGraph g : graphs)
            narr.put(g.toJSON());
        return narr;
    }

    public static <E extends IdxGraph> JSONArray toJSONArray(E[] graphs) {
        JSONArray narr = new JSONArray();
        int cnt = 0;
        for (IdxGraph g : graphs)
            narr.put(g.toJSON());
        return narr;
    }

    /**
     * Save a map of graphs, where each entry has as key a string "labeL, then the graph as value.
     * The label will then be used to tag each graph.
     * @param graphs
     * @return
     */
    public static <E extends IdxGraph> JSONObject toJSON(Map<String, E> graphs) {
        JSONObject json = new JSONObject();
        int cnt = 0;
        for (Map.Entry<String, E> entry : graphs.entrySet())
            json.put(entry.getKey(), entry.getValue().toJSON());
        return json;
    }

    public static <E extends IdxGraph> void saveToJSON(String directory, Collection<E> graphs) throws IOException, ASRException {
        Path filename = Paths.get(directory + "/" + "pogs.json");
        BufferedWriter writer = Files.newBufferedWriter(filename, StandardCharsets.UTF_8);
        IdxGraph.toJSONArray(graphs).write(writer);
        writer.close();
    }

    public static <E extends IdxGraph> void saveToJSON(String directory, Map<String, E> graphs) throws IOException, ASRException {
        StringBuilder sb = new StringBuilder();
        Path filename = Paths.get(directory + "/" + "pogs.json");
        BufferedWriter writer = Files.newBufferedWriter(filename, StandardCharsets.UTF_8);
        IdxGraph.toJSON(graphs).write(writer);
        writer.close();
    }

    public void saveToDOT(String filename) throws IOException {
        FileWriter fwriter=new FileWriter(filename);
        BufferedWriter writer=new BufferedWriter(fwriter);
        writer.write(toDOT());
        writer.close();
        fwriter.close();
    }

    public static void saveToDOT(String directory, IdxGraph... graphs)  throws IOException, ASRException {
        Map<Object, IdxGraph> saveme1 = new HashMap<>();
        for (int idx = 0; idx < graphs.length; idx ++)
            saveme1.put("N" + idx, graphs[idx]);
    }

    public static void saveToDOT(String directory, Map<Object, IdxGraph> graphs) throws IOException, ASRException {
        StringBuilder sb = new StringBuilder();
        FileWriter freadme=new FileWriter(directory + "/README_DOT.txt");
        BufferedWriter readme=new BufferedWriter(freadme);
        int cnt = 0;
        String[] names = new String[graphs.size()];
        for (Map.Entry<Object, IdxGraph> entry : graphs.entrySet())
            names[cnt ++] = entry.getKey().toString();
        Arrays.sort(names);
        for (String name : names) {
            String filename = directory + "/" + toFilename(name) + ".dot";
            sb.append(filename + " ");
            FileWriter fwriter=new FileWriter(filename);
            BufferedWriter writer=new BufferedWriter(fwriter);
            IdxGraph g = graphs.get(name);
            writer.write(g.toDOT());
            writer.newLine();
            writer.close();
            fwriter.close();
        }
        readme.write("Install graphviz\nRun command:\n");
        if (cnt > 1) {
            readme.write("gvpack -u " + sb.toString() + "| dot -Tpdf -o" + directory + "/pogs.pdf");
            readme.newLine();
            readme.write("-- OR --");
            readme.newLine();
            for (String name : names) {

                readme.write("dot -Tpdf " + directory + "/" + toFilename(name) + ".dot -o" + directory + "/" + toFilename(name) + "_pog.pdf");
                readme.newLine();
            }
        } else
            readme.write("dot -Tpdf " + sb.toString() + "-opog.pdf");
        readme.newLine();
        readme.close();
        freadme.close();
    }

    public static final Character[] INVALID_FILENAME_CHARS = {'"', '*', ':', '<', '>', '?', '\\', '|', 0x7F};
    public static final String[] REPLACE_FILENAME_CHARS = {"_dquote_", "_exts_", "_colon_", "_lt_", "_gt_", "_question_", "_backslash_", "_pipe_", "_"};
    public static String toFilename(String proposed) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < proposed.length(); i ++) {
            char c = proposed.charAt(i);
            boolean valid = true;
            for (int j = 0; j < INVALID_FILENAME_CHARS.length; j ++) {
                if (c == INVALID_FILENAME_CHARS[j]) {
                    sb.append(REPLACE_FILENAME_CHARS[j]);
                    valid = false;
                    break;
                }
            }
            if (valid)
                sb.append(c);
        }
        return sb.toString();
    }

    public static void saveToLaTeX(String directory, Map<Object, IdxGraph> graphs) throws IOException, ASRException {
        StringBuilder sb = new StringBuilder();
        FileWriter freadme=new FileWriter(directory + "/README_LaTeX.txt");
        BufferedWriter readme=new BufferedWriter(freadme);
        int cnt = 0;
        String[] names = new String[graphs.size()];
        for (Map.Entry<Object, IdxGraph> entry : graphs.entrySet())
            names[cnt ++] = entry.getKey().toString();
        Arrays.sort(names);
        for (String name : names) {
            String filename = directory + "/" + toFilename(name) + ".tex";
            sb.append(filename + " ");
            FileWriter fwriter=new FileWriter(filename);
            BufferedWriter writer=new BufferedWriter(fwriter);
            IdxGraph g = graphs.get(name);
            writer.write(g.toLaTeXString("E^{" + name + "} = "));
            writer.newLine();
            writer.close();
            fwriter.close();
        }
        readme.write("Install pdflatex\nRun command:\n");
        if (cnt > 1) {
            readme.write("pdflatex -output-directory=" + directory + " -jobname=matrices" + " '\\documentclass[varwidth]{standalone}\\pagestyle{empty}\\begin{document}");
            for (String name : names)
                readme.write("\\input{" + directory + "/" + toFilename(name) + "}\\vspace{0.5in}");
            readme.write("\\end{document}'");
            readme.newLine();
            readme.write("-- OR --");
            readme.newLine();
            for (String name : names) {
                readme.write("pdflatex -output-directory=" + directory + " -jobname=" + toFilename(name) + "_mat '\\documentclass[varwidth]{standalone}\\pagestyle{empty}\\begin{document}\\input{" + directory + "/" + toFilename(name) + "}\\end{document}'");
                readme.newLine();
            }
        } else
            readme.write("pdflatex  -output-directory=" + directory + " -jobname=" + toFilename(names[0]) + "_mat '\\documentclass[varwidth]{standalone}\\pagestyle{empty}\\begin{document}\\input{" + directory + "/" + toFilename(names[0]) + "}\\end{document}'");
        readme.newLine();
        readme.close();
        freadme.close();
    }


    public void saveToMatrix(String filename) throws IOException {
        FileWriter fwriter=new FileWriter(filename);
        BufferedWriter writer=new BufferedWriter(fwriter);
        writer.write(toMatrixString());
        writer.close();
        fwriter.close();
    }

    public static void saveToMatrix(String directory, Map<Object, IdxGraph> graphs) throws IOException, ASRException {
        StringBuilder sb = new StringBuilder();
        FileWriter fwriter=new FileWriter(directory + "/ancestors.m");
        int cnt = 0;
        for (Map.Entry<Object, IdxGraph> entry : graphs.entrySet()) {
            sb.append(entry.getKey() + " = [\n");
            sb.append(entry.getValue().toMatrixString());
            sb.append("];\n");
            cnt += 1;
        }
        BufferedWriter writer=new BufferedWriter(fwriter);
        writer.write(sb.toString());
        writer.close();
        fwriter.close();
    }

/*    public static void saveToMatrix(String directory, IdxGraph... graphs) throws IOException, ASRException {
        StringBuilder sb = new StringBuilder();
        FileWriter fwriter=new FileWriter(directory + "/ancestors.m");
        int cnt = 0;
        for (IdxGraph g : graphs) {
            String name = "N" + Integer.toString(cnt);
            sb.append(name + " = [\n");
            sb.append(g.toMatrixString());
            sb.append("];\n");
            cnt += 1;
        }
        BufferedWriter writer=new BufferedWriter(fwriter);
        writer.write(sb.toString());
        writer.close();
        fwriter.close();
    } */
}
