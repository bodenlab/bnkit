package dat.pog;

import asr.ASRException;
import asr.ASRRuntimeException;
import asr.GRASP;
import bn.prob.EnumDistrib;
import dat.file.Utils;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;
import json.JSONTokener;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

/**
 * Definition of partial-order graph, which is essentially a directed index graph, with virtual start and end nodes.
 * Note that edges can be, but do not have to be instantiated.
 * {@link POGraph.StatusEdge} is defined internally.
 */
public class POGraph extends IdxEdgeGraph<POGraph.StatusEdge> {

    /**
     * Create partial order graph, which is a directed and terminated IdxGraph with some specific properties.
     * @param nNodes     number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
     */
    public POGraph(int nNodes) {
        super(nNodes, false, true);
    }

    public static POGraph createFromAdjacency(int[] startnodes, int[] endnodes, int[][] adjacency) {
        if (startnodes == null || endnodes == null || adjacency == null)
            throw new ASRRuntimeException("Invalid POG adjacency specification");
        if (startnodes.length < 1 && endnodes.length < 1)
            throw new ASRRuntimeException("Start and end indices for nodes are not provided");
        if (adjacency.length < 1)
            throw new ASRRuntimeException("Empty adjacency matrix provided");
        int N = adjacency.length;
        for (int[] a : adjacency) {
            if (a.length != N)
                throw new ASRRuntimeException("Invalid adjacency matrix provided: rows are uneven");
        }
        POGraph pog = new POGraph(N);
        for (int i : startnodes) {
            pog.addNode(i, new Node());
            pog.addEdge(-1, i);
        }
        for (int i : endnodes) {
            if (!pog.isNode(i))
                pog.addNode(i, new Node());
            pog.addEdge(i, N);
        }
        for (int from = 0; from < adjacency.length; from ++) {
            for (int to = 0; to < adjacency[from].length; to ++) {
                if (from != to) {
                    int from1 = from < to ? from : to;
                    int to1 = from < to ? to : from;
                    if (adjacency[from][to] > 0) { // edge present
                        if (!pog.isNode(from1))
                            pog.addNode(from1, new Node());
                        if (!pog.isNode(to1))
                            pog.addNode(to1, new Node());
                        StatusEdge e = pog.getEdge(from1, to1);
                        if (e == null)
                            pog.addEdge(from1, to1);
                        else
                            e.reciprocated = true;
                    }
                }
            }
        }
        return pog;
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    /**
     * Get string representation of instance
     * @return
     */
    public String toString() {
        return super.toString();
    }

    /**
     * Create a JSON representation of the instance
     * @return JSON object of this instance
     */
    public JSONObject toJSON() {
        JSONObject json = super.toJSON();
        // consider if anything needs to be added
        return json;
    }

    /**
     * Decode JSON object into an instance of POGraph
     * @param json
     * @return
     * @throws RuntimeException if there is a formatting problem
     */
    public static POGraph fromJSON(JSONObject json) {
        try {
            String version = json.optString("GRASP_version", null);
            Class datatype = null; //(Class)json.get("Datatype");
            if (version != null)
                if (!(version.equals(GRASP.VERSION) || version.equals("28-Aug-2022")))
                    throw new ASRRuntimeException("Invalid version: " + version);
            if (datatype != null)
                if (!(dat.pog.POGraph.class.isAssignableFrom(datatype)))
                    throw new ASRRuntimeException("Invalid datatype: " + datatype);
            int n = json.getInt("Size");
            boolean terminated = json.getBoolean("Terminated");
            boolean directed = json.getBoolean("Directed");
            if (terminated && directed) {
                POGraph g = new POGraph(n);
                g.useJSON(json); // from IdxGraph
                // POGraph specific fields...
                String edgetype = json.optString("Edgetype", null); // FIXME: could try to cast as Class first
                if (edgetype != null) {
                    JSONArray iarr = json.getJSONArray("Edgeindices");
                    int[][] idxs = new int[iarr.length()][];
                    for (int i = 0; i < iarr.length(); i++) {
                        JSONArray from_to = iarr.getJSONArray(i);
                        idxs[i] = new int[]{from_to.getInt(0), from_to.getInt(1)};
                    }
                    JSONArray earr = json.getJSONArray("Edges");
                    for (int i = 0; i < earr.length(); i++) {
                        JSONObject edge = earr.getJSONObject(i);
                        if (edgetype.equals(dat.pog.POGraph.StatusEdge.class.toString()))
                            g.addEdge(idxs[i][0], idxs[i][1], POGraph.StatusEdge.fromJSON(edge)); // StatusEdge goes here
                        else if (edgetype.equals(dat.pog.POGraph.BidirEdge.class.toString()))
                            g.addEdge(idxs[i][0], idxs[i][1], POGraph.BidirEdge.fromJSON(edge)); // BidirEdge goes here
                        else
                            g.addEdge(idxs[i][0], idxs[i][1]); // Other edges go here
                    }
                }
                return g;
            } else
                throw new ASRRuntimeException("POGraph format is wrong in JSON object: is not terminated or not directed");
        } catch (JSONException e) {
            throw new ASRRuntimeException("POGraph format is wrong in JSON object: " + e.getMessage());
        }
    }

    public static List<POGraph> loadFromJSON(String directory) throws IOException {
        Path filename = Paths.get(directory);
        BufferedReader reader = Files.newBufferedReader(filename, StandardCharsets.UTF_8);
        JSONTokener jtok = new JSONTokener(reader);
        /*
        String line = reader.readLine();
        StringBuilder sb = new StringBuilder();
        while (line != null) {
            String myline = line.trim();
            sb.append(myline);
            line = reader.readLine();
        }
        */
        JSONArray arr = new JSONArray(jtok);
        // expect an array of JSON objects
//        JSONArray arr = new JSONArray(sb.toString());
        List<POGraph> ret = new ArrayList<>();
        for (int i = 0; i < arr.length(); i ++) {
            JSONObject obj = arr.getJSONObject(i);
            ret.add(POGraph.fromJSON(obj));
        }
        reader.close();
        return ret;
    }


    /**
     * Determine the node indices that are ahead (forward direction). Note that it terminates at last node, so will not return the end marker.
     * @param idx current node, use -1 for retrieving the start nodes
     * @return the indices for nodes accessible directly forward from the current node; array is empty if no indices
     */
    public int[] getForward(int idx) {
        int[] next = getNodeIndices(idx, true);
        return next;
    }

    /**
     * Determine the indices for start nodes (same as providing -1 as index in {@link POGraph#getForward(int)})
     * @return the indices for nodes that act as start nodes of the POG; the array is empty if no start nodes
     */
    public int[] getForward() {
        return getForward(-1);
    }

    /**
     * Determine the node indices that are behind (backward direction). Note that it completes at start node, and so does not return the start marker (-1).
     * @param idx current node, use N for retrieving the end nodes, where N is the max number of positions in the POG
     * @return the indices for nodes accessible directly forward from the current node; the array is empty if there are no indices
     */
    public int[] getBackward(int idx) {
        return getNodeIndices(idx, false);
    }

    /**
     * Determine the indices for end nodes (same as providing N as index in {@link POGraph#getBackward(int)} where N is the max number of positions in the POG)
     * @return the indices for nodes that act as end nodes of the POG; the array is empty if no end nodes
     */
    public int[] getBackward() {
        return getBackward(this.nNodes);
    }

    /**
     * Determine a topological ordering of the nodes; note that more than one topological order is likely to exist for bi- and multi-furcating POGs
     *
     * @return an array with the indices of the POG in a topological order
     * TODO: confirm that this implementation can be inherited from IdxGraph
     */
    public int[] getTopologicalOrder() {
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
            int[] next = this.getForward(current);
            for (int idx : next) {
                if (!added[idx]) {
                    inedges[idx] = getCardinality(idx, false) + (isStartNode(idx) ? 1 : 0);
                    added[idx] = true;
                }
                if (inedges[idx] == 1)
                    ret[add_ptr ++] = idx;
                if (inedges[idx] > 1)
                    inedges[idx] --;
            }
        }
        return ret;
    }


    /**
     * Retrieve all edges in the POG, including edges leading to start nodes, and those leading from end nodes to the terminus.
     * @return an array of pairs of int
     */
    public int[][] getPairs() {
        int[][] pairs = new int[getEdgeCount()][];
        int cnt = 0;
        for (int from = -1; from < nNodes; from++) {
            for (int to : getNodeIndices(from, true)) // all forward directed edges
                pairs[cnt++] = new int[]{from, to};
            if (isEndNode(from)) // check if node is also terminal
                pairs[cnt++] = new int[]{from, maxsize()};
        }
        return pairs;
    }

    /**
     * Establish if there is a path from idx1 to idx2.
     * Uses forward recursion.
     * @param idx1 source node index (can use -1)
     * @param idx2 target node index (can use N, where N is the number of nodes in the POG)
     * @return true if there is a path, false otherwise
     */
    public boolean isPath(int idx1, int idx2) {
        if (idx1 == idx2) // by def, a node has a path to itself
            return true;
        if (!isNode(idx1) && idx1 != -1) // the source node must either be the start terminal or a valid node (not just an index)
            return false;
        if (!isNode(idx2) && idx2 != nNodes) // the target node must either be the end terminal or a valid node (not just an index)
            return false;
        int[] forw = getForward(idx1);
        for (int to : forw) {
            if (isPath(to, idx2))
                return true;
        }
        return idx2 == nNodes ? isEndNode(idx1) : false; // if we are searching for end marker
    }

    /**
     * Is this POG contiguous?
     * @return
     */
    public boolean isContiguous() {
        return isPath(-1, nNodes);
    }
    /**
     * Get complete paths with reciprocated edges, and non-reciprocated edges used only if required to complete a path.
     * @return
     */
    public Set<POGEdge> getReciprocated() {
        Set<POGEdge> ret = new HashSet<>();
        throw new RuntimeException("Not implemented");
//        return ret;
    }

    public static int SUPPORTED_PATH_DIJKSTRA = 0;
    public static int SUPPORTED_PATH_ASTAR = 1;
    public static int SUPPORTED_PATH_DEFAULT = SUPPORTED_PATH_DIJKSTRA; // can be overridden in asr.GRASP

    /**
     * Perform search on the POG based on currently set edge weights and status
     * @return indices to the nodes that make up the "consensus" path
     */
    public int[] getMostSupported() {
        return getMostSupported(SUPPORTED_PATH_DEFAULT);
    }

    /**
     * Perform search on the POG based on currently set edge weights and status
     * @return indices to the nodes that make up the "consensus" path
     */
    public int[] getMostSupported(int MODE) {
        if (MODE == SUPPORTED_PATH_ASTAR) {
            AStarSearch search = new AStarSearch(this);
            return search.getOnePath();
        } else { // SUPPORTED_PATH_DIJKSTRA
            DijkstraSearch search = new DijkstraSearch(this);
            return search.getOnePath();
        }
    }

    /** Variable to hold search depths for getDepths */
    private int[] searchdepths;

    /**
     * Helper function to determine search depths of nodes in POG
     * @param toEnd if true, "ground zero" is the end, else the start of the POG
     * @param currdepth
     * @param queue
     * @return
     */
    private List<Integer> findDepthsOf(boolean toEnd, int currdepth, List<Integer> queue) {
        List<Integer> next = new ArrayList<>();
        for (int i : queue) {
            if (searchdepths[i] == 0 || searchdepths[i] > currdepth) {
                searchdepths[i] = currdepth;
                if (toEnd) {
                    for (int j : this.getBackward(i))
                        next.add(j);
                } else { // forward-looking
                    for (int j : this.getForward(i))
                        next.add(j);
                }
            }
        }
        return next;
    }

    /**
     * Determine the search depths from the end or the start (based on given param).
     * @param toEnd if true, "ground zero" is the end, else the start of the POG
     * @return an array with indices for all nodes in the POG, excluding terminal nodes, where each element indicating the search depth
     */
    public int[] getDepths(boolean toEnd) {
        searchdepths = new int[this.maxsize()];
        List<Integer> queue = new ArrayList<>();
        if (toEnd) {
            for (int j : this.getBackward())
                queue.add(j);
        } else {
            for (int j : this.getForward())
                queue.add(j);
        }
        int currdepth = 1;
        while (queue.size() > 0) {
            queue = findDepthsOf(toEnd, currdepth, queue);
            currdepth += 1;
        }
        return searchdepths;
    }

    /**
     * Retrieve all pairs of jump edges, i.e. edges that skip at least one position
     * @return
     */
    public Set<int[]> getIndels() {
        Set<int[]> ivset = new HashSet<>();
        for (int from = -1; from < nNodes; from++) {
            for (int to : getNodeIndices(from, true))
                if (Math.abs(to - from) > 1)
                    ivset.add(new int[] {from, to});
            if (isEndNode(from))
                if (maxsize() - from > 1)
                    ivset.add(new int[] {from, maxsize()});
        }
        return ivset;
    }

    static public Boolean GAP_STATUS_PRESENT = true;
    static public Boolean GAP_STATUS_ABSENT = false;
    static public Boolean GAP_STATUS_UNKNOWN = null;

    /**
     * Check for this POG if the query indel is present, if it is absent, or if it is permissible.
     * @param from start index of indel
     * @param to end index of indel
     * @return
     */
    public Boolean getSimpleGapCode(int from, int to) {
        if (isNode(from) || from == -1) {
            // check if this POG has a jump that matches exactly
            if (this.isEdge(from, to))
                return GAP_STATUS_PRESENT;
            // now we know it is a "jump", so check containment status
            int[] ahead = this.getForward(from); // inspect jumps forward relative the query start position
//            if ((ahead.length > 0 || from == -1) && to - from == 1) // query is not a "jump"" at all, but intersects with one, so must be ...
//                return GAP_STATUS_ABSENT;
            for (int i = 0; i < ahead.length; i++) {
                if (ahead[i] < to) // the gap in this POG ends BEFORE that of the query
                    return GAP_STATUS_ABSENT;
            }
            return GAP_STATUS_UNKNOWN; // the gap ends AFTER that of the query
        } else if (isNode(to) || to == maxsize()) { // the gap starts BEFORE "from"
            int[] before = this.getBackward(to);
//            if ((before.length > 0 || to == maxsize()) && to - from == 1) // query is not a "jump" at all, but intersects with one, so must be ...
//                return GAP_STATUS_ABSENT;
            for (int i = 0; i < before.length; i++) {
                if (before[i] > from) // the gap in this POG begins AFTER that of the query
                    return GAP_STATUS_ABSENT;
            }
            return GAP_STATUS_UNKNOWN;
 //       } else if (to - from == 1) { // not a "jump"" at all, and just internal content
 //           return GAP_STATUS_ABSENT;
        } else { // the query gap started BEFORE and completed AFTER, but there could be internal content
            for (int i = from + 1; i < to - 1; i ++) {
                if (isNode(i))
                    return GAP_STATUS_ABSENT;
            }
            return GAP_STATUS_UNKNOWN;
        }
    }

    /**
     * Create a POG from an EdgeMap
     * @param nPos
     * @param emap
     * @return
     */
    public static POGraph createFromEdgeMap(int nPos, EdgeMap emap) {
        POGraph pog = null;
        pog = new POGraph(nPos);
        for (POGEdge edge : emap.getEdges()) { // go through all edges, doubling up if bi-directional
            synchronized (pog) {
                if (edge.idx1 != -1 && !pog.isNode(edge.idx1))
                    pog.addNode(edge.idx1, new Node());
                if (edge.idx2 != nPos && !pog.isNode(edge.idx2))
                    pog.addNode(edge.idx2, new Node());
                pog.addEdge(edge.idx1, edge.idx2, new StatusEdge(emap.isReciprocated(edge)));
            }
        }
        return pog;
    }

    /**
     * Create a POG from a Directed
     * @param nPos
     * @param emap
     * @return
     */
    public static POGraph createFromEdgeMap(int nPos, EdgeMap.Directed emap) {
        POGraph pog = null;
        pog = new POGraph(nPos);
        for (POGEdge edge : emap.getEdges()) { // go through all edges, doubling up if bi-directional
            synchronized (pog) {
                if (edge.idx1 != -1 && !pog.isNode(edge.idx1))
                    pog.addNode(edge.idx1, new Node());
                if (edge.idx2 != nPos && !pog.isNode(edge.idx2))
                    pog.addNode(edge.idx2, new Node());
                pog.addEdge(edge.idx1, edge.idx2, new BidirEdge(emap.isForward(edge), emap.isBackward(edge)));
            }
        }
        return pog;
    }

    /**
     * Check if node sticks out in one direction
     * @param index node
     * @param direction_start if true, sticking out without an edge that could lead from start; else, sticking out without an edge that could lead to the end
     * @return true if the node sticks out
     */
    private synchronized boolean nibble(int index, boolean direction_start) {
        if (direction_start && index > -1) { // all nodes (except start node) must have backward links
            if (!isStartNode(index)) {
                int[] backw = getBackward(index);
                return (backw.length == 0);
            }
        }
        if (!direction_start && index < nNodes) { // all nodes (except terminal node) must have forward links
            if (!isEndNode(index)) {
                int[] forw = getForward(index);
                return (forw.length == 0);
            }
        }
        return false;
    }

    /**
     * Remove nodes that "stick-out" in either direction.
     * TODO: consider option to not calculate topological order (may not be necessary in most cases)
     * @return number of nodes removed in total
     */
    public int nibble() {
        int cnt = 0; // how many we removed
        // first work from the start
        int[] idxs = getTopologicalOrder();
        for (int i = 0; i < idxs.length; i ++) {
            int idx = idxs[i];
            if (idx >= 0 && idx < nNodes) {
                if (nodes[idx] != null) {
                    if (nibble(idx, true)) {
                        cnt += 1;
                        removeNode(idx);
                    }
                }
            }
        }
        // next, work from the end
        for (int i = idxs.length - 1; i >= 0; i --) {
            int idx = idxs[i];
            if (idx >= 0 && idx < nNodes) {
                if (nodes[idx] != null) {
                    if (nibble(idx, false)) {
                        cnt += 1;
                        removeNode(idx);
                    }
                }
            }
        }
        return cnt;
    }

    /**
     * Retrieve indices in POG that "stick-out" i.e. nodes that protrude (in either direction) forming anomalous start and end nodes.
     * Indices that have direction start-to-end are 1, 2, ... N; Indices from end-to-start are 0, -1, -2, ..., -(N-1).
     * TODO: consider option to not calculate topological order (may not be necessary in most cases)
     * @return columns
     */
    public Set<Integer> getProtruding() {
        // first work from the start
        int[] idxs = getTopologicalOrder();
        Set<Integer> cols = new HashSet<>();
        for (int i = 0; i < idxs.length; i ++) {
            int idx = idxs[i];
            if (idx >= 0 && idx < nNodes) {
                if (nodes[idx] != null) {
                    if (nibble(idx, false))
                        cols.add(idx);
                }
            }
        }
        // next, work from the end
        for (int i = idxs.length - 1; i >= 0; i --) {
            int idx = idxs[i];
            if (idx >= 0 && idx < nNodes) {
                if (nodes[idx] != null) {
                    if (nibble(idx, true))
                        cols.add(-idx); // note that index is negative to mark direction end-->start
                }
            }
        }
        return cols;
    }



    /**
     * Create a POG from an array of edge indices.
     * This is not uncomplicated as each INDEL needs to be considered in a broader context to form a valid, start-to-end POG.
     * @param nNodes number of nodes in new POG
     * @param optional edges that CAN be TRUE, but are precluded by "definitive"
     * @param definitive edges that MUST be TRUE
     * @return new POG
     */
    public static POGraph createFromEdgeIndicesWithoutDeadends(int nNodes, int[][] optional, int[][] definitive) {
        POGraph pog = new POGraph(nNodes);
        for (int[] edge : optional) {
            if (edge.length != 2)
                throw new RuntimeException("Invalid edge specification");
            for (int idx : edge)
                if (!(idx == -1 || idx == pog.nNodes || pog.isNode(idx)))
                    pog.addNode(idx, new StatusNode());
                pog.addEdge(edge[0], edge[1], new StatusEdge(true));

        }
        try {
            pog.saveToDOT(File.createTempFile("POG_", ".dot").getAbsolutePath());
        } catch (IOException e) {
            e.printStackTrace();
        }
        // find and cripple precluded paths, for each definite edge
        for (int[] edge : definitive) {

            POGEdge[] ends = pog.findEndsOfPrecludedEdges(edge[0], edge[1]);
            for (POGEdge end : ends)
                pog.removeEdge(end.idx1, end.idx2);

        }

        // remove all dead-ends
        for (int idx : pog.getForward())
            visitForward(pog, idx);

        for (int idx : pog.getBackward())
            visitBackward(pog, idx);



        for (int idx = 0; idx < pog.nNodes; idx ++) {
            if (pog.isNode(idx)) {
                StatusNode snode = (StatusNode) pog.getNode(idx);
                if (!snode.isStatus(StatusNode.FORWARD)) { // was not reached from start going forward
                    pog.removeNode(idx);
                } else if (!snode.isStatus(StatusNode.BACKWARD)) { // was not reached from end going backward
                    pog.removeNode(idx);
                }
            }
        }


        return pog;
    }

    /**
     * Create a POG from an array of edge indices.
     * This is not uncomplicated as each INDEL needs to be considered in a broader context to form a valid, start-to-end POG.
     * @param nNodes number of nodes in new POG
     * @param optional edges that CAN be TRUE, but are precluded by "definitive"
     * @param definitive edges that MUST be TRUE
     * @return new POG
     */
    public static POGraph createFromEdgeIndicesWithoutDeadends(int nNodes, Set<int[]> optional, Set<int[]> definitive) {
        POGraph pog = new POGraph(nNodes);
        for (int[] edge : optional) {
            if (edge.length != 2)
                throw new RuntimeException("Invalid edge specification");
            for (int idx : edge)
                if (!(idx == -1 || idx == pog.nNodes || pog.isNode(idx)))
                    pog.addNode(idx, new StatusNode());
            pog.addEdge(edge[0], edge[1], new StatusEdge(false));
        }
        for (int[] edge : definitive) {
            if (edge.length != 2)
                throw new RuntimeException("Invalid edge specification");
            for (int idx : edge)
                if (!(idx == -1 || idx == pog.nNodes || pog.isNode(idx)))
                    pog.addNode(idx, new StatusNode());
            pog.addEdge(edge[0], edge[1], new StatusEdge(true));
        }
        // find and cripple precluded paths, for each definite edge
        for (int[] edge : definitive) {
            POGEdge[] ends = pog.findEndsOfPrecludedEdges(edge[0], edge[1]);
            for (POGEdge end : ends)
                pog.removeEdge(end.idx1, end.idx2);
        }
        // remove all dead-ends
        for (int idx : pog.getForward())
            visitForward(pog, idx);
        for (int idx : pog.getBackward())
            visitBackward(pog, idx);
        for (int idx = 0; idx < pog.nNodes; idx ++) {
            if (pog.isNode(idx)) {
                StatusNode snode = (StatusNode) pog.getNode(idx);
                if (!snode.isStatus(StatusNode.FORWARD)) { // was not reached from start going forward
                    pog.removeNode(idx);
                } else if (!snode.isStatus(StatusNode.BACKWARD)) { // was not reached from end going backward
                    pog.removeNode(idx);
                }
            }
        }
        return pog;
    }

    public POGEdge[] findEndsOfPrecludedEdges(int unique_idx, int dest_idx) {
        Set<POGEdge> visited = new HashSet<>();
        Set<Integer> visitedNode = new HashSet();
        int[] next = getForward(unique_idx);
        if (isEndNode(unique_idx)) { // special case...
            int[] nnext = Arrays.copyOf(next, next.length + 1);
            nnext[next.length] = nNodes;
            next = nnext;
        }
        for (int pos : next) {
            if (pos <= dest_idx)
                visited.add(new POGEdge(unique_idx, pos));
        }
        for (int pos : next) {
            if (pos <= dest_idx && !(visitedNode.contains(pos)))
                carryUniqueForward(this, unique_idx, dest_idx, pos, visited, visitedNode);
        }
        Set<POGEdge> preclude = new HashSet<>();
        for (POGEdge pair : visited) {
            if (pair.idx2 == dest_idx && pair.idx1 != unique_idx)
                preclude.add(pair);
        }
        POGEdge[] precarr = new POGEdge[preclude.size()];
        preclude.toArray(precarr);
        return precarr;
    }

    private static void carryUniqueForward(POGraph pog, int unique_idx, int dest_idx, int pos, Set<POGEdge> visited, Set<Integer> visitedNode) {
        boolean all_unique = true;
        if (pog.isNode(pos)) {
            for (int previdx : pog.getBackward(pos)) {
                // make sure {previdx, pos} is in there; if not, {pos, nextidx} is not a unique path to dest_idx
                boolean edge_checked = false;

                POGEdge backwardsEdge = new POGEdge(previdx, pos);

                if (visited.contains(backwardsEdge)){
                    edge_checked = true;
                }
                if (!edge_checked)
                    all_unique = false;
            }
            int[] nextidxs = pog.getForward(pos);
            if (pog.isEndNode(pos)) { // special case...
                int[] nnextidxs = Arrays.copyOf(nextidxs, nextidxs.length + 1);
                nnextidxs[nextidxs.length] = pog.nNodes;
                nextidxs = nnextidxs;
            }
            for (int nextidx : nextidxs) {

                if (nextidx <= dest_idx && all_unique) {
                    visited.add(new POGEdge(pos, nextidx));
                    visitedNode.add(pos);
                }
            }
            for (int nextidx : nextidxs) {
                if (nextidx <= dest_idx && !(visitedNode.contains(nextidx)))
                    carryUniqueForward(pog, unique_idx, dest_idx, nextidx, visited, visitedNode);
            }
        }
    }

    private static void visitForward(POGraph pog, int pos) {
        if (pog.isNode(pos)) {
            StatusNode snode = (StatusNode)pog.getNode(pos);
            if (!snode.isStatus(StatusNode.FORWARD)) { // not visited before going FORWARD
                snode.flipStatus(StatusNode.FORWARD);  // mark as visited going FORWARD
                int[] next = pog.getForward(pos);
                for (int idx : next)
                    visitForward(pog, idx);
            }
        } else if (pos == -1) {
            int[] next = pog.getForward();
            for (int idx : next)
                visitForward(pog, idx);
        }
    }

    private static void visitBackward(POGraph pog, int pos) {
        if (pog.isNode(pos)) {
            StatusNode snode = (StatusNode)pog.getNode(pos);
            if (!snode.isStatus(StatusNode.BACKWARD)) {
                snode.flipStatus(StatusNode.BACKWARD);
                int[] next = pog.getBackward(pos);
                for (int idx : next)
                    visitBackward(pog, idx);
            }
        } else if (pos == pog.nNodes) {
            int[] next = pog.getBackward();
            for (int idx : next)
                visitBackward(pog, idx);
        }
    }


    public void decorateNodes(Object[] states) {
        if (states == null)
            throw new ASRRuntimeException("Attempting to assign uninstantiated states to POG");
        if (states.length != this.nNodes)
            throw new ASRRuntimeException("Attempting to assign invalid states to POG");
        nodes = SymNode.toArray(states);
    }

    public void decorateNodes(EnumDistrib[] distribs) {
        if (distribs == null)
            throw new ASRRuntimeException("Attempting to assign uninstantiated distributions to POG");
        if (distribs.length != this.nNodes)
            throw new ASRRuntimeException("Attempting to assign invalid distributions to POG");
        nodes = EnumNode.toArray(distribs);
    }

    /**
     * Define the standard type of node for POGs that allow for tracing POG consistency
     */
    public static class StatusNode extends Node {
        public static int FORWARD = 0x01;
        public static int BACKWARD = 0x02;
        int status = 0;
        int[] mindist = null;
        public StatusNode() {
        }
        public int getStatus() {
            return status;
        }
        public void setStatus(int status) {
            this.status = status;
        }
        public void flipStatus(int status) {
            this.status |= status;
        }
        public boolean isStatus(int status) {
            return (this.status & status) > 0;
        }
        public void setDistance(boolean toEnd, int dist) {
            if (mindist == null) {
                mindist = new int[2];
                mindist[!toEnd ? 0 : 1] = -1;
            }
            mindist[toEnd ? 0 : 1] = dist;
        }
        public int getDistance(boolean toEnd) {
            if (mindist == null)
                return -1;
            return mindist[toEnd ? 0 : 1];
        }
    }

    /**
     * Define the standard type of edge for POGs as we use them in GRASP, with weight and support for "bi-directionality" status
     */
    public static class StatusEdge extends Edge implements WeightedEdge {
        private double weight = 0;  // default
        private boolean reciprocated = false;
        public StatusEdge(boolean reciprocated) {
            this.reciprocated = reciprocated;
        }
        public StatusEdge(boolean reciprocated, int weight) {
            this.reciprocated = reciprocated;
            this.weight = weight;
        }
        public StatusEdge(boolean reciprocated, double weight) {
            this.reciprocated = reciprocated;
            this.weight = weight;
        }
        public boolean getReciprocated() {
            return reciprocated;
        }
        public void setWeight(int weight) {
            this.weight = weight;
        }
        public void setWeight(double weight) {
            this.weight = weight;
        }
        @Override
        public double getWeight() {
            return weight;
        }
        @Override
        public String toString() {
            return "<" + (getReciprocated()?"R":":") + String.format("%3.1f", getWeight()) + ">";
        }
        @Override
        public String toDOT() {
            return  ((label == null) ? "" : ("label=\"" + label + "\",")) +
                    ("penwidth=" + (getWeight() + 1) + ",") +
                    ((fontname == null) ? "" : ("fontname=\"" + fontname + "\",")) +
                    (reciprocated ? "" : ("style=\"dashed\""));
        }

        /**
         * Convert instance to JSON
         * @return
         */
        @Override
        public JSONObject toJSON() {
            JSONObject edge = super.toJSON();
            edge.put("Recip", getReciprocated());
            edge.put("Weight", getWeight());
            return edge;
        }

        /**
         * Convert JSON to instance
         * @param json
         * @return
         */
        public static StatusEdge fromJSON(JSONObject json) {
            String label = json.optString("Label", null);
            Boolean recip = json.optBoolean("Recip", false);
            Double weight = json.optDouble("Weight", 0);
            StatusEdge e = new StatusEdge(recip, weight);
            if (label != null)
                e.setLabel(label);
            return e;
        }

        @Override
        public Object[] toObjectArray() {
            if (label != null)
                return new Object[] {getReciprocated(), getWeight(), getLabel()};
            else
                return new Object[] {getReciprocated(), getWeight()};
        }
        public static StatusEdge fromObjectArray(Object[] input) {
            if (input.length >= 1) {
                try {
                    StatusEdge e = new StatusEdge((Boolean) input[0]);
                    if (input.length > 1) {
                        try {
                            e.setWeight((Double) input[1]);
                        } catch (ClassCastException e1) {
                            e.setWeight((Integer) input[1]);
                        }
                    }
                    if (input.length > 2)
                        e.setLabel(input[2].toString());
                    return e;
                } catch (ClassCastException e2) {
                    throw new RuntimeException("Invalid format for StatusEdge: " + e2.getMessage());
                }
            }
            return null;
        }
    }

    /**
     * Define the standard type of edge for POGs as we use them in GRASP, with weight and support for "bi-directionality" status
     */
    public static class BidirEdge extends StatusEdge {
        private boolean forward, backward;
        public BidirEdge(boolean forward, boolean backward) {
            super(forward && backward);
            this.forward = forward;
            this.backward = backward;
        }

        public boolean isForward() {
            return forward;
        }

        public boolean isBackward() {
            return backward;
        }

        /**
         * Convert instance to JSON
         * @return
         */
        @Override
        public JSONObject toJSON() {
            JSONObject edge = super.toJSON();
            edge.put("Forward", isForward());
            edge.put("Backward", isBackward());
            return edge;
        }

        /**
         * Convert JSON to instance
         * @param json
         * @return
         */
        public static StatusEdge fromJSON(JSONObject json) {
            String label = json.optString("Label", null);
            Boolean forward = json.optBoolean("Forward", false);
            Boolean backward = json.optBoolean("Backward", false);
            Double weight = json.optDouble("Weight", 0);
            BidirEdge e = new BidirEdge(forward, backward);
            e.setWeight(weight);
            if (label != null)
                e.setLabel(label);
            return e;
        }

        @Override
        public Object[] toObjectArray() {
            if (label != null)
                return new Object[] {forward, backward, getWeight(), getLabel()};
            else
                return new Object[] {forward, backward, getWeight()};
        }
        public static BidirEdge fromObjectArray(Object[] input) {
            if (input.length >= 2) {
                try {
                    BidirEdge e = new BidirEdge((Boolean) input[0], (Boolean) input[1]);
                    if (input.length > 2) {
                        try {
                            e.setWeight((Double) input[2]);
                        } catch (ClassCastException e1) {
                            e.setWeight((Integer) input[2]);
                        }
                    }
                    if (input.length > 3)
                        e.setLabel(input[3].toString());
                    return e;
                } catch (ClassCastException e2) {
                    throw new RuntimeException("Invalid format for BidirEdge: " + e2.getMessage());
                }
            }
            return null;
        }

    }

}
