package dat.pog;

import asr.ASRRuntimeException;
import bn.prob.EnumDistrib;

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

    /**
     * Get string representation of instance
     * @return
     */
    public String toString() {
        return super.toString();
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
     * @return an array with the indices of the POG in a topological order
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
            for (int to : getNodeIndices(from, true))
                pairs[cnt++] = new int[]{from, to};
            if (isEndNode(from))
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
        if (idx1 == idx2)
            return true;
        if (!isNode(idx1) && idx1 != -1)
            return false;
        if (!isNode(idx2) && idx2 != size())
            return false;
        int[] forw = getForward(idx1);
        for (int to : forw) {
            if (isPath(to, idx2))
                return true;
        }
        return idx2 == nNodes ? isEndNode(idx1) : false; // if we are searching for end marker
    }


    /**
     * Get complete paths with reciprocated edges, and non-reciprocated edges used only if required to complete a path.
     * @return
     */
    public Set<POGEdge> getReciprocated() {
        Set<POGEdge> ret = new HashSet<>();

        return ret;
    }

    /**
     * Perform A* search on the POG
     * @return
     */
    public GraphSearch getMostSupported() {
        GraphSearch astar = new GraphSearch(this);
        //astar.runAStar();
        return astar;
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
            // check containment status
            int[] ahead = this.getForward(from); // inspect jumps forward relative the query start position
            for (int i = 0; i < ahead.length; i++) {
                if (ahead[i] < to) // the gap in this POG ends BEFORE that of the query
                    return GAP_STATUS_ABSENT;
            }
            return GAP_STATUS_UNKNOWN; // the gap ends AFTER that of the query
        } else if (isNode(to) || to == maxsize()) { // the gap starts BEFORE "from"
            int[] before = this.getBackward(to);
            for (int i = 0; i < before.length; i++) {
                if (before[i] > from) // the gap in this POG begins AFTER that of the query
                    return GAP_STATUS_ABSENT;
            }
            return GAP_STATUS_UNKNOWN;
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
        boolean tidy = true; // treat reciprocated edges as definitive, and non-reciprocated as optional
        if (tidy) {
            pog = new POGraph(nPos);
            Set<POGEdge> optional = new HashSet<>();
            for (POGEdge edge : emap.getRecip()) {
                if (edge.idx1 != -1 && !pog.isNode(edge.idx1))
                    pog.addNode(edge.idx1, new Node());
                if (edge.idx2 != nPos && !pog.isNode(edge.idx2))
                    pog.addNode(edge.idx2, new Node());
                pog.addEdge(edge.idx1, edge.idx2, new StatusEdge(true));
            }
            for (POGEdge edge : emap.getEdges()) {
                if (!emap.isReciprocated(edge)) {
                    if (pog.isNode(edge.idx1) || pog.isNode(edge.idx2)) {
                        if (!pog.isPath(edge.idx1, edge.idx2))
                            optional.add(edge);
                    }
                }
            }
            for (POGEdge edge : optional) {
                if (edge.idx1 != -1 && !pog.isNode(edge.idx1))
                    pog.addNode(edge.idx1, new Node());
                if (edge.idx2 != nPos && !pog.isNode(edge.idx2))
                    pog.addNode(edge.idx2, new Node());
                pog.addEdge(edge.idx1, edge.idx2, new StatusEdge(false));
            }
            // TODO: remove deadends
        } else {
            pog = new POGraph(nPos);
            for (POGEdge edge : emap.getEdges()) {
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
     * Create a POG from an array of edge indices.
     * This is not uncomplicated as each INDEL needs to be considered in a broader context to form a valid, start-to-end POG.
     * @param nNodes number of nodes in new POG
     * @param optional edges that CAN be TRUE
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
            pog.addEdge(edge[0], edge[1]);
        }
        // find and cripple precluded paths, for each definite edge
        for (int[] edge : definitive) {
            int[][] ends = pog.findEndsOfPrecludedEdges(edge[0], edge[1]);
            for (int[] end : ends)
                pog.removeEdge(end[0], end[1]);
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
     * @param optional edges that CAN be TRUE
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
            int[][] ends = pog.findEndsOfPrecludedEdges(edge[0], edge[1]);
            for (int[] end : ends)
                pog.removeEdge(end[0], end[1]);
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

    public int[][] findEndsOfPrecludedEdges(int unique_idx, int dest_idx) {
        Set<int[]> visited = new HashSet<>();
        int[] next = getForward(unique_idx);
        if (isEndNode(unique_idx)) { // special case...
            int[] nnext = Arrays.copyOf(next, next.length + 1);
            nnext[next.length] = nNodes;
            next = nnext;
        }
        for (int pos : next) {
            if (pos <= dest_idx)
                visited.add(new int[] {unique_idx, pos});
        }
        for (int pos : next) {
            if (pos <= dest_idx)
                carryUniqueForward(this, unique_idx, dest_idx, pos, visited);
        }
        Set<int[]> preclude = new HashSet<>();
        for (int[] pair : visited) {
            if (pair[1] == dest_idx && pair[0] != unique_idx)
                preclude.add(pair);
        }
        int[][] precarr = new int[preclude.size()][];
        preclude.toArray(precarr);
        return precarr;
    }

    private static void carryUniqueForward(POGraph pog, int unique_idx, int dest_idx, int pos, Set<int[]> visited) {
        boolean all_unique = true;
        if (pog.isNode(pos)) {
            for (int previdx : pog.getBackward(pos)) {
                // make sure {previdx, pos} is in there; if not, {pos, nextidx} is not a unique path to dest_idx
                boolean edge_checked = false;
                for (int[] edge : visited) {
                    if (edge[0] == previdx && edge[1] == pos) {
                        edge_checked = true;
                        break;
                    }
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
                if (nextidx <= dest_idx && all_unique)
                    visited.add(new int[]{pos, nextidx});
            }
            for (int nextidx : nextidxs) {
                if (nextidx <= dest_idx)
                    carryUniqueForward(pog, unique_idx, dest_idx, nextidx, visited);
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
            throw new ASRRuntimeException("Retrieving states for sequence that has not been inferred");
        nodes = SymNode.toArray(states);
    }

    public void decorateNodes(EnumDistrib[] distribs) {
        if (distribs == null)
            throw new ASRRuntimeException("Retrieving distributions for sequence that has not been inferred");
        nodes = EnumNode.toArray(distribs);
    }

    public static class StatusNode extends Node {
        public static int FORWARD = 0x01;
        public static int BACKWARD = 0x02;
        int status = 0;
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
    }


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
    }

}
