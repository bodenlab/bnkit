package dat.pog;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class POGraph extends IdxGraph {

    /**
     * Create partial order graph, which is a directed and terminated IdxGraph with some specific properties.
     * @param nNodes     number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
     */
    public POGraph(int nNodes) {
        super(nNodes, false, true);
    }

    public String toString() {
        return super.toString();
    }

    /**
     * Determine the node indices that are ahead (forward direction)
     * @param idx current node, use -1 for retrieving the start nodes
     * @return the indices for nodes accessible directly forward from the current node
     */
    public int[] getForward(int idx) {
        int[] next = getNodeIndices(idx, true);
        return next;
    }

    /**
     * Determine the indices for start nodes
     * @return the indices for nodes that act as start nodes of the POG
     */
    public int[] getForward() {
        return getForward(-1);
    }

    /**
     * Determine the node indices that are behind (backward direction)
     * @param idx current node, use N for retrieving the start nodes, where N is the nu
     * @return the indices for nodes accessible directly forward from the current node
     */
    public int[] getBackward(int idx) {
        return getNodeIndices(idx, false);
    }

    public int[] getBackward() {
        return getBackward(this.nNodes);
    }

    public int[][] getEdges() {
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
            return GAP_STATUS_UNKNOWN; // the gap end AFTER that of the query
        } else if (isNode(to) || to == maxsize()) { // the gap has started BEFORE "from"
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
     * Create a POG from an array of edge indices.
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

}
