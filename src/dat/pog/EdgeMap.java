package dat.pog;

import java.util.*;

/**
 * Edge container class to build a POGraph.
 * Each member is identified by its source and target indices (by reference to the original alignment).
 */
public class EdgeMap {

    private final Map<Integer, Set<Integer>> links = new HashMap<>(); // map that will match every idx with those that it connects with
    private final Set<POGEdge> edges = new HashSet<>();
    private final Set<POGEdge> recip = new HashSet<>(); // reciprocated edges only

    public EdgeMap() {
    }

    public Set<POGEdge> getEdges() {
        return edges;
    }

    public Set<POGEdge> getRecip() { return recip; }

    public boolean isReciprocated(POGEdge edge) {
        return recip.contains(edge);
    }

    /**
     * Add an edge to the set of edges (idx1 --> idx2) associated with the POG.
     * When the edge is added a second time, it is automatically labeled reciprocated.
     * @param idx1 first index
     * @param idx2 second index
     */
    public POGEdge add(int idx1, int idx2) {
        POGEdge e = new POGEdge(idx1, idx2);
        if (edges.contains(e))
            recip.add(e);
        else
            edges.add(e);
        Set<Integer> mylinks = links.get(idx1);
        if (mylinks == null) {
            mylinks = new HashSet<>();
            links.put(idx1, mylinks);
        }
        mylinks.add(idx2);
        mylinks = links.get(idx2);
        if (mylinks == null) {
            mylinks = new HashSet<>();
            links.put(idx2, mylinks);
        }
        mylinks.add(idx1);
        return e;
    }

    public void remove(POGEdge pe) {
        edges.remove(pe);
        recip.remove(pe);
    }
    /**
     * Return the total number of edges.
     * @return number of edges (reciprocated or not)
     */
    public int size() {
        return edges.size();
    }

    /**
     * Return the total number of edges, or only reciprocated.
     * @param reciprocated_only if true, determine number of reciprocated edges, else, all edges
     * @return number of edges, reciprocated if qualified by parameter
     */
    public int size(boolean reciprocated_only) {
        if (reciprocated_only)
            return recip.size();
        else
            return edges.size();
    }

    /**
     * Return number of edges that link with specified index
     * @param idx index
     * @return number of edges (regardless of order)
     */
    public int size(int idx) {
        Set<Integer> mylinks = links.get(idx);
        if (mylinks == null)
            return 0;
        else
            return mylinks.size();
    }

    /**
     * Retrieve all indices that are linked to specified index
     * @param idx
     * @return set with all indices (forward or backward)
     */
    public Set<Integer> getLinkedIndices(int idx) {
        Set<Integer> mylinks = links.get(idx);
        if (mylinks == null)
            return Collections.EMPTY_SET;
        return mylinks;
    }

    /**
     * Check if the edge between two indices idx1 --> idx2 exists
     * @param idx1 first index
     * @param idx2 second index (order is relevant)
     * @param reciprocated_only only check reciprocated edges
     * @return true if the edge exists, false otherwise
     */
    public boolean isEdge(int idx1, int idx2, boolean reciprocated_only) {
        if (reciprocated_only)
            return recip.contains(new POGEdge(idx1, idx2));
        else
            return edges.contains(new POGEdge(idx1, idx2));
    }

    /**
     * Check if the edge between two indices exists, regardless of whether the edge is reciprocated
     * @param idx1 first index
     * @param idx2 second index (order is relevant)
     * @return true if the edge exists, false otherwise
     */
    public boolean isEdge(int idx1, int idx2) {
        return isEdge(idx1, idx2, false);
    }

    public static class Directed extends EdgeMap {
        Set<POGEdge> forward = new HashSet<>();
        Set<POGEdge> backward = new HashSet<>();
        /**
         * Add an edge to the set of edges (idx1 --> idx2) associated with the POG.
         * When the edge is added a second time, it is automatically labeled reciprocated.
         * @param idx1 first index
         * @param idx2 second index
         * @param directionForward set to true if the edge is to be flagged Forward-looking, set to false for Backward-looking
         */
        public POGEdge add(int idx1, int idx2, boolean directionForward) {
            POGEdge e = add(idx1, idx2);
            if (directionForward)
                forward.add(e);
            else
                backward.add(e);
            return e;
        }

        public boolean isForward(POGEdge e) {
            return forward.contains(e);
        }

        public boolean isBackward(POGEdge e) {
            return backward.contains(e);
        }

    }
}

