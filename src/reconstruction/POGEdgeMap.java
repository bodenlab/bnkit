package reconstruction;

import java.util.*;

/**
 * Edge data structure for a POGraph.
 * Each member is identified by its source and target indices (reference to the original alignment).
 * An edge is either reciprocated (bi-directionally supported) or not.
 */
public class POGEdgeMap {

    public class POGEdge {
        final int idx1, idx2; // direction idx1 --> idx2

        /**
         * Create an edge from idx1 to idx2
         * @param idx1 first index
         * @param idx2 second index
         */
        public POGEdge(int idx1, int idx2) {
            this.idx1 = idx1;
            this.idx2 = idx2;
        }
        @Override
        public String toString() {
            return idx1 + "-->" + idx2;
        }
        @Override
        public boolean equals(Object o) {
            try {
                POGEdge e = (POGEdge)o;
                return (e.idx1 == this.idx1 && e.idx2 == this.idx2);
            } catch (ClassCastException e) {
                e.printStackTrace();
                return false;
            }
        }
        @Override
        public int hashCode() {
            return Objects.hash(idx1, idx2);
        }
        public int[] getIndices() {
            return new int[] {idx1, idx2};
        }
    }

    private final Map<Integer, Set<Integer>> links = new HashMap<>();
    private final Set<POGEdge> edges = new HashSet<>();
    private final Set<POGEdge> recip = new HashSet<>(); // reciprocated edges only
    private boolean STRICT_DIRECTIONALITY = true; // only allow indices to be linked ONE way

    public POGEdgeMap() {
    }

    public Set<POGEdge> getEdges() {
        return edges;
    }

    public boolean isReciprocated(POGEdge edge) {
        return recip.contains(edge);
    }

    /**
     * Add an edge to the set of edges (idx1 --> idx2) associated with the POG.
     * When the edge is added a second time, it is automatically labeled reciprocated.
     * @param idx1 first index
     * @param idx2 second index
     */
    public void add(int idx1, int idx2) {
        POGEdge e = new POGEdge(idx1, idx2);
//        if (STRICT_DIRECTIONALITY) {
//            POGEdge r = new POGEdge(idx2, idx1);
//            if (edges.contains(r)) {
//                System.err.println("Attempting to add invalid edge");
//                return;
//            }
//        }
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
}

