package dat.pog;

import asr.ASRRuntimeException;

import java.util.*;

/**
 * Container class to perform a search on a graph.
 * The idea is to allow sub-classes to be specific to some class of graphs (sub-classing IdxEdgeGraph)
 * and operate on some edge type and node constraints.
 */
public class GraphSearch<E extends Edge & WeightedEdge> {

    final IdxEdgeGraph<E> g;    // the graph that is being searched
    final int N;                // total number of nodes (excluding terminators)

    /**
     * Create a graph search instance and perform the search
     */
    public GraphSearch(IdxEdgeGraph<E> g) {
        this.g = g;
        N = g.maxsize();
    }

    /**
     * Implements the value with which edges are (typically) ranked.
     * Sub-classes should override.
     * @param edge
     * @return
     */
    public double getPriority(E edge) {
        return edge.getWeight();
    };

    /**
     * Class that defines a search node, which is used in a PriorityQueue by A*
     */
    public class SearchNode implements Comparable<SearchNode> {

        final int from, to;
        final E edge;
        final double f;

        /**
         * Create a search node from the two indices that are joined by the edge, which needs to be explored.
         * A cost is assigned based on the actual cost to reach the source node, and the weight on the edge, factoring-in
         * whether the edge is reciprocated or not
         * @param from source node index
         * @param to target node index
         */
        public SearchNode(int from, int to) {
            this.from = from;
            this.to = to;
            this.edge = g.getEdge(from, to);
            this.f = getPriority(this.edge);
        }

        /**
         * The search node can be compared to any other node, ensuring that the priority queue polls
         * the node with the lowest estimated distance
         * @param o other node
         * @return -1 if this node is lower, +1 if this node is greater, or 0 if they are equal
         */
        @Override
        public int compareTo(SearchNode o) {
            double f1 = this.f;
            double f2 = o.f;
            return f1 < f2 ? -1 : f1 > f2 ? +1 : 0;
        }
    }

}

/**
 * Find the paths from start to end that has the fewest number of non-reciprocated edges
 */
class DijkstraSearch<E extends POGraph.StatusEdge> extends GraphSearch<E> {

    final Set[] closed;     // indices linking back, ultimately to -1 (start node)
    final double[] actual;  // actual costs/distance when determined and applicable; each element is -1 before determined
    final int start, goal;  // start and goal nodes for search

    /**
     * Create a search instance and perform the search
     */
    public DijkstraSearch(POGraph g) {
        this((IdxEdgeGraph<E>) g, -1, g.maxsize(), PRIORITY_RECIP_WEIGHT);
    }

    public DijkstraSearch(POGraph g, int mode) {
        this((IdxEdgeGraph<E>) g, -1, g.maxsize(), mode);
    }

    /**
     * Create a search instance and perform the search
     */

    public DijkstraSearch(IdxEdgeGraph<E> g, int start, int goal) {
        this(g, start, goal, PRIORITY_WEIGHT);
    }

    /**
     * Create a search instance and perform the search
     */
    public DijkstraSearch(IdxEdgeGraph<E> g, int start, int goal, int mode) {
        super(g);
        this.setPriorityMode(mode);
        actual = new double[N + 1]; // space to store actual costs
        Arrays.fill(actual, Double.POSITIVE_INFINITY);
        closed = new Set[N + 1]; // space to store linkages leading-to a node index
        this.start = start;
        this.goal = goal;
        if (start >= 0) { // terminated graphs can't refer to index for initial terminal "-1"
            actual[start] = 0; // the cost is by definition zero for the start node
            closed[start] = new HashSet();
        }
        int current = start;
        Set<Integer> visited = new HashSet<>();
        while (current != goal) {
            // track visited nodes (by index) and cheapest path to each node; this is in "closed"
            // find the lowest-cost path from "start" (-1) to "anywhere"; this is in "actual"
            int[] nexts;
            nexts = g.getNodeIndices(current); // all out-going edges; forward-looking if directed
            if (goal == N && g.isEndNode(current)) { // if goal is the terminus, need also check if the current node is an end node
                // add N to next
                int[] xnexts = new int[nexts.length + 1];
                for (int i = 0; i < nexts.length; i++)
                    xnexts[i] = nexts[i];
                xnexts[nexts.length] = N;
                nexts = xnexts;
            }
            visited.add(current);
            for (int next : nexts) { // add all nodes coming-up
                double cost_sofar = (current == start ? 0 : actual[current]);
                E edge = g.getEdge(current, next);
                double weight = getPriority(edge);
                if (actual[next] > cost_sofar + weight) { // if better...
                    actual[next] = cost_sofar + weight;
                    closed[next] = new HashSet();
                    closed[next].add(current);
                } else if (actual[next] == cost_sofar + weight && closed[next] != null) { // else if just equal...
                    closed[next].add(current);
                }
            }
            double cheapest_cost = Double.POSITIVE_INFINITY;
            // set "current" to lowest-cost node that can be reached from "start", perhaps via any of the other "visited"
            for (int i = 0; i < actual.length; i ++) {
                if (!visited.contains(i)) {
                    if (cheapest_cost > actual[i]) {
                        current = i;
                        cheapest_cost = actual[i];
                    }
                }
            }
            if (cheapest_cost == Double.POSITIVE_INFINITY)
                break;
        }
    }

    public static final int PRIORITY_RECIPROCATED = 0x01;
    public static final int PRIORITY_WEIGHT = 0x02;
    public static final int PRIORITY_RECIP_WEIGHT = 0x04;

    private int PRIORITY_MODE = 0x02;
    public void setPriorityMode(int mode) {
        PRIORITY_MODE = mode;
    }

    @Override
    public double getPriority(E edge) {
        switch (PRIORITY_MODE) {
            case PRIORITY_RECIPROCATED:
                return edge.getReciprocated() ? 0 : 1;
            case PRIORITY_WEIGHT:
                return edge.getWeight();
            case PRIORITY_RECIP_WEIGHT:
                return edge.getReciprocated() ? edge.getWeight() : 1000 * edge.getWeight();
            default:
                return 0;
        }
    }

    /**
     * Retrieve one optimal path
     * @return an array with an optimal path as specified by the node indices, in order; null if no path was found.
     */
    public int[] getOnePath() {
        List<Integer> path = new ArrayList<>();
        int current = goal;
        if (closed[goal] == null)
            return null;
        //path.add(current);
        while (true) {
            current = (Integer)closed[current].iterator().next();
            if (current == start)
                break;
            path.add(current);
        }
        int[] arr = new int[path.size()];
        for (int i = 0; i < arr.length; i ++)
            arr[i] = path.get(arr.length - i - 1);
        return arr;
    }

    /**
     * Retrieve the edges that make up all optimal paths (can be one, multiple or none).
     * @return the set of edges that make up all optimal paths; the set is empty if no optimal paths exist
     */
    public Set<POGEdge> getOptimal() {
        // Find all edges that can form part of an optimal path (can be many paths)
        if (closed[goal] == null)
            return Collections.EMPTY_SET;
        Stack<Integer> stack = new Stack<>();
        stack.push(goal);
        Set<POGEdge> optimal = new HashSet<>();
        Set<Integer> visited = new HashSet<>();
        while (!stack.empty()) {
            int idx2 = stack.pop();
            if (!visited.contains(idx2)) {
                for (int idx1 : (Set<Integer>)closed[idx2]) {
                    optimal.add(new POGEdge(idx1, idx2));
                    if (idx1 != -1)
                        stack.push(idx1);
                }
            }
            visited.add(idx2);
        }
        return optimal;
    }

    /**
     * Retrieve the actual, optimal cost of visiting the node with the specified index
     * @param idx the node index
     * @return the cost
     */
    public double getCost(int idx) {
        if (idx == -1)
            return 0;
        if (actual[idx] >= 0)
            return actual[idx];
        else
            throw new ASRRuntimeException("Invalid A* search: cannot determine cost onwards index " + idx);
    }

    /**
     * Retrieve the cost of reaching the end of the search
     * @return the cost, lower is better
     */
    public double getCost() {
        return actual[this.goal];
    }

}

class AStarSearch extends GraphSearch {

    final Set[] closed;     // indices linking back, ultimately to -1 (start node)
    final double[] actual;  // actual costs when determined and applicable; each element is -1 before determined
    final double minW;      // smallest weight, if non-negative, else set to 1
    final int[] toposort;   // nodes in topological sort
    final int[] mapidx;     // offset by +1, so that index -1 works (at 0), and also N (at N + 1)

    /**
     * Create an A* search instance and perform the search
     */
    public AStarSearch(POGraph pog) {
        super(pog);
        mapidx = new int[N + 2];
        toposort = pog.getTopologicalOrder();
        for (int i = 0; i < toposort.length; i++)
            mapidx[toposort[i] + 1] = i;
        mapidx[N + 1] = N;
        actual = new double[N + 1]; // space to store actual costs
        Arrays.fill(actual, -1);
        closed = new Set[N + 1]; // space to store linkages leading-to a node index
        double saveW = Double.MAX_VALUE;
        for (POGraph.StatusEdge e : pog.getEdges().values()) {
            if (e.getWeight() < saveW)
                saveW = e.getWeight();
        }
        minW = saveW <= 0 ? 1 : saveW;
        PriorityQueue<SearchNode> pq = new PriorityQueue<>();
        int current = -1;
        double cost = 0; // the cost to reach the start node is 0
        while (true) {
            // first, check if this is an end node, and if so add path to terminate
            if (pog.isEndNode(current)) {
                // create a new path to terminate graph, in turn calculate costs
                SearchNode node = new SearchNode(current, N);
                pq.add(node); // this will need eventual exploration
            }
            // second, check what the next frontier of nodes are, and add their paths
            int[] nexts = pog.getForward(current);
            for (int next : nexts) { // add all nodes coming-up
                // create a new path to a search node, in turn calculate costs
                SearchNode node = new SearchNode(current, next);
                pq.add(node); // this will need eventual exploration
            }
            // third, poll the most promising node...
            SearchNode promise = pq.poll();
            if (promise == null) // no more nodes to explore, so terminate
                break;
            current = promise.to;
            //cost = promise.g;
            if (cost < 0)
                System.out.println("\t" + cost);
            if (actual[current] < 0 || actual[current] >= cost) { // not yet set, or just better
                if (actual[current] > cost || closed[current] == null)
                    closed[current] = new HashSet();
                closed[current].add(promise.from);
                actual[current] = cost;
            }
            if (current == N)
                break;
        }

    }

    @Override
    public double getPriority(Edge edge) {
        double f = 0; // FIXME getCost(from) + ((mapidx[to + 1] - mapidx[from + 1]) * edge.getWeight()) * (edge.getReciprocated() ? 1 : N) + getH(this.to);
        return f;
    }

    /**
     * Retrieve one optimal path
     * @return an array with an optimal path as specified by the node indices, in order; null if no path was found.
     */
    public int[] getOnePath() {
        List<Integer> path = new ArrayList<>();
        int current = N;
        if (closed[N] == null)
            return null;
        while (true) {
            current = (Integer)closed[current].iterator().next();
            if (current < 0)
                break;
            path.add(current);
        }
        int[] arr = new int[path.size()];
        for (int i = 0; i < arr.length; i ++)
            arr[i] = path.get(arr.length - i - 1);
        return arr;
    }

    /**
     * Retrieve the edges that make up all optimal paths (can be one, multiple or none).
     * @return the set of edges that make up all optimal paths; the set is empty if no optimal paths exist
     */
    public Set<POGEdge> getOptimal() {
        // Find all edges that can form part of an optimal path (can be many paths)
        if (closed[N] == null)
            return Collections.EMPTY_SET;
        Stack<Integer> stack = new Stack<>();
        Set<POGEdge> optimal = new HashSet<>();
        Set<Integer> visited = new HashSet<>();
        while (!stack.empty()) {
            int idx2 = stack.pop();
            if (!visited.contains(idx2)) {
                for (int idx1 : (Set<Integer>)closed[idx2]) {
                    stack.push(idx1);
                    optimal.add(new POGEdge(idx1, idx2));
                }
            }
            visited.add(idx2);
        }
        return optimal;
    }

    /**
     * Retrieve the actual, optimal cost of visiting the node with the specified index
     * @param idx the node index
     * @return the cost
     */
    public double getCost(int idx) {
        if (idx == -1)
            return 0;
        if (actual[idx] >= 0)
            return actual[idx];
        else
            throw new ASRRuntimeException("Invalid A* search: cannot determine cost onwards index " + idx);
    }

    /**
     * Retrieve the cost of reaching the end of the POGraph
     * @return the cost, lower is better
     */
    public double getCost() {
        return actual[N];
    }

    /**
     * The admissible heuristic function used for finding the optimal path.
     * The default function uses the topological order of the node relative to the end, multiplied by the minimum weight of the POG
     * @param idx the node index
     * @return the (optimistic) estimate of the remaining distance to reach the end
     */
    public double getH(int idx) {
        return (N - mapidx[idx + 1]) * minW; // weight is minimally 1
    }


}