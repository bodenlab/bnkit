package reconstruction;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Stack;
import json.JSONArray;
import json.JSONObject;
import vis.Defines;
import vis.POAGJson;

/**
 * A class that performs two functions:
 *      1. Parsing the JSON object
 *      2. Developing the consensus from this
 */
public class Consensus {

    HashMap<Integer, Node> nodeMap;
    ArrayList<Edge> edges;
    Node initialNode;
    Node finalNode;

    public Consensus(POAGJson poagJson) {

    }

    /**
     * This class initialiser takes the unformatted JSON object from the controller in the GRASP
     * application.
     * Need to match up the nodes and the edges. We do
     * @param unformattedJson
     */
    public Consensus(JSONObject unformattedJson) {
        JSONArray jsonNodes = unformattedJson.getJSONArray("nodes");
        JSONArray jsonEdges = unformattedJson.getJSONArray("edges");

        edges = new ArrayList<>();
        nodeMap = new HashMap<>();

        int initialId = 10000;
        int endId = 0;

        // Iterate through each node and add it to a map on ID
        for (int n = 0; n < jsonNodes.length(); n ++) {
            JSONArray node = (JSONArray) jsonNodes.get(n);
            int nodeId = node.getInt(Defines.G_ID);
            Character base = (char) (int) node.get(Defines.G_CHAR);
            nodeMap.put(nodeId, new Node(base, nodeId));
            if (nodeId < initialId) {
                initialId = nodeId;
            }
            if (nodeId > endId) {
                endId = nodeId;
            }
        }

        // Iterate through each edge and add it to the edge array
        for (int e = 0; e < jsonEdges.length(); e ++) {
            JSONArray edgeJson = (JSONArray) jsonEdges.get(e);
            int fromId = edgeJson.getInt(Defines.E_FROM);
            int toId = edgeJson.getInt(Defines.E_TO);
            boolean reciprocated = edgeJson.getInt(Defines.E_RECIPROCATED) == Defines.TRUE ? true : false;
            double weight = edgeJson.getDouble(Defines.E_WEIGHT);
            Edge edge = new Edge(fromId, toId, weight, reciprocated);
            edges.add(edge);
            // Add the edge to the node map
            nodeMap.get(fromId).addOutEdge(edge);
        }

        // Run consensus gen
        initialNode = nodeMap.get(initialId);
        finalNode = nodeMap.get(endId);
    }


    /**
     * Method of adding in heuristics to the AStar search algorithm.
     * The aim is to get the path that has the most sequences in agreement.
     * However we also want the longest sequence (as this will contain the
     * fewest gaps).
     *
     * To do this we use the edge.getSequences().size() and multiply it by a
     * factor that penalises long gappy regions i.e. (1 + (1/gap_size))
     *
     * ToDo: optimise the heuristic function.
     *
     * @return
     */
    private double heuristicCostEstimate(Edge edge, Node from, Node to, boolean isBidirectional) {
        int multiplier = 1;
        if (!isBidirectional) {
            multiplier = 1000;
        }
        int positionDiff = java.lang.Math.abs(to.getId() - from.getId());
        positionDiff = (positionDiff > 0) ? positionDiff : 1;

        // Edge weight is out of 100
        double val =  multiplier * ((100 - edge.getWeight() + 1) * positionDiff);
        if (val < 0) {
            val = Double.MAX_VALUE;
        }
        return Math.abs(val);

    }


    /**
     * Reconstructs the consensus sequence based on the A* search. We need to
     * reverse the path and add in the gaps.
     *
     * @param cameFrom
     * @param current
     * @param gappy
     * @return
     */
    private String reconstructPath(HashMap<Node, Path> cameFrom, Node current, boolean gappy) {
        Stack<Character> sequence = new Stack<>();
        String sequenceString = "";
        while (cameFrom.keySet().contains(current)) {
            Path prevPath = cameFrom.get(current);
            Node prevNode = prevPath.getNode();
            // Set the edge to have a true consensus flag
            prevPath.getEdge().setConsensus(true);
            // If we have a character we want to add it
            if (current.getBase() != null && current != initialNode && current != finalNode) {
                sequence.push(current.getBase());
            }
            // Set to be the consensus path
            current.setConsensus(true);
            // If we have a gappy sequence we need to add in the gaps
            if (gappy) {
                int cameFromPosition = current.getId();
                int nextPosition = prevNode.getId();
                int numGaps = -1 * ((nextPosition - cameFromPosition) + 1);
                if (numGaps > 0) {
                    for (int g = 0; g < numGaps; g++) {
                        sequence.push('-');
                    }
                }
            }
            cameFrom.remove(current);
            current = prevNode;
        }
        // Reverse and create a string
        while (!sequence.empty()) {
            sequenceString += sequence.pop();
        }
        return sequenceString;
    }


    /**
     * Gets the consensus sequences using an A star search algorithm.
     *
     * @param gappy
     * @return
     */
    public String getSupportedSequence(boolean gappy) {
        // Intanciate the comparator class
        Comparator<Node> comparator = new NodeComparator();
        // Already visited nodes
        ArrayList<Node> closedSet = new ArrayList<>();
        // Unvisted nodes keep track of the best options
        PriorityQueue<Node> openSet = new PriorityQueue<>(10, comparator);
        // Add the initial node to the open set
        openSet.add(initialNode);
        // Storing the previous node
        HashMap<Node, Path> cameFrom = new HashMap<>();
        // Map with heuristics
        HashMap<Node, Double> cost = new HashMap<>();
        // Add the initial node cost
        cost.put(initialNode, new Double(0));
        while (!openSet.isEmpty()) {
            Node current = openSet.poll();
            if (current.equals(finalNode)) {
                // Reconstruct the path
                return reconstructPath(cameFrom, current, gappy);
            }
            // Otherwise add this to the closedSet
            closedSet.add(current);
            for (int n = 0; n < current.getOutEdges().size(); n++) {
                Edge next = current.getOutEdges().get(n);
                Node neighbor = nodeMap.get(next.getToId());
                double thisCost = heuristicCostEstimate(next, current, neighbor, current.getOutEdges().get(n).reciprocated);
                if (closedSet.contains(neighbor)) {
                    continue; // ignore as it has already been visited
                }
                // Otherwise we set the cost to this node
                double tentativeCost = cost.get(current) + thisCost;

                // Check if we have discovered a new node
                if (!openSet.contains(neighbor)) {
                    // Assign the cost to the node
                    neighbor.setCost(tentativeCost);
                    openSet.add(neighbor);
                } else if (tentativeCost >= cost.get(neighbor)) {
                    continue; // This isn't a better path
                }
                // If we have made it here this is the best path so let's
                cameFrom.put(neighbor, new Path(current, next));
                cost.put(neighbor, tentativeCost);
            }
        }
        return null;
    }


    /**
     * Comparator class for ensuring the nodes with the largest cost are
     * kept at the front of the queue.
     *
     * Here, we want to maximise the cost as we are trying to get either:
     * 		a. the longest sequence, or
     * 		b. the path which the greatest number of sequences agree on.
     */
    public class NodeComparator implements Comparator<Node>
    {
        @Override
        public int compare(Node x, Node y)
        {
            if (x.getCost() < y.getCost())
            {
                return -1;
            }
            if (x.getCost() > y.getCost())
            {
                return 1;
            }
            return 0;
        }
    }

    /**
     * Helper class to store the path. Keeps track of the edge and the node
     * that lead to a particular path. We need to keep track of the edge to
     * be able to set it as 'consensus' and the node to be able to get the
     * character.
     */
    public class Path {
        private Node node;
        private Edge edge;

        public Path(Node node, Edge edge) {
            this.edge = edge;
            this.node = node;
        }

        public Node getNode() { return this.node; }

        public Edge getEdge() { return this.edge; }
    }

    public class Edge {
        private int fromId;
        private int toId;
        private double weight;
        private boolean reciprocated;
        private boolean consensus;

        public Edge(int fromId, int toId, double weight, boolean reciprocated) {
            this.fromId = fromId;
            this.toId = toId;
            this.weight = weight;
            this.reciprocated = reciprocated;
        }

        public void setConsensus(boolean consensus) { this.consensus = consensus; }

        public boolean getReciprocated () { return this.reciprocated; }

        public int getFromId() { return this.fromId; }

        public int getToId() { return this.toId; }

        public double getWeight() {return this.weight; }

    }

    public class Node {
        private char base;
        private int id;
        private ArrayList<Edge> outEdges;
        private double cost = Double.MAX_VALUE;
        private boolean consensus = false;

        public Node(char base, int id) {
            this.base = base;
            this.id = id;
            this.outEdges = new ArrayList<>();
        }

        public ArrayList<Edge> getOutEdges() {
            return this.outEdges;
        }

        public double getCost() { return this.cost; }

        public void setCost(double cost) {
           this.cost = cost;
        }

        public void addOutEdge(Edge edge) {
            this.outEdges.add(edge);
        }

        public void setConsensus(boolean consensus) { this.consensus = consensus;}

        public boolean getConsensus() {
            return consensus;
        }

        public Character getBase() { return this.base; }


        public int getId()  { return this.id; }
    }
}
