package reconstruction;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;
import json.JSONArray;
import json.JSONObject;
import vis.Defines;

public class ConsensusObject {
    /**
     * A class that performs two functions:
     *      1. Parsing the JSON object
     *      2. Developing the consensus from this
     */
    Map<Integer,  Node> nodeMap;
    ArrayList< Edge> edges;
    int[] consensusByIndex = null; // the indices that make up the most supported/consensus sequence

    // These are the nodes that are dummy start and end nodes in the original version.
    // We keep them as they are used on the front end for signifying the start and end terminals.
    Node initialNode;
    Node finalNode;
    int numberPositionsInSequence;

    // Map with heuristics
    Map<Integer, Double> cost;
    Map<Integer, Map<Integer, Integer>> edgeCounts;
    int[] countArray;

    int numberSeqsUnderNode;
    ArrayList<Integer> possibleInitialIds;
    ArrayList<Integer> possibleFinalIds;

    Map<Integer, Edge> initialAndFinalEdges;
    Map<Integer, Node> initialAndFinalNodeMap;

    // These are our best (or optimal) initial nodes. We use these to determine the consensus seq.
    Node bestInitialNode;
    Node bestFinalNode;

    boolean updateEdgeCounts; // Used to toggel whether or not we want to update the cost of the edges based on the
    // edge count map that is passed in.

    public ConsensusObject(Map<Integer, Map<Integer, Integer>> edgeCounts, int numberSeqsUnderNode) {
        this.numberSeqsUnderNode = numberSeqsUnderNode;
        this.edgeCounts = edgeCounts;
        // For now let's just automatically set to update the edge counts cost
        this.updateEdgeCounts = true;
        cost = new HashMap<>();
    }

    private int getBestCount(ArrayList<Integer> possibleIds) {
        // Here we want to get the best count in the node dictionary.
        int bestCount = 0;
        for (Integer initialId: possibleIds) {
            int thisCount = countArray[initialId];
            if (thisCount > bestCount) {
                bestCount = thisCount;
            }
        }
        return bestCount;
    }


    /**
     *
     * @param updateEdgeCost
     */
    public void setUpdateEdgeCounts(boolean updateEdgeCost) {
        this.updateEdgeCounts = updateEdgeCost;
    }

    /**
     * Build the mapping of the edges.
     *
     * @return
     */
    public TreeObject buildTreeFromNewick(String treeAsNewick) {
        // Get the reconstructed tree from the database
        // Parse the tree
        if (treeAsNewick != null) {
            return new TreeObject(treeAsNewick);
        }
        System.out.println("No tree string passed.");
        return null;
    }


    /**
     *
     * @param nwkfilepathOrString
     * @param load
     * @param edgeCounts
     * @return
     */
    public TreeObject getEdgeMapping(String nwkfilepathOrString, boolean load, Map<Integer, Map<Integer, Map<Integer, Integer>>> edgeCounts) {
        // Get the reconstructed tree from the database
        TreeObject tree;
        if (load) {
            tree = new TreeObject(nwkfilepathOrString, true);
        } else {
            tree = new TreeObject(nwkfilepathOrString);
        }
        TreeNodeObject root = tree.getRoot();
        // Now we want to build the edge count map up recursively.
        root.buildEdgeCountMap(edgeCounts, tree.getSeqIdMap());
        // Return the tree containing the mapping (i.e. the counts).
        return tree;
    }


    /**
     *
     * @param nwkfilepathOrString
     * @param load
     * @param nodeLabel
     * @param edgeCounts
     * @return
     */
    public TreeNodeObject getEdgeMappingForNode(String nwkfilepathOrString, boolean load, String nodeLabel, Map<Integer, Map<Integer, Map<Integer, Integer>>> edgeCounts) {
        TreeObject tree;
        if (load) {
            tree = new TreeObject(nwkfilepathOrString, true);
        } else {
            tree = new TreeObject(nwkfilepathOrString);
        }
        TreeNodeObject root = tree.getRoot();
        // Now we want to build the edge count map up recursively.
        root.buildEdgeCountMap(edgeCounts, tree.getSeqIdMap());
        TreeNodeObject node = null;
        try {
            node = tree.getNodeByLabel(nodeLabel.split("_")[0]);
        } catch (Exception e) {
            try {
                node = tree.getNodeByLabel(nodeLabel);
            } catch (Exception e1) {
                System.out.print(e.getStackTrace() + "" + e1.getStackTrace());
                return null;
            }
        }

        // Return the tree containing the mapping (i.e. the counts).
        return node;
    }



    /**
     To ensure that if we have a tie we actually choose the best one (lets for now say the
     best one is the one which has a higher sequence content in that column (very unlikely
     that this would be the same)
     Set the best first and last node ID's
     * @param ids
     * @param tagToLookFor
     * @return
     */
    public int getBestId(ArrayList<Integer> ids, Character tagToLookFor) {
        ArrayList<Integer> bestIds = new ArrayList<>();
        double errorMargin = 0.05 * numberSeqsUnderNode;
        int bestCount = getBestCount(ids);
        for (Integer i: ids) {
            // Check that the best IDS are with the error margin deemed to be 5% of the sequence count
            if (countArray[i] > bestCount - errorMargin) {
                bestIds.add(i);
            }
        }

        if (bestIds.size() == 0) {
            System.out.println("ERROR!");
            return ids.get(0);
        }
        if (bestIds.size() == 1) {
            return bestIds.get(0);
        }

        ArrayList<Integer> bestIdsWithMet = new ArrayList<>();
        if (tagToLookFor != null) {
            for (Integer i : bestIds) {
                if (nodeMap.get(i).getBase() == tagToLookFor) {
                    bestIdsWithMet.add(i);
                }
            }

            if (bestIdsWithMet.size() == 1) {
                return bestIdsWithMet.get(0);
            }
        }
        // If none had a met start tag we want to change it back to the original bestIds
        if (bestIdsWithMet.size() == 0) {
            bestIdsWithMet = bestIds;
        }

        double bestWeightVal = 0.0;
        int bestFinalId = 0;
        ArrayList<Integer> bestBiDir = new ArrayList<>();
        for (Integer i : bestIdsWithMet) {
            if (nodeMap.get(i).getBiDir()) {
                bestBiDir.add(i);
            }
        }
        if (bestBiDir.size() == 1) {
            return bestBiDir.get(0);
        }

        if (bestBiDir.size() == 0) {
            bestBiDir = bestIdsWithMet;
        }
        for (Integer i: bestBiDir) {
            double weightVal = countArray[i];
            if (weightVal > bestWeightVal) {
                bestWeightVal = weightVal;
                bestFinalId = i;
            }
        }
        return bestFinalId;
    }

    /**
     * This is used so we can create a consensus object and add the unformatted json after,
     * used to make it more efficient when we build the map at the start.
     * @param unformattedJson
     */
    public void setJsonObject(JSONObject unformattedJson) {
        formatJSON(unformattedJson);
    }

    /**
     * Format the JSON object into Java objects for processing.
     * @param unformattedJson
     */
    private void formatJSON(JSONObject unformattedJson) {
        JSONArray jsonNodes = unformattedJson.getJSONArray("nodes");
        JSONArray jsonEdges = unformattedJson.getJSONArray("edges");
        initialAndFinalEdges = new HashMap<>();
        edges = new ArrayList<>();
        nodeMap = new HashMap<>();
        initialAndFinalNodeMap = new HashMap<>();

        int initialId = 10000;
        int endId = 0;

        // Iterate through each node and add it to a map on ID
        for (int n = 0; n < jsonNodes.length(); n++) {
            JSONArray node = (JSONArray) jsonNodes.get(n);
            int nodeId = node.getInt(Defines.G_ID);
            Character base = (char) (int) node.get(Defines.G_CHAR);
            nodeMap.put(nodeId, new Node(base, nodeId, node, n));
            if (nodeId < initialId) {
                initialId = nodeId;
            }
            if (nodeId > endId) {
                endId = nodeId;
            }
        }

        possibleInitialIds = new ArrayList<>();
        possibleFinalIds = new ArrayList<>();

        // Iterate through each edge and add it to the edge array
        for (int e = 0; e < jsonEdges.length(); e++) {
            JSONArray edgeJson = (JSONArray) jsonEdges.get(e);
            int fromId = edgeJson.getInt(Defines.E_FROM);
            int toId = edgeJson.getInt(Defines.E_TO);
            double weight = edgeJson.getDouble(Defines.E_WEIGHT);
            if (updateEdgeCounts) {
                try {
                    weight = (edgeCounts.get(fromId).get(toId)/(double)numberSeqsUnderNode) * 100;
                    edgeJson.put(Defines.E_WEIGHT, weight);
                } catch (Exception excep) {
                    // Skip this edge as this is not possible, we want to remove it.
                    continue;
                }
            }

            boolean reciprocated = edgeJson.getInt(Defines.E_RECIPROCATED) == Defines.TRUE;
            Edge edge = new Edge(fromId, toId, weight, reciprocated, edgeJson, e);
            edges.add(edge);
            // Add the edge to the node map
            nodeMap.get(fromId).addOutEdge(edge);
            // If the from ID is the initial node id then we want to keep track of the possible
            // node IDs that were part of this.
            if (fromId == initialId) {
                possibleInitialIds.add(toId);
                initialAndFinalEdges.put(toId, edge);
                // Want to set that this node is bi-dir
                Node node = nodeMap.get(toId);
                node.setBiDir(edge.getReciprocated());
                initialAndFinalNodeMap.put(toId, node);

            }
            // Similarly, if the toId is the final node ID we want to keep track of this
            if (toId == endId) {
                possibleFinalIds.add(fromId);
                initialAndFinalEdges.put(fromId, edge);
                // set whether it is bidirectional or not
                Node node = nodeMap.get(fromId);
                node.setBiDir(edge.getReciprocated());
                initialAndFinalNodeMap.put(fromId, nodeMap.get(fromId));

            }
        }
        this.initialNode = nodeMap.get(initialId);
        this.finalNode = nodeMap.get(endId);
        // Calculate the best initial and final IDs
        // Set the number of nodes to be the end ID
        numberPositionsInSequence = nodeMap.size();
    }



    /**
     * Converts the consensus object into the JSON representation required by the grasp app.
     * @return
     */
    public JSONObject getAsJson() {
        JSONArray nodesJSON = new JSONArray();
        JSONArray edgesJSON = new JSONArray();
        JSONObject jsonMap = new JSONObject();
        for (Edge e: edges) {
            edgesJSON.put(e.arrayPos, e.edgeAsJson);
        }
        for (Integer n: nodeMap.keySet()) {
            nodesJSON.put( nodeMap.get(n).arrayPos, nodeMap.get(n).nodeAsJson);
        }
        jsonMap.put("nodes", nodesJSON);
        jsonMap.put("edges", edgesJSON);
        return jsonMap;
    }


    /**
     *
     * We want the smallest cost estimate to get from the start to the end.
     *
     * We have penalties for things such as a non-bidirectional edge.
     *      1. weight penalty = numberPositionsInSequence * 2
     * We also have the number of sequences that have taken an edge, we convert this into a % by
     * dividing by the number of sequences under this node
     *      2. 1 - edge weight (we do 1 - since we want the largest possible amount)
     * Method of adding in heuristics to the AStar search algorithm.
     * The aim is to get the path that has the most sequences in agreement.
     * However we also want the longest sequence (as this will contain the
     * fewest gaps). Hence we multiply by the position change.
     *
     *
     * numberPositionsInSequence * 2 is the maximum value we could get for a bidirectional
     * path - hence we use this as the cost.
     * ToDo: optimise the heuristic function.
     *
     * @return
     */
    private double heuristicCostEstimate(Edge edge, Node from, Node to, boolean isBidirectional) {

        int multiplier = 1;
        int positionDiff = java.lang.Math.abs(to.getId() - from.getId());
        positionDiff = (positionDiff > 0) ? positionDiff : 1;

        // This has a 0 cost
        Double weight = 0.0;
        try {
            // Note we don't use edge.getWeight here, as sometimes the person won't have
            // wanted to override this.
            weight = edgeCounts.get(from.getId()).get(to.getId())/(double)(this.numberSeqsUnderNode);

        } catch (Exception e) {
            // This mean this edge doesn't exist, so we can't go down here (for this - lets return the max value)
            // ToDo: Do we just want to add a penalty? Saying no seqs go down here... another option.
            System.out.println("NO WEIGHTS " + from.getId() + "->" + to.getId());
            multiplier = 2 * numberPositionsInSequence;
        }

        // Ensure we don't add our multiplier to start and end nodes
        if (!isBidirectional && to.getId() != finalNode.getId() && from.getId() != initialNode.getId()) {
            multiplier = 2 * numberPositionsInSequence;
        }
        if (!isBidirectional) // give bi-directional edges a tiny advantage, regardless
            multiplier *= 1.001;

        double val = 1 - weight;
        if (val < 0) {
            // That is very strange and we need to return an error
            System.out.println("ERROR: RUNNING: " + from.getId() + "->" + to.getId() + "VAL < 0: " + val);
            val = numberPositionsInSequence * 2;
        }

        // We add 1 on here since just incase we have a 0 cost (i.e. all the sequences travel here,
        // we still want to assign this a cost of "moving". and ensure this is multiplied by the
        // number of poistions we travel as this is additive.
        val =  ((multiplier * val) + 1) * positionDiff;

        return val;
    }


    /**
     * Reconstructs the consensus sequence based on the A* search. We need to
     * reverse the path and add in the gaps.
     *
     * @param cameFrom
     * @param current
     * @param gappy
     * @return
     * @deprecated just kept until function of replacement is verified
     */
    private String reconstructPathLegacy(Map< Node,  Path> cameFrom,  Node current, boolean gappy) {
        Stack<Character> sequence = new Stack<>();
        // Add the initial Base that we decided on during the pre-processing stage
        String sequenceString = "";
        // Set the initial and final edges

        while (cameFrom.keySet().contains(current)) {
            Path prevPath = cameFrom.get(current);
            prevPath.getEdge().setConsensus(true);
            Node prevNode = prevPath.getNode();
            // Set the edge to have a true consensus flag
            prevPath.getEdge().setConsensus(true);
            // If we have a character we want to add it
            if (current.getBase() != null && current != finalNode) {
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
     * Reconstructs the consensus sequence based on the A* search. We need to
     * reverse the path and add in the gaps.
     *
     * @param cameFrom
     * @param current
     * @return
     */
    public int[] updateSupportedPath(Map< Node,  Path> cameFrom, Node current) {
        Stack<Integer> sequence = new Stack<>();
        while (cameFrom.keySet().contains(current)) {
            Path prevPath = cameFrom.get(current);
            prevPath.getEdge().setConsensus(true);
            Node prevNode = prevPath.getNode();
            // Set the edge to have a true consensus flag
            prevPath.getEdge().setConsensus(true);
            // If we have a character we want to add it
            if (current.getBase() != null && current != finalNode)
                sequence.push(current.getId());
            // Set to be the consensus path
            current.setConsensus(true);
            cameFrom.remove(current);
            current = prevNode;
        }

        // Reverse and create a string
        int[] ret = new int[sequence.size()];
        int i = 0;
        while (!sequence.empty()) {
            ret[i]= sequence.pop();
            i ++;
        }
        if (consensusByIndex == null)
            consensusByIndex = ret;
        return ret;
    }

    public Node getLowestCostNode(ArrayList<Node> openSet) {
        double minCost = Double.MAX_VALUE;
        ArrayList<Node> bests = new ArrayList<>();
        for (Node n: openSet) {
            if (cost.get(n.getId()) < minCost) {
                minCost = cost.get(n.getId());
            }
        }
        if (minCost == Double.MAX_VALUE) {
            int i = 0;
        }
        // If there are multiple bests we want to tie break on the one with the smaller nodeId
        for (Node n: openSet) {
            if (cost.get(n.getId()) == minCost) {
                bests.add(n);
            }
        }

        // If the length of the array is > 1 we want to return the node with the  lowest ID
        Node best = null;
        int lowestId = 1000000000;
        for (Node n: bests) {
            if (n.getId() < lowestId) {
                lowestId = n.getId();
                best = n;
            }
        }
        // Remove the best node from the openset
        openSet.remove(best);

        return best;
    }

    /**
     * Get the consensus sequence using an A star search algorithm.
     * Determines the indices that form part of the most supported sequence.
     * The logic here comes from (Ariane's old) getSupportedSequence, and an array
     * is stored internal to this class to save compute time...
     * @return
     */
    public int[] getSupportedIndices() {
        // Instantiate the comparator class
        Comparator< Node> comparator = new NodeComparator();
        // Already visited nodes
        ArrayList< Node> closedSet = new ArrayList<>();
        // Unvisited nodes keep track of the best options
        // ToDo Work out why priority queue isn't working
        //PriorityQueue< Node> openSet = new PriorityQueue<>(1000, comparator);
        ArrayList<Node> openSet = new ArrayList<>();
        // Add the initial node to the open set
        //for (Integer i: possibleInitialIds) {
        //   Node n = nodeMap.get(i);
        openSet.add(initialNode);
        //if (i == this.initialNode.getId()) {
        cost.put(initialNode.getId(), 0.0);

        // Storing the previous node
        Map< Node,  Path> cameFrom = new HashMap<>();

        boolean printout = false;
        while (!openSet.isEmpty()) {
            Node current = getLowestCostNode(openSet); //openSet.poll();
            if (current.id == 1446)
                System.out.println(current.id);
            if (current == null) {
                current = openSet.get(0);
                openSet.remove(current);
            }
            if (current.equals(finalNode)) {
                // Reconstruct the path
                int[] idxs = updateSupportedPath(cameFrom, current);
                return idxs;
            }
            // Otherwise add this to the closedSet
            closedSet.add(current);

            if (printout) {
                System.out.println("Looking at edges from: " + current.getId());
            }
            for (int n = 0; n < current.getOutEdges().size(); n++) {
                Edge next = current.getOutEdges().get(n);
                Node neighbor = nodeMap.get(next.getToId());
                double thisCost = heuristicCostEstimate(next, current, neighbor, current.getOutEdges().get(n).reciprocated);
                if (closedSet.contains(neighbor)) {
                    // Check if this path is better and update the path to get to the neighbour
                    if (cost.get(neighbor.getId()) > thisCost) {
                        //cameFrom.put(neighbor, new Path(current, next));
                        cost.put(neighbor.getId(), thisCost);
                    }
                    continue; // ignore as it has already been visited
                }
                // Otherwise we set the cost to this node
                double tentativeCost = cost.get(current.getId()) + thisCost;

                // Check if we have discovered a new node
                if (!openSet.contains(neighbor)) {
                    // Assign the cost to the node
                    neighbor.setCost(tentativeCost);
                    cost.put(neighbor.getId(), tentativeCost);
                    openSet.add(neighbor);
                } else if (tentativeCost > cost.get(neighbor.getId())) {
                    if (printout) {
                        System.out.println("WORSE : " + current.getBase() + "-" + neighbor.getBase() + ": " + neighbor.getId() + " , " + neighbor.getBase() + ":" + thisCost + ", " + tentativeCost + " vs." + neighbor.getCost());
                    }
                    continue; // This isn't a better path
                }
                cost.put(neighbor.getId(), tentativeCost);
                neighbor.setCost(tentativeCost);

                if (printout) {
                    System.out.println("BETTER : " + current.getBase() + "-" + neighbor.getBase() + ": " + neighbor.getId() + " , " + neighbor.getBase() + ":" + thisCost + ", " + tentativeCost + " vs." + cost.get(neighbor.getId()));
                }
                // Check if we already have this in the camefrom path, if so remove
                // ToDo: Check if overriding is causing issues?
                if (cameFrom.get(neighbor) != null) {
                    //System.out.println("ALREADY HAD PATH, BEING OVERRIDDEN, " + cameFrom.get(neighbor).edge.fromId + "->" + cameFrom.get(neighbor).edge.toId + ", " + cameFrom.get(neighbor).node.base + " to " + current.base + " path:" + next.fromId + " ->" + next.toId);
                }
                // If we have made it here this is the best path so let's
                cameFrom.put(neighbor, new Path(current, next));
            }
        }
        return null;
    }

    /**
     * Gets the consensus sequences using an A star search algorithm.
     *
     * @param gappy
     * @return
     */
    public String getSupportedSequence(boolean gappy) {
        if (consensusByIndex == null)
            consensusByIndex = getSupportedIndices();
        int[] idxs = consensusByIndex;
        StringBuilder sb = new StringBuilder();
        int sidx = 0; // index in sequence
        for (int i = 0; i < idxs.length; i ++) {
            int pidx = idxs[i]; // index in POG
            if (gappy) {
                for (; sidx < pidx; sidx ++)
                    sb.append('-');
            }
            sb.append(nodeMap.get(pidx).getBase());
            sidx ++;
        }
        if (gappy) {
            for (; sidx < this.finalNode.getId(); sidx ++)
                sb.append('-');
        }
        return sb.toString();
    }

//        // Instantiate the comparator class
//        Comparator< Node> comparator = new NodeComparator();
//        // Already visited nodes
//        ArrayList< Node> closedSet = new ArrayList<>();
//        // Unvisited nodes keep track of the best options
//        // ToDo Work out why priority queue isn't working
//        //PriorityQueue< Node> openSet = new PriorityQueue<>(1000, comparator);
//        ArrayList<Node> openSet = new ArrayList<>();
//        // Add the initial node to the open set
//        //for (Integer i: possibleInitialIds) {
//        //   Node n = nodeMap.get(i);
//        openSet.add(initialNode);
//        //if (i == this.initialNode.getId()) {
//        cost.put(initialNode.getId(), 0.0);
//
//        // Storing the previous node
//        Map< Node,  Path> cameFrom = new HashMap<>();
//
//        boolean printout = false;
//        while (!openSet.isEmpty()) {
//            Node current = getLowestCostNode(openSet); //openSet.poll();
//            if (current == null) {
//                current = openSet.get(0);
//                openSet.remove(current);
//            }
//            if (current.equals(finalNode)) {
//                // Reconstruct the path
//                return reconstructPath(cameFrom, current, gappy);
//            }
//            // Otherwise add this to the closedSet
//            closedSet.add(current);
//
//            if (printout) {
//                System.out.println("Looking at edges from: " + current.getId());
//            }
//            for (int n = 0; n < current.getOutEdges().size(); n++) {
//                Edge next = current.getOutEdges().get(n);
//                Node neighbor = nodeMap.get(next.getToId());
//                double thisCost = heuristicCostEstimate(next, current, neighbor, current.getOutEdges().get(n).reciprocated);
//                if (closedSet.contains(neighbor)) {
//                    // Check if this path is better and update the path to get to the neighbour
//                    if (cost.get(neighbor.getId()) > thisCost) {
//                        //cameFrom.put(neighbor, new Path(current, next));
//                        cost.put(neighbor.getId(), thisCost);
//                    }
//                    continue; // ignore as it has already been visited
//                }
//                // Otherwise we set the cost to this node
//                double tentativeCost = cost.get(current.getId()) + thisCost;
//
//                // Check if we have discovered a new node
//                if (!openSet.contains(neighbor)) {
//                    // Assign the cost to the node
//                    neighbor.setCost(tentativeCost);
//                    cost.put(neighbor.getId(), tentativeCost);
//                    openSet.add(neighbor);
//                } else if (tentativeCost > cost.get(neighbor.getId())) {
//                    if (printout) {
//                        System.out.println("WORSE : " + current.getBase() + "-" + neighbor.getBase() + ": " + neighbor.getId() + " , " + neighbor.getBase() + ":" + thisCost + ", " + tentativeCost + " vs." + neighbor.getCost());
//                    }
//                    continue; // This isn't a better path
//                }
//                cost.put(neighbor.getId(), tentativeCost);
//                neighbor.setCost(tentativeCost);
//
//                if (printout) {
//                    System.out.println("BETTER : " + current.getBase() + "-" + neighbor.getBase() + ": " + neighbor.getId() + " , " + neighbor.getBase() + ":" + thisCost + ", " + tentativeCost + " vs." + cost.get(neighbor.getId()));
//                }
//                // Check if we already have this in the camefrom path, if so remove
//                // ToDo: Check if overriding is causing issues?
//                if (cameFrom.get(neighbor) != null) {
//                    //System.out.println("ALREADY HAD PATH, BEING OVERRIDDEN, " + cameFrom.get(neighbor).edge.fromId + "->" + cameFrom.get(neighbor).edge.toId + ", " + cameFrom.get(neighbor).node.base + " to " + current.base + " path:" + next.fromId + " ->" + next.toId);
//                }
//                // If we have made it here this is the best path so let's
//                cameFrom.put(neighbor, new Path(current, next));
//            }
//        }
//        return null;
//    }


    /**
     * Comparator class for ensuring the nodes with the largest cost are
     * kept at the front of the queue.
     *
     * Here, we want to maximise the cost as we are trying to get either:
     * 		a. the longest sequence, or
     * 		b. the path which the greatest number of sequences agree on.
     */
    public class NodeComparator implements Comparator< Node>
    {
        @Override
        public int compare( Node x,  Node y)
        {
            if (cost.get(x.getId()) < cost.get(y.getId()))
            {
                return -1;
            }
            if (cost.get(x.getId()) > cost.get(y.getId()))
            {
                return 1;
            }
            if (x.getId() < y.getId()) {
                return -1;
            }
            if (x.getId() > y.getId()) {
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
        private  Node node;
        private  Edge edge;

        public Path( Node node,  Edge edge) {
            this.edge = edge;
            this.node = node;
        }

        public  Node getNode() { return this.node; }

        public  Edge getEdge() { return this.edge; }
    }

    public class Edge {
        private int fromId;
        private int toId;
        private double weight;
        private boolean reciprocated;
        private boolean consensus;
        private JSONArray edgeAsJson;
        private String id;
        private int arrayPos;

        public Edge(int fromId, int toId, double weight, boolean reciprocated, JSONArray edgeAsJson, int arrayPos) {
            this.arrayPos = arrayPos;
            this.fromId = fromId;
            this.toId = toId;
            this.weight = weight;
            this.reciprocated = reciprocated;
            this.edgeAsJson = edgeAsJson;
            this.edgeAsJson.put(Defines.E_CONSENSUS, Defines.FALSE);
            this.id = fromId + "-" + toId;
        }

        public void setConsensus(boolean consensus) {
            this.consensus = consensus;
            edgeAsJson.put(Defines.E_CONSENSUS, Defines.TRUE);
        }

        public boolean getReciprocated () { return this.reciprocated; }

        public int getFromId() { return this.fromId; }

        public int getToId() { return this.toId; }

        public double getWeight() { return this.weight; }

        public void setWeight(double weight) { this.weight = weight; }

    }

    public class Node {
        private char base;
        private int id;
        private ArrayList< Edge> outEdges;
        private double cost = Double.MAX_VALUE;
        private boolean consensus = false;
        private JSONArray nodeAsJson;
        int arrayPos = 0;
        private boolean isBiDir = false;

        public Node(char base, int id, JSONArray nodeAsJson, int arrayPos) {
            this.base = base;
            this.arrayPos = arrayPos;
            this.id = id;
            this.outEdges = new ArrayList<>();
            this.nodeAsJson = nodeAsJson;
            this.nodeAsJson.put(Defines.G_CONSENSUS, Defines.FALSE);
        }

        public void setBiDir(boolean isBiDir) {
            this.isBiDir = isBiDir;
        }

        public boolean getBiDir() {
            return this.isBiDir;
        }

        public ArrayList< Edge> getOutEdges() {
            return this.outEdges;
        }

        public double getCost() { return this.cost; }

        public void setCost(double cost) {
            this.cost = cost;
        }

        public void addOutEdge( Edge edge) {
            this.outEdges.add(edge);
        }

        public void setConsensus(boolean consensus) {
            this.nodeAsJson.put(Defines.G_CONSENSUS, Defines.TRUE);
            this.consensus = consensus;
        }

        public boolean getConsensus() {
            return consensus;
        }

        public Character getBase() { return this.base; }

        public int getId()  { return this.id; }
    }
}