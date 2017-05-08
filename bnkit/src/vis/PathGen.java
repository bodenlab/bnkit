/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vis;

import api.PartialOrderGraph;

import java.util.*;

/**
 *
 * @author ariane
 */
public class PathGen {

    /**
     * A queue holding the states to be searched.
     */
    PriorityQueue<QueueTiny<Object>> closedSet;
    PriorityQueue<QueueTiny<Object>> openSet;
    int depth;
    int costPenalty;

    PartialOrderGraph poag;

    HashMap<Integer, Integer> cameFrom;
    HashMap<Integer, Double> stepScoreMap;
    HashMap<Integer, Double> totalScoreMap;

    Integer startNode;
    Integer goalNode;

    public Queue<Integer[]> gaps;
    private HashMap<Integer, Integer> searchedList = null;

    private final int size;

    Integer[] sorted;
    public int goalID;
    public int startID;
    HashMap<Integer, Node> nodes;

    ArrayList<Integer> totalNodeList;
    /**
     * The queue entry currently being processed.
     */
    protected QueueTiny<Object> currentEntry;
    QueueComparator queueComparator;

    public PathGen(PartialOrderGraph poag) {
        depth = 0;
        this.poag = poag;
        size = poag.getNodeIDs().length;
        nodes = new HashMap<>();
        searchedList = initSearchList();
        gaps = new LinkedList<>();
        // Get path depths fills everything out
        HashMap<Integer, List<Integer>> pathDepths = getPathDepths();
        // Here we want to add in any nodes which we may have missed at the 
        // start of the POAG. There may be multiple starting nodes.
    }

    public Integer choose_starting_node() {
        // Performing a topological sort sometimes causes the ordering between
        // the MSA and the inferred graph to differ. Instead assume that it is
        // already sorted.
        sorted = poag.getNodeIDs();
        ArrayList<Integer> startingNodeIds = new ArrayList<>();
        Double maxOutWeight = 0.0;
        Integer maxOutWeightID = sorted[0];
        Integer currNodeId = sorted[0];
        Map<Integer, Double> nodeOutWeights;
        for (int i = 0; i < size; i++) {
            Integer[] prevs = poag.getPreviousNodeIDs(sorted[i]);
            if (prevs[0] == null) {
                nodeOutWeights = poag.getOutEdgeWeights(sorted[i]);
                // Add it to our starting node ids
                startingNodeIds.add(sorted[i]);
                for (Map.Entry<Integer, Double> out : nodeOutWeights.entrySet()) {
                    Double outWeight = out.getValue();
                    if (outWeight > maxOutWeight) {
                        currNodeId = sorted[i];
                        maxOutWeightID = currNodeId; // set the id to be that 
                        maxOutWeight = outWeight; // with the largest out weight

                    }
                }
            }
        }
        // Add the nodes which won't be part of the main path with a depth of 1
        int depth = 1;
        for (Integer id : startingNodeIds) {
            if (id != maxOutWeightID) {
                Node n = new Node(id, id, depth, poag.getCharacterDistribution(id), poag.getOutEdgeWeights(id), poag.getSeqChars(id));
                nodes.put(id, n);
            }
        }
        // Return the node with the highest out weight
        return currNodeId;
    }

    public void initAStarSearch(Integer startNode, Integer goalNode) {
        int queueSize = 10000; //Arbitary should change.
        this.goalNode = goalNode;
        this.startNode = startNode;

        costPenalty = 101; // The max cost is 100 (if the chance is 100% we get 
        // a cost of a move = 1, as every move should have some penalty 
        // associated (NEED TO ASK)

        // The comparator currently just compares on the cost + heuristic
        queueComparator = new QueueComparator();
        closedSet = new PriorityQueue<>(queueSize, queueComparator);

        // Open nodes
        openSet = new PriorityQueue<>(queueSize, queueComparator);

        // Keeping track of how we got somewhere
        cameFrom = new HashMap<>();

        // Need to set the initial heuristics
        initScoreMaps();

        // Add the start node
        totalScoreMap.put(startNode, getHeuristic(startNode, goalNode));
    }

    /**
     * Initialises the score maps setting a max heuristic.
     */
    private void initScoreMaps() {
        // At the start the score is heuristic
        totalScoreMap = new HashMap<>();
        Integer[] poagIDs = poag.getNodeIDs();
        for (Integer poagID : poagIDs) {
            totalScoreMap.put(poagID, getHeuristic(0, goalNode));
        }
    }

    /**
     * Gets the path with the least cost associated from an initial to goal
     * node.
     *
     * @return
     */
    public List<Integer> getMainPath() {
        // Add the first position to the open set
        openSet.add(new QueueTiny<>(startNode, 0, getHeuristic(0, 0), startNode, null));
        totalScoreMap.put(startNode, getHeuristic(0, 0));
        cameFrom.put(0, startNode);
        while (!openSet.isEmpty()) {
            currentEntry = openSet.remove();
            closedSet.add(currentEntry);
            if (Objects.equals(currentEntry.getCurrentPosition(), goalNode)) {
                // Check to make sure we didn't get straight there from the 
                // initial node, if we did return null and don't re add the gap
                if (Objects.equals(cameFrom.get(goalNode), startNode)) {
                    return null;
                }
                return (getPath());
            }
            searchNeighbours();
        }
        return null;
    }

    /**
     * Searches the neighboring nodes and updates the cost to these nodes if it
     * is less from this current node.
     */
    private void searchNeighbours() {
        int currentID = currentEntry.getCurrentPosition();
        Map<Integer, Double> successors = poag.getOutEdgeWeights(currentID);
        Integer[] successorIds = poag.getNextNodeIDs(currentID);
        double heuristic = 100000; // An arbitarilty large initial heuristic.
        double currentCost = totalScoreMap.get(currentID);

        for (Integer neighborID : successorIds) {

            // This is if we are searching between gaps - we already know the
            // direct path is best so we want to skip the optimal path and find
            // the secondary one.
            if (neighborID == goalNode && currentID == startNode) {
                continue;
            }
            int pos = isNodeSearched(neighborID);
            // If there has been no weight associated it will error out so we 
            // need to give it the maxinum cost of 0 (which is what the label says)
            double posCost = costPenalty;
            try {
                posCost = costPenalty - (100 * successors.get(neighborID)); //Because we want it to have some cost                                
            } catch (Exception e) {
                posCost = 200; // This has no weight associated need to ask!
                System.err.println("The ids didnt match the length of out "
                        + " weights: " + Arrays.toString(successorIds) + "\n" + successors.toString());
            }
            heuristic = getHeuristic(0, 0);
            // Means it hasn't been searched yet otherwise we ignore it
            if (pos > 0) {
                // Get the total score to the 
                double totalMoveCost = currentCost + posCost;
                openSet.add(new QueueTiny<>(neighborID, totalMoveCost, heuristic, neighborID, currentID));
                // If its less than this is a better move than the current score
                // to get to that position - best path so far
                try {
                    if (totalMoveCost < totalScoreMap.get(neighborID)) {
                        cameFrom.put(neighborID, currentID);
                        totalScoreMap.put(neighborID, totalMoveCost);
                    }
                } catch (Exception e) {
                    //System.err.println("This node wasn't searchable? + " + neighborID);
                }
            }
        }
        // Mark this node as searched
        setNodeSearched(currentID);
    }

    /**
     * Gets the heuristic between two nodes.
     */
    private Double getHeuristic(Integer startNode, Integer goalNode) {
        return 100.0 * (goalNode - startNode);
    }

    /**
     * Gets the path - steps back through the path with the least cost.
     *
     * @return most likely path.
     */
    private List<Integer> getPath() {
        int currId = currentEntry.getCurrentPosition();
        String totalPath = "" + currId;
        List<Integer> pathReverse = new ArrayList<>();
        pathReverse.add(currId);
        int prevId = 0;
        while (currId != startNode) {
            prevId = currId;
            try {
                currId = cameFrom.get(prevId);
            } catch (Exception e) {
                System.err.println(e);
                return null;
            }
            totalPath += "-" + currId;
            pathReverse.add(currId);
            // Means there is a gap so add it to the gaps
            if (prevId - currId != 1) {
                // If the gap size is two we know we're just going to want to 
                // add the node one out from the two (add one to the depth)
                if (prevId - currId == 2) {
                    int gapNodeId = prevId + 1;
                    int prevY = depth + 1;
                    Node n = new Node(gapNodeId, gapNodeId, prevY, poag.getCharacterDistribution(gapNodeId), poag.getOutEdgeWeights(gapNodeId), poag.getSeqChars(gapNodeId));
                    nodes.put(gapNodeId, n);
                } else {
                    addGaps(currId, prevId);
                }
            }
        }
        if (prevId - startNode != 1) {
            if (prevId - currId == 2) {
                int gapNodeId = prevId + 1;
                int prevY = depth + 1;
                Node n = new Node(gapNodeId, gapNodeId, prevY, poag.getCharacterDistribution(gapNodeId), poag.getOutEdgeWeights(gapNodeId), poag.getSeqChars(gapNodeId));
                nodes.put(gapNodeId, n);
            } else {
                addGaps(currId, prevId);
            }
        }
        return pathReverse;
    }

    /**
     * Gets the search indication for a node
     *
     * @param id
     * @return 0 for not searched and -1 if it has been searched
     */
    public int isNodeSearched(int id) {
        return searchedList.get(id);
    }

    /**
     * Adds gaps to the stored gap list associated with the poag.
     *
     * @param start
     * @param end
     */
    public void addGaps(int start, int end) {
        Integer[] gap = {start, end};
        gaps.add(gap);
    }

    /**
     * Gets the search indication for a node
     *
     * @param id
     * @return 0 for not searched and -1 if it has been searched
     */
    public int setNodeSearched(int id) {
        return searchedList.put(id, -1);
    }

    /**
     * The nodes need to be reset to be not searched once we have the main path
     * and are looking for subsequent paths.
     */
    public void resetSearchList() {
        Set<Integer> keys = searchedList.keySet();
        keys.stream().forEach((Integer k) -> {
            searchedList.put(k, 1);
        });
    }

    /**
     * Initialise the search list to indicate none of the nodes have been
     * searched yet (set all to 0)
     */
    private HashMap<Integer, Integer> initSearchList() {
        HashMap<Integer, Integer> search = new HashMap<>();
        // Performing a topological sort sometimes causes the ordering between
        // the MSA and the inferred graph to differ. Instead assume that it is
        // already sorted.
        // sorted = poag.sort();
        sorted = poag.getNodeIDs();
        goalID = sorted[sorted.length - 1];
        // Choose the starting node and add those which won't be covered by the
        // search path to the list
        startID = choose_starting_node();
        for (int i = 0; i < size; i++) {
            search.put(sorted[i], 1);
        }
        return search;
    }

    /**
     * The comparator for the queue. This stores then with the lowest cost as
     * the highest priority.
     */
    public class QueueComparator implements Comparator<QueueTiny<Object>> {

        @Override
        public int compare(QueueTiny<Object> o1, QueueTiny<Object> o2) {
            return Double.compare(o1.totalCost + o1.heuristicEstimate,
                    o2.totalCost + o2.heuristicEstimate);
        }
    }

    /**
     * Returns a map of paths and the depth at which these occur
     *
     * @return paths
     */
    public HashMap<Integer, List<Integer>> getPathDepths() {
        Integer[] nodeys = poag.getNodeIDs();
        int numNodes = nodeys.length;

        HashMap<Integer, List<Integer>> paths = new HashMap<>();
        // want there to be a distinct x position for each node
        // This gets the main path which will be centered
        resetSearchList();
        initAStarSearch(startID, goalID);
        List<Integer> path = getMainPath();
        paths.put(depth, path);
        addNodes(path, depth);
        // Now we need to get each of the subsequent paths
        // A depth will be associated for each itteration - this will be
        // visualised as further away from the central line.
        depth++;

        int prevGapStart = 10000000; // Something large (note this needs to be done better)
        int gapEnd = 0;
        int gapStart = 0;
        while (!gaps.isEmpty()) {
            // Get the next gap from the FIFO gap queue set up in the POAG
            Integer[] gap = gaps.remove();
            gapStart = gap[0];
            gapEnd = gap[1];

            // Set up the AStar search environment with the start and end 
            // This also clears the searched nodes in the search map
            resetSearchList();
            initAStarSearch(gapStart, gapEnd);
            path = getMainPath();

            // If the path was null we don't want to add it - see code below
            // for alternative method (ugly with nested while loops to step
            // through the gap.
            if (path == null) {
                //System.err.println("depth: " + depth + " Path: was none for gap: " + Arrays.toString(gap));
                continue;
            }
            //System.err.println("Depth: " + depth + ", Path: " + path);  
            paths.put(depth, path);
            addNodes(path, depth);
            // Check if we have reached the end of the first itteration of gaps
            // We can tell this because the gaps will start again from the 
            // End of the path.
            if (gapStart > prevGapStart) {
                depth++;
            }

            prevGapStart = gapStart;
        }

        // Reset the nodes to have correct x coords
        return paths;
    }

    private void addNodes(List<Integer> path, int depth) {
        int id;

        for (int i : path) {
            id = i;
            Node n = nodes.get(id);
            if (n == null) {
                // if we don't have the node already then we want to add it
                n = new Node(id, id, depth, poag.getCharacterDistribution(id), poag.getOutEdgeWeights(id), poag.getSeqChars(id));
            }
            // Otherwise, check if the depth is smaller for this one as we always
            // want the smallest depth, keep the origional x
            if (n.getY() > depth) {
                n.setY(depth);
            }
            nodes.put(id, n);
        }
    }
}