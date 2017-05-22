package api;

import dat.POGraph;

import java.util.HashMap;
import java.util.Map;

/**
 * Partial order graph data structure.
 *
 * Created by marnie on 28/3/17.
 */
public class PartialOrderGraph {

    private POGraph graph = null;

    /**
     * Constructors.
     */
    public PartialOrderGraph() {
        graph = new POGraph();
    }

    public PartialOrderGraph(String dotStructure) {
        graph = new POGraph(dotStructure);
    }

    /**
     * Get the sequence with the most edge support through the graph. Support is based on maximizing immediate edge
     * weights.
     *
     * @return      most supported sequence based on edge weights
     */
    public String getConsensusSequence() {
        return graph.getSupportedSequence();
    }

    /**
     * Get all node IDs in the graph.
     *
     * @return      IDs of nodes in the graph
     */
    public Integer[] getNodeIDs() {
        Integer[] ids = new Integer[graph.getNumNodes()];
        graph.getNodeIDs().toArray(ids);
        return ids;
    }

    /**
     * Get the IDs of the next nodes of the node with the given ID.
     *
     * @param id    ID of current node
     * @return      IDs of 'next' nodes
     */
    public Integer[] getNextNodeIDs(int id) {
        graph.setCurrent(id);
        Integer[] ids = new Integer[graph.getNextIDs().size()];
        graph.getNextIDs().toArray(ids);
        return ids;
    }

    /**
     * Get the weights of the out edges of the node with the given ID.
     *
     * @param id    ID of current node
     * @return      Map between next node ID and edge weight
     */
    public Map<Integer, Double> getOutEdgeWeights(int id) {
        graph.setCurrent(id);
        return new HashMap<>(graph.getEdgeWeights());
    }

    /**
     * Get the character state of the node with the given ID. Returns null if the state is unknown.
     *
     * @param id    ID of current node
     * @return      Character state of node, or null
     */
    public String getLabel(int id) {
        graph.setCurrent(id);
        return graph.getCurrentLabel();
    }

    /**
     * Get the probability distribution of characters in the node with the given ID.
     *
     * @param id    ID of current node
     * @return      Map between character and probability
     */
    public Map<Character, Double> getCharacterDistribution(int id) {
        graph.setCurrent(id);
        return new HashMap<>(graph.getCharacterDistribution());
    }

    /**
     * Perform a topological sort and get an array of sorted node IDs.
     *
     * @return      Sorted node IDs
     */
    public Integer[] sort() {
        Integer[] ids = new Integer[graph.getNumNodes()];
        graph.topologicalSort().toArray(ids);
        return ids;
    }
}
