package api;

import dat.POGraph;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
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

    public PartialOrderGraph(POGraph graph) { this.graph = graph; }

    public PartialOrderGraph(String dotStructure) throws IOException {
        graph = new POGraph(dotStructure);
    }

    /**
     * Get the sequence with the most edge support through the graph. Support is based on maximizing immediate edge
     * weights.
     *
     * @return      most supported sequence based on edge weights
     */
    public String getConsensusSequence() {
        return graph.getSupportedSequence(false);
    }

    /**
     * Get the sequence with the most edge support through the graph. Support is based on maximizing immediate edge
     * weights. Represents 'gaps' (i.e. jumps in the partial order graph)
     *
     * @return      most supported gappy sequence based on edge weights
     */
    public String getConsensusGappySequence() {
        return graph.getSupportedSequence(true);
    }


    /**
     * Get indication of consensus membership of the node with the provided ID.
     *
     * @param id    Node ID
     * @return      Flag indicating consensus membership of the node with the provided ID
     */
    public Boolean getConsensusMembership(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        return graph.getCurrentConsensusFlag();
    }

    /**
     * Get the ID of the next node in the consensus path.
     *
     * @param id    ID of the current node
     * @return      ID of the next node in the consensus path.
     */
    public Integer getNextConsensusID(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        return graph.getNextConsensusID();
    }

    /**
     * Get the set of next node IDs that have reciprocated edges.
     *
     * @param id    ID of the current node
     * @return      Set of next node IDs that have reciprocated edges
     */
    public Integer[] getReciprocatedNextIDs(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        Integer[] ids = new Integer[graph.getReciprocatedNextIDs().size()];
        graph.getReciprocatedNextIDs().toArray(ids);
        return ids;
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
    public Integer[] getNextNodeIDs(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        Integer[] ids = new Integer[graph.getNextIDs().size()];
        graph.getNextIDs().toArray(ids);
        return ids;
    }

    public Integer getFinalNodeID(){
        return graph.getFinalNodeID();
    }

    public Integer getInitialNodeID() { return graph.getInitialNodeID(); }

    /**
     * Get the IDs of the previous nodes of the node with the given ID.
     *
     * @param id    ID of current node
     * @return      IDs of 'previous' nodes
     */
    public Integer[] getPreviousNodeIDs(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        List<Integer> prevIds = graph.getPreviousIDs();
        if (prevIds == null)
            return null;
        Integer[] ids = new Integer[graph.getPreviousIDs().size()];
        graph.getPreviousIDs().toArray(ids);
        return ids;
    }

    /**
     * Get the characters in the node with the given ID.
     *
     * @param id    ID of current node
     * @return      Map of sequence characters
     */
    public Map<Integer, Character> getSeqChars(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        return new HashMap<>(graph.getSequenceCharacterMapping());
    }

    /**
     * Return the number of sequences in the graph.
     *
     * @return  number of sequences
     */
    public int getNumberSequences() {
        return graph.getSequences().size();
    }


    /**
     * Get the weights of the out edges of the node with the given ID.
     *
     * @param id    ID of current node
     * @return      Map between next node ID and edge weight
     */
    public Map<Integer, Double> getOutEdgeWeights(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        return new HashMap<>(graph.getNextEdgeWeights());
    }

    /**
     * Get the character state of the node with the given ID. Returns null if the state is unknown.
     *
     * @param id    ID of current node
     * @return      Character state of node, or null
     */
    public String getLabel(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        return graph.getCurrentLabel();
    }

    /**
     * Get the probability distribution of characters in the node with the given ID.
     *
     * @param id    ID of current node
     * @return      Map between character and probability
     */
    public Map<Character, Double> getCharacterDistribution(Integer id) {
        if (!graph.setCurrent(id))
            return null;
        Map<Character, Double> dist = graph.getCharacterDistribution();
        return (dist == null) ? null :  new HashMap<>(graph.getCharacterDistribution());
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

    /**
     * String representation of partial order graph.
     *
     * @return      String representation
     */
    public String toString() {
        return graph.toString();
    }
}
