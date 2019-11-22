package api;

import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.POGraph;
import reconstruction.ConsensusObject;
import vis.POAGJson;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Partial order graph data structure.
 *
 * Created by marnie on 28/3/17.
 * Patched by mikael 10/10/2019
 * FIXME: clean out a lot of deprecated functions relating to the old definition of consensus/most supported sequence.
 */
public class PartialOrderGraph {

    private POGraph graph = null;
    private ConsensusObject consensus = null;

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

    public EnumSeq getMostSupported(boolean GAPPY) {
        if (consensus == null) {
            consensus = new ConsensusObject(graph.getEdgeCountsNode(), graph.getNumSeqsUnderNode());
            // Here we could use the POAGJson object itself rather than the actual JSON object
            consensus.setJsonObject(new POAGJson(this, GAPPY).toJSON());
            consensus.getSupportedIndices();
        }
        char[] supportedSeq = consensus.getSupportedSequence(GAPPY).toCharArray();
        Character[] arr = new Character[supportedSeq.length];
        for (int j = 0; j < arr.length; j ++)
            arr[j] = supportedSeq[j];
        EnumSeq ancseq = GAPPY ? new EnumSeq.Gappy(Enumerable.aacid_ext) : new EnumSeq(Enumerable.aacid_ext);
        ancseq.set(arr);
        return ancseq;
    }

    public EnumDistrib[] getDistribMostSupported(boolean GAPPY) {
        if (consensus == null) {
            consensus = new ConsensusObject(graph.getEdgeCountsNode(), graph.getNumSeqsUnderNode());
            // Here we could use the POAGJson object itself rather than the actual JSON object
            consensus.setJsonObject(new POAGJson(this, GAPPY).toJSON());
            consensus.getSupportedIndices();
        }
        int[] idxs = consensus.getSupportedIndices();
        EnumDistrib[] ret = new EnumDistrib[GAPPY ? getFinalNodeID() : idxs.length];
        int j = 0; // index for walking through the nodes in the POG
        for (int i = 0; i < ret.length; i ++) {
            int pogidx = idxs[j];
            if (GAPPY) {
                if (i == pogidx) {
                    ret[i] = graph.getNode(pogidx).getDistrib();
                    j++;
                } else
                    ret[i] = null;
            } else {
                ret[i] = graph.getNode(pogidx).getDistrib();
                j++;
            }
        }
        return ret;
    }

    public int[] getIndicesMostSupported(boolean GAPPY) {
        if (consensus == null) {
            consensus = new ConsensusObject(graph.getEdgeCountsNode(), graph.getNumSeqsUnderNode());
            // Here we could use the POAGJson object itself rather than the actual JSON object
            consensus.setJsonObject(new POAGJson(this, GAPPY).toJSON());
            consensus.getSupportedIndices();
        }
        int[] idxs = consensus.getSupportedIndices();
        if (!GAPPY)
            return idxs;
        int[] ret = new int[getFinalNodeID()];
        for (int i = 0; i < ret.length; i ++)
            ret[i] = i;
        return ret;
    }

    /**
     * Get the sequence with the most edge support through the graph. Support is based on maximizing immediate edge
     * weights.
     *
     * @return      most supported sequence based on edge weights
     * @deprecated
     */
    public String getConsensusSequence() {
        return graph.getSupportedSequence(false);
    }

    /**
     * Get the sequence with the most edge support through the graph. Support is based on maximizing immediate edge
     * weights.
     *
     * @return      most supported sequence based on edge weights
     * @deprecated
     */
    public String getConsensusSequence(boolean gappy) {
        return graph.getSupportedSequence(gappy);
    }

    /**
     * Get the sequence with the most edge support through the graph. Support is based on maximizing immediate edge
     * weights. Represents 'gaps' (i.e. jumps in the partial order graph)
     *
     * @return      most supported gappy sequence based on edge weights
     * @deprecated
     */
    public String getConsensusGappySequence() {
        return graph.getSupportedSequence(true);
    }


    /**
     * Get indication of consensus membership of the node with the provided ID.
     *
     * @param id    Node ID
     * @return      Flag indicating consensus membership of the node with the provided ID
     * @deprecated
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
     * @deprecated
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
