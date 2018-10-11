/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vis;

import api.PartialOrderGraph;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.HashMap;
import java.util.Map;
import json.JSONArray;
import json.JSONObject;

/**
 *
 * @author ariane
 */
public class POAGJson {

    PartialOrderGraph poag;
    JSONObject jsonMap;
    Integer[] nodeIds;

    public POAGJson(PartialOrderGraph poag) {
        this.poag = poag;
        poag.getConsensusSequence(); // populate consensus flags
        nodeIds = poag.getNodeIDs();
        jsonMap = new JSONObject();
    }

    /**
     * Helper function to round to 2 decimal places for the vis.
     * @param value
     * @param places
     * @return
     */
    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    /**
     * Returns a JSON representation of the histogram,
     *
     * @param map
     * @return
     */
    private JSONArray map2JSON(Map<Character, Double> map, int multiplier) {
        JSONArray arr = new JSONArray();
        for (Map.Entry<Character, Double> list : map.entrySet()) {
            JSONArray arrEntry = new JSONArray();
            if (round(list.getValue(), 4) != 0) {
                arrEntry.put(Defines.G_LABEL, list.getKey().charValue());
                arrEntry.put(Defines.G_VALUE, round(list.getValue() * multiplier, 2));
                arr.put(arrEntry);
            }
        }
        return arr;
    }

    /**
     * Returns a JSON representation of the histogram,
     *
     * @param map
     * @return
     */
    private JSONArray seq2JSON(Map<Character, Double> map) {
        JSONArray arr = new JSONArray();
        return arr.put(Defines.G_CHAR, map2JSON(map, 1));
    }

    /**
     * Gets the outgoing edges of a node.
     * @param nodeId
     * @return
     */
    private Map<Integer, Double> getOutEdges (int nodeId) {
        return poag.getOutEdgeWeights(nodeId);
    }

    /**
     * Gets the graph for a specific Node Id.
     * @param nodeId
     * @return
     */
    private Map<Character, Double> getGraph (int nodeId) {
        return poag.getCharacterDistribution(nodeId);
    }

    /**
     * Makes a distribution of the sequence characters so that we can easily
     * display it in a pie chart when the amino acid hasn't been chosen yet.
     * @param nodeId
     * @return
     */
    private Map<Character, Double> getSeq (int nodeId, Map<Character, Double> graph) {
        Map<Character, Double> seqChars = new HashMap<>();
        if (poag.getSeqChars(nodeId) == null) {
            return seqChars;
        }
        int numSeqs = poag.getSeqChars(nodeId).size();
        for (Character base : graph.keySet()) {
            if (graph.get(base) > 0.01) {
                seqChars.put(base, graph.get(base) * numSeqs);
            }
        }
        return seqChars;
        }


    /**
     * Converts a node object to a JSON array representation.
     * @param nodeId
     * @return
     */
    private JSONArray nodeToJsonArray (Integer nodeId) {
        JSONArray jsonNode = new JSONArray();
        Map<Character, Double> graph = getGraph(nodeId);
        JSONArray seq = seq2JSON(getSeq(nodeId, graph));
        jsonNode.put(Defines.G_ID, nodeId);
        jsonNode.put(Defines.G_LABEL, nodeId.equals(poag.getFinalNodeID()) ? 'f' : nodeId.equals(poag.getInitialNodeID()) ? 'o' : poag.getLabel(nodeId).toCharArray()[0]);
        jsonNode.put(Defines.G_X, nodeId);
        jsonNode.put(Defines.G_GRAPH, graph == null ? null : map2JSON(graph, 100));
        jsonNode.put(Defines.G_SEQ, seq);
        jsonNode.put(Defines.G_MUTANTS, seq);
        jsonNode.put(Defines.G_CONSENSUS, poag.getConsensusMembership(nodeId));
        return jsonNode;
    }


    /**
     * Converts an Edge to a JSON array
     * @param nodeFromId
     * @param nodeToId
     * @param weight
     * @return
     */
    private JSONArray edgeToJsonArray (int nodeFromId, int nodeToId, double weight) {
        JSONArray jsonEdge = new JSONArray();
        try {
            jsonEdge.put(Defines.E_FROM, nodeFromId);
            jsonEdge.put(Defines.E_TO, nodeToId);
            jsonEdge.put(Defines.E_WEIGHT, (int) (100 * weight));
            jsonEdge.put(Defines.E_CONSENSUS, (poag.getConsensusMembership(nodeFromId)
                    && nodeToId == poag.getNextConsensusID(nodeFromId)) ? Defines.TRUE : Defines.FALSE);
            Integer[] reciprocatedIds = poag.getReciprocatedNextIDs(nodeFromId);
            int reciprocated = Defines.FALSE;
            for (int r = 0; r < reciprocatedIds.length; r++) {
                if (reciprocatedIds[r] == nodeToId) {
                    reciprocated = Defines.TRUE;
                    break;
                }
            }
            jsonEdge.put(Defines.E_RECIPROCATED, reciprocated);
            jsonEdge.put(Defines.E_SINGLE, (weight * poag.getNumberSequences() <= 1.5) ? Defines.TRUE : Defines.FALSE);
        } catch (Exception e) {
            System.err.println("Error with edge: " + nodeFromId + e.getMessage());
        }
        return jsonEdge;
    }


    /**
     * Converts an array of Nodes to a JSON Object.
     * @return
     */
    public JSONObject toJSON() {
        JSONArray nodesJSON = new JSONArray();
        JSONArray edgesJSON = new JSONArray();
        for (Integer nodeId : nodeIds) {
            nodesJSON.put(nodeToJsonArray(nodeId));
            Map<Integer, Double> outEdges = getOutEdges(nodeId);
            for (Map.Entry<Integer, Double> outEdge : outEdges.entrySet()) {
                edgesJSON.put(edgeToJsonArray(nodeId, outEdge.getKey(), outEdge.getValue()));
            }
        }
        jsonMap.put("nodes", nodesJSON);
        jsonMap.put("edges", edgesJSON);
        return jsonMap;
    }


}
