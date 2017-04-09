/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package poag;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import json.JSONArray;
import json.JSONObject;
import api.PartialOrderGraph;

/**
 *
 * @author ariane
 */
public class POAGJson {

    PartialOrderGraph poag;
    HashMap<Integer, Node> nodes;
    HashMap<Integer, Integer> xcoords;
    JSONObject jsonMap;

    public POAGJson(PartialOrderGraph poag, HashMap<Integer, Node> nodes) {
        this.poag = poag;
        this.nodes = nodes;
        jsonMap = new JSONObject();
        xcoords = new HashMap<>();

    }

    /**
     * For the POAG's where the state hasn't been decided we pass the label
     * in a format so that a pie chart can easily be rendered on the
     * JavaScript side
     * @param map
     * @return 
     */
    private JSONObject seq2JSON(Map<Character, Integer> map) {
        JSONObject chars = new JSONObject();
        JSONArray bars = new JSONArray();
        int numChars = map.size();
        for (Map.Entry<Character, Integer> list : map.entrySet()) {

                JSONObject bar = new JSONObject();
                bar.put("label", list.getKey().toString());
                bar.put("value", list.getValue());
                bars.put(bar);
        }
        chars.put("chars", bars);
        return chars;
    }

    /**
     * Returns a JSON representation of the histogram,
     * @param map
     * @return 
     */
    private JSONObject map2JSON(Map<Character, Double> map) {
        JSONObject graph = new JSONObject();
        JSONArray bars = new JSONArray();
        for (Map.Entry<Character, Double> list : map.entrySet()) {
            if (list.getValue() * 100 > 1) {
                JSONObject bar = new JSONObject();
                bar.put("x_label", list.getKey().toString());
                bar.put("value", (int) (list.getValue() * 100));
                bars.put(bar);
            }
        }
        graph.put("bars", bars);
        return graph;
    }

    public void toJSON() {
        JSONObject nodesJSON = new JSONObject();
        int x;
        int y;
        JSONObject reactions = new JSONObject();
        Node n;
        Set<Integer> ns = nodes.keySet();
        for (int i : ns) {
            JSONObject thisNode = new JSONObject();
            //           JSONObject nodesReactions = new JSONObject();
            n = nodes.get(i);
            x = n.getX();
            y = n.getY();
            int id = n.getID();
            String nid = "node-" + id;

            thisNode.put("label", poag.getLabel(id));
            thisNode.put("x", x);
            thisNode.put("y", y);
            thisNode.put("graph", map2JSON(n.getGraph()));
            thisNode.put("seq", seq2JSON(n.getSeq()));
            nodesJSON.put(nid, thisNode);
            // Need to get the out egdes into JSON reaction object
            Map<Integer, Double> outedges = n.getOutedges();
            for (Map.Entry<Integer, Double> outNodes : outedges.entrySet()) {
                try {
                    JSONObject thisReaction = new JSONObject();
                    int n2id = outNodes.getKey();
                    Node tempNode = nodes.get(n2id);
                    String rid = "edges_" + id + ":" + n2id;
                    int x2 = tempNode.getX();
                    int y2 = tempNode.getY();
                    int weight = (int) (100 * outNodes.getValue());
                    thisReaction.put("from", id);
                    thisReaction.put("to", n2id);
                    thisReaction.put("x1", x);
                    thisReaction.put("x2", x2);
                    thisReaction.put("y1", y);
                    thisReaction.put("y2", y2);
                    thisReaction.put("weight", weight);
                    reactions.put(rid, thisReaction);
                } catch (Exception e) {
                    System.err.println("Error with reaction: " + nid);
                }

            }
            // Want to add each of the edge weights to the reaction nodes
        }
        jsonMap.put("nodes", nodesJSON);
        jsonMap.put("edges", reactions);

        System.err.println(jsonMap);
    }

}
