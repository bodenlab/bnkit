/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vis;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import json.JSONArray;
import json.JSONObject;
import api.PartialOrderGraph;

/**
 * A class which assigns x and y coordinates to a POAG so that it can 
 * be visualised on the front end in a more readable manner.
 * Takes in a POAG and returns a JSON object representative of that POAG.
 * @author ariane
 */
public class POAGJson {

    PartialOrderGraph poag;
    HashMap<Integer, Node> nodes;
    HashMap<Integer, Integer> xcoords;
    JSONObject jsonMap;
    PathGen pathGen;

    public POAGJson(PartialOrderGraph poag) {
        this.poag = poag;
        // PathGen generates the path from the POAG
        pathGen = new PathGen(poag);
        this.nodes = pathGen.nodes;
        jsonMap = new JSONObject();
        xcoords = new HashMap<>();
    }

    /**
     * For the POAG's where the state hasn't been decided we pass the label in a
     * format so that a pie chart can easily be rendered on the JavaScript side
     *
     * @param map
     * @return
     */
    private JSONObject seq2JSON(Map<Character, Integer> map) {
        JSONObject chars = new JSONObject();
        JSONArray bars = new JSONArray();
        for (Map.Entry<Character, Integer> list : map.entrySet()) {
            JSONObject bar = new JSONObject();
            bar.put("label", list.getKey().toString());
            bar.put("value", list.getValue());
            bars.put(bar);
            break;
        }
        chars.put("chars", bars);
        return chars;
    }

    /**
     * Returns a JSON representation of the histogram.
     * @param sequences with the probabilities assigned
     * @return histogram as JSON
     */
    private JSONObject map2JSON(Map<Character, Double> map) {
        JSONObject graph = new JSONObject();
        JSONArray bars = new JSONArray();
        for (Map.Entry<Character, Double> list : map.entrySet()) {
            JSONObject bar = new JSONObject();
            bar.put("x_label", list.getKey().toString());
            bar.put("value", (list.getValue() * 100));
            bars.put(bar);
        }
        graph.put("bars", bars);
        return graph;
    }

    
    /**
     * Creates a JSON object with the nodes and edges of the POAG as two 
     * sub objects.
     * Nodes are a JSON object with each node stored with its x coordinate as 
     * the ID.
     * Edges are similarly stored, the edge ID is a concatenation of the node
     * to and node from ID's which this edge joins.
     * @return JSON representation of the nodes and edges for the POAG
     */
    public JSONObject toJSON() {
        JSONArray nodesJSON = new JSONArray();
        int x;
        int y;
        int max_depth = 0;
        JSONObject reactions = new JSONObject();
        Node n;
        Set<Integer> ns = nodes.keySet();
        for (int i : ns) {
            JSONObject thisNode = new JSONObject();
            n = nodes.get(i);
            x = n.getX();
            y = n.getY();
            int id = n.getID();
            String nid = "node-" + id;
            // Extra things for the multi view poag
            thisNode.put("id", id);
            // Ones needed for tthe actual poag
            thisNode.put("label", poag.getLabel(id));
            thisNode.put("x", x);
            thisNode.put("y", y);
            thisNode.put("graph", map2JSON(n.getGraph()));
            thisNode.put("seq", seq2JSON(n.getSeq()));
            nodesJSON.put(thisNode);
            //Check if it is the max depth
            if (y > max_depth) {
                max_depth = y;
            }
            // Need to get the out egdes into JSON reaction object
            Map<Integer, Double> outedges = n.getOutedges();
            // Add the number of out edges so we can know whether this node
            // is of interest.
            thisNode.put("num_out_edges", outedges.size());
            System.err.println(outedges.size());
            for (Map.Entry<Integer, Double> outNodes : outedges.entrySet()) {
                try {
                    JSONObject thisEdge = new JSONObject();
                    int n2id = outNodes.getKey();
                    Node tempNode = nodes.get(n2id);
                    String rid = "edges_" + id + ":" + n2id;
                    int x2 = tempNode.getX();
                    int y2 = tempNode.getY();
                    // Can reduce the specificity by casting to an integer
                    double weight = 100 * outNodes.getValue();
                    thisEdge.put("from", id);
                    thisEdge.put("to", n2id);
                    thisEdge.put("x1", x);
                    thisEdge.put("x2", x2);
                    thisEdge.put("y1", y);
                    thisEdge.put("y2", y2);
                    thisEdge.put("weight", weight);
                    reactions.put(rid, thisEdge);
                } catch (Exception e) {
                    System.err.println("Error with reaction: " + nid);
                }

            }
            // Want to add each of the edge weights to the reaction nodes
        }
        JSONObject metadata = new JSONObject();
        jsonMap.put("nodes", nodesJSON);
        jsonMap.put("edges", reactions);
        jsonMap.put("max_depth", max_depth);
        return jsonMap;
    }

}
