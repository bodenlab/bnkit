/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vis;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author ariane
 */
public class Node {
    private int x;
    private int y;
    private final Map<Character, Double> graph;
    private Map<Character, Double> seqChars;
    private final int id;
    private Map<Integer, Double> outedges;
    
    public Node(int id, int x, int y, Map<Character, Double> graph, 
                        Map<Integer, Double> outedges, Map<Integer, Character> seqChars) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.graph = graph;
        this.outedges = outedges;
        this.seqChars = makeDistSeqChars(seqChars);
    }
    
    /**
     * Makes a distribution of the sequence characters so that we can easily
     * display it in a pie chart when the amino acid hasn't been chosen yet.
     * @param rawSeq
     * @return 
     */
    private Map<Character, Double> makeDistSeqChars(Map<Integer, Character> rawSeq) {
        seqChars = new HashMap<>();

        int numSequences = rawSeq.size();
        for (Character base : graph.keySet())
            seqChars.put(base, graph.get(base)*numSequences);
        /*for (Map.Entry<Integer, Character> list : rawSeq.entrySet()) {
            
            // Try and see if we have already stored the character otherwise add
            // it
            try {
               int count = seqChars.get(list.getValue());
               // If we could get it we increment the count of seqChars
               count ++;
               seqChars.put(list.getValue(), count);
            } catch(Exception e) {
                seqChars.put(list.getValue(), 1); // The initial count
            }
        }*/
        return seqChars;
    }
    
    public Map<Character, Double> getSeq() {
        return seqChars;
    }
           
    public int getX() {
        return x;
    }
    
    public int getID() {
        return id;
    }
    
        
    public void setY(int y) {
        this.y =  y;
    }
    
    public int getY() {
        return y;
    }
    
    public Map<Integer, Double> getOutedges() {
        return outedges;
    }
    
    public String toStr() {
        return ("id:" + id + ", x:" + x + ", y:" + y);
    }
    public Map<Character, Double> getGraph() {
        return graph;
    }
}
