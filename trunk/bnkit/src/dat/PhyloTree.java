/*
 * binfkit -- software for bioinformatics
 * Copyright (C) 2014  M. Boden et al.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package dat;

import dat.EnumSeq.Alignment;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Class to represent a phylogenetic tree.
 * Rooted, bifurcating tree for representing phylogenetic relationships.
 * Functionality includes labeling and traversing nodes; reading and writing to Newick format;
 * Programmers should note that almost all functionality is implemented through recursion.
 * @author mikael
 */
public class PhyloTree {
    
    final private Node root; // the root of the tree
    static int count;
    
    /**
     * Private constructor for tree from a root with all nodes connected off that.
     * Use factory methods to construct trees.
     * @param root root node
     */
    private PhyloTree(Node root) {
        this.root = root;
    }

    /**
     * String representation in the Newick format.
     * @return string representation of tree
     */
    @Override
    public String toString() {
        return root.toString();
    }
    
    public String printValues() {
        return root.printValue();
    }
    
    public Node[] toNodesBreadthFirst() {
        List<Node> done = new ArrayList<>();
        List<Node> queue = new ArrayList<>();
        queue.add(root);
        done.add(root);
        expandNodes(done, queue);
        Node[] arr = new Node[done.size()];
        done.toArray(arr);
        return arr;
    }
    
    private void expandNodes(List<Node> done, List<Node> queue) {
        if (queue.isEmpty())
            return;
        Node head = queue.remove(0);
        for (Node child : head.getChildren()) {
            done.add(child);
            queue.add(child);
        }
        expandNodes(done, queue);
    }
    
    public String[] toStringsBreadthFirst() {
        Node[] nodes = toNodesBreadthFirst();
        String[] arr = new String[nodes.length];
        for (int i = 0; i < nodes.length; i ++) 
            arr[i] = nodes[i].label.toString();
        return arr;
    }
    
    /**
     * Get root node of tree.
     * @return the root of the tree
     */
    public Node getRoot() {
        return root;
    }
    
    /**
     * Find the node with the specified label.
     * @param content label or label
     * @return matching node, or null if not found
     */
    public Node find(Object content) {
        return root.find(content);
    }
    
    /**
     * Set sequence of nodes that are matched (by name) by the sequences in the alignment.
     * @param aln sequence alignment
     */
    public void setAlignment(EnumSeq.Alignment aln) {
        int nseqs = aln.getHeight();
        for (int i = 0; i < nseqs; i ++) {
            EnumSeq seq = aln.getEnumSeq(i);
            Node node = find(seq.getName());
            if (node != null) {
                node.setSequence(seq);
            }
        }
    }
    
    /**
     * Find index of first comma at the current level (non-embedded commas are ignored) or end of string.
     * @param str a Newick string
     * @return index of the first comma or end-of-string
     */
    private static int getComma(String str) {
        if (str.length() == 0)
            return -1;
        int mylevel = 0;
        char[] chararr = str.toCharArray();
        for (int i = 0; i < chararr.length; i ++) {
            if (chararr[i] == '(') mylevel += 1;
            else if (chararr[i] == ')') mylevel -= 1;
            else if (chararr[i] == ',' && mylevel == 0) return i;
        }
        return str.length();
    }
    
    /**
     * Utility method for recursively parse an embedded string on the Newick format.
     * @param str text on Newick format
     * @return the root node of tree
     */
    private static Node parseNewick(String str) {
        Node node = null;
        int start_index = str.indexOf('('); // start parenthesis
        int end_index = str.lastIndexOf(')'); // end parenthesis
        if (start_index == -1 && end_index == -1) { // we are at leaf (no parentheses)
            int split_index = str.indexOf(':'); // check if a distance is specified
            if (split_index == -1) // no distance
                node = new Node(str);
            else { // there's a distance
                node = new Node(str.substring(0, split_index));
                double dist = Double.parseDouble(str.substring(split_index + 1, str.length()));
                node.setDistance(dist);
            }
        } else if (start_index >= 0 && end_index >= 0) { // balanced parentheses
            //end_index = str.length() - end_index - 1; // correct index to refer from start instead of end of string
            String embed = str.substring(start_index + 1, end_index);
            String tail = str.substring(end_index + 1, str.length());
            int split_index = tail.indexOf(':'); // check if a distance is specified
            if (split_index == -1) { // no distance
                node = new Node("N" + count + "_" + tail.substring(0, tail.length() - 1));
                count += 1;
            } else { // there's a distance
                node = new Node("N" + count + "_" + tail.substring(0, split_index));
                double dist = Double.parseDouble(tail.substring(split_index + 1, tail.length()));
                node.setDistance(dist);
                count +=1;
            }
            // find where the commas are, and create children of node
            int comma = getComma(embed);
            while (comma != -1) {
                String process_me = embed.substring(0, comma);
                node.addChild(parseNewick(process_me));
                if (comma + 1 > embed.length())
                    break;
                embed = embed.substring(comma + 1);
                comma = getComma(embed);
            }
        }
        return node;
    }
    
    /**
     * Factory method to create a tree instance from a Newick formatted file.
     * @param filename name of file
     * @return instance of tree
     */
    public static PhyloTree loadNewick(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        StringBuilder sb = new StringBuilder();
        String line = null;
        int cnt = 1;
        while ((line = reader.readLine()) != null)
            sb.append(line.trim());
        String newick = sb.toString();
        Node root = parseNewick(newick);
        PhyloTree t = new PhyloTree(root);
        reader.close();
        return t;
    }
    
    public Object setContentByParsimony(String[] names, Object[] symbols) {
        Map<String, Object> map = new HashMap<>();
        for (int i = 0; i < names.length; i ++)
            map.put(names[i], symbols[i]);
        return setContentByParsimony(map);
    }
    
    public Object setContentByParsimony(Map<String, Object> map) {
        Set<Object> values = new HashSet<>(map.values());
        Object[] unique = values.toArray();
        Node root = this.getRoot();
        root.forwardParsimony(map, unique);
        return root.backwardParsimony(unique);
    }
    
    /**
     * Class for nodes that make up tree.
     * Note recursive definition. Supports any branching factor.
     */
    public static class Node {
        final private List<Node> children; // the children of this node
        final private Object label; // arbitrary label of node
        private Object value; // arbitrary value of node
        private EnumSeq sequence = null; // sequence
        private double[] scores = null; 
        private int[][] traceback = null; 
        private Double dist = null; // optional distance (from this node to its parent)
        
        /**
         * Construct node from label/label.
         * @param label label/label
         */
        public Node(Object label) {
            this.label = label;
            this.children = new ArrayList<>();
        }
        
        public Object getLabel() {
            return label;
        }
        
        public Object getValue() {
            return value;
        }

        public Object setValue(Object value) {
            this.value = value;
            return value;
        }
        
        public void setSequence(EnumSeq seq) {
            this.sequence = seq;
        }
        
        public EnumSeq getSequence() {
            return sequence;
        }
        
        /**
         * String representation of the node and its children (recursively) that uses the Newick format.
         * @return string representation
         */
        public String toString() {
            StringBuilder sb = new StringBuilder();
            String dstr = null;
            int nchildren = children.size();
            int cnt = 0;
            for (Node child : children) {
                sb.append(child.toString());
                if (++cnt < nchildren)
                    sb.append(",");
            }
            if (dist != null)
                dstr = ":" + dist.toString(); 
            if (nchildren < 1) 
                return label.toString() + ((dist != null) ? (dstr) : (""));
            else 
                return "(" + sb.toString() + ")" + label.toString() + ((dist != null) ? (dstr) : (""));
        }
        
        public String printValue() {
            StringBuilder sb = new StringBuilder();
            String dstr = null;
            int nchildren = children.size();
            int cnt = 0;
            for (Node child : children) {
                sb.append(child.printValue());
                if (++cnt < nchildren)
                    sb.append(",");
            }
            if (dist != null)
                dstr = ":" + dist.toString(); 
            if (nchildren < 1) 
                return label.toString() + "_" + value.toString() + ((dist != null) ? (dstr) : (""));
            else 
                return "(" + sb.toString() + ")" + label.toString() + "_" + value.toString() + ((dist != null) ? (dstr) : (""));
        }
        
        /**
         * Add child to node.
         * @param child 
         */
        public void addChild(Node child) {
            children.add(child);
        }
        
        /**
         * Retrieve all the children of the node.
         * @return 
         */
        public Collection<Node> getChildren() {
            return children;
        }
        
        /**
         * Find node by label/label. 
         * Searches the tree recursively using the current node as root.
         * @param content 
         * @return the node that contains the specified label, or null if not found
         */
        public Node find(Object label) {
            if (this.label.equals(label)) 
                return this;
            else {
                for (Node child : children) {
                    Node node = child.find(label);
                    if (node != null)
                        return node;
                }
                return null;
            }
        }
        
        /**
         * Set the distance for this node (from this node to its parent)
         * @param dist the distance
         */
        public void setDistance(double dist) {
            this.dist = dist;
        }
        
        /**
         * Retrieve the distance of this node (from this node to its parent)
         * @return the distance
         */ 
        public double getDistance() {
            if (this.dist == null)
                throw new RuntimeException("Node " + this + " with content " + label + " does not have a distance");
            return this.dist;
        }

        /**
         * Internal function that operates recursively to first initialise each node (forward),
         * stopping only once a sequence has been assigned to the node, 
         * then to propagate scores from sequence assigned nodes to root (backward).
         * @param assign map with assignments (named nodes and corresponding values)
         * @param unique all possible values, in order
         */
        protected double[] forwardParsimony(Map<String, Object> assign, Object[] unique) {
            this.scores = new double[unique.length]; 
            Object sym = assign.get(label);
            if (sym != null) { // this node is instantiated
                int index;
                for (index = 0; index < unique.length; index ++) 
                    if (sym.equals(unique[index]))
                        break;
                setValue(unique[index]);
                Arrays.fill(this.scores, Double.POSITIVE_INFINITY);
                this.scores[index] = 0; // the actual symbol is scored 0, all others impossible do positive infinity
                return this.scores;
            } else { // this node is NOT instantiated
                if (this.children == null) { // no children, ouch...
                    throw new RuntimeException("Leaf " + this + " has not been assigned a value");
                } else { // recurse into children nodes...
                    double[][] cscores = new double[children.size()][];
                    this.traceback = new int[unique.length][children.size()];
                    for (int c = 0; c < children.size(); c ++) {
                        Node child = children.get(c);
                        cscores[c] = child.forwardParsimony(assign, unique);
                    }
                    
                    // now we have children scores (cscores)
                    for (int i = 0; i < scores.length; i ++) { // loop through each possible parent assignment, let's see what symbol in each child that best supports this
                        for (int c = 0;  c < children.size(); c ++) {
                            int best_index = 0;
                            double best_score = Double.POSITIVE_INFINITY;
                            for (int j = 0; j < cscores[c].length; j ++) { // loop through each possible child value to score parent value
                                if (cscores[c][j] + (i == j ? 0 : 1) < best_score) {
                                    best_score = cscores[c][j] + (i == j ? 0 : 1);
                                    best_index = j;
                                }
                            }
                            scores[i] += best_score; // the best we can do with this child
                            traceback[i][c] = best_index; // ... and that score is based on this child symbol
                        }
                    }
                    return this.scores;
                }
            }
        }
        
        protected Object backwardParsimony(Object[] unique) {
            return backwardParsimony(null, unique);
        }
        
        protected Object backwardParsimony(Object parent_symbol, Object[] unique) {
            int parent_index = 0;
            if (parent_symbol == null) { // root
                for (int i = 1; i < scores.length; i ++) {
                    if (scores[i] < scores[parent_index])
                        parent_index = i;
                }
            } else {
                for (parent_index = 0; parent_index < unique.length; parent_index ++) {
                    if (parent_symbol.equals(unique[parent_index]))
                        break;
                }
            }
            // now we know the index of the parent
            if (this.children != null) { // recurse into children nodes...
                for (int c = 0; c < children.size(); c ++) {
                    int child_index = traceback[parent_index][c];
                    Node child = children.get(c);
                    child.backwardParsimony(unique[child_index], unique);
                }
            }
            return setValue(unique[parent_index]);
        }
        
    }
    
    public static void main(String[] args) {
        Node root = parseNewick("((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2");
        System.out.println(root);
        root = parseNewick("(((E:3.9,F:4.5,A,B,C)ef:2.5,G:0.3)efg:7,x,z,q,w,e,r,t)");
        System.out.println(root);
        try {
            PhyloTree cyp3 = PhyloTree.loadNewick("/Users/mikael/simhome/ASR/gap_example.nwk");
            System.out.println(cyp3);
            Alignment aln = new Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/gap_example.aln", Enumerable.aacid));
            cyp3.setAlignment(aln);
            cyp3.setContentByParsimony(aln.getNames(), aln.getGapColumn(2));
            System.out.println(cyp3.printValues());
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
}
