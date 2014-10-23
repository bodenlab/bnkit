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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
    public String toString() {
        return root.toString();
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
            if (split_index == -1) // no distance
                node = new Node(tail);
            else { // there's a distance
                node = new Node(tail.substring(0, split_index));
                double dist = Double.parseDouble(tail.substring(split_index + 1, tail.length()));
                node.setDistance(dist);
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
    public static PhyloTree loadNewick(String filename) {
        try {
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
        } catch (IOException e) {
            return null;
        }
    }
    
    /**
     * Class for nodes that make up tree.
     * Note recursive definition. Supports any branching factor.
     */
    public static class Node {
        final private Set<Node> children; // the children of this node
        final private Object content; // arbitrary content/label of node
        private Double dist = null; // optional distance (from this node to its parent)
        
        /**
         * Construct node from content/label.
         * @param content content/label
         */
        public Node(Object content) {
            this.content = content;
            this.children = new HashSet<>();
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
                return content.toString() + ((dist != null) ? (dstr) : (""));
            else 
                return "(" + sb.toString() + ")" + content.toString() + ((dist != null) ? (dstr) : (""));
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
        public Set<Node> getChildren() {
            return children;
        }
        
        /**
         * Find node by content/label. 
         * Searches the tree recursively using the current node as root.
         * @param content 
         * @return the node that contains the specified content, or null if not found
         */
        public Node find(Object content) {
            if (this.content.equals(content)) 
                return this;
            else {
                for (Node child : children) {
                    Node node = child.find(content);
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
            return this.dist;
        }
    }
    
    public static void main(String[] args) {
        Node root = parseNewick("((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2");
        System.out.println(root);
        root = parseNewick("(((E:3.9,F:4.5,A,B,C)ef:2.5,G:0.3)efg:7,x,z,q,w,e,r,t)");
        System.out.println(root);
        PhyloTree cyp3 = PhyloTree.loadNewick("/Users/mikael/workspace/binf/data/cyp3.newick");
        System.out.println(cyp3);
    }
}
