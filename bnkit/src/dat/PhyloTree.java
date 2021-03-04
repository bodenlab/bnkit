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
import java.util.*;

/**
 * Class to represent a phylogenetic tree.
 * Rooted, bifurcating tree for representing phylogenetic relationships.
 * Functionality includes labeling and traversing nodes; reading and writing to Newick format;
 * Programmers should note that almost all functionality is implemented through recursion.
 * @author mikael
 * @deprecated
 */
public class PhyloTree {

    private Node root = null; // the root of the tree
    private ArrayList<Node> nodeList = new ArrayList<>();

    /**
     * Private constructor for tree from a root with all nodes connected off that.
     * Use factory methods to construct trees.
     */
    private PhyloTree(ArrayList<Node> nodeList) {
        this.nodeList = nodeList;
    }

    public PhyloTree() {

    }

    public static PhyloTree load(String filename, String format) throws IOException {
        if (format.equalsIgnoreCase("newick") || format.equalsIgnoreCase("nwk"))
            return new PhyloTree().loadNewick(filename);
        else
            throw new RuntimeException("Unknown format: " + format);
    }

    /**
     * String representation in the Newick format.
     *
     * @return string representation of tree
     */
    @Override
    public String toString() {
        return root.toString();
    }

    /**
     * Print the tree as a string with values associated with each node.
     *
     * @return a text string
     */
    public String printValues() {
        return root.printValue();
    }

    /**
     * Convert nodes to an array, breadth-first
     *
     * @return
     */
    public Node[] toNodesBreadthFirst() {
        Node[] arr = new Node[nodeList.size()];
        nodeList.toArray(arr);
        return arr;
    }

    public String[] toStringsBreadthFirst() {
        Node[] nodes = toNodesBreadthFirst();
        String[] arr = new String[nodes.length];
        for (int i = 0; i < nodes.length; i++)
            arr[i] = nodes[i].label.toString();
        return arr;
    }

    /**
     * Get root node of tree.
     *
     * @return the root of the tree
     */
    public Node getRoot() {
        return root;
    }

    /**
     * Find the node with the specified label.
     *
     * @param content label or label
     * @return matching node, or null if not found
     */
    public Node find(Object content) {
        return root.find(content);
    }


    public void removeInternalLabels() {
        root.removeInternalLabels();
    }

    /**
     * Set sequence of nodes that are matched (by name) by the sequences in the alignment.
     *
     * @param aln sequence alignment
     */
    public void setAlignment(EnumSeq.Alignment aln) {
        int nseqs = aln.getHeight();
        for (int i = 0; i < nseqs; i++) {
            EnumSeq seq = aln.getEnumSeq(i);
            Node node = find(seq.getName());
            if (node != null) {
                node.setSequence(seq);
            }
        }
    }

    /**
     * Find index of first comma at the current level (non-embedded commas are ignored) or end of string.
     *
     * @param str a Newick string
     * @return index of the first comma or end-of-string
     */
    private static int getComma(String str) {
        if (str.length() == 0)
            return -1;
        int mylevel = 0;
        char[] chararr = str.toCharArray();
        for (int i = 0; i < chararr.length; i++) {
            if (chararr[i] == '(') mylevel += 1;
            else if (chararr[i] == ')') mylevel -= 1;
            else if (chararr[i] == ',' && mylevel == 0) return i;
        }
        return str.length();
    }

    /**
     * Utility method to parse an embedded string on the Newick format.
     *
     * @param parent the parent of the current node
     * @return the root node of tree
     */
    private Node parseNewick(String newickStr, Node parent) {
        root = parseNewick(newickStr, parent, new ArrayList<>(), 0);
        return root;
    }

    /**
     * Utility method to parse an embedded string on the Newick format.
     *
     * @return the tree
     */
    public PhyloTree parseNewick(String newickStr) {
        root = parseNewick(newickStr, null, new ArrayList<>(), 0);
        return this;
    }

    /**
     * Helper function to parse a leaf (extent sequence) in a Newick file.
     *
     * @param str    The Newick String
     * @param parent Parent Node
     * @return
     */
    private Node parseLeafNewick(String str, Node parent) {
        Node node;
        String label;
        int splitIdx = str.indexOf(':'); // check if a distance is specified
        if (splitIdx == -1) {// no distance
            return new Node(str, parent, null);
        } else { // there's a distance
            label = str.substring(0, splitIdx).trim();
            try {
                double dist = Double.parseDouble(str.substring(splitIdx + 1, str.length()));
                if (dist == 0.0) {
                    dist = 0.00001;
                    //System.err.println("Distance 0.0 parsed in tree file (child=" + label + "): Representing distance as " + Double.toString(dist));
                }
                node = new Node(label, parent, dist);
                if (root == null) {
                    root = node;
                }
                return node;
            } catch (NumberFormatException ex) {
                throw new RuntimeException("Error: A distance value in your Newick file couldn't be parsed as a number  \n \nThe value was - " + str.substring(splitIdx + 1, str.length()));
            }
        }
    }


    /**
     * Helper function to parse an internal node (i.e. the template for an ancestor) in the
     * Newick file.
     *
     * @param embed   Part of Newick String containing the ancestral node
     * @param tail    End of the String
     * @param parent  Parent of the Node
     * @param nodeIds List of traversed NodeIds
     * @param count   Number of nodeIds visited
     * @return
     */
    private Node parseInternalNewick(String embed, String tail, Node parent, ArrayList<Integer> nodeIds, int count) {
        String label;
        Node node;
        int splitIdx = tail.indexOf(':'); // check if a distance is specified
        if (splitIdx == -1) { // no distance
            if (!tail.isEmpty() && tail.substring(0, tail.length() - 1) != null && !tail.substring(0, tail.length() - 1).isEmpty()) {
                label = tail.substring(splitIdx + 1, tail.length()).replace(";", "");
                node = new Node("N" + count + "_" + label, parent, null);
            } else {
                node = new Node("N" + count, parent, null);
            }
        } else { // there's a distance
            if (tail.substring(0, splitIdx) != null && !tail.substring(0, splitIdx).isEmpty()) {
                label = "N" + count + "_" + tail.substring(0, splitIdx);
            } else {
                label = "N" + count;
            }
            try {
                double dist = Double.parseDouble(tail.substring(splitIdx + 1, tail.length()).replace(";", ""));
                if (dist == 0.0) {
                    dist = 0.00001;
                    //System.err.println("Distance 0.0 parsed in tree file (child=" + label + "): Representing distance as " + Double.toString(dist));
                }
                node = new Node(label, parent, dist);
                if (root == null) {
                    root = node;
                }
            } catch (NumberFormatException ex) {
                throw new RuntimeException("Error: A distance value in your Newick file couldn't be parsed as a number  \n \nThe value was - " + tail.substring(splitIdx + 1, tail.length()).replace(";", ""));
            }
        }
        nodeIds.add(count);
        // find where the commas are, and create children of node
        int comma = getComma(embed);
        String toProcess;
        while (comma != -1) {
            toProcess = embed.substring(0, comma);
            //GOING TO HAVE TO PASS PARENT NODE WITH RECURSION TO RECORD IT
            // get unique ID to pass through
            while (nodeIds.contains(count)) {
                count++;
            }
            node.addChild(parseNewick(toProcess, node, nodeIds, count));
            if (comma + 1 > embed.length()) {
                break;
            }
            embed = embed.substring(comma + 1);
            comma = getComma(embed);
        }
        return node;
    }

    /**
     * Utility method for recursively parse an embedded string on the Newick format.
     * MB-Fix: fixed a bug that meant that labels were missing the last character.
     * (Only last node or any node if distance is not given.)
     *
     * @param parent the parent of the current node
     * @return the root node of tree
     */
    private Node parseNewick(String str, Node parent, ArrayList<Integer> nodeIds, int count) {
        Node node = null;
        str = str.replace("\t", "");
        int startIdx = str.indexOf('('); // start parenthesis
        int endIdx = str.lastIndexOf(')'); // end parenthesis
        if (startIdx == -1 && endIdx == -1) { // we are at leaf (no parentheses)
            node = parseLeafNewick(str, parent);
        } else if (startIdx >= 0 && endIdx >= 0) { // balanced parentheses
            String embed = str.substring(startIdx + 1, endIdx);
            String tail = str.substring(endIdx + 1, str.length());
            node = parseInternalNewick(embed, tail, parent, nodeIds, count);
        }
        if (!nodeList.contains(node)) {
            nodeList.add(node);
        }
        return node;
    }

    /**
     * Factory method to create a tree instance from a Newick formatted file.
     *
     * @param filename name of file
     * @return instance of tree
     */
    public PhyloTree loadNewick(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        StringBuilder sb = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null)
            sb.append(line.trim());

        root = parseNewick(sb.toString(), null); //null parent for root
        reader.close();
        return this;
    }


    /**
     * Perform parsimony calcs on this tree.
     * Uses a simple 0/1 substitution penalty based on equivalence or not of states
     *
     * @param names   the names of the leaves
     * @param symbols the symbols/states/values assigned to the leaves (same order as names)
     */
    public void setContentByParsimony(String[] names, Object[] symbols) {
        Map<String, Object> map = new HashMap<>();
        for (int i = 0; i < names.length; i++)
            map.put(names[i], symbols[i]);
        setContentByParsimony(map);
    }

    /**
     * Perform parsimony calcs on this tree.
     * Uses a simple 0/1 substitution penalty based on equivalence or not of states
     *
     * @param map contains key/value pairs for leaf names and states
     */
    public void setContentByParsimony(Map<String, Object> map) {
        // reset node values (otherwise stores previous uses)
        reset_values(this.getRoot());
        Set<Object> values = new HashSet<>(map.values());
        Object[] unique = values.toArray();
        Node root = this.getRoot();
        root.forwardParsimony(map, unique);
        root.backwardParsimony(unique);
    }

    /**
     * Perform parsimony calcs on this tree.
     * Uses a simple 0/1 substitution penalty based on equivalence or not of states
     *
     * @param map    contains key/value pairs for leaf names and states
     * @param unique array with the unique states in order of preference
     *               (which decides what states that are preferred when multiple states contribute to optimal solution)
     */
    public void setContentByParsimony(Map<String, Object> map, Object[] unique) {
        // reset node values (otherwise stores previous uses)
        reset_values(this.getRoot());
        Node root = this.getRoot();
        root.forwardParsimony(map, unique);
        root.backwardParsimony(unique);
    }

    private void reset_values(Node node) {
        if (node.getChildren().isEmpty())
            return;
        node.setValue(null);
        for (Node child : node.getChildren())
            reset_values(child);
    }

    public static void main(String[] args) {

        try {
            PhyloTree mytree = new PhyloTree().loadNewick("/Users/mikael/cloudstor/DHAD_Jan2019/dhad_clustal_08012019/r_7500_9112_dhad_01012018.nwk");
            mytree.removeInternalLabels();
            System.out.println(mytree.toString());
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
        //null parent for root
        PhyloTree phyloTree = new PhyloTree();
        Node root = phyloTree.parseNewick("((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2", null);
        System.out.println(root);
        //null parent for root
        root = phyloTree.parseNewick("(((E:3.9,F:4.5,A,B,C)ef:2.5,G:0.3)efg:7,x,z,q,w,e,r,t)", null);
        System.out.println(root);
        try {
            PhyloTree edge1 = phyloTree.loadNewick("/Users/mikael/simhome/ASR/edge1.nwk");
            System.out.println(edge1);
            Alignment aln = new Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/gap1.aln", Enumerable.aacid));
            edge1.setAlignment(aln);
            edge1.setContentByParsimony(aln.getNames(), aln.getGapColumn(1));
            System.out.println(edge1.printValues());
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Class for nodes that make up tree.
     * Note recursive definition. Supports any branching factor.
     */
    public static class Node {
        final private List<Node> children;  // the children of this node
        private Object label;               // arbitrary label of node
        // MB-Fix: now a list of values, previously just a single "value"
        private List<Object> values = null; // values of node
        private EnumSeq sequence = null;    // sequence
        private double[] scores = null;     // the optimal score for each parent value
        // MB-Fix: added third index, to allow multiple child symbols to contribute to the same parent score
        private int[][][] traceback = null; // [parent value][child branch] = [value in child that gives optimal value, next value in same child ...]
        private Double dist = null;         // optional distance (from this node to its parent)
        private List<Double> modelProb = new ArrayList<>(); // Every node has a single parent (bar root) so it can carry
        // the value for the edge
        private Node parent;

        /**
         * Construct node from label/label.
         *
         * @param label label/label
         */
        public Node(String label, Node parent, Double dist) {
            this.label = label;
            this.children = new ArrayList<>();
            this.parent = parent;
            this.dist = dist;
            this.values = new ArrayList<>();
        }


        public Object getLabel() {
            return label;
        }

        public void setLabel(String label) {
            this.label = label;
        }

        // MB-Fix: re-directs to a new function, defaults to return the first symbol if multiple available
        public Object getValue() {
            return getValue(0);
        }

        // MB-Fix: new function signature, using index to retrieve specific value
        public Object getValue(int index) {
            if (this.values == null)
                return null;
            else if (this.values.size() < index + 1)
                return null;
            return values.get(index);
        }

        public double[] getScores() {
            return scores;
        }

        // MB-Fix: new function, returns all values; should be run after parsimony
        public List<Object> getValues() {
            return values;
        }

        // MB-Fix: same signature, handles multiple calls to save each internally.
        // Note it no longer returns value, instead it returns true if new value was set, false, if value was already set (used to optimise traversal below)
        public boolean setValue(Object value) {
            if (value == null) {
                this.values = null;
                return true;
            }
            if (this.values == null)
                this.values = new ArrayList<>();
            else if (this.values.contains(value))
                return false;
            this.values.add(value);
            return true;
        }

        public List<Double> getModelProb() {
            return modelProb;
        }

        public Double getModelProb(int c) {
            return modelProb.get(c);
        }

        /**
         * When adding modelProb values you must iterate over the alignment in order so each value is added in
         * the position corresponding to the column it represents
         *
         * @param probability
         */
        public void addModelProb(Double probability) {
            this.modelProb.add(probability);
        }

        /**
         * Can only be used after the initial assignment of all modelProb values (using addLikelihood()) otherwise the list will be unpopulated
         * and you will get an indexOutOfBoundsException
         *
         * @param modelProb
         * @param column
         */
        public void setModelProb(Double modelProb, int column) {
            try {
                this.modelProb.set(column, modelProb);
            } catch (IndexOutOfBoundsException iob) {
                System.out.println("Model probability list must be initialised using addModelProb() prior to setting specific columns");
            }
        }

        public void setParent(Node parent) {
            this.parent = parent;
        }

        public Node getParent() {
            return parent;
        }

        public void setSequence(EnumSeq seq) {
            this.sequence = seq;
        }

        public EnumSeq getSequence() {
            return sequence;
        }

        /**
         * String representation of the node and its children (recursively) that uses the Newick format.
         *
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

        // MB-Fix: added to simplify printing of multiple values
        private static String concat(List<Object> values) {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < values.size(); i++) {
                Object y = values.get(i);
                sb.append(y + (i == values.size() - 1 ? "" : "/"));
            }
            return sb.toString();
        }

        // MB-Fix: added printing of multiple values
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
                return label.toString() + "=" + concat(values) + ((dist != null) ? (dstr) : (""));
            else
                return "(" + sb.toString() + ")" + label.toString() + "=" + concat(values) + ((dist != null) ? (dstr) : (""));
        }

        /**
         * Add child to node.
         *
         * @param child
         */
        public void addChild(Node child) {
            children.add(child);
        }


        /**
         * Remove child from node.
         *
         * @param child
         */
        public void removeChild(Node child) {
            children.remove(child);
        }


        /**
         * Retrieve all the children of the node.
         *
         * @return
         */
        public Collection<Node> getChildren() {
            return children;
        }

        /**
         * Find node by label/label.
         * Searches the tree recursively using the current node as root.
         *
         * @param label
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

        private void removeInternalLabels() {
            if (children.size() > 0) {
                this.setLabel("");
                for (Node child : children) {
                    child.removeInternalLabels();
                }
            }
        }


        /**
         * Set the distance for this node (from this node to its parent)
         *
         * @param dist the distance
         */
        public void setDistance(double dist) {
            this.dist = dist;
        }

        /**
         * Retrieve the distance of this node (from this node to its parent)
         *
         * @return the distance
         */
        public double getDistance() {
            if (this.dist == null)
                throw new RuntimeException("Node " + this + " with content " + label + " does not have a distance");
            return this.dist;
        }

        /**
         * Internal function that operates recursively to first initialise each node (forward),
         * stopping only once a value has been assigned to the node,
         * then to propagate scores from assigned nodes to root (backward).
         * Extended to deal with multiple values contributing to optimal scores
         *
         * @param assign map with assignments (named nodes and corresponding values)
         * @param unique all possible values, in order; order matters if SET_ONE_TARGET_PARSIMONY is true
         * @return the scores of the unique values at the root
         */
        protected double[] forwardParsimony(Map<String, Object> assign, Object[] unique) {
            this.scores = new double[unique.length]; // a score for each possible value; pushed from leaves to this node; starts at 0
            this.values = null;
            Object sym = assign.get(this.label);     // check if this node has a definitive state (symbol) assigned to it (via its name/label)
            if (sym != null) {                       // yes, this node is instantiated
                int index;                           // index of the assigned state
                for (index = 0; index < unique.length; index++)
                    if (sym.equals(unique[index]))   // found it
                        break;
                setValue(unique[index]);             // definitive state is assigned to node
                Arrays.fill(this.scores, Double.POSITIVE_INFINITY); // the default score is Infinity
                this.scores[index] = 0; // the actual assignment is scored 0, all others impossible, so positive infinity
                return this.scores;     // done with this leaf node
            } else { // this node is NOT instantiated
                if (this.children == null) { // leaf node, no children, ouch... problem
                    throw new RuntimeException("Leaf " + this + " has not been assigned a value");
                } else { // recurse into children nodes...
                    // determine scores contributed by each child (cscores) BEFORE substitution penalties
                    double[][] cscores = new double[children.size()][];
                    for (int c = 0; c < children.size(); c++) {
                        Node child = children.get(c);
                        cscores[c] = child.forwardParsimony(assign, unique); // one score from child c for each symbol
                    }
                    // traceback array needs to hold all child symbol indices that contribute to (indexed) parent symbol score via (indexed) child
                    this.traceback = new int[unique.length][children.size()][];
                    double best_parent_score = Double.POSITIVE_INFINITY; // need to work out best parent score for "scores"; they all start at 0 (init at allocation above)
                    for (int c = 0; c < children.size(); c++) {
                        // loop through each possible parent assignment, record what symbol in each child that best supports this (adding substitution penalties as we go)
                        for (int i = 0; i < scores.length; i++) {
                            double best_score = Double.POSITIVE_INFINITY;
                            int best_cnt = 0; // this is how many symbols in child that need to be recorded (multiple are possible)
                            // next, loop through each possible value in this child to find best parent score (may be the score from multiple origins)
                            for (int j = 0; j < cscores[c].length; j++) {
                                int subst_penalty = (i == j ? 0 : 1);
                                double parent_score = cscores[c][j] + subst_penalty;
                                if (parent_score < best_score) { // new best, reset count to 1
                                    best_score = parent_score;
                                    best_cnt = 1;
                                } else if (parent_score == best_score) { // new equal best, add +1 to count
                                    best_cnt++;
                                } // else, let this parent score slip
                            }
                            // now we know what the best parent score is; work out all assignments in child c that give it (could be multiple)
                            traceback[i][c] = new int[best_cnt]; // allocate space to hold their indices
                            int k = 0; // index for holding possible origins
                            for (int j = 0; j < cscores[c].length; j++) { // loop through each possible child symbol, again adding substitution penalties
                                int subst_penalty = (i == j ? 0 : 1);
                                double parent_score = cscores[c][j] + subst_penalty;
                                if (parent_score == best_score)
                                    traceback[i][c][k++] = j;
                            }
                            scores[i] += best_score; // the best we can do with parent symbol i, from child c
                        } // finished the score for a parent i, for one child c
                    } // finished all children here
                    return this.scores; // done, return parent scores (recursing them as child scores up the tree)
                }
            }
        }

        public double getParsimonyScore() {
            return getParsimonyScore(null);
        }

        private double getParsimonyScore(Object parent_state) {
            Object my_state = this.values.get(0);
            if (my_state == null)
                return 0;
            double score = 0;
            if (parent_state != null)
                score = (parent_state == my_state ? 0 : 1);
            for (Node child : children) {
                score += child.getParsimonyScore(my_state);
            }
            return score;
        }
    /*        private double getParsimonyScore(Object parent_state, int parent_state_index) {
                Object my_state = this.values.get(0);
                if (my_state == null)
                    return 0;
                double score = 0;
                if (parent_state != null)
                    score = (parent_state == my_state ? 0 : 1);
                for (int c = 0; c < children.size(); c ++) {
                    Node child = children.get(c);
                    int y = random.nextInt(traceback[child_state_index][c].length);
                    int traceback[child_state_index][c][y];
                    score += child.getParsimonyScore(my_state);
                }
                return score;
            }*/

        /**
         * If one or multiple solutions should be identified
         */
        static public boolean SET_ONE_TARGET_PARSIMONY = false;
        /**
         * If solutions are listed in order or if they should be randomised
         */
        static public boolean SET_RANDOM_PARSIMONY = false;

        private static Random random = new Random(System.currentTimeMillis());

        /**
         * Shuffles the elements in the array randomly.
         * Code based on methods in java.util.Collections.shuffle();
         */
        protected static void shuffle(int[] array) {
            int count = array.length;
            for (int i = count; i > 1; i--) {
                swap(array, i - 1, random.nextInt(i));
            }
        }

        /**
         * Helper function to shuffle.
         *
         * @param array
         * @param i
         * @param j
         */
        private static void swap(int[] array, int i, int j) {
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
        }

        /**
         * Generates an array with values from 0 to specified n.
         * The order is either ascending or shuffled, depending on parameter shuffled.
         *
         * @param n        number of elements, populating the array with 0 up to n - 1
         * @param shuffled if true, the array will be returned with elements in random order
         * @return the array
         */
        protected static int[] range(int n, boolean shuffled) {
            int[] ret = new int[n];
            for (int i = 0; i < n; i++)
                ret[i] = i;
            shuffle(ret);
            return ret;
        }

        /**
         * Two backwardParsimony functions handle the "root" and internal nodes, respectively.
         * This handles the node as if it was root; goes by scores assigned to states/symbols.
         *
         * @param unique array with the states/symbols that can be assigned, in the same order as provided to forwardParsimony.
         */
        protected void backwardParsimony(Object[] unique) {
            int best_index = 0; // find one index with the best score (could be many but one is enough)
            for (int i = 1; i < scores.length; i++) {
                if (scores[i] < scores[best_index])
                    best_index = i;
            }
            // Go through each score and when it is "best", recurse into each child, propagating the state for the score
            // This iteration could be randomised (SET_RANDOM_PARSIMONY=true) so that the states are assigned in different order
            boolean butt_out = false;
            for (int parent_index : range(scores.length, SET_RANDOM_PARSIMONY)) {
                if (scores[parent_index] == scores[best_index]) {
                    // now we know the index of (one of) the optimal parent state/s
                    if (setValue(unique[parent_index])) { // will check so that the state is assigned only once... setValue returns false if already set
                        if (this.children != null) {      // recurse into children nodes...
                            for (int c = 0; c < children.size(); c++) {
                                Node child = children.get(c);
                                // Go through each optimal child state
                                for (int child_index : range(traceback[parent_index][c].length, SET_RANDOM_PARSIMONY)) {
                                    int best_index_in_child = traceback[parent_index][c][child_index]; // index of optimal child state
                                    child.backwardParsimony(best_index_in_child, unique);
                                    // if we are interested in only one optimal solution, we can butt out...
                                    if (SET_ONE_TARGET_PARSIMONY) {
                                        butt_out = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if (butt_out)
                    break;
            }
        }

        /**
         * Two backwardParsimony functions handle the "root" and internal nodes, respectively.
         * This function handles internal nodes and should only be called by the formes; goes by state assigned to parent.
         *
         * @param assign_me the state assigned to this node (by its parent)
         * @param unique    array with the states/symbols that can be assigned, in the same order as provided to forwardParsimony.
         */
        private void backwardParsimony(Object assign_me, Object[] unique) {
            int best_index;
            for (best_index = 0; best_index < unique.length; best_index++) {
                if ((assign_me == null && unique[best_index] == null) || (assign_me != null & assign_me.equals(unique[best_index])))
                    break;
            }
            backwardParsimony(best_index, unique);
        }

        /**
         * Two backwardParsimony functions handle the "root" and internal nodes, respectively.
         * This function handles internal nodes and should only be called by the formes; goes by state assigned to parent.
         *
         * @param value_index the index of the state assigned to this node
         * @param unique      array with the states/symbols that can be assigned, in the same order as provided to forwardParsimony.
         */
        private void backwardParsimony(int value_index, Object[] unique) {
            // now we know the index of the parent
            if (setValue(unique[value_index])) { // will return false if already set; if so, no point in recursing into children
                if (this.children != null) {
                    for (int c = 0; c < children.size(); c++) { // iterate over children...
                        // For each child: choose the first optimal state, or go through each. Ordered randomly, or dictated by traceback/parent.
                        for (int child_state_index : range(traceback[value_index][c].length, SET_RANDOM_PARSIMONY)) {
                            int best_index_in_child = traceback[value_index][c][child_state_index];
                            Node child = children.get(c);
                            child.backwardParsimony(best_index_in_child, unique);
                            if (SET_ONE_TARGET_PARSIMONY)
                                break;
                        }
                    }
                }
            }
        }
    }
}
