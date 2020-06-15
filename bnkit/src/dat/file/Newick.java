package dat.file;

import dat.phylo.Tree;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Collection of utility methods for operating on data relevant to Newick files, and serving the phylo.Tree class.
 */
public class Newick {

    /**
     * Prevent instantiation of the class.
     */
    private Newick() {
    }

    /**
     * Find index of first comma at the current level (non-embedded commas are ignored) or end of string.
     * @param str a Newick string, e.g. "(((A:0.3,B:0.4):0.1,C:0.4):0.2,D:0.5);"
     * @return index of the first comma or end-of-string
     */
    public static int getComma(String str) {
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
     * Helper function to parse a leaf in a Newick file.
     *
     * @param str    The Newick String
     * @param parent Parent Node
     * @return
     */
    private static Tree.BranchPoint parseLeaf(String str, Tree.BranchPoint parent) {
        String label;
        int splitIdx = str.indexOf(':'); // check if a distance is specified
        if (splitIdx == -1) { // no distance
            return new Tree.BranchPoint(str, parent, null);
        } else { // there's a distance
            label = str.substring(0, splitIdx).trim();
            try {
                double dist = Double.parseDouble(str.substring(splitIdx + 1));
                if (dist == 0.0) {
                    dist = 0.00001;
                }
                return new Tree.BranchPoint(label, parent, dist);
            } catch (NumberFormatException ex) {
                throw new RuntimeException("A distance value couldn't be parsed as a number. The value is \"" + str.substring(splitIdx + 1)+ "\" about here: " + str);
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
    private static Tree.BranchPoint parseInternal(String embed, String tail, Tree.BranchPoint parent, ArrayList<Integer> nodeIds, int count) {
        String label;
        Tree.BranchPoint branchPoint;
        int splitIdx = tail.indexOf(':'); // check if a distance is specified
        if (splitIdx == -1) { // no distance
            if (!tail.isEmpty() && tail.substring(0, tail.length() - 1) != null && !tail.substring(0, tail.length() - 1).isEmpty()) {
                label = tail.substring(splitIdx + 1).replace(";", "");
                branchPoint = new Tree.BranchPoint(label, parent, null);
            } else {
                branchPoint = new Tree.BranchPoint("N" + count, parent, null);
            }
        } else { // there's a distance
            if (tail.substring(0, splitIdx) != null && !tail.substring(0, splitIdx).isEmpty()) {
                label = tail.substring(0, splitIdx);
            } else {
                label = "N" + count;
            }
            try {
                double dist = Double.parseDouble(tail.substring(splitIdx + 1).replace(";", ""));
                if (dist == 0.0) {
                    dist = 0.00001;
                }
                branchPoint = new Tree.BranchPoint(label, parent, dist);
            } catch (NumberFormatException ex) {
                throw new RuntimeException("A distance value couldn't be parsed as a number. The value is \"" + tail.substring(splitIdx + 1).replace(";", "") + "\" about here: " + tail);
            }
        }
        branchPoint.setAncestor(count);
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
            branchPoint.addChild(parse(toProcess, branchPoint, nodeIds, count));
            if (comma + 1 > embed.length()) {
                break;
            }
            embed = embed.substring(comma + 1);
            comma = getComma(embed);
        }
        return branchPoint;
    }

    /**
     * Utility method for recursively parse an embedded string on the Newick format.
     * MB-Fix: fixed a bug that meant that labels were missing the last character.
     * (Only last node or any node if distance is not given.)
     *
     * @param parent the parent of the current node
     * @return the root node of tree
     */
    private static Tree.BranchPoint parse(String str, Tree.BranchPoint parent, ArrayList<Integer> nodeIds, int count) {
        Tree.BranchPoint branchPoint;
        // str = str.replace("\t", "");
        int startIdx = str.indexOf('('); // start parenthesis
        int endIdx = str.lastIndexOf(')'); // end parenthesis
        if (startIdx == -1 && endIdx == -1) { // we are at leaf (no parentheses)
            branchPoint = Newick.parseLeaf(str, parent);
        } else if (startIdx >= 0 && endIdx >= startIdx) { // balanced parentheses
            String embed = str.substring(startIdx + 1, endIdx);
            String tail = str.substring(endIdx + 1);
            branchPoint = parseInternal(embed, tail, parent, nodeIds, count);
        } else {
            if (startIdx >=0)
                throw new RuntimeException("Missing \")\" in Newick string, before here: " + str.substring(startIdx));
            else
                throw new RuntimeException("Missing \"(\" in Newick string, matching this: " + str.substring(0, endIdx));
        }
        return branchPoint;
    }

    /**
     * Utility method to parse an embedded string on the Newick format, and attach to a given parent.
     * @param parent the parent of the current node (can be null if no parent)
     * @return the root node of sub-tree
     */
    private static Tree.BranchPoint parse(String newickStr, Tree.BranchPoint parent) {
        Tree.BranchPoint subroot = Newick.parse(newickStr, parent, new ArrayList<>(), 0);
        if (parent != null)
            parent.addChild(subroot);
        return subroot;
    }

    /**
     * Factory method to parse an embedded string on the Newick format and make a tree out of it.
     * @return the tree
     */
    public static Tree parse(String str) {
        // Prepare to extract (quoted) labels, and/or annotations to internal nodes.
        StringBuilder sb = new StringBuilder(); // building a new, clean Newick string
        Map<String, String> matchLabel = new HashMap<>(); // map to match managed label with actual label
        int ptr = 0; // where we are in the input string
        int startIdx1 = str.indexOf('\''); // start quote
        int endIdx1 = startIdx1 >= 0 ? str.substring(startIdx1 + 1).indexOf('\'') : -1; // end quote
        int startIdx2 = str.indexOf('['); // start block
        int endIdx2 = startIdx2 >= 0 ? str.substring(startIdx2 + 1).indexOf(']') : -1; // end block
        while ((startIdx1 >= 0 && endIdx1 >= 0) || (startIdx2 >=0 && endIdx2 >= 0)) {
            if (startIdx1 < startIdx2 && startIdx1 >= 0 || startIdx2 == -1) { // single quotes are first, may actually contain block
                sb.append(str, ptr, ptr + startIdx1); // the bit before the first quote
                String tok = str.substring(ptr + startIdx1, ptr + startIdx1 + 1 + endIdx1 + 1);
                sb.append(tok);
                matchLabel.put(tok, tok.substring(1, tok.length() - 1));
                ptr += startIdx1 + 1 + endIdx1 + 1;
            } else { // block
                sb.append(str, ptr, ptr + startIdx2); // the bit before the first quote
                String tok = str.substring(ptr + startIdx2, ptr + startIdx2 + 1 + endIdx2 + 1);
                sb.append(tok);
                int idx = tok.indexOf("label=\"");
                if (idx >= 0) {
                    int idxx = tok.substring(idx + 7).indexOf("\"");
                    String extracted = tok.substring(idx + 7, idx + 7 + idxx);
                    matchLabel.put(tok, extracted);
                } else
                    matchLabel.put(tok, tok.substring(1, tok.length() - 1));
                ptr += startIdx2 + 1 + endIdx2 + 1;
            }
            startIdx1 = str.substring(ptr).indexOf('\''); // start quote
            endIdx1 = startIdx1 >= 0 ? str.substring(ptr + startIdx1 + 1).indexOf('\'') : -1; // end quote
            startIdx2 = str.substring(ptr).indexOf('['); // start quote
            endIdx2 = startIdx2 >= 0 ? str.substring(ptr + startIdx2 + 1).indexOf(']') : -1; // end quote
        }
        if (startIdx1 >= 0)
            throw new RuntimeException("Non-matched quotes, about here: " + str.substring(ptr + startIdx1));
        if (startIdx2 >= 0)
            throw new RuntimeException("Non-matched block, about here: " + str.substring(ptr + startIdx2));
        sb.append(str, ptr, str.length());
        Tree.BranchPoint root = Newick.parse(sb.toString(), null);
        Integer count = 0;
        List<Tree.BranchPoint> all = root.getSubtree();
        for (Tree.BranchPoint bp : all) {
            String preferred = matchLabel.get(bp.getLabel());
            if (preferred != null)
                bp.setLabel(preferred);
            if (bp.getID() == null) // if leaf, and label is not set, OR if ancestor and ancestor ID is not set
                bp.setAncestor(count ++);
        }
        Tree t = new Tree(root);
        return t;
    }

    /**
     * Factory method to create a tree instance from a Newick formatted file.
     *
     * @param filename name of file
     * @return instance of tree
     */
    public static Tree load(String filename) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        StringBuilder sb = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null)
            sb.append(line.trim());
        reader.close();
        return Newick.parse(sb.toString());
    }

    private static char[] forbidden = new char[] {':','(',')',','}; // characters that must be tracked before Newick parsing
    private static char[] embed = new char[] {':','(',')',',',' ','@','#','[',']','&','$','{','}','=','!','*','+'}; // characters that will be embedded with single quotes
    private static String wrap(String label) {
        boolean found = false;
        for (char ch : embed) {
            if (label.indexOf(ch) != -1) {
                found = true;
                break;
            }
        }
        if (found)
            return '\'' + label + '\'';
        else
            return label;
    }
    private static boolean isInvalid(String label) {
        for (char ch : forbidden) {
            if (label.indexOf(ch) != -1)
                return true;
        }
        return false;
    }

    public static int MODE_DEFAULT  = 0;
    public static int MODE_ANCESTOR = 1;
    public static int MODE_STRIPPED = 2;

    public static String sprint(Tree.BranchPoint node, int MODE) {
        StringBuilder sb = new StringBuilder();
        String dstr = null;
        List<Tree.BranchPoint> children = node.getChildren();
        int cnt = 0;
        for (Tree.BranchPoint child : children) {
            sb.append(sprint(child, MODE));
            if (++cnt < children.size())
                sb.append(",");
        }
        String label;
        if (!node.isLeaf() && (MODE == MODE_ANCESTOR || MODE == MODE_STRIPPED))
            label = MODE == MODE_STRIPPED ? "" : "N" + node.getAncestor();
        else
            label = MODE == MODE_STRIPPED ? node.getLabel().toString() : wrap(node.getLabel().toString());
        if (MODE == MODE_STRIPPED)
            if (isInvalid(label))
                throw new RuntimeException("Invalid label: \"" + label + "\". Please use alternative format or DEFAULT mode.");
        try {
            double dist = node.getDistance();
            dstr = ":" + dist;
            if (node.isLeaf())
                return label + dstr;
            else
                return "(" + sb.toString() + ")" + label + dstr;
        } catch (RuntimeException e) { // distance is not set
            if (node.isLeaf())
                return label;
            else
                return "(" + sb.toString() + ")" + label;
        }
    }

    public static void save(Tree tree, String filename) throws IOException {
        save(tree, filename, MODE_DEFAULT);
    }

    public static void save(Tree tree, String filename, int MODE) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
        String s = sprint(tree.getRoot(), MODE);
        bw.write(s + ";\n");
        bw.close();
    }

}
