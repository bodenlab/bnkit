package reconstruction;

import bn.BNet;
import bn.BNode;
import bn.ctmc.PhyloBNet;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JTT;
import dat.EnumSeq;
import dat.Enumerable;
import dat.PhyloTree;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.*;

/**
 * Created by aesseb on 10/6/2015.
 * @deprecated
 */
public class Analysis {

    private PhyloTree tree = new PhyloTree();
    private PhyloBNet pbn; //only to be used for navigating branches
    private EnumSeq.Alignment<Enumerable> aln;
    private Map<Object, Map<Object, Integer>> transitions = new HashMap<>();

    /**
     * Constructor following on immediately from performing a reconstruction
     * Currently runs the methods to model transitions across all columns in the alignment
     *
     * @param reconstruction
     */
    public Analysis(ASR reconstruction) {
        tree = reconstruction.getTree();
        pbn = reconstruction.getPbn();
        aln = reconstruction.getAln();

        evolEdges();
        getDistanceMatrix();

        //Write the output
        try {
            PrintWriter writer = new PrintWriter("newick_transformed.txt", "UTF-8");
            for (int col = 0; col < aln.getWidth(); col++) {
                List<PhyloTree.Node> visited = new ArrayList<>();
                Collection<PhyloTree.Node> children = tree.getRoot().getChildren();
                Collection<PhyloTree.Node> newChildren = new ArrayList<>(children);
                String output = "(";
                String newick = constructNewick(col, tree.getRoot(), newChildren, visited, output);
//                writer.println(col);
                writer.println(newick);
            }
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
    }

    /**
     * Constructor which assumes the reconstructed sequence and reconstructed tree have been saved and
     * can be loaded
     *
     * @param file_tree
     * @param file_aln
     */
    public Analysis(String file_tree, String file_aln) {
        try {
            tree = tree.loadNewick(file_tree); //load tree - tree not in Newick
            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive

            List<EnumSeq.Gappy<Enumerable>> seqs = EnumSeq.Gappy.loadClustal(file_aln, Enumerable.aacid);
            aln = new EnumSeq.Alignment<>(seqs);
            pbn = PhyloBNet.create(tree, new JTT());

        } catch (IOException ex) {
            ex.printStackTrace();
        }

//        BNet[] models = createNtrainBN();
//        inferEdgeValues(models); //Stored within each node

    }

    /**
     * Constructor which takes the output from ASR
     * FIXME need to work on this
     *
     * @param result_file
     */
    public Analysis(String result_file, int type) {
        switch (type) {
            case 1:
                //JSON output
                break;
            case 2:
                //XML output
                break;
        }

    }

    public void evolEdges() {

        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //Explore all nodes to get all branches
        for (int n = 0; n < nodes.length; n++) {
            PhyloTree.Node child = nodes[n];
            PhyloTree.Node parent = child.getParent(); //works for bifurcating trees as a child has a single parent
            if (parent == null)
                continue; //node is root so no branch to record
            BNet model = pbn.getBN();
            for (int c = 0; c < aln.getWidth(); c++) {
                Object childState = child.getSequence().get()[c];
                Object parentState = parent.getSequence().get()[c];
                SubstNode aachd = (SubstNode) model.getNode(child.getLabel().toString());
                double prob = 0.0;
                if (childState == null && parentState == null) {
                    //FIXME - issues with gap character
                    System.out.println();
                } else if (childState == null) { //childState is gap and has to be handled
                    double[][] probs = aachd.getModel().getProbs(aachd.getTime()); //get probability matrix for node
                    int index = aachd.getModel().getDomain().getIndex(parentState); //get index of parent state
                    double[] colProb = probs[index];
                    for (int i = 0; i < colProb.length; i++) {
                        if (i == index)
                            continue;
                        prob += colProb[i];
                    }
                } else if (parentState == null){
                    //FIXME - issues with gap character
                    System.out.println();
                } else {
                    Object[] ap = {parentState};
                    prob = aachd.get(ap, childState);
                }
                child.addModelProb(prob); //Child stores the likelihoods for the edge connecting it to its parent
            }
        }
    }


    public String constructNewick(int col, PhyloTree.Node node, Collection<PhyloTree.Node> children, List<PhyloTree.Node> visited, String output) {

        if (children.iterator().hasNext()) {
            PhyloTree.Node newNode = children.iterator().next(); //identify child of interest
            children.remove(newNode);
            visited.add(newNode); //add child node to visited list
            if (newNode.getChildren().size() > 0) { //this child has children, recurse
                Collection<PhyloTree.Node> nextChildren = newNode.getChildren(); //get next set of children to explore
                Collection<PhyloTree.Node> newChildren = new ArrayList<>(nextChildren);
                output += "(";
                return constructNewick(col, newNode, newChildren, visited, output);
            } else { //we have reached the end of this branch
                if (children.iterator().hasNext()) { //there is another child so we need a comma to separate them
                    output += newNode.getLabel() + ":" + newNode.getModelProb(col) + ",";
                } else { //this is the last child so close brackets
                    output += newNode.getLabel() + ":" + newNode.getModelProb(col) + ")";
                }
                return constructNewick(col, node, children, visited, output); //recurse with same node to explore any other children
            }
        } else { //we are back tracking up the tree
            if (node.getParent() != null) { //to handle root
                PhyloTree.Node parent = node.getParent(); //identify parent from bnet in phylo tree
                Collection<PhyloTree.Node> back = parent.getChildren(); //get the children of the parent of current node
                visited.add(node); //add node to visited list
                Collection<PhyloTree.Node> newChildren = new ArrayList<>(back);
                newChildren.removeAll(visited);
                if (newChildren.iterator().hasNext()) {
                    output += node.getLabel() + ":" + node.getModelProb(col) + ",";
                } else {
                    output += node.getLabel() + ":" + node.getModelProb(col) + ")";
                }
                return constructNewick(col, parent, newChildren, visited, output); //the parent goes back to being the node of interest
                //potentially incorrect use of recursion...but it works
            } else { //explored all branches so finish
                output += node.getLabel();
                output += ";";
                return output;
            }
        }
    }

    public void getDistanceMatrix(){

        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst();
        for (int i = 0; i < aln.getWidth(); i++) {
            try {
                PrintWriter writer = new PrintWriter("distanceMatrix_" + i + ".txt", "UTF-8");
                Double[][] matrix = new Double[nodes.length][nodes.length];//include row name
                int l1 = 0;
                for (PhyloTree.Node n1 : nodes) {
                    int l2 = 0;
                    for (PhyloTree.Node n2 : nodes) {
                        double distance = getDistance(n1, n2, i);
                        matrix[l1][l2] = distance;
                        matrix[l2][l1] = distance;
                        l2++;
                    }
                    l1++;
                }
                writer.write("col_" + i + "\t");
                for (int a = 0; a < nodes.length; a++) { //write header
                    if (a == nodes.length - 1) {
                        writer.write((String)nodes[a].getLabel() + "\n");
                    } else {
                        writer.write((String) nodes[a].getLabel() + "\t");
                    }
                }
                for (int j = 0; j < matrix.length; j++) {
                    writer.write((String)nodes[j].getLabel() + "\t");
                    for (int k = 0; k < matrix[j].length; k++) {
                        if (k == matrix[j].length - 1) {
                            writer.write(matrix[j][k] + "\n");
                        } else {
                            writer.write(matrix[j][k] + "\t");
                        }
                    }
                }
                writer.close();
            } catch (FileNotFoundException fnf) {
                System.out.println(fnf.getStackTrace());
            } catch (UnsupportedEncodingException use) {
                System.out.println(use.getStackTrace());
            }
        }
    }

    public Double getDistance(PhyloTree.Node n1, PhyloTree.Node n2, int col) {

        if (n1.equals(n2))
            return 0.0;
        if (n1 == tree.getRoot()) {
            PhyloTree.Node temp = n1;
            n1 = n2;
            n2 = temp;
        }

        String n1Lab = (String) n1.getLabel();
        BNode node1 = pbn.getBN().getNode(n1Lab);
        List<BNode> n1Ancs = pbn.getBN().getAncestors(node1);
        List<BNode> n1Decs = pbn.getBN().getDescendants(node1);

        String n2Lab = (String) n2.getLabel();
        BNode node2 = pbn.getBN().getNode(n2Lab);
        List<BNode> n2Ancs = pbn.getBN().getAncestors(node2);

        double distance = 0.0;
        if (n1Ancs.contains(node2)) {
            distance += n1.getModelProb().get(col);
            for (BNode n1a : n1Ancs) {
                if (n1a.equals(node2)) { //Hit node of interest so don't add any more to distance
                    break;
                }
                distance += tree.find(n1a.getName()).getModelProb().get(col);
            }
        } else if (n1Decs.contains(node2)) {
            //Need to look up from node 2 because of way likelihood is stored + you will have issues with parent having
            //two children and picking the 'right' distance to add
            distance += n2.getModelProb().get(col);
            for (BNode n2a : n2Ancs) {
                if (n2a.equals(node1)) { //Hit node of interest so don't add any more to distance
                    break;
                }
                distance += tree.find(n2a.getName()).getModelProb().get(col);
            }
        } else {
            //Node pair is on other side of tree so sum all ancestors of both nodes
            distance += n1.getModelProb().get(col);
            for (BNode n1a : n1Ancs) {
                if (!tree.find(n1a.getName()).equals(tree.getRoot())) {
                    distance += tree.find(n1a.getName()).getModelProb().get(col);
                }
            }
            distance += n2.getModelProb().get(col);
            for (BNode n2a : n1Ancs) {
                if (!tree.find(n2a.getName()).equals(tree.getRoot())) {
                    distance +=tree.find(n2a.getName()).getModelProb().get(col);
                }
            }
        }
        return distance;
    }


    /**
     * A method to explore all branches of a phylogenetic tree based on
     * EXTANT nodes.
     *
     * @return map of branches keyed by leaf node
     */
    public Map<Object, List<BNode>> exploreBranches() {

        //Structure to store the set of branches in the tree
        Map<Object, List<BNode>> result = new HashMap<>();

        //A branch is generated for each extant node in the tree
        PhyloTree.Node[] extNodes = getExtantNodes();
        for (int n = 0; n < extNodes.length; n++) {
            String nodeLab = extNodes[n].getLabel().toString();
            BNode node = pbn.getBN().getNode(nodeLab);
            List<BNode> ancs = pbn.getBN().getAncestors(node);
            result.put(nodeLab, ancs);
        }
        return result;
    }

    /**
     * A method to explore all branches of a phylogenetic tree based on
     * PROVIDED nodes.
     *
     * @param nodes - list of nodes which represent the 'leaves' of the branches to be explored
     * @return map of branches keyed by leaf node
     */
    public Map<Object, List<BNode>> exploreBranches(PhyloTree.Node[] nodes) {

        //Structure to store the set of branches in the tree
        Map<Object, List<BNode>> result = new HashMap<>();

        for (int n = 0; n < nodes.length; n++) {
            String nodeLab = (String) nodes[n].getLabel();
            BNode node = pbn.getBN().getNode(nodeLab);
            List<BNode> ancs = pbn.getBN().getAncestors(node);
            result.put(nodeLab, ancs);
        }
        return result;
    }


    /**
     * Based on the branches of interest from exploreBranches(), generate an alignment for each branch including
     * ancestral sequence. Then calculate the entropy for each column in the alignment
     *
     * @param branchMap - keyed by leaf node of branch
     * @return an array of entropy scores for each branch across the alignment
     */
    public Object[][] calculateBranchEntropy(Map<Object, List<Object>> branchMap) {
        //For each branch in tree, get the sequence for each node and create an alignment
        int r = 0; //count rows
        Object[][] store = new Object[branchMap.size()][];
        for (Map.Entry<Object, List<Object>> b : branchMap.entrySet()) {
            List<EnumSeq.Gappy<Enumerable>> asrs = new ArrayList<>();
            Object ext = b.getKey();
            List<Object> brNodes = b.getValue();
            PhyloTree.Node[] tNodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
            for (PhyloTree.Node tNode : tNodes) {
                if (brNodes.contains(tNode.getLabel())) {
                    asrs.add((EnumSeq.Gappy<Enumerable>) tNode.getSequence());
                }
            }
            EnumSeq.Alignment aln_asr = new EnumSeq.Alignment(asrs);
            //For each column in branch specific alignment (based on ancestors), calculate an entropy score
            Object[] row = new Object[aln_asr.getWidth() + 1]; //hard code branch name into array
            row[0] = ext; //Need to know which set of scores belongs to which branch
            for (int i = 1; i < aln_asr.getWidth() + 1; i++) {
                Object[] col = aln_asr.getColumn(i);
                double entropy = getShannonEntropy(col);
                row[i] = entropy;
            }
            store[r] = row;
            r++;
        }
        return store;
    }

    /**
     * Take the 2D array generated by one of the branch processing methods and write it to a file
     *
     * @param branchValues
     * @param file_name    - include only the base of the name - no file extensions
     */
    public void writeBranches(Object[][] branchValues, String file_name) {
        //save the results to a file
        try {
            PrintWriter writer = new PrintWriter(file_name + ".txt", "UTF-8");
            for (int b = 0; b < branchValues.length; b++) {
                for (int d = 1; d < branchValues[b].length; d++) { //ignore the label in column 0
                    if (d < branchValues[b].length - 1) {
                        writer.print(branchValues[b][d] + "\t");
                    } else {
                        writer.print(branchValues[b][d]);
                    }
                }
                writer.print("\n");
            }
            writer.close();
            PrintWriter writerLab = new PrintWriter(file_name + "lab.txt", "UTF-8");
            for (int b = 0; b < branchValues.length; b++) {
                for (int d = 0; d < branchValues[b].length; d++) { //include label in column 0
                    if (d < branchValues[b].length - 1) {
                        writerLab.print(branchValues[b][d] + "\t");
                    } else {
                        writerLab.print(branchValues[b][d]);
                    }
                }
                writerLab.print("\n");
            }
            writerLab.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
    }

    /**
     * Calculate the Shannon Entropy given a column
     * from an alignment
     *
     * @param column
     * @return entropy
     */
    public double getShannonEntropy(Object[] column) {
        Object[] aacid = Enumerable.aacid.getValues();
        List<Object> col = Arrays.asList(column);
        double entropy = 0;
        for (int k = 0; k < aacid.length; k++) {
            Object aa = aacid[k];
            double count = (double) Collections.frequency(col, aa);
            double prob = count / column.length;
            if (prob > 0.0) {
                entropy = entropy + (prob * java.lang.Math.log(prob));
            }
        }
        if (entropy != 0.0) {
            entropy = -entropy;
        }
        return entropy;
    }

    private PhyloTree.Node[] getInternalNodes(){
        String rootname = tree.getRoot().toString();
        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
        List<PhyloTree.Node> intNodes = new ArrayList<>();
        for (PhyloTree.Node node : nodes) {
            String nodename = node.toString();
            if (rootname.equals(nodename))
                continue;
            //Not root and not leaf
            if (node.getChildren().toArray().length > 0){
                intNodes.add(node);
            }
        }
        PhyloTree.Node[] intNodesA = new PhyloTree.Node[nodes.length - intNodes.size() - 1];
        return intNodes.toArray(intNodesA);
    }

    private PhyloTree.Node[] getExtantNodes(){
        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
        List<PhyloTree.Node> extNodes = new ArrayList<>();
        for (PhyloTree.Node node : nodes) {
            //Only leaf nodes have no children
            if (node.getChildren().toArray().length == 0){
                extNodes.add(node);
            }
        }
        PhyloTree.Node[] intNodesA = new PhyloTree.Node[nodes.length - extNodes.size() - 1];
        return extNodes.toArray(intNodesA);
    }

    private static String replacePunct(String str) {
        return str.replace('.', '_');
    }

}
