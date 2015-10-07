package bn.reconstruction;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import bn.alg.EM;
import bn.alg.LearningAlg;
import bn.alg.VarElim;
import bn.ctmc.PhyloBNet;
import bn.ctmc.matrix.JTT;
import bn.node.CPT;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.PhyloTree;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.*;

/**
 * Created by aesseb on 10/6/2015.
 */
public class Analysis {

    private PhyloTree tree;
    private PhyloBNet pbn; //only to be used for navigating branches
    private EnumSeq.Alignment<Enumerable> aln;
    private Map<String, String> mapForNodes;

    private Map<String, PhyloTree.Node> branchParent;
    private Map<String, PhyloTree.Node> branchChild;

//    private List<EnumSeq.Gappy<Enumerable>> seqs;
//    private PhyloBNet[] pbnets;
//    private double[] R; //Rates at positions in alignment
//    private EnumDistrib[] margin_distribs; //Marginal distributions for nodes
//    private boolean use_sampled_rate = false;
//    private String asr_root; //Reconstructed sequence
//    private GammaDistrib gd; //Calculated gamma distribution
//    private BNet[] models; //model that will be trained for each column in alignment
//    private List<String> indexForNodes;

    public Analysis(ASR reconstruction) {
        tree = reconstruction.getTree();
        pbn = reconstruction.getPbn();
        aln = reconstruction.getAln();
        mapForNodes = reconstruction.getMapForNodes();

        labelEdges();
        BNet[] models = createNtrainBN();
        inferEdgeValues(models);
    }

    public Analysis(String file_tree, String file_aln){
        try {
            tree = PhyloTree.loadNewick(file_tree); //load tree - tree not in Newick
            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
            mapForNodes = new HashMap<>(); // Shortname --> Newick string for subtree
            for (PhyloTree.Node n : nodes) {
                mapForNodes.put(n.getLabel().toString(), replacePunct(n.toString()));
            }

            List<EnumSeq.Gappy<Enumerable>> seqs = EnumSeq.Gappy.loadClustal(file_aln, Enumerable.aacid);
            aln = new EnumSeq.Alignment<>(seqs);
            pbn = PhyloBNet.create(tree, new JTT());

        } catch (IOException ex) {
            ex.printStackTrace();
        }

        labelEdges();
        BNet[] models = createNtrainBN();
        inferEdgeValues(models);

    }

    /**
     * FIXME need to work on this
     * @param result_file
     */
    public Analysis(String result_file, int type){
        switch(type){
            case 1:
                //JSON output
                break;
            case 2:
                //XML output
                break;
        }

    }


    /**
     * A method to explore all branches of a phylogenetic tree based on
     * EXTANT nodes. An alignment is generated for each branch and its ancestors
     * and an entropy score is calculated for each column in the branch specific
     * alignment.
     * The scores are then written to a file which can be read by an R script to
     * make the output look pretty
     */
    public void exploreBranches() {

        //Structure to store the set of branches in the tree
        Map<Object, List<BNode>> result = new HashMap<>();

        //A branch is generated for each extant node in the tree
        PhyloTree.Node[] extNodes = getExtantNodes();
        for (int n = 0; n < extNodes.length; n++) {
            String nodeLab = (String)extNodes[n].getLabel();
            String nodeName = mapForNodes.get(nodeLab);
            BNode node = pbn.getBN().getNode(nodeName);
            List<BNode> ancs = pbn.getBN().getAncestors(node);
            result.put(nodeLab, ancs);
        }

        //This part is up for debate - currently each node in the branch is recorded
        //using the BNet node name not the phyloTree label
        //This converts them all to phyloTree labels
        //This makes it easier to extract sequence stored in the phyloTree structure
        Map<Object, List<Object>> modified = new HashMap<>();
        List<String> extant = new ArrayList<>();
        for (int i = 0 ; i < extNodes.length; i++) {
            String n = replacePunct(extNodes[i].toString());
            extant.add(n);
        }

        //Convert from BNode to phylo tree label
        for (Map.Entry<Object, List<BNode>> res : result.entrySet()) {
            Object key = res.getKey();
            List<BNode> branch = res.getValue();
            List<Object> branchNodes = new ArrayList<>();
            for (BNode node : branch) {
                String nodeTwo = node.toString();
                String nodeName = null;
                if (nodeTwo.contains("|")) {
                    nodeName = nodeTwo.substring(3, nodeTwo.indexOf("|"));
                } else {
                    nodeName = nodeTwo.substring(3, nodeTwo.length() - 1);
                }
                String nodeLabel = null;
                //FIXME there must be a better way to extract to node label!
                for (Map.Entry<String, String> e : mapForNodes.entrySet()) {
                    String value = e.getValue();
                    if (value.equals(nodeName)){
                        nodeLabel = e.getKey();
                        continue;
                    }
                }
                branchNodes.add(nodeLabel);
            }
            //Where the newly edited node labels are stored
            modified.put(key, branchNodes);
        }

        //For each branch in tree, get the sequence for each node and create an alignment
        List<Object> branchOrder = new ArrayList<>();
        int r = 0; //count rows
        Object[][] store = new Object[modified.size()][];
        for (Map.Entry<Object, List<Object>> b : modified.entrySet()) {
            List<EnumSeq.Gappy<Enumerable>> asrs = new ArrayList<>();
            Object ext = b.getKey();
            branchOrder.add(ext);//Important to know which branch is associated with which set of values
            List<Object> brNodes = b.getValue();
            PhyloTree.Node[] tNodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
            for (PhyloTree.Node tNode: tNodes) {
                if (brNodes.contains(tNode.getLabel())) {
                    asrs.add((dat.EnumSeq.Gappy<dat.Enumerable>)tNode.getSequence());
                }
            }
            EnumSeq.Alignment aln_asr = new EnumSeq.Alignment(asrs);
            //For each column in branch specific alignment (based on ancestors), calculate an entropy score
            Object[] row = new Object[aln_asr.getWidth()];
            for (int i = 0; i < aln_asr.getWidth(); i++) {
                Object[] col = aln_asr.getColumn(i);
                double entropy = getShannonEntropy(col);
                row[i] = entropy;
            }
            store[r] = row;
            r++;
        }

        //save the results to a file
        try {
            PrintWriter writer = new PrintWriter("test_matrix.txt", "UTF-8");
            for (int b = 0; b < store.length; b++) {
                for (int d = 0; d < store[b].length; d++) {
                    if (d < store[b].length -1) {
                        writer.print(store[b][d] + "\t");
                    } else {
                        writer.print(store[b][d]);
                    }
                }
                writer.print("\n");
            }
            writer.close();
            PrintWriter writerLab = new PrintWriter("test_matrix_lab.txt", "UTF-8");
            for (int b = 0; b < store.length; b++) {
                String label = (String)branchOrder.get(b);
                writerLab.print(label+"\t");
                for (int d = 0; d < store[b].length; d++) {
                    if (d < store[b].length -1) {
                        writerLab.print(store[b][d] + "\t");
                    } else {
                        writerLab.print(store[b][d]);
                    }
                }
                writerLab.print("\n");
            }
            writerLab.close();
        }
        catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
        System.out.println();
    }

    /**
     * Calculate the Shannon Entropy given a column
     * from an alignment
     * @param column
     * @return entropy
     */
    public double getShannonEntropy(Object[] column) {
        Object[] aacid = Enumerable.aacid.getValues();
        List<Object> col = Arrays.asList(column);
        double entropy = 0;
        for (int k = 0; k < aacid.length; k++) {
            Object aa = aacid[k];
            double count = (double)Collections.frequency(col, aa);
            double prob = count/column.length;
            if (prob > 0.0){
                entropy = entropy + (prob*java.lang.Math.log(prob));
            }
        }
        if (entropy != 0.0) {
            entropy = -entropy;
        }
        return entropy;
    }

    public BNet[] createNtrainBN() {
        BNet[] columnModels = new BNet[aln.getWidth()];
        Map<Object, Map<Object,List<Object>>> trainingData = new HashMap<>();
        int edges = tree.toNodesBreadthFirst().length - 1; //edges = total nodes - 1
        for (int a = 0; a < aln.getWidth(); a++) { //for every column in alignment
            BNet curNet = network();//create fresh BN for modelling
            PhyloTree.Node root = tree.getRoot();
            Collection<PhyloTree.Node> children = root.getChildren();
            Map<Object, List<Object>> store = new HashMap<>();
            Set<PhyloTree.Node> visited = new HashSet<>();
            Collection<PhyloTree.Node> newChildren = new HashSet<>(children);
            getTransitions(a, root, newChildren, store, root, visited);
            trainingData.put(a, store);

            Object[][] columnTransitions = new Object[edges][];
            int e = 0;
            for (Map.Entry<Object, List<Object>> trans : store.entrySet()) {
                Object aaParent = trans.getKey();
                for (Object aaChild : trans.getValue()) {
                    Object[] tran = {null, null, aaParent, aaChild};//nodes in curNet are ordered parent, aaparent, child, aachild - data must reflect this
                    columnTransitions[e] = tran;
                    e++;
                }
            }

            List<BNode> bNodes = Arrays.asList(curNet.getNode("Parent"), curNet.getNode("Child"), curNet.getNode("AAparent"), curNet.getNode("AAchild"));
            LearningAlg em = new EM(curNet);
            em.train(columnTransitions, bNodes);
            columnModels[a] = curNet;
//            BNBuf.save(curNet, "col_" + a + "_network.out");
            System.out.println();
        }
        return columnModels;
    }

    public void inferEdgeValues(BNet[] models) {

        Map<Integer, Map<String,Double>> store = new HashMap<>();
        for (int c = 0; c < models.length; c++) {
            Map<String, Double> branchRes = new HashMap<>();
            for (Map.Entry<String, PhyloTree.Node> branch : branchParent.entrySet()) {
                String bName = branch.getKey();
                Object parent = branch.getValue().getSequence().get()[c];
                Object child = branchChild.get(bName).getSequence().get()[c];

                BNet model = models[c];
                model.getNode("AAparent").setInstance(parent);
                model.getNode("AAchild").setInstance(child);
                VarElim ve = new VarElim();
                ve.instantiate(model);
                double llh = ve.logLikelihood();
                branchRes.put(bName, llh);
            }
            store.put(c, branchRes);
        }
    }

    private void labelEdges() {
        PhyloTree.Node root = tree.getRoot();
        Collection<PhyloTree.Node> children = root.getChildren();
        branchParent = new HashMap<>();
        branchChild = new HashMap<>();
        Collection<PhyloTree.Node> newChildren = new HashSet<>(children);
        int indx = 0;
        Set<PhyloTree.Node> visited = new HashSet<>();
        expandNode(root, newChildren, indx, branchParent, branchChild, visited);
        System.out.println();

    }
//
    private void expandNode(PhyloTree.Node node, Collection<PhyloTree.Node> children, int indx, Map<String, PhyloTree.Node> branchParent,Map<String, PhyloTree.Node> branchChild, Set<PhyloTree.Node> visited){
        //node = current node = parent
        if (children.iterator().hasNext()) {
            PhyloTree.Node newNode = children.iterator().next();
            //newNode = child
            children.remove(newNode); //we are now looking at this child so remove it from set
            visited.add(newNode); //add child node to visited list
            Object test = newNode.getChildren();
            if (newNode.getChildren().size() > 0) { //this child has children, recurse
                String bName = "b" + indx;
                indx++;
                branchParent.put(bName, node);
                branchChild.put(bName, newNode);
                Collection<PhyloTree.Node> nextChildren = newNode.getChildren(); //get next set of children to explore
                Collection<PhyloTree.Node> newChildren = new HashSet<>(nextChildren);
                expandNode(newNode, newChildren, indx, branchParent, branchChild, visited);
            }
            else { //we have reached the end of this branch
                String bName = "b" + indx;
                indx++;
                branchParent.put(bName, node);
                branchChild.put(bName, newNode);
                expandNode(node, children, indx, branchParent, branchChild, visited);
            }
        } else { //we are back tracking up the tree
            if (node != tree.getRoot()) { //to handle root

                //need to identify parent of current node - tricky when you're looking at a phylo tree not phylo BNet
                String nodeLab = (String)node.getLabel();
                String nodeName = mapForNodes.get(nodeLab); //get identifier of current node for bnet
                BNode bNode = pbn.getBN().getNode(nodeName); //find current node in bnet of phylo tree
                String bParent = bNode.getParents().get(0).getName(); //only ever single parent in phylo tree

                //FIXME there must be a better way to extract to node label!
                //get the node label of the parent of current node to enable back tracking
                String parentLab = null;
                for (Map.Entry<String, String> e : mapForNodes.entrySet()) {
                    String value = e.getValue();
                    if (value.equals(bParent)){
                        parentLab = e.getKey();
                        break;
                    }
                }
                PhyloTree.Node newNode = tree.find(parentLab); //identify parent from bnet in phylo tree
                Collection<PhyloTree.Node> back = newNode.getChildren(); //get the children of the parent of current node
                visited.add(node); //add node to visited list
                Collection<PhyloTree.Node> newChildren = new HashSet<>(back);
                newChildren.removeAll(visited);
                expandNode(newNode, newChildren, indx, branchParent, branchChild, visited);; //the parent goes back to being the node of interest
                //potentially incorrect use of recursion...but it works
            } else { //explored all branches so finish
                return;
            }
        }
    }

    private BNet network() {
        String[] stringVars = {"1","2","3","4","5"};
        EnumVariable P = Predef.Nominal(stringVars, "Parent");
        EnumVariable C = Predef.Nominal(stringVars, "Child");
        EnumVariable AAP = Predef.AminoAcidExt("AAparent");
        EnumVariable AAC = Predef.AminoAcidExt("AAchild");

        CPT p = new CPT(P);
        CPT c = new CPT(C, P);
        CPT aap = new CPT(AAP, P);
        CPT aac = new CPT(AAC, C);

//        aap.tieTo(aac);

        BNet bn = new BNet();
        bn.add(p,c,aap,aac);
        return bn;
    }


    /**
     * A recursive method to traverse the tree depth first and document all transitions from parent to child.
     * This method carries all variables through recursion
     * @param col column in alignment
     * @param node current node of interest
     * @param children children left to explore
     * @param store where results are kept
     * @param root root of the tree - as base case
     * @param visited record of which nodes have been visited
     * @return a map with parent as key and all possible children store in a list as the value
     */
    private void getTransitions(int col, PhyloTree.Node node, Collection<PhyloTree.Node> children, Map<Object, List<Object>> store, PhyloTree.Node root, Set<PhyloTree.Node> visited) {
        //node = parent
        //newNode = child
        char parentState = node.getSequence().toString().charAt(col); //get state of parent

        if (children.iterator().hasNext()) {
            PhyloTree.Node newNode = children.iterator().next(); //identify child of interest
            children.remove(newNode);
            char childState = newNode.getSequence().toString().charAt(col); //get state of child of interest
            visited.add(newNode); //add child node to visited list
            if (newNode.getChildren().size() > 0) { //this child has children, recurse
                Collection<PhyloTree.Node> nextChildren = newNode.getChildren(); //get next set of children to explore
                Collection<PhyloTree.Node> newChildren = new HashSet<>(nextChildren);
                if (store.containsKey(parentState)) { //you've seen the parent state before
                    store.get(parentState).add(childState); //record parent -> child relationship
                } else { //this is a new parent state to record
                    List<Object> states = new ArrayList<>();
                    states.add(childState);
                    store.put(parentState, states); //record parent -> child relationship
                }
                getTransitions(col, newNode, newChildren, store, root, visited); //recurse with child of interest as new parent
            }
            else { //we have reached the end of this branch
                store.get(parentState).add(childState);
                getTransitions(col, node, children, store, root, visited); //recurse with child of interest as new parent
            }
        } else { //we are back tracking up the tree
            if (node != root) { //to handle root

                //need to identify parent of current node - tricky when you're looking at a phylo tree not phylo BNet
                String nodeLab = (String)node.getLabel();
                String nodeName = mapForNodes.get(nodeLab); //get identifier of current node for bnet
                BNode bNode = pbn.getBN().getNode(nodeName); //find current node in bnet of phylo tree
                String bParent = bNode.getParents().get(0).getName(); //only ever single parent in phylo tree

                //FIXME there must be a better way to extract to node label!
                //get the node label of the parent of current node to enable back tracking
                String parentLab = null;
                for (Map.Entry<String, String> e : mapForNodes.entrySet()) {
                    String value = e.getValue();
                    if (value.equals(bParent)){
                        parentLab = e.getKey();
                        break;
                    }
                }
                PhyloTree.Node parent = tree.find(parentLab); //identify parent from bnet in phylo tree
                Collection<PhyloTree.Node> back = parent.getChildren(); //get the children of the parent of current node
                visited.add(node); //add node to visited list
                Collection<PhyloTree.Node> newChildren = new HashSet<>(back);
                newChildren.removeAll(visited);
                getTransitions(col, parent, newChildren, store, root, visited); //the parent goes back to being the node of interest
                //potentially incorrect use of recursion...but it works
            } else { //explored all branches so finish
                return;
            }
        }
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
