package bn.reconstruction;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.Predef;
import bn.alg.EM;
import bn.alg.LearningAlg;
import bn.alg.VarElim;
import bn.ctmc.PhyloBNet;
import bn.ctmc.matrix.JTT;
import bn.file.BNBuf;
import bn.node.CPT;
import bn.prob.EnumDistrib;
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

        BNet trainedNet = learnParameters();
        BNet[] models = createNtrainBN(trainedNet);
//        inferEdgeValues(models); //Stored within each node
        recordEdgeValues(models);
        transformLikelihood();

//        Map<Object, List<Object>> branches = exploreBranches();
//        recordBranchValues(branches);

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
            tree = PhyloTree.loadNewick(file_tree); //load tree - tree not in Newick
            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive

            List<EnumSeq.Gappy<Enumerable>> seqs = EnumSeq.Gappy.loadClustal(file_aln, Enumerable.aacid);
            aln = new EnumSeq.Alignment<>(seqs);
            pbn = PhyloBNet.create(tree, new JTT());

        } catch (IOException ex) {
            ex.printStackTrace();
        }

        BNet[] models = createNtrainBN();
        inferEdgeValues(models); //Stored within each node

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

        if (n1.equals(n2)) {
            return 0.0;
        }

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

        boolean high;
        boolean low;
        boolean side;
        double distance = 0.0;
        if (n1Ancs.contains(node2)) {
            low = true;
            distance += n1.getLikelihood().get(col);
            for (BNode n1a : n1Ancs) {
                if (n1a.equals(node2)) { //Hit node of interest so don't add any more to distance
                    break;
                }
                distance += tree.find(n1a.getName()).getLikelihood().get(col);
            }
        } else if (n1Decs.contains(node2)) {
            high = true;
            //Need to look up from node 2 because of way likelihood is stored + you will have issues with parent having
            //two children and picking the 'right' distance to add
            distance += n2.getLikelihood().get(col);
            for (BNode n2a : n2Ancs) {
                if (n2a.equals(node1)) { //Hit node of interest so don't add any more to distance
                    break;
                }
                distance += tree.find(n2a.getName()).getLikelihood().get(col);
            }
        } else {
            side = true;
            //Node pair is on other side of tree so sum all ancestors of both nodes
            distance += n1.getLikelihood().get(col);
            for (BNode n1a : n1Ancs) {
                if (!tree.find(n1a.getName()).equals(tree.getRoot())) {
                    distance += tree.find(n1a.getName()).getLikelihood().get(col);
                }
            }
            distance += n2.getLikelihood().get(col);
            for (BNode n2a : n1Ancs) {
                if (!tree.find(n2a.getName()).equals(tree.getRoot())) {
                    distance +=tree.find(n2a.getName()).getLikelihood().get(col);
                }
            }
        }
        return distance;
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

    public BNet learnParameters() {
        BNet net = network();

        Map<Object, Map<Object, List<Object>>> trainingData = new HashMap<>();
        int edges = tree.toNodesBreadthFirst().length - 1; //edges = total nodes - 1
        //Generate all transitions across all columns in alignment
        for (int a = 0; a < aln.getWidth(); a++) { //for every column in alignment
            BNet curNet = network();//create fresh BN for modelling
            PhyloTree.Node root = tree.getRoot();
            Collection<PhyloTree.Node> children = root.getChildren();
            Map<Object, List<Object>> store = new HashMap<>();
            List<PhyloTree.Node> visited = new ArrayList<>();
            Collection<PhyloTree.Node> newChildren = new ArrayList<>(children);
            getTransitions(a, root, newChildren, store, root, visited);
            trainingData.put(a, store);
        }

        Object[][] columnTransitions = new Object[edges * aln.getWidth()][];
        try {
            PrintWriter transitionWriter = new PrintWriter("raw_transitions.txt", "UTF-8");
            //Transform this into a 2d array to be passed to the network storing ALL transitions

            int e = 0;
            for (Map.Entry<Object, Map<Object, List<Object>>> colTran : trainingData.entrySet()) {
                Map<Object, List<Object>> store = colTran.getValue();
                for (Map.Entry<Object, List<Object>> trans : store.entrySet()) {
                    Object aaParent = trans.getKey();
                    for (Object aaChild : trans.getValue()) {
                        Object[] tran = {null, null, aaParent, aaChild};//data must reflect order of nodes passed to em.train
                        columnTransitions[e] = tran;
                        e++;
                        transitionWriter.println(aaParent + "\t" + aaChild);

                        //Count transition
                        if (transitions.containsKey(aaParent)) {
                            if (transitions.get(aaParent).containsKey(aaChild)) {
                                transitions.get(aaParent).put(aaChild, transitions.get(aaParent).get(aaChild) + 1);
                            } else {
                                transitions.get(aaParent).put(aaChild, 1);
                            }

                        } else {
                            transitions.put(aaParent, new HashMap<Object, Integer>());
                            transitions.get(aaParent).put(aaChild, 1);
                        }
                    }
                }
            }
            transitionWriter.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }



        //Train the network with full set of transitions
        List<BNode> bNodes = Arrays.asList(net.getNode("Parent"), net.getNode("Child"), net.getNode("AAparent"), net.getNode("AAchild"));
        LearningAlg em = new EM(net);
        em.train(columnTransitions, bNodes);
        BNBuf.save(net, "full_network.out");

        BNode aapar = net.getNode("AAparent");
        BNode aachild = net.getNode("AAchild");

        try {
            PrintWriter writer = new PrintWriter("transition_llhs.txt", "UTF-8");
            for (Object aap : Enumerable.aacid_ext.getValues()) {
                aapar.setInstance(aap);
                for (Object aac : Enumerable.aacid_ext.getValues()) {
                    aachild.setInstance(aac);
                    VarElim ve = new VarElim();
                    ve.instantiate(net);
                    double llh = ve.logLikelihood();
                    writer.println(aap + "\t" + aac + "\t" + llh);
                }
            }
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
        return net;
    }

    /**
     * Creates a Bayesian network for each column in the alignment. This network structure is specified in network().
     * When passing in a trainedNet, the distribution for AAparent and AAchild will be set and NOT TRAINED
     * (FIXME create method which has the network passed in?) Issue with copying the network
     * Each network is trained based on all transitions seen in phylogenetic tree (including ancestral sequence) for a
     * specific column in the alignment. Transitions are recorded using getTransitions().
     *
     * @param trainedNet    a network that has parameters for aaParent and aaChild
     * @return an array of trained Bayesian networks
     */
    public BNet[] createNtrainBN(BNet trainedNet) {
        BNet[] columnModels = new BNet[aln.getWidth()];
        Map<Object, Map<Object, List<Object>>> trainingData = new HashMap<>();
        int edges = tree.toNodesBreadthFirst().length - 1; //edges = total nodes - 1

        //Extract parameters for AAparent and AAchild
        BNode aaPar = trainedNet.getNode("AAparent");
        BNode aaChld = trainedNet.getNode("AAchild");
        BNode parent = trainedNet.getNode("Parent");
        BNode child = trainedNet.getNode("Child");

        for (int a = 0; a < aln.getWidth(); a++) { //for every column in alignment
            BNet curNet = network(); //generate a fresh network to store the new values
            BNode aaParT = curNet.getNode("AAparent");
            BNode aaChldT = curNet.getNode("AAchild");

            for (String paramPar : parent.getVariable().getParams().split(";")) {
                parent.setInstance(paramPar);
                Object[] parKey = {parent.getInstance()};
                Distrib aaParDistrib = aaPar.getDistrib(parKey);
                aaParT.put(aaParDistrib, parKey);

                child.setInstance(paramPar);
                Object[] chldKey = {child.getInstance()};
                aaChldT.put(aaParDistrib, chldKey);
            }

            aaParT.setTrainable(false);
            aaChldT.setTrainable(false);

            PhyloTree.Node root = tree.getRoot();
            Collection<PhyloTree.Node> children = root.getChildren();
            Map<Object, List<Object>> store = new HashMap<>();
            List<PhyloTree.Node> visited = new ArrayList<>();
            Collection<PhyloTree.Node> newChildren = new ArrayList<>(children);
            getTransitions(a, root, newChildren, store, root, visited);
            trainingData.put(a, store);

            Object[][] columnTransitions = new Object[edges][];
            int e = 0;
            for (Map.Entry<Object, List<Object>> trans : store.entrySet()) {
                Object aaParent = trans.getKey();
                for (Object aaChild : trans.getValue()) {
                    Object[] tran = {null, null, aaParent, aaChild};//data must reflect order of nodes passed to em.train
                    columnTransitions[e] = tran;
                    e++;
                }
            }

            System.out.println(a);

            List<BNode> bNodes = Arrays.asList(curNet.getNode("Parent"), curNet.getNode("Child"), curNet.getNode("AAparent"), curNet.getNode("AAchild"));
            LearningAlg em = new EM(curNet);
            em.train(columnTransitions, bNodes);
            columnModels[a] = curNet;
//            BNBuf.save(curNet, "network_col_" + a + ".out");

            BNode aapar = curNet.getNode("AAparent");
            BNode aachild = curNet.getNode("AAchild");

//            try {
//                PrintWriter writer = new PrintWriter("col_"+a+"_transition_llhs.txt", "UTF-8");
//                for (Object aap : Enumerable.aacid_ext.getValues()) {
//                    aapar.setInstance(aap);
//                    for (Object aac : Enumerable.aacid_ext.getValues()) {
//                        aachild.setInstance(aac);
//                        VarElim ve = new VarElim();
//                        ve.instantiate(curNet);
//                        double llh = ve.logLikelihood();
//                        writer.println(aap + "\t" + aac + "\t" + llh);
//                    }
//                }
//                writer.close();
//            } catch (FileNotFoundException fnf) {
//                System.out.println(fnf.getStackTrace());
//            } catch (UnsupportedEncodingException use) {
//                System.out.println(use.getStackTrace());
//            }
        }
        return columnModels;
    }

    /**
     * Creates a Bayesian network for each column in the alignment. This network is specified in network().
     * (FIXME create method which has the network passed in?) Issue with copying the network
     * Each network is trained based on all transitions seen in phylogenetic tree (including ancestral sequence) for a
     * specific column in the alignment. Transitions are recorded using getTransitions().
     *
     * @return an array of trained Bayesian networks
     */
    public BNet[] createNtrainBN() {
        BNet[] columnModels = new BNet[aln.getWidth()];
        Map<Object, Map<Object, List<Object>>> trainingData = new HashMap<>();
        int edges = tree.toNodesBreadthFirst().length - 1; //edges = total nodes - 1
        for (int a = 0; a < aln.getWidth(); a++) { //for every column in alignment
            BNet curNet = network();//create fresh BN for modelling
            PhyloTree.Node root = tree.getRoot();
            Collection<PhyloTree.Node> children = root.getChildren();
            Map<Object, List<Object>> store = new HashMap<>();
            List<PhyloTree.Node> visited = new ArrayList<>();
            Collection<PhyloTree.Node> newChildren = new ArrayList<>(children);
            getTransitions(a, root, newChildren, store, root, visited);
            trainingData.put(a, store);

            Object[][] columnTransitions = new Object[edges][];
            int e = 0;
            for (Map.Entry<Object, List<Object>> trans : store.entrySet()) {
                Object aaParent = trans.getKey();
                for (Object aaChild : trans.getValue()) {
                    Object[] tran = {null, null, aaParent, aaChild};//data must reflect order of nodes passed to em.train
                    columnTransitions[e] = tran;
                    e++;
                }
            }

            System.out.println(a);

            List<BNode> bNodes = Arrays.asList(curNet.getNode("Parent"), curNet.getNode("Child"), curNet.getNode("AAparent"), curNet.getNode("AAchild"));
            LearningAlg em = new EM(curNet);
            em.train(columnTransitions, bNodes);
            columnModels[a] = curNet;
//            BNBuf.save(curNet, "network_col_" + a + ".out");
            System.out.println();
        }
        return columnModels;
    }

    /**
     * Based on the trained network for each column, instantiate the parent and child node and get the likelihood of
     * the model in that state. Record this value for every edge, across every column
     *
     * @param models
     */
    public void inferEdgeValues(BNet[] models) {

        Map<Integer, Map<String, Double>> store = new HashMap<>();
        for (int c = 0; c < models.length; c++) {
            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //Explore all nodes to get all branches
            for (int n = 0; n < nodes.length; n++) {
                if (nodes[n].getParent() == null) {
                    continue; //node is root so no branch to record
                } else {
                    PhyloTree.Node child = nodes[n];
                    PhyloTree.Node parent = child.getParent(); //works for bifurcating trees as a child has a single parent

                    Object childState = child.getSequence().get()[c];
                    Object parentState = parent.getSequence().get()[c];
                    char gap = '-';
                    if (childState == null)
                        childState = gap;
                    if (parentState == null)
                        parentState = gap;
                    BNet model = models[c];
                    model.getNode("AAparent").setInstance(parentState);
                    model.getNode("AAchild").setInstance(childState);
                    VarElim ve = new VarElim();
                    ve.instantiate(model);
                    double llh = ve.logLikelihood();
                    child.addLikelihood(llh); //Child stores the likelihoods for the edge connecting it to its
                }
            }
        }
    }

    public void recordEdgeValues(BNet[] models) {

        try {
            PrintWriter writer = new PrintWriter("col_llhs.txt", "UTF-8");

            Map<Integer, Map<String, Double>> store = new HashMap<>();

            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //Explore all nodes to get all branches
            for (int n = 0; n < nodes.length; n++) {
                PhyloTree.Node child = nodes[n];
                PhyloTree.Node parent = child.getParent(); //works for bifurcating trees as a child has a single parent
                for (int c = 0; c < models.length; c++) {
                    if (nodes[n].getParent() == null) {
                        continue; //node is root so no branch to record
                    } else {
                        Object childState = child.getSequence().get()[c];
                        Object parentState = parent.getSequence().get()[c];
                        char gap = '-';
                        if (childState == null)
                            childState = gap;
                        if (parentState == null)
                            parentState = gap;
                        BNet model = models[c];
                        model.getNode("AAparent").setInstance(parentState);
                        model.getNode("AAchild").setInstance(childState);
                        VarElim ve = new VarElim();
                        ve.instantiate(model);
                        double llh = ve.logLikelihood();
                        child.addLikelihood(llh); //Child stores the likelihoods for the edge connecting it to its
                        if (c == 0) {
                            writer.write(nodes[n].getLabel() + "\t" + llh);
                        } else {
                            writer.write("\t" + llh);
                        }
                    }
                }
                if (nodes[n].getParent() == null)
                    continue; //node is root so no branch to record
                writer.write("\n");
            }
            writer.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
        }
    }

    /**
     * To display likelihood values they must be transformed so they can be placed in a Newick string as
     * a branch length
     */
    public void transformLikelihood() {

        List<Double> original = new ArrayList<>();
        List<Double> transformed = new ArrayList<>();
        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst();
        for (int c = 0; c < aln.getWidth(); c++) { //for every column in alignment
            List<Double> colVals = new ArrayList<>();
            for (int n = 0; n < nodes.length; n++) { //for every node in network get the value for the column
                if (nodes[n].getLikelihood().size() > 0) //if node is not root
                    colVals.add(nodes[n].getLikelihood(c));
            }
            double max = Collections.max(colVals);
            double min = Collections.min(colVals);
            for (int n = 0; n < nodes.length; n++) { //for every node in network replace the llh value with transformed value
                if (nodes[n].getLikelihood().size() > 0) { //if node is not root
                    double cl = nodes[n].getLikelihood(c);
                    original.add(cl);
                    double transform = (max - cl) / (max - min); //treat max as new '0' value
                    if (Double.isNaN(transform)) { //all edges have equal likelihoods -> division by 0
                        transform = 0.01; //FIXME what is a good neutral branch size here?
                    }
                    else {
                        transform += 0.01; //FIXME - can't have branch size equal to 0.0 either?
                    }
                    transformed.add(transform);
                    nodes[n].setLikelihood(transform, c);
                }
            }
        }

        try {
            PrintWriter writer = new PrintWriter("original.txt", "UTF-8");
            for (double d : original) {
                writer.print(d + "\n");
            }
            writer.close();

            PrintWriter writerT = new PrintWriter("transform.txt", "UTF-8");
            for (double d : transformed) {
                writerT.print(d + "\n");
            }
            writer.close();
            writerT.close();
        } catch (FileNotFoundException fnf) {
            System.out.println(fnf.getStackTrace());
        } catch (UnsupportedEncodingException use) {
            System.out.println(use.getStackTrace());
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
                    output += newNode.getLabel() + ":" + newNode.getLikelihood(col) + ",";
                } else { //this is the last child so close brackets
                    output += newNode.getLabel() + ":" + newNode.getLikelihood(col) + ")";
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
                    output += node.getLabel() + ":" + node.getLikelihood(col) + ",";
                } else {
                    output += node.getLabel() + ":" + node.getLikelihood(col) + ")";
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


    /**
     * Create the network for modelling transitions across the alignment.
     * If aaParDistrib and aaChldDistrib are not null, set the distributions
     * for the nodes and set nodes to not trainable
     * @return
     */
    private BNet network() {
//        String[] stringVars = {"1","2","3", "4", "5", "6", "7"};
        String[] stringVars = {"1","2","3"};
        EnumVariable P = Predef.Nominal(stringVars, "Parent");
        EnumVariable C = Predef.Nominal(stringVars, "Child");
        EnumVariable AAP = Predef.AminoAcidExt("AAparent");
        EnumVariable AAC = Predef.AminoAcidExt("AAchild");

        CPT p = new CPT(P);
        CPT c = new CPT(C, P);
        CPT aap = new CPT(AAP, P);
        CPT aac = new CPT(AAC, C);

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
    private void getTransitions(int col, PhyloTree.Node node, Collection<PhyloTree.Node> children, Map<Object, List<Object>> store, PhyloTree.Node root, List<PhyloTree.Node> visited) {
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
                Collection<PhyloTree.Node> newChildren = new ArrayList<>(nextChildren);
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
                if (store.containsKey(parentState)) { //you've seen the parent state before
                    store.get(parentState).add(childState); //record parent -> child relationship
                } else { //this is a new parent state to record
                    List<Object> states = new ArrayList<>();
                    states.add(childState);
                    store.put(parentState, states); //record parent -> child relationship
                }
                getTransitions(col, node, children, store, root, visited); //recurse with same node to explore any other children
            }
        } else { //we are back tracking up the tree
            if (node.getParent() != null) { //to handle root
                PhyloTree.Node parent = node.getParent(); //identify parent from bnet in phylo tree
                Collection<PhyloTree.Node> back = parent.getChildren(); //get the children of the parent of current node
                visited.add(node); //add node to visited list
                Collection<PhyloTree.Node> newChildren = new ArrayList<>(back);
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
