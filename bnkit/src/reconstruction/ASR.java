/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package reconstruction;

import bn.BNet;
import bn.BNode;
import bn.alg.*;
import bn.ctmc.PhyloBNet;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.Gap;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.PhyloTree;
import dat.Variable;

import java.io.*;
import java.util.*;

import dat.file.FastaWriter;
import json.JSONObject;


/**
 * @author Mikael Boden
 * @author Alex Essebier
 * @deprecated
 */
public class ASR {
    
    private PhyloTree tree = new PhyloTree();
    private PhyloBNet pbn; //only to be used for navigating branches
    private List<EnumSeq.Gappy<Enumerable>> seqs; //Ignores gaps in sequence
    private EnumSeq.Alignment<Enumerable> aln; //Store gaps in sequence
//    private PhyloBNet[] pbnets;
//    private PhyloBNet[] pbgapnets;
    private double[] R; //Rates at positions in alignment
    private EnumDistrib[] margin_distribs; //Marginal distributions for nodes
    private boolean use_sampled_rate = false;
    
    private List<String> indexForNodes;
    //Joint reconstruction of tree
    private Object[][] asr_matrix; //storage of all reconstructed sequences
    private double sampled_rate = // sampled rate, copy from a previous 1.0-rate run
	0.15599004226404184;
    
    public ASR(String file_tree, String file_aln, String inference) {
        loadData(file_tree, file_aln);
        PhyloBNet[] pbnets = createNetworks();
        if (inference.equals("Joint")) {
            queryNetsJoint(pbnets);
            getSequences();
        } else if (inference.equals("Marginal")) {
            System.out.println("*Information*\nNo node specification: returning marginal distribution of root node");
            queryNetsMarg(pbnets);
        } else {
            throw new RuntimeException("Error: Inference must be either 'Joint' or 'Marginal'");
//            System.exit(1);
        }
    }

    /**
     * For using the new gap model - incomplete and untested
     * @param file_tree
     * @param file_aln
     * @param inference
     * @param Q
     * @param F
     */
    public ASR(String file_tree, String file_aln, String inference, double[][] Q, double[] F) {
        loadData(file_tree, file_aln);
//        PhyloBNet[] pbnets = createNetworks(Q, F);
        PhyloBNet[] pbnets = createNetworks();
        if (inference.equals("Joint")) {
            queryNetsJoint(pbnets);
            getSequences();
        } else if (inference.equals("Marginal")) {
            System.out.println("*Information*\nNo node specification: returning marginal distribution of root node");
            queryNetsMarg(pbnets);
        } else {
            throw new RuntimeException("Error: Inference must be either 'Joint' or 'Marginal'");
//            System.exit(1);
        }

        Character gap = "-".charAt(0);
        List<EnumSeq.Gappy<Enumerable>> knownSeqs = null;
        try {
//            knownSeqs = EnumSeq.Gappy.loadFasta("bnkit/data/cyp2f_aln_full.fa", Enumerable.aacid, gap);
            knownSeqs = EnumSeq.Gappy.loadFasta("bnkit/data/simprot/cyp2f_r3.aln", Enumerable.aacid, gap);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
        double mismatches = 0.0;
        double comparisons = 0.0;

        for (EnumSeq.Gappy ks : knownSeqs) {
            String kName = ks.getName();
            EnumSeq ts = tree.find(kName).getSequence();
            String tName = ts.getName();
            if (kName.equals(tName)) {
                for (int i = 0; i < ts.length(); i++) {
                    comparisons++;
                    if (!(ks.get()[i] == ts.get()[i])){
                        mismatches ++;
                    }
                }
            }
        }
        double results = 1-(mismatches/comparisons);
        System.out.println(results);
    }

    public ASR(String file_tree, String file_aln, String inference, String nodeLabel) {
        loadData(file_tree, file_aln);
        PhyloBNet[] pbnets = createNetworks();
        if (inference.equals("Joint")) {
            System.out.println("*Information*\nUsing joint probability so node specification will be ignored");
            queryNetsJoint(pbnets);
            getSequences();
        } else if (inference.equals("Marginal")) {
            if (tree.find(nodeLabel) == null) {
                throw new RuntimeException("Error: Invalid node label" + nodeLabel + " - exiting");
//                System.exit(1);
            }
            queryNetsMarg(pbnets, tree.find(nodeLabel));
        } else {
            throw new RuntimeException("Error: Inference must be either 'Joint' or 'Marginal' - exiting");
//            System.exit(1);
        }
    }

    /**
  * Load the supplied tree and alignment files
  * Create phylogenetic tree, create node index, create node map
  * Store sequences and alignment
  * @param file_tree
  * @param file_aln 
  */
    public void loadData(String file_tree, String file_aln) {
        try {
            tree = tree.loadNewick(file_tree); //load tree - tree not in Newick
            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
            indexForNodes = new ArrayList<>(); // Newick string for subtree
            for (PhyloTree.Node n : nodes) {
                indexForNodes.add(n.getLabel().toString()); //link node to row number (index) in various matrices
            }

            BufferedReader aln_file = new BufferedReader(new FileReader(file_aln));

            String line = aln_file.readLine();
            if (line.startsWith("CLUSTAL")) {
                seqs = EnumSeq.Gappy.loadClustal(file_aln, Enumerable.aacid);
            } else if (line.startsWith(">")) {
                Character gap = "-".charAt(0);
                seqs = EnumSeq.Gappy.loadFasta(file_aln, Enumerable.aacid, gap);
            } else {
                throw new RuntimeException("Incorrect sequence or alignment format (requires FASTA or Clustal format .aln, .fa or .fasta)");
            }

            //Duplicate extant node names not allowed - will influence reconstruction outcomes
            PhyloTree.Node[] extNodes = getExtantNodes();
            Set<String> eNodes = new HashSet<>();
            for (PhyloTree.Node en : extNodes) {
                if (!eNodes.add(en.getLabel().toString())) {
                    throw new RuntimeException("Extant node names must be unique - " + en.getLabel().toString() + " is duplicated");
                }
            }
            //Duplicate sequence names not allowed - will influence reconstruction outcomes
            Set<String> seqNames = new HashSet<>();
            for (EnumSeq s : seqs) {
                if (!seqNames.add(s.getName())) {
                    throw new RuntimeException("Sequence names must be unique - " + s.getName() + " is duplicated");
                }
            }
            //Check whether the provided sequences match up to the provided tree - seqs.size() is used to monitor extant
            //sequences
            if (!eNodes.equals(seqNames)) {
                throw new RuntimeException("The sequence names in the provided alignment must all have a match" +
                        " in the provided tree");
            }

            aln = new EnumSeq.Alignment<>(seqs);
            pbn = PhyloBNet.create(tree, new JTT());
            
            } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
    /**
     * Using the stored alignment, create a PhyloBNet representing each column
     * in alignment. These networks will be queried to reconstruct ancestral
     * sequences. Populate pbnets array with networks.
     * pbnets array networks are instantiated with sequence from aln
     * pbgapnets array networks are instantiated with reduced gap/character alphabet
     */
    public PhyloBNet[] createNetworks(){
        //Each column/position in alignment/sequence gets a network
        PhyloBNet[] pbnets = new PhyloBNet[aln.getWidth()];
        
        String[] names = aln.getNames(); //seq_name_id - names of sequences
        List<String> labels = getLabels(names);

        //Create network for each column in alignment
        //Instantiate nodes
        for (int col = 0; col < aln.getWidth(); col ++) {
            Object[] gaps = aln.getGapColumn(col); // array with true for gap, false for symbol
            Object[] column = aln.getColumn(col);  // array for symbols, null for gaps
            tree.setContentByParsimony(names, gaps);
            PhyloBNet pbn;
            PhyloBNet pbnGap;
            if (use_sampled_rate) {
                pbn = PhyloBNet.create(tree, new JTT(), sampled_rate);
                pbnGap = PhyloBNet.createGap(tree, new Gap(), sampled_rate);
            } else {
                //creates BNet beginning with root then recursively
                //traversing subtrees
                pbn = PhyloBNet.create(tree, new JTT());
                pbnGap = PhyloBNet.createGap(tree, new Gap());
            }
            pbnets[col] = pbn;

            for (int i = 0; i < labels.size(); i ++) {
                BNode bnode = pbn.getBN().getNode(labels.get(i));// set node according to alignment
                BNode bnodeG = pbnGap.getBN().getNode(labels.get(i));// set node according to gap or character
                if (column[i] == null) {
                    bnode.setInstance(column[i]);// set node according to alignment
                    SubstNode sbnode = (SubstNode) bnode;
                    sbnode.setGap(true);
                    bnodeG.setInstance('G');// set node as gap
                } else {
                    bnode.setInstance(column[i]);// set node according to alignment
                    bnodeG.setInstance('C');// set node as character
                }
            }

            // Determine which internal nodes should be gaps based on gap evolutionary model
            Variable.Assignment[] states = queryGapNetJoint(pbnGap);
            for (Variable.Assignment a0 : states) {
                EnumVariable asr_var = (EnumVariable)a0.var;
                Object asr_val = a0.val;
                SubstNode bnode = (SubstNode)pbn.getBN().getNode(asr_var.getName());
                if (bnode == null)
                    System.out.println();
                if (asr_val.equals('G'))
                    bnode.setGap(true);
                else
                    bnode.setGap(false);
            }
        }
        return(pbnets);
    }

//    public PhyloBNet[] createNetworks(double[][] Q, double[] F){
//        //Each column/position in alignment/sequence gets a network
//        PhyloBNet[] pbnets = new PhyloBNet[aln.getWidth()];
//
//        String[] names = aln.getNames(); //seq_name_id - names of sequences
//        List<String> labels = getLabels(names);
//
//        //Create network for each column in alignment
//        //Instantiate nodes
//        for (int col = 0; col < aln.getWidth(); col ++) {
////            System.out.println(col);
//            Object[] gaps = aln.getGapColumn(col); // array with true for gap, false for symbol
//            Object[] column = aln.getColumn(col);  // array for symbols, null for gaps
//            tree.setContentByParsimony(names, gaps);
//            PhyloBNet pbn;
//            PhyloBNet pbnGap;
//            if (use_sampled_rate) {
//                pbn = PhyloBNet.create(tree, new JTT(), sampled_rate);
//                pbnGap = PhyloBNet.createGap(tree, new gap(Q, F), sampled_rate);
//            } else {
//                //creates BNet beginning with root then recursively
//                //traversing subtrees
//                pbn = PhyloBNet.create(tree, new JTT());
//                pbnGap = PhyloBNet.createGap(tree, new gap(Q, F));
//            }
//            pbnets[col] = pbn;
//
//            for (int i = 0; i < labels.size(); i ++) {
//                BNode bnode = pbn.getBN().getNode(labels.get(i));// set node according to alignment
//                BNode bnodeG = pbnGap.getBN().getNode(labels.get(i));// set node according to gap or character
//                if (column[i] == null) {
//                    bnode.setInstance(column[i]);// set node according to alignment
//                    SubstNode sbnode = (SubstNode) bnode;
//                    sbnode.setGap(true);
//                    bnodeG.setInstance('G');// set node as gap
//                } else {
//                    bnode.setInstance(column[i]);// set node according to alignment
//                    bnodeG.setInstance('C');// set node as character
//                }
//            }
//
//            // Determine which internal nodes should be gaps based on gap evolutionary model
//            Variable.Assignment[] states = queryGapNetJoint(pbnGap);
//            for (Variable.Assignment a0 : states) {
//                EnumVariable asr_var = (EnumVariable)a0.var;
//                Object asr_val = a0.val;
//                SubstNode bnode = (SubstNode)pbn.getBN().getNode(asr_var.getName());
//                if (bnode == null)
//                    System.out.println();
//                if (asr_val.equals('G'))
//                    bnode.setGap(true);
//                else
//                    bnode.setGap(false);
//            }
//        }
//        return(pbnets);
//    }
    
    /**
     * Query all networks in pbnets array using marginal probability
     * Populate margin_distribs for each column in alignment
     */
    public void queryNetsMarg(PhyloBNet[] pbnets) {
        this.margin_distribs = new EnumDistrib[aln.getWidth()];
        BNode root = null;
        for (int col = 0; col < aln.getWidth(); col ++) {
            PhyloBNet pbn = pbnets[col];
            BNet bn = pbn.getBN();
            root = pbn.getRoot();
            VarElim ve = new VarElim();
            ve.instantiate(bn);

            pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
            pbn.collapseSingles();
            
            margin_distribs[col] = getMarginalDistrib(ve, root.getVariable());
        }        
    }

    /**
     * Query all networks in pbnets array using marginal probability and a
     * specified node of interest.
     * Populate margin_distribs for each column in alignment
     * @param node to be queried
     */
    public void queryNetsMarg(PhyloBNet[] pbnets, PhyloTree.Node node) {
        this.margin_distribs = new EnumDistrib[aln.getWidth()];
        String nodeName = node.getLabel().toString();
//        BNode bnode = pbn.getBN().getNode(nodeName);
        for (int col = 0; col < aln.getWidth(); col ++) {
            PhyloBNet pbn = pbnets[col];
            BNet bn = pbn.getBN();
            BNode bnode = pbn.getBN().getNode(nodeName);
            VarElim ve = new VarElim();
            ve.instantiate(bn);

            pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
            pbn.collapseSingles();

            margin_distribs[col] = getMarginalDistrib(ve, bnode.getVariable());
        }
    }
    
    /**
     * Query all networks in pbnets array using MPE
     * Populate asr_matrix with reconstructed sequences
     * Populate rate matrix with calculated rate for each position in alignment
     *
     */
    public void queryNetsJoint(PhyloBNet[] pbnets) {
        // joint reconstruction for tree
        this.asr_matrix = new Object[indexForNodes.size()][aln.getWidth()];
        this.R = new double[aln.getWidth()]; //Rate matrix
        for (int col = 0; col < aln.getWidth(); col ++) {
            PhyloBNet pbn = pbnets[col];
            BNet bn = pbn.getBN();
            VarElim ve = new VarElim();
            ve.instantiate(bn);

            //introduce 'gap' tag to detect locations we think are truly gaps
            pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
            pbn.collapseSingles();
            
            for (Variable.Assignment a0 : getJointAssignment(ve)) {
                EnumVariable asr_var = (EnumVariable)a0.var;
                Object asr_val = a0.val;
                int index = indexForNodes.indexOf(asr_var.getName());
                if (index >= 0) 
                    //index = current node
                    //col = position in alignment
                    asr_matrix[index][col] = asr_val;
                BNode node = bn.getNode(asr_var);
                node.setInstance(asr_val);
            }
            R[col] = pbn.getRate(); //calculates rate based on evidence provided
            //All nodes instantiated with MPE - rate across entire tree?
        }        
    }
    
    /**
     * Query a specific column in the alignment using marginal probability
     * Update marginal distribution for column
     * @param col
     * @return EnumDistrib describing the marginal distribution for the column given the root node
     */
    public EnumDistrib queryNetMarg(PhyloBNet[] pbnets, int col) {
        BNode root = null;
        PhyloBNet pbn = pbnets[col];
        BNet bn = pbn.getBN();
        root = pbn.getRoot();
        VarElim ve = new VarElim();
        ve.instantiate(bn);

        pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
        pbn.collapseSingles();

        return(getMarginalDistrib(ve, root.getVariable()));
    }
    
    /**
     * Query a specific column in the alignment using MPE
     * Update ASR for column
     * Update rate for column
     *
     */
    public Variable.Assignment[] queryGapNetJoint(PhyloBNet pbn) {
        BNet bn = pbn.getBN();
        VarElim ve = new VarElim();
        ve.instantiate(bn);
            
        return(getJointAssignment(ve));
    }
    
    private EnumDistrib getMarginalDistrib(VarElim ve, Variable queryNode) {
        EnumDistrib d_marg = null;
        try {
            Query q_marg = ve.makeQuery(queryNode);
            CGTable r_marg = (CGTable)ve.infer(q_marg);
            d_marg = (EnumDistrib)r_marg.query(queryNode);
        } catch (NullPointerException npe) { //When node of interest has been removed from network of interest
            if (npe.toString().contains("Invalid query")) {
                double[] empty = new double[Enumerable.aacid.size()];
                for (int d = 0; d < Enumerable.aacid.size(); d++) {
                    empty[d] = 0.0;
                }
                d_marg = new EnumDistrib(Enumerable.aacid, empty);
            } else {
                npe.printStackTrace();
            }

        }
        return d_marg;
    }
    
    private Variable.Assignment[] getJointAssignment(VarElim ve) {
        Query q_joint = ve.makeMPE();
        CGTable r_joint = (CGTable)ve.infer(q_joint);
        Variable.Assignment[] a = r_joint.getMPE();
        return a;
    }
    
    private Variable.Assignment[] getJointAssignment(VarElim ve, Variable[] queryNode) {
        Query q_joint = ve.makeMPE(queryNode);
        CGTable r_joint = (CGTable)ve.infer(q_joint);
        Variable.Assignment[] a = r_joint.getMPE();
        return a;
    }
    
    /**
     * Extract sequences for all nodes from the asr matrix
     */
    public void getSequences(){
        //Set extant nodes
        for (int j = 0; j < seqs.size(); j++) {
            EnumSeq.Gappy<Enumerable> seq = seqs.get(j);
            tree.find(seq.getName()).setSequence(seq);
        }

        //Retrieve and store reconstructions for each node
        List<EnumSeq.Gappy<Enumerable>> asrs = new ArrayList<>();
        List<PhyloTree.Node> extant = Arrays.asList(getExtantNodes());
        //asr_matrix stores joint reconstruction - MPE assignment of each internal node in each network
        for (int row = 0; row < asr_matrix.length; row ++) { 
            Object[] asr_obj = asr_matrix[row]; //retrieve MP sequence for node/row
            if (extant.contains(tree.find(indexForNodes.get(row)))) { //extant node so try to manually retrieve sequence
                EnumSeq.Gappy<Enumerable> extSeq = (EnumSeq.Gappy<Enumerable>) tree.find(indexForNodes.get(row)).getSequence();
                asrs.add(extSeq);
            } else { //reconstructed node so get sequence from matrix
                EnumSeq.Gappy<Enumerable> myasr = new EnumSeq.Gappy<>(Enumerable.aacid_alt);
                myasr.set(asr_obj);
                myasr.setName(indexForNodes.get(row));
                asrs.add(myasr);
            }
        }

        String rootname = tree.getRoot().toString();
        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
        //Create a new alignment from the reconstructed sequences
        //Print reconstruction for each internal node
        EnumSeq.Alignment aln_asr = new EnumSeq.Alignment(asrs);
        for (int i = 0; i < aln_asr.getHeight(); i ++) {
            EnumSeq.Gappy<Enumerable> asr_seq = aln_asr.getEnumSeq(i);
            String nodename = asr_seq.getName();
            if (rootname.equals(nodename))
                nodes[i].setSequence(asr_seq);
                
            //Not root and not leaf
            if (nodes[i].getChildren().toArray().length > 0){
                nodes[i].setSequence(asr_seq); //set sequence for internal nodes
            }
        }
        aln = aln_asr;
    }

    /**
     * Estimates parameters of gamma distribution then creates a gamma distribution
     */
    public GammaDistrib calcGammaDistrib(){
        //estimates parameters of gamma distribution
        double alpha = GammaDistrib.getAlpha(R);
        double beta = 1 / alpha;
        //Creates a gamma distribution
        return new GammaDistrib(alpha, 1/beta);
    }
    
    public boolean save(String id, boolean infSequence) {
        if (infSequence) {
//            saveJSON(id + "_JSON_output.txt");
            saveALN(id + "_aln_full.fa");
            saveTree(id + "_new_tree.nwk");
            saveRate(id + "_rates.txt");
        } else {
            saveTree(id + "_new_tree.nwk");
            saveDistrib(id + "_distribution.txt");
        }
        return true;
    }

    public boolean saveJSON(String filename) {
        try{
            Writer writer = new PrintWriter(filename, "UTF-8");

            //Create the top level object - reconstruction
            JSONObject recon = new JSONObject();

            //Create and populate the root node object
            JSONObject root = new JSONObject();
            root.put("Sequence",tree.getRoot().getSequence().toString());
            root.put("SeqName", tree.getRoot().getLabel());
            root.put("NewickRep", tree.getRoot().toString());
            JSONObject rates = new JSONObject();
            for (int j = 0; j < R.length; j++){
                rates.put(Integer.toString(j), R[j]);
            }
            root.put("Rates", rates);

            //Add the root information to the reconstruction object
            recon.put("Root", root);

            //Create and populate the internal nodes object
            PhyloTree.Node[] nodes = getInternalNodes();
            for (PhyloTree.Node n: nodes) {
                JSONObject node = new JSONObject();
                node.put("SeqName",n.getLabel());
                node.put("Sequence", n.getSequence().toString());
//                node.put("NewickRep",nodes[i].toString());
                //Using append rather than put automatically creates an array
                recon.append("ReconstructedNodes", node);
            }

            //Create and populate the extant nodes object
            PhyloTree.Node[] exNodes = getExtantNodes();
            for (PhyloTree.Node n: exNodes) {
                JSONObject node = new JSONObject();
                node.put("SeqName",n.getLabel());
                node.put("Sequence", n.getSequence().toString());
//                node.put("NewickRep",nodes[i].toString());
                recon.append("ExtantNodes", node);
            }

            JSONObject fin = new JSONObject();
            fin.put("Reconstruction", recon);
            fin.write(writer);
            writer.close();

        } catch (UnsupportedEncodingException uee) {
            uee.printStackTrace();
            return false;
        } catch (FileNotFoundException fnf) {
            fnf.printStackTrace();
            return false;
        } catch (IOException ioe) {
            ioe.printStackTrace();
            return false;
        }
        return true;
    }

    public void saveALN(String filename) {
        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst();
        EnumSeq.Gappy<Enumerable>[] allSeqs = new EnumSeq.Gappy[nodes.length];
        for (int n = 0; n < nodes.length; n++) {
            //Update sequence name at this point
            EnumSeq.Gappy<Enumerable> seq = (EnumSeq.Gappy) nodes[n].getSequence();

            String seqName = nodes[n].toString();
            String seqLab = seq.getName();
            if (seqLab != null) {
                seq.setName(seqLab + " " + seqName + ";"); //Newick strings require a ';' to indicate completion
            }
            allSeqs[n] = seq;
        }
        try {
            FastaWriter fw = new FastaWriter(filename);
            fw.save(allSeqs);
            fw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public void saveTree(String filename) {
        try {
            Writer writer = new PrintWriter(filename, "UTF-8");
            String newick = tree.getRoot().toString();
            writer.write(newick);
            writer.write(";\n");
            writer.close();
        } catch (UnsupportedEncodingException uee) {
            uee.printStackTrace();
        } catch (FileNotFoundException fnf) {
            fnf.printStackTrace();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public void saveDistrib(String filename){
        try {
            Writer writer = new PrintWriter(filename, "UTF-8");
            Object[] aacid = Enumerable.aacid_ext.getValues();
            Object[][] margMatrix = new Object[aacid.length][margin_distribs.length + 1];
            writer.write("columns\t");
            for (int i = 1; i < margin_distribs.length; i++) { //write header
                if (i == margin_distribs.length - 1)
                    writer.write(i + "\n");
                else
                    writer.write(i + "\t");
            }
            for (int j = 0; j < aacid.length; j++) { //fill in row names
                margMatrix[j][0] = aacid[j];
            }
            for (int k = 1; k < margin_distribs.length + 1; k++) {
                EnumDistrib distr = margin_distribs[k-1];
                for (int a = 0; a < aacid.length; a++) {
                    if (aacid[a].equals('-')) {
                        Object na = "NA";
                        margMatrix[a][k] = na;
                    } else {
                        double val = distr.get(aacid[a]);
                        margMatrix[a][k] = val;
                    }
                }
            }
            for (int j = 0; j < aacid.length; j++) {
                for (int k = 0; k < margMatrix[j].length; k++) {
                    if (k == margMatrix[j].length - 1) {
                        writer.write(margMatrix[j][k] + "\n");
                    } else {
                        writer.write(margMatrix[j][k] + "\t");
                    }
                }
            }
            writer.close();
        } catch (UnsupportedEncodingException uee) {
            uee.printStackTrace();
        } catch (FileNotFoundException fnf) {
            fnf.printStackTrace();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    public void saveRate(String filename) {
        try {
            Writer writer = new PrintWriter(filename, "UTF-8");
            String rates = R.toString();
            writer.write(rates);
            writer.close();
        } catch (UnsupportedEncodingException uee) {
            uee.printStackTrace();
        } catch (FileNotFoundException fnf) {
            fnf.printStackTrace();
        } catch (IOException ioe) {
            ioe.printStackTrace();
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
        PhyloTree.Node[] intNodesA = new PhyloTree.Node[nodes.length - seqs.size() - 1];
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

        PhyloTree.Node[] intNodesA = new PhyloTree.Node[nodes.length - seqs.size() - 1];
        return extNodes.toArray(intNodesA);
    }

    public PhyloTree getTree() {
        return tree;
    }

    public PhyloBNet getPbn() {return pbn;}
    
    public List<EnumSeq.Gappy<Enumerable>> getSeqs(){
        return seqs;
    }

    public EnumSeq.Alignment<Enumerable> getAln() {
        return aln;
    }
    
//    //FIXME: Can you update the pbnets array without modifying the aln. Where
//    //do the two interact?
//    /**
//     * Get the array of pbnets representing the alignment
//     * @return array of pbnets
//     */
//    public PhyloBNet[] getPbnets() {
//        return pbnets;
//    }
//    /**
//     * Unlikely to be used
//     * Set the array of pbnets representing the alignment
//     * @param pbnets
//     */
//    public void setPbnets(PhyloBNet[] pbnets) {
//        this.pbnets = pbnets;
//    }
//    /**
//     * Update a particular pbnet in the pbnet array
//     * Use: when modifying single column in alignment, update array of networks
//     * @param index
//     * @param pbnet
//     */
//    public void setPbnet(int index, PhyloBNet pbnet) {
//        pbnets[index] = pbnet;
//    }
//    /**
//     * Remove a specific pbnet from the array
//     * Use: when removing a column from an alignment, update array of networks
//     * @param index
//     */
//    public void removePbnet(int index) {
//        PhyloBNet[] pbnetsNew = new PhyloBNet[pbnets.length - 1];
//        List<PhyloBNet> list = new ArrayList<PhyloBNet>(Arrays.asList(pbnets));
//        list.remove(index);
//        pbnets = list.toArray(pbnetsNew);
//    }
        
    public double[] getRates(){
        return R;
    }
    
    public EnumDistrib[] getMarginDistribs(){
        return margin_distribs;
    }
    
    
    /**
     * Checks names of sequences from alignment for amendments and records 
     * appropriate substring label
     * @param names - array of names from alignment file
     * @return labels - list of labels with string modifications
     */
    private static List<String> getLabels(String[] names) {
        //seq_name_id - names of sequences
        List<String> labels = new ArrayList<>();
           
        for (int i = 0; i < names.length; i ++) {
            //check if any names amended and modify if they are
            int index = names[i].indexOf("/"); // in this aln file, names have been amended
            if (index > 0)
                labels.add(names[i].substring(0, index));
            else
                labels.add(names[i]);
        }
        return labels;
    }

}
