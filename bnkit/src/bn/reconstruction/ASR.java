/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bn.reconstruction;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import bn.alg.*;
import bn.ctmc.PhyloBNet;
import bn.ctmc.matrix.JTT;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.PhyloTree;
import dat.Variable;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.*;

import json.JSONObject;


/**
 * @author Mikael Boden
 * @author Alex Essebier
 */
public class ASR {
    
    private PhyloTree tree;
    private PhyloBNet pbn; //only to be used for navigating branches
    private List<EnumSeq.Gappy<Enumerable>> seqs;
    private EnumSeq.Alignment<Enumerable> aln;
    private PhyloBNet[] pbnets;
    private double[] R; //Rates at positions in alignment
    private EnumDistrib[] margin_distribs; //Marginal distributions for nodes
    private boolean use_sampled_rate = false;
    private String asr_root; //Reconstructed sequence 
    private GammaDistrib gd; //Calculated gamma distribution
    private BNet[] models; //model that will be trained for each column in alignment
    
    
    private List<String> indexForNodes;
    private Map<String, String> mapForNodes;
    //Joint reconstruction of tree
    private Object[][] asr_matrix; //storage of all reconstructed sequences
    private double sampled_rate = // sampled rate, copy from a previous 1.0-rate run
	0.15599004226404184;
    
    public ASR(String file_tree, String file_aln) {
        loadData(file_tree, file_aln);
        createNetworks();
        queryNetsMarg();
        queryNetsJoint();
        getSequences();
        calcGammaDistrib();
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
            tree = PhyloTree.loadNewick(file_tree); //load tree - tree not in Newick
            PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
            indexForNodes = new ArrayList<>(); // Newick string for subtree
            mapForNodes = new HashMap<>(); // Shortname --> Newick string for subtree
            for (PhyloTree.Node n : nodes) {
                //n string format internal 'N0_'
                //n string format extant (no children) 'seq_name_id'
                //n.toString() recursive Newick representation node and children
                //creates subtrees?
                indexForNodes.add(replacePunct(n.toString())); 
                mapForNodes.put(n.getLabel().toString(), replacePunct(n.toString()));
            }
            
            seqs = EnumSeq.Gappy.loadClustal(file_aln, Enumerable.aacid);
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
     */
    public void createNetworks(){
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
            if (use_sampled_rate)
                pbn = PhyloBNet.create(tree, new JTT(), sampled_rate);
            else
                //creates BNet beginning with root then recursively
                //traversing subtrees
                pbn = PhyloBNet.create(tree, new JTT());
            pbnets[col] = pbn;

            // set variables according to alignment
            for (int i = 0; i < labels.size(); i ++) {
                String shortname = labels.get(i);
                String longname = mapForNodes.get(shortname);
                if (longname != null) {
                    BNode bnode = pbn.getBN().getNode(longname);
                    bnode.setInstance(column[i]);
                }
            }
        }
        this.pbnets = pbnets;
    }
    
    /**
     * Query all networks in pbnets array using marginal probability
     * Populate margin_distribs for each column in alignment
     */
    public void queryNetsMarg() {
        this.margin_distribs = new EnumDistrib[aln.getWidth()];
        BNode root = null;
        for (int col = 0; col < aln.getWidth(); col ++) {
            PhyloBNet pbn = pbnets[col];
            BNet bn = pbn.getBN();
            //FIXME
            root = pbn.getRoot(); //Possibly in wrong location??
            //Root can change in purge and collapse steps?
            VarElim ve = new VarElim();
            ve.instantiate(bn);

            int purged_leaves = pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
            int collapsed_nodes = pbn.collapseSingles();

//            Variable testNode = bn.getAlphabetical().get(1).getVariable();
//            EnumDistrib test = getMarginalDistrib(ve, testNode);
            
            margin_distribs[col] = getMarginalDistrib(ve, root.getVariable());
        }        
    }
    
    /**
     * Query all networks in pbnets array using MPE
     * Populate asr_matrix with reconstructed sequences
     * Populate rate matrix with calculated rate for each position in alignment
     *
     */
    public void queryNetsJoint() {
        // joint reconstruction for tree
        this.asr_matrix = new Object[indexForNodes.size()][aln.getWidth()];
        this.R = new double[aln.getWidth()]; //Rate matrix
        for (int col = 0; col < aln.getWidth(); col ++) {
            PhyloBNet pbn = pbnets[col];
            BNet bn = pbn.getBN();
            VarElim ve = new VarElim();
            ve.instantiate(bn);

            int purged_leaves = pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
            int collapsed_nodes = pbn.collapseSingles();
            
            for (Variable.Assignment a0 : getJointAssignment(ve)) {
                EnumVariable asr_var = (EnumVariable)a0.var;
                Object asr_val = a0.val;
                int index = indexForNodes.indexOf(replacePunct(asr_var.getName()));
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
     */
    public void queryNetMarg(int col) {
        BNode root = null;
        PhyloBNet pbn = pbnets[col];
        BNet bn = pbn.getBN();
        //FIXME
        root = pbn.getRoot(); //Possibly in wrong location??
        //Root can change in purge and collapse steps?
        VarElim ve = new VarElim();
        ve.instantiate(bn);

        int purged_leaves = pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
        int collapsed_nodes = pbn.collapseSingles();

        margin_distribs[col] = getMarginalDistrib(ve, root.getVariable());
    }
    
    /**
     * Query a specific column in the alignment using MPE
     * Update ASR for column
     * Update rate for column
     * @param col - column in alignment
     */
    public void queryNetJoint(int col) {
        PhyloBNet pbn = pbnets[col];
        BNet bn = pbn.getBN();
        VarElim ve = new VarElim();
        ve.instantiate(bn);

        int purged_leaves = pbn.purgeGaps(); //Remove leaves with gap (i.e. uninstantiated)
        int collapsed_nodes = pbn.collapseSingles();
            
        for (Variable.Assignment a0 : getJointAssignment(ve)) {
            EnumVariable asr_var = (EnumVariable)a0.var;
            Object asr_val = a0.val;
            int index = indexForNodes.indexOf(replacePunct(asr_var.getName()));
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
    
    private EnumDistrib getMarginalDistrib(VarElim ve, Variable queryNode) {
        Query q_marg = ve.makeQuery(queryNode);
        CGTable r_marg = (CGTable)ve.infer(q_marg);
        EnumDistrib d_marg = (EnumDistrib)r_marg.query(queryNode);
        return d_marg;
    }
    
    private Variable.Assignment[] getJointAssignment(VarElim ve) {
        Query q_joint = ve.makeMPE();
        CGTable r_joint = (CGTable)ve.infer(q_joint);
        Variable.Assignment[] a = r_joint.getMPE();
        return a;
    }
    
    private Variable.Assignment[] getJointAssignment(VarElim ve, Variable queryNode) {
        Query q_joint = ve.makeMPE(queryNode);
        CGTable r_joint = (CGTable)ve.infer(q_joint);
        Variable.Assignment[] a = r_joint.getMPE();
        return a;
    }
    
    /**
     * Extract sequences for internal nodes from the asr matrix
     * Identify sequence of root node and print internal node results
     */
    public void getSequences(){
        //Retrieve and store reconstructions for each latent node
        List<EnumSeq.Gappy<Enumerable>> asrs = new ArrayList<>();
        //asr_matrix stores joint reconstruction - MPE assignment of each node in each network
        for (int row = 0; row < asr_matrix.length; row ++) { 
            Object[] asr_obj = asr_matrix[row]; //retrieve MP sequence for node/row
            EnumSeq.Gappy<Enumerable> myasr = new EnumSeq.Gappy<>(Enumerable.aacid_alt);
            myasr.set(asr_obj);
            myasr.setName(indexForNodes.get(row));
            asrs.add(myasr);
        }
        
        String rootname = replacePunct(tree.getRoot().toString());
        PhyloTree.Node[] nodes = tree.toNodesBreadthFirst(); //tree to nodes - recursive
        //Create a new alignment from the reconstructed sequences
        //Print reconstruction for each internal node
        EnumSeq.Alignment aln_asr = new EnumSeq.Alignment(asrs);
        for (int i = 0; i < aln_asr.getHeight(); i ++) {
            EnumSeq.Gappy<Enumerable> asr_seq = aln_asr.getEnumSeq(i);
            String nodename = asr_seq.getName();
            if (rootname.equals(nodename))
                this.asr_root = asr_seq.toString();
                nodes[i].setSequence(asr_seq);
                //System.out.println(asr_seq.getName() + "\t" + asr_seq.toString());
                
            //Not root and not leaf
            if (nodes[i].getChildren().toArray().length > 0){
                nodes[i].setSequence(asr_seq);
            } else {
                nodes[i].setSequence(asr_seq); //set sequence for internal nodes
            }
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

    /**
     * Estimates parameters of gamma distribution then creates a gamma distribution
     */
    public void calcGammaDistrib(){
        //estimates parameters of gamma distribution
        double alpha = GammaDistrib.getAlpha(R);
        double beta = 1 / alpha;
//        System.out.println("Gamma alpha = " + alpha + " beta = " + beta);
        //Creates a gamma distribution
        this.gd = new GammaDistrib(alpha, 1/beta);
    }
    
    public boolean save(String filename) {
        
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
//            JSONObject margDistribs = new JSONObject();
//            Object[] aacid = Enumerable.aacid.getValues();
//            for (int k = 0; k < margin_distribs.length; k++) {
//                EnumDistrib distr = margin_distribs[k];
//                JSONObject position = new JSONObject();
//                for (int a = 0; a < aacid.length; a++) {
//                    double val = distr.get(aacid[a]);
//                    position.put(aacid[a].toString(), val);
//                }
//                margDistribs.put(Integer.toString(k), position);                
//            }
//            root.put("MarginalDistribs", margDistribs);
            
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
            
//            //Create and populate the gamma information
//            JSONObject gammaAB = new JSONObject();
//            gammaAB.put("alpha", gd.getAlpha());
//            gammaAB.put("beta", gd.getBeta());
//            recon.put("gammaDistrib", gammaAB);
            
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
        for (int j = 0; j < seqs.size(); j++) {
            EnumSeq.Gappy<Enumerable> seq = seqs.get(j);
            tree.find(seq.getName()).setSequence(seq);
        }
        PhyloTree.Node[] intNodesA = new PhyloTree.Node[nodes.length - seqs.size() - 1];
        return extNodes.toArray(intNodesA);
    }
    
    public String getAsrSeq(){
        return asr_root;
    }
    
    public PhyloTree getTree() {
        return tree;
    }
    
    public List<EnumSeq.Gappy<Enumerable>> getSeqs(){
        return seqs;
    }

    public EnumSeq.Alignment<Enumerable> getAln() {
        return aln;
    }
    
    //FIXME: Can you update the pbnets array without modifying the aln. Where
    //do the two interact?
    /**
     * Get the array of pbnets representing the alignment
     * @return array of pbnets
     */
    public PhyloBNet[] getPbnets() {
        return pbnets;
    }
    /**
     * Unlikely to be used
     * Set the array of pbnets representing the alignment
     * @param pbnets 
     */
    public void setPbnets(PhyloBNet[] pbnets) {
        this.pbnets = pbnets;
    }
    /**
     * Update a particular pbnet in the pbnet array
     * Use: when modifying single column in alignment, update array of networks
     * @param index
     * @param pbnet 
     */
    public void setPbnet(int index, PhyloBNet pbnet) {
        pbnets[index] = pbnet;
    }
    /**
     * Remove a specific pbnet from the array
     * Use: when removing a column from an alignment, update array of networks
     * @param index 
     */
    public void removePbnet(int index) {
        PhyloBNet[] pbnetsNew = new PhyloBNet[pbnets.length - 1];
        List<PhyloBNet> list = new ArrayList<PhyloBNet>(Arrays.asList(pbnets));
        list.remove(index);
        pbnets = list.toArray(pbnetsNew);
    }
        
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
    
    private static String replacePunct(String str) {
        return str.replace('.', '_');
    }    
    
}
