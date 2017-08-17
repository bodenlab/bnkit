package reconstruction;


import alignment.MSA;
import api.PartialOrderGraph;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.PhyloBNet;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JTT;
import bn.prob.EnumDistrib;
import dat.*;
import dat.file.AlnWriter;
import dat.file.FastaWriter;

import java.io.*;
import java.util.*;

/**
 * Reconstruct ancestral sequences using partial order graphs to represent indels. Each node of the resulting phylogenetic tree
 * contains a partial order alignment graph indicative of the ancestral sequence alignment at that position.
 *
 * @author Mikael Boden
 * @author Alex Essebier
 * @author Marnie Lamprecht
 *
 */
public class ASRPOG {

	private PhyloTree phyloTree; 										// Phylogenetic tree structure
	private List<EnumSeq.Gappy<Enumerable>> extantSequences;			// List of sequences (label,bases)
	private List<String> ancestralSeqLabels;							// Ancestral sequences labels (internal nodes of phylogenetic tree structure)
	private POGraph pogAlignment;										// partial order alignment graph structure template
	private EnumDistrib[] marginalDistributions; 						// Marginal distributions for nodes if doing a marginal reconstruction
	private String marginalNode = null;									// Label of node to perform marginal reconstruction (if applicable)
	private Map<String, List<Inference>> ancestralInferences;			// stores updates to the POGStructure for the ancestral node <node label, changes>
	private Double[] rates = null; 										// Rates at positions in alignment

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile		filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 */
	public ASRPOG(String alignmentFile, String treeFile, boolean jointInference, boolean performMSA) throws IOException {
		performASR("", treeFile, alignmentFile, jointInference, performMSA);
	}

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile	filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param marginalNode	node label for maginal inference
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 */
	public ASRPOG(String alignmentFile, String treeFile, String marginalNode, boolean performMSA) throws IOException {
		this.marginalNode = marginalNode;
		performASR("", treeFile, alignmentFile, false, performMSA);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean performMSA) throws IOException {
		performASR(pog, treeFile, sequenceFile, jointInference, performMSA);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 * This one
	 * @param msa				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 */
	public ASRPOG(POGraph msa, String treeFile, String sequenceFile, boolean jointInference, boolean performMSA) throws IOException {
		performASR(msa, treeFile, sequenceFile, jointInference, performMSA);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog			POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile	filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param marginalNode	node label for maginal inference
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, String marginalNode, boolean performMSA) throws IOException {
		this.marginalNode = marginalNode;
		performASR(pog, treeFile, sequenceFile, false, performMSA);
	}

	/**
	 * Constructs partial order alignment graphs for each of the internal nodes of the phylogenetic tree.
	 * 
	 * @param filepath	filepath to store the internal POAG alignment information to
	 */
	public void saveGraph(String filepath) {
		for (String phyloNodeLabel : ancestralSeqLabels)
			saveGraph(filepath, phyloNodeLabel);
	}
	
	/**
	 * Constructs partial order alignment graphs for the given internal node of the phylogenetic tree.
	 * 
	 * @param filepath	filepath to store the POAG alignment information to
	 * @param nodeLabel	label of the internal node to generate file for (internal labels are specified in the phylogenetic tree .nwk file)
	 */
	public void saveGraph(String filepath, String nodeLabel) {
		// save partial order graph alignment to text file
		POGraph ancestor = getAncestor(nodeLabel);
		ancestor.saveToDot(filepath + nodeLabel);
	}

	/**
	 * Get the partial order alignment graph for the given internal node of the phylogenetic tree.
	 *
	 * @param nodeLabel	Ancestral node
	 * @return			Partial order graph of the given ancestral node
	 */
	public PartialOrderGraph getGraph(String nodeLabel) {
		return new PartialOrderGraph(getAncestor(nodeLabel));
	}

	/**
	 * Get the multiple sequence alignment partial order graph.
	 *
	 * @return	Partial order graph representing sequence alignment
	 */
	public PartialOrderGraph getMSAGraph() {
		return new PartialOrderGraph(this.pogAlignment);
	}
	
	/**
	 * Save multiple sequence alignment partial order alignment graph as a dot file in the given output filepath.
	 *
	 * @param filepath	Output filepath
	 */
	public void saveMSAGraph(String filepath) {
		pogAlignment.saveToDot(filepath + "MSA");
	}

	/**
	 * Save the reconstructed sequences that have the most support through the partial order graph in FASTA format.
	 * Saved in output path as "reconstructed_sequences.fasta"
	 *
	 * @param filepath	Output filepath
	 */
	public void saveSupportedAncestors(String filepath) throws IOException {
		Map<String, String> ancestralSeqs = new HashMap<>();
		if (marginalNode == null)
			for (String phyloNodeLabel : ancestralSeqLabels)
				ancestralSeqs.put(phyloNodeLabel, phyloTree.find(phyloNodeLabel).getSequence().toString());
		else
			ancestralSeqs.put(marginalNode, phyloTree.find(marginalNode).getSequence().toString());

		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath + "_recon.fa", false));
		for (String node : ancestralSeqs.keySet()) {
			bw.write(">" + node);
			bw.newLine();
			bw.write(ancestralSeqs.get(node));
			bw.newLine();
			bw.newLine();
		}
		bw.close();
	}

	/**
	 * Save the reconstructed sequences that have the most support through the partial order graph in FASTA format.
	 * Saved in output path as "reconstructed_sequences.fasta"
	 *
	 * @param filepath	Output filepath
	 * @param label		label of ancestor to save (tree node label)
	 */
	public void saveSupportedAncestor(String filepath, String label) throws IOException {
		Map<String, String> ancestralSeqs = new HashMap<>();
		if (marginalNode == null)
			ancestralSeqs.put(label, phyloTree.find(label).getSequence().toString());
		else
			ancestralSeqs.put(marginalNode, phyloTree.find(marginalNode).getSequence().toString());

		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath + "_recon.fa", false));
		for (String node : ancestralSeqs.keySet()) {
			bw.write(">" + node);
			bw.newLine();
			bw.write(ancestralSeqs.get(node));
			bw.newLine();
			bw.newLine();
		}
		bw.close();
	}

	/**
	 * Save ASR information; ALN, tree, rates, marginal distribution (if applicable).
	 *
	 * @param filepath		filename to save ASR
	 * @param infSequence   sequences inferred (true) or marginal (false)
	 * @param format		format to save alignment to, "fasta" or "clustal"
	 */
	public void save(String filepath, boolean infSequence, String format) throws IOException {
		if (infSequence) {
			saveALN(filepath + "_aln_full", format);
			saveTree(filepath + "_new_tree.txt");
			saveRate(filepath + "_rates.txt");
		} else {
			saveTree(filepath + "_new_tree.txt");
			saveDistrib(filepath + "_distribution.txt");
		}
	}

	/**
	 * Save ALN
	 *
	 * @param filepath	filepath to save ALN
	 * @param format	format to save ALN, "clustal" or "fasta", default: fasta
	 */
	public void saveALN(String filepath, String format) throws IOException {
		PhyloTree.Node[] nodes = phyloTree.toNodesBreadthFirst();
		EnumSeq.Gappy<Enumerable>[] allSeqs = new EnumSeq.Gappy[nodes.length];
		for (int n = 0; n < nodes.length; n++) {
			EnumSeq.Gappy<Enumerable> seq = (EnumSeq.Gappy) nodes[n].getSequence();
			String seqName = nodes[n].toString();
			String seqLab = seq.getName();
			if (seqLab != null)
				seq.setName(seqLab + " " + seqName + ";"); //Newick strings require a ';' to indicate completion
			allSeqs[n] = seq;
		}
		if (format.equalsIgnoreCase("clustal")) {
			AlnWriter aw = new AlnWriter(filepath + ".aln");
			aw.save(allSeqs);
			aw.close();
		} else {
			FastaWriter fw = new FastaWriter(filepath + ".fa");
			fw.save(allSeqs);
			fw.close();
		}
	}

	/**
	 * Save rates
	 *
	 * @param filename filepath to save rates
	 */
	public void saveRate(String filename) throws IOException {
		if (rates == null)
			rates = new Double[]{};
		Writer writer = new PrintWriter(filename, "UTF-8");
		for (Integer nodeId : pogAlignment.getNodeIDs()) {
			pogAlignment.setCurrent(nodeId);
			writer.write(pogAlignment.getCurrentId() + ":");
			for (Character base : pogAlignment.getCurrentBases())
				writer.write(base + ",");
			writer.write(" " + (rates[nodeId] == null ? "NA" : Double.toString(rates[nodeId])) + "\n");
		}
		writer.close();
	}

	/**
	 * Save marginal distribution
	 *
	 * @param filename	filepath to save distribution
	 */
	public void saveDistrib(String filename) throws IOException {
		if (marginalNode == null) {
			System.err.println("Marginal reconstruction has not been performed. Cannot save distribution.");
			return;
		}
		POGraph ancestor = getAncestor(marginalNode);
		Writer writer = new PrintWriter(filename + "_dist.tsv", "UTF-8");
		Object[] aacid = Enumerable.aacid.getValues();
		writer.write("ID\t");
		// Header is amino acid characters
		for (int i = 0; i < aacid.length; i++)  //write header
			if (i == aacid.length - 1)
				writer.write(aacid[i] + "\n");
			else
				writer.write(aacid[i] + "\t");
		// Each row is a graph node ID and character distribution
		for (int k = 0; k < ancestor.getNumNodes(); k++) {
			writer.write(Integer.toString(ancestor.getNodeIDs().get(k)) + "\t");
			ancestor.setCurrent(ancestor.getNodeIDs().get(k));
			Map<Character, Double> distribution = ancestor.getCharacterDistribution();
			for (int a = 0; a < aacid.length; a++) {
				if (!distribution.containsKey(aacid[a]))
					writer.write("NA");
				else
					writer.write(Double.toString(distribution.get(aacid[a])));
				if (a == aacid.length - 1)
					writer.write("\n");
				else
					writer.write("\t");
			}
		}
		writer.close();
	}

	/**
	 * Save distribution of characters for MSA
	 *
	 * @param filename	filepath to save distribution
	 */
	public void saveMSADistrib(String filename) throws IOException {
		Writer writer = new PrintWriter(filename + "_dist.tsv", "UTF-8");
		Object[] aacid = Enumerable.aacid.getValues();
		writer.write("ID\t");
		// Header is amino acid characters
		for (int i = 0; i < aacid.length; i++)  //write header
			if (i == aacid.length - 1)
				writer.write(aacid[i] + "\n");
			else
				writer.write(aacid[i] + "\t");
		// Each row is a graph node ID and character distribution
		for (int k = 0; k < pogAlignment.getNumNodes(); k++) {
			writer.write(Integer.toString(pogAlignment.getNodeIDs().get(k)) + "\t");
			pogAlignment.setCurrent(pogAlignment.getNodeIDs().get(k));
			Map<Character, Double> distribution = pogAlignment.getCharacterDistribution();
			for (int a = 0; a < aacid.length; a++) {
				if (!distribution.containsKey(aacid[a]))
					writer.write("NA");
				else
					writer.write(Double.toString(distribution.get(aacid[a])));
				if (a == aacid.length - 1)
					writer.write("\n");
				else
					writer.write("\t");
			}
		}
		writer.close();
	}

	/**
	 * Save the phylogenetic tree with the ancestral nodes entered.
	 *
	 * @param filepath	filepath to save phylogenetic tree (.nwk)
	 */
	public void saveTree(String filepath) throws IOException {
		Writer writer = new PrintWriter(filepath, "UTF-8");
		String newick = phyloTree.getRoot().toString();
		writer.write(newick);
		writer.write(";\n");
		writer.close();
	}

	/* ****************************************************************************************************************************************************
	 * 																PRIVATE METHODS
	 * ****************************************************************************************************************************************************/

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param pog				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 */
	private void performASR(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean performMSA) throws RuntimeException, IOException {
		loadData(treeFile, sequenceFile);
		if (pog == null || pog.equals(""))	// load graph structure from alignment file
			pog = sequenceFile;
		if (performMSA) {
			MSA alignment = new MSA(pog);
			pogAlignment = alignment.getMSAGraph();
		} else
			pogAlignment = new POGraph(pog, sequenceFile);
		pog = null;

		// perform inference
		if (jointInference) {
			marginalDistributions = null;
			queryBNJoint();
		} else if (marginalNode != null && phyloTree.find(marginalNode) != null) {
			queryBNMarginal(marginalNode);
		} else {
			if (marginalNode == null)
				System.out.println("No node was specified for the marginal inference: inferring the root node");
			else
				throw new RuntimeException("Incorrect internal node label provided for marginal reconstruction: " + marginalNode + " tree: " + phyloTree.toString());
			marginalNode = phyloTree.getRoot().getLabel().toString();
			queryBNMarginal(phyloTree.getRoot().getLabel().toString());
		}
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param msa				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param parsimony			flag to identify gaps in reconstruction using parsimony (true) or maximum likelihood (false)
	 */
	private void performASR(POGraph msa, String treeFile, String sequenceFile, boolean jointInference, boolean parsimony) throws RuntimeException, IOException {
		loadData(treeFile, sequenceFile);

		pogAlignment = msa;

		// perform inference
		if (jointInference) {
			marginalDistributions = null;
			queryBNJoint();
		} else if (marginalNode != null && phyloTree.find(marginalNode) != null) {
			queryBNMarginal(marginalNode);
		} else {
			if (marginalNode == null)
				System.out.println("No node was specified for the marginal inference: inferring the root node");
			else
				throw new RuntimeException("Incorrect internal node label provided for marginal reconstruction: " + marginalNode + " tree: " + phyloTree.toString());
			marginalNode = phyloTree.getRoot().getLabel().toString();
			queryBNMarginal(phyloTree.getRoot().getLabel().toString());
		}
	}

	/**
	 * Set the sequence of the given ancestral node.
	 *
	 * @param phyloNodeLabel	label of ancestral node
	 */
	private void populateTreeNodeSeq(String phyloNodeLabel) {
		EnumSeq.Gappy<Enumerable> seq = new EnumSeq.Gappy<>(Enumerable.aacid_ext);
		seq.setName(phyloNodeLabel);
		String s = getAncestor(phyloNodeLabel).getSupportedSequence();
		seq.setInfo(s);
		Object[] chars = new Object[s.length()];
		for (int c = 0; c < s.length(); c++)
			chars[c] = s.toCharArray()[c];
		seq.length();
		seq.set(chars);
		phyloTree.find(phyloNodeLabel).setSequence(seq);
	}

	/**
	 * Generates the ancestor by applying the inference changes stored in the ancestralInferences log.
	 *
	 * @param label	Ancestral sequence label
	 * @return		POGStructure of ancestor
	 */
	private POGraph getAncestor(String label) {
		POGraph ancestor = new POGraph(pogAlignment);
		if (label.equalsIgnoreCase("root"))
			label = (String)phyloTree.getRoot().getLabel();
		List<Inference> inferences = ancestralInferences.get(label);
		for (Inference inferredBase : inferences)
			if (ancestor.setCurrent(inferredBase.pogId)) {
				// set inferred base, or remove node if inferred to be removed
				if (inferredBase.pogId != -1)
					if (inferredBase.base == '-')
						// remove node from ancestor
						ancestor.removeNode();
					else
						ancestor.setBase(inferredBase.base);
				// if node is still there, remove transitions that are not inferred using maximum parsimony
				if (ancestor.getCurrentId() == inferredBase.pogId)
					for (Integer nextId : ancestor.getNextIDs())
						if (!inferredBase.transitions.contains(nextId))
							ancestor.removeNextTransition(nextId);
			}

		// if marginal, update each node with character distribution
		if (marginalNode != null)
			for (Integer nodeId : ancestor.getNodeIDs()) {
				double[] dist = marginalDistributions[nodeId].get();
				HashMap<Character, Double> distribution = new HashMap<>();
				for (int ind = 0; ind < dist.length; ind++)
					distribution.put((char) marginalDistributions[nodeId].getDomain().get(ind), dist[ind]);
				ancestor.setCurrent(nodeId);
				ancestor.setCharacterDistribution(distribution);
			}
		return ancestor;
	}

	/**
	 * Loads the phylogenetic tree into the BN tree structure, performs MSA using the input sequences and generates the partial order alignment graph of the MSA
	 * 
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile	filepath to the sequences (expected extension .aln, .fa or .fasta)
	 */
	private void loadData(String treeFile, String sequenceFile) throws IOException {
		// load extant sequences
		BufferedReader aln_file = new BufferedReader(new FileReader(sequenceFile));
		String line = aln_file.readLine();
		if (line.startsWith("CLUSTAL")) {
			extantSequences = EnumSeq.Gappy.loadClustal(sequenceFile, Enumerable.aacid_ext);
		} else if (line.startsWith(">")) {
			extantSequences = EnumSeq.Gappy.loadFasta(sequenceFile, Enumerable.aacid_ext, '-');
		} else {
			throw new RuntimeException("Alignment should be in Clustal or Fasta format");
		}
		aln_file.close();

		// initialise ancestral sequence labels and list for tracking changes to the POG structure for the ancestral nodes
		ancestralSeqLabels = new ArrayList<>();
		ancestralInferences = new HashMap<>();

		// create phylogenetic tree structure
		phyloTree = PhyloTree.loadNewick(treeFile);

		// Check if there are duplicate extant node names in the phylogenetic tree
		// Duplicate extant node names not allowed - will influence reconstruction outcomes
		PhyloTree.Node[] nodes = phyloTree.toNodesBreadthFirst(); //tree to nodes - recursive
		List<PhyloTree.Node> extNodes = new ArrayList<>();
		for (PhyloTree.Node node : nodes)
			if (node.getChildren().toArray().length == 0)
				extNodes.add(node);
			else
				ancestralSeqLabels.add((String) node.getLabel());

		Set<String> eNodes = new HashSet<>();
		for (PhyloTree.Node en : extNodes)
			try {
				if (!eNodes.add(en.getLabel().toString()))
					throw new RuntimeException("Extant node names must be unique - " + en.getLabel().toString() + " is duplicated");
			} catch (Exception e) {
				e.printStackTrace();
			}

		// Check if there are duplicate sequence names in the extant sequences
		// Duplicate sequence names not allowed - will influence reconstruction outcomes
		Set<String> seqNames = new HashSet<>();
		for (EnumSeq seq : extantSequences)
			if (!seqNames.add(seq.getName()))
				throw new RuntimeException("Sequence names must be unique - " + seq.getName() + " is duplicated");

		// Check if the provided extant sequences match up to the provided tree
		if (!eNodes.equals(seqNames)) {
		    // find labels that don't match
            String seqLabels = "";
            for (String seqLabel : seqNames)
                if (!eNodes.contains(seqLabel))
                    seqLabels += " " + seqLabel;
            String eLabels = "";
            for (String eLabel : eNodes)
                if (!seqNames.contains(eLabel))
                    eLabels += " " + eLabel;
            throw new RuntimeException("The sequence names in the provided alignment must all have a match" +
                    " in the provided tree. Unique labels in the alignment: " + seqLabels + ": unique labels in the tree: " + eLabels);
        }

		// save sequence information in internal nodes of the phylogenetic tree
		for (EnumSeq.Gappy<Enumerable> extant : extantSequences)
			phyloTree.find(extant.getName()).setSequence(extant);
	}
	
	/**
	 * Create Bayesian network for node in the partial order alignment graph that contains multiple characters. 
	 * 
	 * @return	Bayesian networks for node position in pogAlignment
	 */
	private PhyloBNet createCharacterNetwork(){
		// create a bayesian network with the phylogenetic tree structure and the JTT substitution model for amino acids
		PhyloBNet phyloBN = PhyloBNet.create(phyloTree, new JTT());														
		Map<Integer, Character> sequenceCharacterMapping = pogAlignment.getSequenceCharacterMapping();
		
		// for all extant sequences, if sequence is in this alignment, find where the location is in the phylogenetic tree and assign base to that position in the bayesian network
		// for all sequences not in this alignment, find where the location is in the phylogenetic tree and assign a gap to that position in the bayesian network
		for (int extantSeq = 0; extantSeq < extantSequences.size(); extantSeq++)
			if (sequenceCharacterMapping.containsKey(extantSeq))
				// find where the sequence is in the BN and set the base character
				phyloBN.getBN().getNode(extantSequences.get(extantSeq).getName()).setInstance(sequenceCharacterMapping.get(extantSeq));	
			else {
				// sequence is not part of character inference, check for gap
				SubstNode snode = (SubstNode)phyloBN.getBN().getNode(extantSequences.get(extantSeq).getName());
				snode.setGap(true);
			}

		// remove all nodes that are 'gaps'
		phyloBN.purgeGaps();
		
		return phyloBN;
	}

	private Map<String, Integer[]> getPhyloTransitions() {

		Integer nodeId = pogAlignment.getCurrentId();

		Map<String, Integer[]> phyloTransition = new HashMap<>();

		// populate tree for transitional inference using max parsimony
		// construct string[] of extant names
		// construct object[] of 'next' nodes (for transitions)
		String[] extants = new String[extantSequences.size()];
		Integer[] nextNodes = new Integer[extantSequences.size()];
		Map<Integer, List<Integer>> nodeSeqs = pogAlignment.getSequenceNodeMapping();
		for (int seqId = 0; seqId < extantSequences.size(); seqId++) {
			extants[seqId] = extantSequences.get(seqId).getName();
			if (nodeId == -1) // initial node, set next node as the starting node
				nextNodes[seqId] = nodeSeqs.get(seqId).get(0);
			else {
				int ind = nodeSeqs.get(seqId).indexOf(nodeId);
				if (ind != -1 && ind + 1 < nodeSeqs.get(seqId).size())
					nextNodes[seqId] = nodeSeqs.get(seqId).get(ind + 1);
			}
		}

		phyloTree.setContentByParsimony(extants, nextNodes);

		for (String phyloNode : ancestralSeqLabels) {
			List<Object> values = phyloTree.find(phyloNode).getValues();
			if (values == null) {
				values = new ArrayList<>();
				values.add(null);
			}
			Integer[] ids = new Integer[values.size()];
			for (int i = 0; i < values.size(); i++)
				ids[i] = (Integer) values.get(i);
			phyloTransition.put(phyloNode, ids);
		}

		return phyloTransition;
	}

	/**
	 * Infer gap/base character of each partial order alignment graph structure at each internal node of the phylogenetic tree using joint inference.
	 */
	private void queryBNJoint(){

		// infer base/gap of each aligned node
		List<Integer> nodeIDs = new ArrayList<>();
		nodeIDs.add(-1);
		nodeIDs.addAll(pogAlignment.getNodeIDs());
		rates = new Double[nodeIDs.get(nodeIDs.size()-1) + 1]; //Rate matrix
		for (Integer nodeId : nodeIDs) {
			pogAlignment.setCurrent(nodeId);

			//if ((pogAlignment.getCurrentBase() != null || pogAlignment.getCurrentId() == -1) && pogAlignment.getNumEdgesOut() == 1)
			//	continue;
			
			Variable.Assignment[] charAssignments = new Variable.Assignment[]{};

			// get gap or character inference at phylogenetic nodes
			Map<String, Integer[]> phyloTransition = getPhyloTransitions();

			// perform character inference if not the dummy initial node
			if (nodeId != -1) {
				VarElim ve = new VarElim();
				PhyloBNet charNet = createCharacterNetwork();
				rates[nodeId] = charNet.getRate();
				ve.instantiate(charNet.getBN());
				charAssignments = getJointAssignment(ve);
			}

			// for each node in the phylogenetic tree, if character is inferred at position, set inferred base,
			// otherwise gap is inferred, remove from alignment and set transitions of previous nodes to the
			// inferred transition
			for (String phyloNode : ancestralSeqLabels) {
				Character base = '-';
				// check for inferred base character
				for (Variable.Assignment varassign : charAssignments)
					if (phyloNode.equals(varassign.var.getName())) {
						base = (char) varassign.val;
						break;
					}

				// store inferred transitions
				List<Integer> transitionIds = new ArrayList<>();
				Integer[] transitions = phyloTransition.get(phyloNode);
				for (int i = 0; i < transitions.length; i++)
					transitionIds.add(transitions[i]);

				if (!ancestralInferences.containsKey(phyloNode))
					ancestralInferences.put(phyloNode, new ArrayList<>());
				ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), base, transitionIds));
			}
		}

		// save sequence information in internal nodes of the phylogenetic tree
		for (String phyloNodeLabel : ancestralSeqLabels)
			populateTreeNodeSeq(phyloNodeLabel);
	}

	/**
	 * Infer gap/base character of each partial order alignment graph structure at each internal node of the phylogenetic tree using marginal inference.
	 *
	 * @param phyloNode		ancestral node to perform marginal reconstruction of
	 */
	private void queryBNMarginal(String phyloNode) {
		marginalDistributions = new EnumDistrib[pogAlignment.getNumNodes()];

		// infer base/gap of each aligned node

		List<Integer> nodeIDs = new ArrayList<>();
		nodeIDs.add(-1);
		nodeIDs.addAll(pogAlignment.getNodeIDs());
		for (Integer nodeId : nodeIDs) {
			pogAlignment.setCurrent(nodeId);

			// get gap or character inference at phylogenetic nodes
			Map<String, Integer[]> phyloTransition = getPhyloTransitions();

			// perform character inference if not the dummy initial node
			Character base = '-';
			if (nodeId != -1) {
				VarElim ve = new VarElim();
				PhyloBNet phyloNet = createCharacterNetwork();
				ve.instantiate(phyloNet.getBN());
				// get character variable representing current phylogenetic tree node using marginal probability
				EnumVariable charNode = null;
				for (EnumVariable internalNode : phyloNet.getInternal())
					if (internalNode.getName().equalsIgnoreCase(marginalNode))
						charNode = internalNode;
				if (charNode != null) {
					marginalDistributions[nodeId] = getMarginalDistrib(ve, charNode);
					base = (char) marginalDistributions[nodeId].getMax();
				}
			}

			// store inferred transitions
			List<Integer> transitionIds = new ArrayList<>();
			Integer[] transitions = phyloTransition.get(phyloNode);
			for (int i = 0; i < transitions.length; i++)
				transitionIds.add(transitions[i]);

			if (!ancestralInferences.containsKey(phyloNode))
				ancestralInferences.put(phyloNode, new ArrayList<>());
			ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), base, transitionIds));
		}

		// save sequence information in internal nodes of the phylogenetic tree
		populateTreeNodeSeq(marginalNode);
		ancestralSeqLabels = new ArrayList<>();
		ancestralSeqLabels.add(marginalNode);
	}

	/**
	 * Gets the joint probability assignment of instantiated bayesian network associated with VarElim ve
	 * 
	 * @param ve	Instantiated variable elimination for bayesian inference of joint probabilities
	 * @return		likelihood assignment of each node from ve
	 */
    private Variable.Assignment[] getJointAssignment(VarElim ve) {
        Query q_joint = ve.makeMPE();
        CGTable r_joint = (CGTable)ve.infer(q_joint);
        return r_joint.getMPE();
    }
    
    /**
     * Gets the marginal distribution of the queryNode using instantiated bayesian network associated with VarElim ve
     * 
     * @param ve			instantiated bayesian network
     * @param queryNode		node to query
     * @return				marginal distribution of query node using ve
     */
    private EnumDistrib getMarginalDistrib(VarElim ve, EnumVariable queryNode) {
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

	/**
	 * Check to see if a node has an inferred character that is being contributed by one but not both children
	 * @param node node to check
	 */

	public void checkBranchIsolation(String node){
		System.out.println("Yo");

	}

    /**
     * Helper class to store changes to an ancestral graph node
     * 
     * Information:
     * 		- POG structure index 
     * 		- Inferred base character: base character that is inferred or '-' to represent a gap (i.e. that the node needs to be deleted when updating the structure)
     */
    private class Inference {
    	Integer pogId;
    	char base;
    	List<Integer> transitions;

    	public Inference(int id, char ch, List<Integer> tr){
    		pogId = id;
    		base = ch;
    		transitions = tr;
    	}
    	public String toString(){
    		return pogId + "->" + base;
    	}
    }
}