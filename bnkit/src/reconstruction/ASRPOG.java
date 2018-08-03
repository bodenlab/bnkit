package reconstruction;


import alignment.MSA;
import api.PartialOrderGraph;
import bn.alg.CGTable;
import bn.alg.Inference;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.PhyloBNet;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.Dayhoff;
import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.LG;
import bn.ctmc.matrix.WAG;
import bn.prob.EnumDistrib;
import com.sun.xml.internal.bind.v2.TODO;
import dat.*;
import dat.file.AlnWriter;
import dat.file.FastaWriter;
import json.JSONArray;
import json.JSONObject;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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

	private PhyloTree phyloTree = null; 								// Phylogenetic tree structure
	private List<EnumSeq.Gappy<Enumerable>> extantSequences = null;		// List of sequences (label,bases)
	private List<String> ancestralSeqLabels = null;						// Ancestral sequences labels (internal nodes of phylogenetic tree structure)
	private POGraph pogAlignment = null;								// partial order alignment graph structure template
	private EnumDistrib[] marginalDistributions = null; 				// Marginal distributions for nodes if doing a marginal reconstruction
	private String marginalNode = null;									// Label of node to perform marginal reconstruction (if applicable)
	private Map<String, List<Inference>> ancestralInferences = null;	// stores updates to the POGStructure for the ancestral node <node label, changes>
	private Double[] rates = null; 										// Rates at positions in alignment
	private String model = "JTT";										// Evolutionary model to use for character inference
	private int threads = 1;											// Number of threads to use for performing the reconstruction
	private boolean performMSA = false;									// Flag to track whether the input sequences required alignment or not

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile		filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 * @param model				evolutionary model to use for inference (e.g. JTT, LG, WAG, Dayhoff)
	 * @param threads			number of threads to use for reconstruction. Default: 1
	 */
	public ASRPOG(String alignmentFile, String treeFile, boolean jointInference, boolean performMSA, String model, int threads) throws IOException, InterruptedException {
		setupASRPOG(model, null, threads);
		performASR("", treeFile, alignmentFile, jointInference, performMSA);
	}

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile	filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param marginalNode	node label for maginal inference
	 * @param performMSA		flag for indicating whether to perform the alignment prior to reconstruction
	 * @param model				evolutionary model to use for inference (e.g. JTT, LG, WAG, Dayhoff)
	 * @param threads			number of threads to use for reconstruction. Default: 1
	 */
	public ASRPOG(String alignmentFile, String treeFile, String marginalNode, boolean performMSA, String model, int threads) throws IOException, InterruptedException {
		setupASRPOG(model, marginalNode, threads);
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
	 * @param model				evolutionary model to use for inference (e.g. JTT, LG, WAG, Dayhoff)
	 * @param threads			number of threads to use for reconstruction. Default: 1
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean performMSA, String model, int threads) throws IOException, InterruptedException {
		setupASRPOG(model, null, threads);
		performASR(pog, treeFile, sequenceFile, jointInference, performMSA);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 * This one
	 * @param msa				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param model				evolutionary model to use for inference (e.g. JTT, LG, WAG, Dayhoff)
	 * @param threads			number of threads to use for reconstruction. Default: 1
	 */
	public ASRPOG(POGraph msa, String treeFile, String sequenceFile, boolean jointInference, String model, int threads) throws IOException, InterruptedException {
		setupASRPOG(model, null, threads);
		performASR(msa, treeFile, sequenceFile, jointInference);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog			POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile	filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param marginalNode	node label for maginal inference
	 * @param performMSA	flag for indicating whether to perform the alignment prior to reconstruction
	 * @param model			evolutionary model to use for inference (e.g. JTT, LG, WAG, Dayhoff)
	 * @param threads		number of threads to use for reconstruction. Default: 1
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, String marginalNode, boolean performMSA, String model, int threads) throws IOException, InterruptedException {
		setupASRPOG(model, marginalNode, threads);
		performASR(pog, treeFile, sequenceFile, false, performMSA);
	}


	public ASRPOG(String model, int threads) {
		setupASRPOG(model, null, threads);
	}

	public ASRPOG(String model, int threads, String marginalNode) {
		setupASRPOG(model, marginalNode, threads);
	}

	public ASRPOG(String model, int threads, JSONObject inferences, List<EnumSeq.Gappy<Enumerable>> sequences, String tree) {
		setupASRPOG(model, null, threads);
		extantSequences = new ArrayList<>(sequences);
		importInferencesFromJSON(inferences);
		phyloTree = PhyloTree.parseNewick(tree);
		pogAlignment = new POGraph(sequences);
	}

	public void runReconstruction(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean performMSA) throws IOException, InterruptedException {
		performASR(pog, treeFile, sequenceFile, jointInference, performMSA);
	}

	public void runReconstruction(POGraph msa, String treeFile, String sequenceFile, boolean jointInference) throws IOException, InterruptedException {
		performASR(msa, treeFile, sequenceFile, jointInference);
	}

	public void runReconstruction(String treeNewick, List<EnumSeq.Gappy<Enumerable>> sequences, boolean jointInference, POGraph msa) throws InterruptedException {

		extantSequences = new ArrayList<>(sequences);

		// create phylogenetic tree structure
		phyloTree = PhyloTree.parseNewick(treeNewick);

		// Check if there are duplicate extant node names in the phylogenetic tree
		// Duplicate extant node names not allowed - will influence reconstruction outcomes
		// Check if there are duplicate sequence names in the extant sequences
		// Duplicate sequence names not allowed - will influence reconstruction outcomes
		// Check if the provided extant sequences match up to the provided tree
		checkData();

		if (msa == null)
			pogAlignment = new POGraph(extantSequences);
		else
			pogAlignment = new POGraph(msa);

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
	public POGraph getMSAGraph() {
		return new POGraph(this.pogAlignment);
	}

	/**
	 * Get the multiple sequence alignment partial order graph.
	 *
	 * @return	Partial order graph representing sequence alignment
	 */
	public PartialOrderGraph getPartialOrderGraph() {
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
	 * Save the ASR Object in JSON format to the given filepath.
	 *
	 * @param filepath		Filepath to save JSON export details
	 * @throws IOException
	 */
	public void saveJSONExport(String filepath) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath + "export.json", false));
		bw.write(exportToJSON().toString());
		bw.close();
	}

	public JSONObject exportInferencesToJSON() {
		JSONObject infJSON = new JSONObject();

		// add ancestral inferences
		JSONArray allInferences = new JSONArray();
		for (String ancestor : ancestralInferences.keySet()) {
			JSONObject ancestorsInferences = new JSONObject();
			JSONArray inferences = new JSONArray();
			for (Inference i : ancestralInferences.get(ancestor)) {
				JSONObject infDetails = new JSONObject();
				infDetails.put("id", i.pogId);
				infDetails.put("base", i.base + "");
				JSONArray transitions = new JSONArray();
				for (Integer transition : i.transitions)
					transitions.put(transition);
				infDetails.put("transitions", transitions);
				inferences.put(infDetails);
			}
			ancestorsInferences.put("label", ancestor);
			ancestorsInferences.put("inferences", inferences);
			allInferences.put(ancestorsInferences);
		}
		infJSON.put("inferences", allInferences);

		return infJSON;
	}

	public void importInferencesFromJSON(JSONObject inferences) {
		ancestralInferences = new HashMap<>();
		ancestralSeqLabels = new ArrayList<>();

		JSONArray allInfArray = inferences.getJSONArray("inferences");
		for (Object ancestor : allInfArray){
			JSONObject anc = (JSONObject)ancestor;
			List<Inference> infs = new ArrayList<>();
			JSONArray ancInfArray = anc.getJSONArray("inferences");
			for (Object inf : ancInfArray) {
				JSONObject infJSON = (JSONObject)inf;
				List<Integer> transitions = new ArrayList<>();
				JSONArray jTransitions = infJSON.getJSONArray("transitions");
				for (int i = 0; i < jTransitions.length(); i++)
					transitions.add(jTransitions.getInt(i));
				Character base = infJSON.getString("base").toCharArray()[0];
				infs.add(new Inference(infJSON.getInt("id"), base, transitions));
			}
			ancestralInferences.put(anc.getString("label"), infs);
			ancestralSeqLabels.add(anc.getString("label"));
		}
	}

	/**
	 * Export the ASR object to JSON format.
	 *
	 * @return	JSON representation of current ASR details
	 */
	public JSONObject exportToJSON() {
		JSONObject asrJSON = new JSONObject();

		// Add extant sequences
		JSONArray extants = new JSONArray();
		for (EnumSeq.Gappy<Enumerable> extantSeq : extantSequences) {
			JSONObject extant = new JSONObject();
			extant.put("label", extantSeq.getName());
			extant.put("sequence", extantSeq.toString());
			extants.put(extant);
		}
		asrJSON.put("extants", extants);

		// add ancestral inferences
		JSONArray allInferences = new JSONArray();
		for (String ancestor : ancestralInferences.keySet()) {
			JSONObject ancestorsInferences = new JSONObject();
			JSONArray inferences = new JSONArray();
			for (Inference i : ancestralInferences.get(ancestor)) {
				JSONObject infDetails = new JSONObject();
				infDetails.put("id", i.pogId);
				infDetails.put("base", i.base + "");
				JSONArray transitions = new JSONArray();
				for (Integer transition : i.transitions)
					transitions.put(transition);
				infDetails.put("transitions", transitions);
				inferences.put(infDetails);
			}
			ancestorsInferences.put("label", ancestor);
			ancestorsInferences.put("inferences", inferences);
			allInferences.put(ancestorsInferences);
		}
		asrJSON.put("inferences", allInferences);

		asrJSON.put("model", model);
		asrJSON.put("threads", threads);
		asrJSON.put("marginal_node", marginalNode);
		asrJSON.put("phylotree", phyloTree.toString());

		return asrJSON;
	}


	/**
	 * Save the reconstructed sequences that have the most support through the partial order graph in FASTA format.
	 * Saved in output path as "reconstructed_sequences.fasta"
	 *
	 * @param filepath	Output filepath
	 * @param gappy		Flag to save gappy sequence (true) or not (false)
	 */
	public void saveSupportedAncestors(String filepath, boolean gappy) throws IOException {
		String[] labels = new String[ancestralSeqLabels.size()];
		ancestralSeqLabels.toArray(labels);
		saveSupportedAncestors(filepath, labels, gappy);
	}

	/**
	 * Save the reconstructed sequences that have the most support through the partial order graph in FASTA format.
	 * Saved in output path as "reconstructed_sequences.fasta"
	 *
	 * @param filepath	Output filepath
	 * @param nodes		Array of ancestral nodes to save
	 * @param gappy		Flag to save gappy sequence (true) or not (false)
	 */
	public void saveSupportedAncestors(String filepath, String[] nodes, boolean gappy) throws IOException {
		Map<String, String> ancestralSeqs = new HashMap<>();
		for (String node : nodes) {
			POGraph ancestor = getAncestor(node);
			ancestralSeqs.put(node, ancestor.getSupportedSequence(gappy));
		}
		File directory = new File(filepath);
		if (! directory.exists())
			directory.mkdir();

		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath + "_recon.fa", false));
		bw.write(">" + nodes[nodes.length - 1]);
		bw.newLine();
		bw.write(ancestralSeqs.get(nodes[nodes.length - 1]));
		bw.newLine();
		bw.newLine();
		for (int i = 0; i < nodes.length - 1; i ++) {
			bw.write(">" + nodes[i]);
			bw.newLine();
			bw.write(ancestralSeqs.get(nodes[i]));
			bw.newLine();
			bw.newLine();
		}
		/*
		for (String node : ancestralSeqs.keySet()) {
			bw.write(">" + node);
			bw.newLine();
			bw.write(ancestralSeqs.get(node));
			bw.newLine();
			bw.newLine();
		}
		*/
		bw.close();
	}


	/**
	 * Save the reconstructed sequences that have the most support through the partial order graph in FASTA format.
	 * Saved in output path as "reconstructed_sequences.fasta"
	 *
	 * @param filepath	Output filepath
	 * @param label		label of ancestor to save (tree node label)
	 * @param gappy		Flag to save gappy sequence (true) or not (false)
	 */
	public void saveSupportedAncestor(String filepath, String label, boolean gappy) throws IOException {
		if (label.equalsIgnoreCase("root"))
			label = (String)phyloTree.getRoot().getLabel();

		POGraph ancestor = getAncestor(label);

		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath + "_recon.fa", false));
		bw.write(">" + label);
		bw.newLine();
		bw.write(ancestor.getSupportedSequence(gappy));
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
		ArrayList<EnumSeq.Gappy<Enumerable>> allSeqs = new ArrayList<>();

		for (String node : ancestralSeqLabels) {
			POGraph ancestor = getAncestor(node);
			EnumSeq.Gappy<Enumerable> seq = new EnumSeq.Gappy<>(Enumerable.aacid_ext);
			String seqName = ancestor.getSupportedSequence(true);
			String seqLab = node;
			seq.setName(seqLab + " " + seqName + ";"); //Newick strings require a ';' to indicate completion
			Object[] s = new Object[seqName.length()];
			for (int n = 0; n < seqName.length(); n++)
				s[n] = seqName.toCharArray()[n];
			seq.set(s);
			allSeqs.add(seq);
		}

		EnumSeq.Gappy<Enumerable>[] seqs = new EnumSeq.Gappy[allSeqs.size()];
		for (int n = 0; n < allSeqs.size(); n++)
			seqs[n] = allSeqs.get(n);
		if (format.equalsIgnoreCase("clustal")) {
			AlnWriter aw = new AlnWriter(filepath + ".aln");
			aw.save(seqs);
			aw.close();
		} else {
			FastaWriter fw = new FastaWriter(filepath + ".fa");
			fw.save(seqs);
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

	/**
	 * Return a collection of the children of a given node
	 *
	 * @param node	node to get children of
	 */
	public Collection<PhyloTree.Node> getChildren(String node) {
		return this.phyloTree.find(node).getChildren();
	}


	public PhyloTree.Node getParent(String node) {
		return this.phyloTree.find(node).getParent();
	}


	public String getReconstructedNewick() {
		return phyloTree.getRoot().toString();
	}

	/**
	 * Return a dictionary mapping ancestral label to inferred ancestral sequence
	 *
	 */

	public Map<String, String> getAncestralDict(){
		Map<String, String> ancestralDict = new HashMap<>();
		for (String label : this.ancestralSeqLabels) {
			List<ASRPOG.Inference> inferences = this.ancestralInferences.get(label);
			String inferenceString = "";
			for (ASRPOG.Inference inference : inferences) {
				inferenceString += inference.base;
			}

			ancestralDict.put(label, inferenceString);
		}
		return ancestralDict;
	}

	public String getRootLabel() {
		return (String)phyloTree.getRoot().getLabel();
	}

	public Map<String, List<Inference>> getAncestralInferences(){
		return this.ancestralInferences;
	}

	public List<String> getAncestralSeqLabels(){
		return this.ancestralSeqLabels;

	}

	/**
	 *
	 * Return the marginal distributions
	 *
	 */
	public EnumDistrib[] getMarginalDistributions(){
		return this.marginalDistributions;
	}


	/**
	 * Get the current node ID of the alignment PO graph.
	 *
	 * @return	ID of current node in alignment PO graph
	 */
	public int getGraphReconNodeId() {
		return pogAlignment.getCurrentId();
	}

	/* ****************************************************************************************************************************************************
	 * 																PRIVATE METHODS
	 * ****************************************************************************************************************************************************/

	/**
	 * Set conditions for reconstruction
	 * @param model		Evolutionary model
	 * @param node		Marginal node (or null)
	 * @param threads	Number of threads for processing
	 */
	private void setupASRPOG(String model, String node, int threads) {
		this.threads = threads;
		if (model != null)
			this.model = model;
		this.marginalNode = node;

		// initialise ancestral sequence labels and list for tracking changes to the POG structure for the ancestral nodes
		ancestralSeqLabels = new ArrayList<>();
		ancestralInferences = new HashMap<>();
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param pog				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 */
	private void performASR(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean performMSA) throws RuntimeException, IOException, InterruptedException {
		this.performMSA = performMSA;
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
	 */
	private void performASR(POGraph msa, String treeFile, String sequenceFile, boolean jointInference) throws RuntimeException, IOException, InterruptedException {
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
	 * Generates the ancestor by applying the inference changes stored in the ancestralInferences log.
	 *
	 * @param label	Ancestral sequence label
	 * @return		POGStructure of ancestor
	 */
	public POGraph getAncestor(String label) {
		POGraph ancestor = new POGraph(pogAlignment);
		if (label.equalsIgnoreCase("root"))
			label = (String)phyloTree.getRoot().getLabel();
		List<Inference> inferences = ancestralInferences.get(label);
		for (Inference inferredBase : inferences)
			if (ancestor.setCurrent(inferredBase.pogId)) {
				// set inferred base, or remove node if inferred to be removed
				if (!inferredBase.pogId.equals(ancestor.getFinalNodeID()) && !inferredBase.pogId.equals(pogAlignment.getInitialNodeID()))
					if (inferredBase.base == '-')
						// remove node from ancestor
						ancestor.removeNode();
					else
						ancestor.setBase(inferredBase.base);
				// if node is still there, check the parsimonious transitions. Because we are doing both backwards and forwards parsimony,
				// get the intersection of the previous/next nodes and all inferred transitions of the current node
				if (ancestor.getCurrentId().equals(ancestor.getFinalNodeID()) && !inferredBase.pogId.equals(ancestor.getFinalNodeID()) || (!ancestor.getCurrentId().equals(ancestor.getFinalNodeID()) && !ancestor.getCurrentId().equals(inferredBase.pogId)))
					continue;


				// find union of next/previous transitions and remove transitions that are not inferred
				// 'next':
				ArrayList<Integer> keepNext = new ArrayList<>();
				for (Integer nextId : ancestor.getNextIDs()) {
					// Identify if edge is reciprocated by backwards parsimony (if so, flag)
					Inference next = null;
					for (Inference i : inferences)
						if ((i.pogId.equals(ancestor.getFinalNodeID()) && (nextId.equals(ancestor.getFinalNodeID()))) || i.pogId.equals(nextId)) {
							next = i;
							break;
						}
					if (inferredBase.transitions.contains(nextId) && next != null && next.transitions.contains(inferredBase.pogId))
						ancestor.setReciprocated(nextId);
					if (next != null && next.transitions.contains(inferredBase.pogId))
						keepNext.add(next.pogId);
				}

				// take the union of the inferred transitions (i.e. use all provided transitions from maximum parsimony)
				for (Integer nextId : ancestor.getNextIDs())
					if (inferredBase.transitions.contains(nextId) && !keepNext.contains(nextId))
						keepNext.add(nextId);
				// check previous transitions of future nodes to see if there has been a parsimonious edge to this node
				// if so, add to keep
				for (Integer nextId : ancestor.getNextIDs())
					if (!keepNext.contains(nextId) && !nextId.equals(pogAlignment.getFinalNodeID()))
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

		try {

			if (line.startsWith("CLUSTAL")) {
				extantSequences = EnumSeq.Gappy.loadClustal(sequenceFile, Enumerable.aacid_ext);
			} else if (line.startsWith(">")) {
				extantSequences = EnumSeq.Gappy.loadFasta(sequenceFile, Enumerable.aacid_ext, '-');
			} else {
				throw new RuntimeException("Incorrect sequence or alignment format (requires FASTA or Clustal format .aln, .fa or .fasta)");
			}
			aln_file.close();

		}

		catch (NullPointerException npe) {
			throw new RuntimeException("Error: Incorrect sequence or alignment format (requires FASTA or Clustal format .aln, .fa or .fasta)");
//			System.exit(1);


		}


		// create phylogenetic tree structure
		phyloTree = PhyloTree.loadNewick(treeFile);

		// Check if there are duplicate extant node names in the phylogenetic tree
		// Duplicate extant node names not allowed - will influence reconstruction outcomes
		// Check if there are duplicate sequence names in the extant sequences
		// Duplicate sequence names not allowed - will influence reconstruction outcomes
		// Check if the provided extant sequences match up to the provided tree
		checkData();

		// save sequence information in internal nodes of the phylogenetic tree
		for (EnumSeq.Gappy<Enumerable> extant : extantSequences)
			phyloTree.find(extant.getName()).setSequence(extant);
	}

	private void checkData(){
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
			if (!eNodes.add(en.getLabel().toString()))
				throw new RuntimeException("Error: Extant node names must be unique.\nDuplicate names are - " + en.getLabel().toString());


		// Check if there are duplicate sequence names in the extant sequences
		// Duplicate sequence names not allowed - will influence reconstruction outcomes
		Set<String> seqNames = new HashSet<>();
		for (EnumSeq seq : extantSequences)
			if (!seqNames.add(seq.getName()))
				throw new RuntimeException("Error: Sequence names must be unique.\nDuplicate names are - " + seq.getName());

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
            throw new RuntimeException("Error: The sequence names in the provided alignment must all have a match" +
                    " in the provided tree.\nUnique labels in the alignment: " + seqLabels + ": unique labels in the tree: " + eLabels);
        }
	}
	
	/**
	 * Create Bayesian network for node in the partial order alignment graph that contains multiple characters. 
	 * 
	 * @return	Bayesian networks for node position in pogAlignment
	 */
	private PhyloBNet createCharacterNetwork(){
		// create a bayesian network with the phylogenetic tree structure and the JTT substitution model for amino acids
		PhyloBNet phyloBN;
		if (this.model.equalsIgnoreCase("Dayhoff"))
			phyloBN = PhyloBNet.create(phyloTree, new Dayhoff());
		else if (this.model.equalsIgnoreCase("LG"))
			phyloBN = PhyloBNet.create(phyloTree, new LG());
		else if (this.model.equalsIgnoreCase("WAG"))
			phyloBN = PhyloBNet.create(phyloTree, new WAG());
		else
			phyloBN = PhyloBNet.create(phyloTree, new JTT());
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

		Map<String, Integer[]> phyloTransition = new HashMap<>();

		// populate tree for transitional inference using max parsimony

		// get ordered list of unique transitions based on MSA and num. seqs on the 'out' edges: forward and backwards
		Map<String, Object> mapNext = pogAlignment.getNextMapping();				// map of extant label and 'next' transition
		Map<String, Object> mapPrevious = pogAlignment.getPrevMapping();			// map of extant label and 'previous' transition
		ArrayList<Integer> orderedNext = pogAlignment.getOrderedNext();
		Object[] uniqueForward = new Object[orderedNext.size()];
		orderedNext.toArray(uniqueForward);
		ArrayList<Integer> orderedPrev = pogAlignment.getOrderedPrev();
		Object[] uniqueBackward = new Object[orderedPrev.size()];
		orderedPrev.toArray(uniqueBackward);

		// 'Next' transitions
		phyloTree.setContentByParsimony(mapNext, uniqueForward);
		for (String phyloNode : ancestralSeqLabels) {
			List<Object> values = phyloTree.find(phyloNode).getValues();
			if (values == null) {
				values = new ArrayList<>();
				values.add(pogAlignment.getFinalNodeID());
			}
			Integer[] ids = new Integer[values.size()];
			for (int i = 0; i < values.size(); i++)
				ids[i] = (Integer) values.get(i);
			phyloTransition.put(phyloNode, ids);
		}

		// 'Previous' transitions
		if (orderedPrev.isEmpty())
			return phyloTransition;

		phyloTree.setContentByParsimony(mapPrevious, uniqueBackward);
		for (String phyloNode : ancestralSeqLabels) {
			List<Object> values = phyloTree.find(phyloNode).getValues();
			if (values == null) {
				values = new ArrayList<>();
				values.add(pogAlignment.getFinalNodeID());
			}
			Integer[] ids = new Integer[values.size() + phyloTransition.get(phyloNode).length];
			for (int i = 0; i < phyloTransition.get(phyloNode).length; i++)
				ids[i] = phyloTransition.get(phyloNode)[i];
			for (int i = 0; i < values.size(); i++)
				ids[i+phyloTransition.get(phyloNode).length] = (Integer) values.get(i);
			phyloTransition.put(phyloNode, ids);
		}

		return phyloTransition;
	}

	// MB's laptop on 2U1 data (145 seqs) queryBNJoint for different number of threads:
	// 0: 30 secs; 1: 29 secs; 2: 22 secs; 3: 21 secs; 4: 20 secs; 5: 19 secs; 6: 19 secs; 7: 19 secs; 8: 20 secs

	/**
	 * Infer gap/base character of each partial order alignment graph structure at each internal node of the phylogenetic tree using joint inference.
	 */
	private void queryBNJoint() throws InterruptedException {

//		long startTime = System.nanoTime();

		// infer base/gap of each aligned node
		List<Integer> nodeIDs = new ArrayList<>();
		nodeIDs.add(pogAlignment.getInitialNodeID());
		nodeIDs.addAll(pogAlignment.getNodeIDs());
		rates = new Double[nodeIDs.get(nodeIDs.size()-1) + 1]; //Rate matrix

		// First, create a queue with which all nodes will be inferred (by NodeID)
		Queue<Integer> nodeQueue = new PriorityQueue<>();
		for (Integer nodeId : nodeIDs)
			nodeQueue.add(nodeId);

		// Create a thread coordinator, set the number of threads it should be using
		JointInferenceExecutor batch = null;
		if (this.threads > 1)
			batch = new JointInferenceExecutor(this.threads);

		// Here's where all results will end up
		Map<Integer, Variable.Assignment[]> result = new HashMap<>(nodeIDs.size());

		// Populate the coordinator with jobs, and run them successively;
		// results are kept for later (associated with IDs) but the actual internals of the inference are not
		while (!nodeQueue.isEmpty()) {
			Integer nodeId = nodeQueue.poll(); // next job
			pogAlignment.setCurrent(nodeId);
			// perform character inference if not the dummy initial node or dummy final node
			if (!nodeId.equals(pogAlignment.getInitialNodeID())) {
				VarElim ve = new VarElim();
				PhyloBNet charNet = createCharacterNetwork();
				rates[nodeId] = charNet.getRate();
				ve.instantiate(charNet.getBN());
				if (this.threads <= 1) {
					Variable.Assignment[] charAssignments = getJointAssignment(ve);
					result.put(nodeId, charAssignments);
				} else
					batch.addJointInference(nodeId, ve);
			}
			if (batch != null && nodeQueue.isEmpty()) { // job list complete OR this was the last, so need to run batch
				Map<Integer, Variable.Assignment[]> rets = batch.run();
				result.putAll(rets);
			}
		}
		// last time through all nodes; this time serially, extracting results
		nodeIDs.add(pogAlignment.getFinalNodeID()); // add dummy final node for identifying backwards transitions
		for (Integer nodeId : nodeIDs) {
			pogAlignment.setCurrent(nodeId);
			Variable.Assignment[] charAssignments = new Variable.Assignment[]{};
			if (!nodeId.equals(pogAlignment.getFinalNodeID()) && !nodeId.equals(pogAlignment.getInitialNodeID()))
				charAssignments = result.get(nodeId);
			// get gap or character inference at phylogenetic nodes
			Map<String, Integer[]> phyloTransition = getPhyloTransitions();
			// for each node in the phylogenetic tree, if character is inferred at position, set inferred base,
			// otherwise gap is inferred, remove from alignment and set transitions of previous nodes to the
			// inferred transition
			for (String phyloNode : ancestralSeqLabels) {
				Character base = '-';
				// check for inferred base character
				for (Variable.Assignment varassign : charAssignments) {
					if (phyloNode.equals(varassign.var.getName())) {
						base = (char) varassign.val;
						break;
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
		}
//		long elapsedTimeNs = System.nanoTime() - startTime;
//		System.out.printf("Elapsed time in secs: %5.3f\n", elapsedTimeNs/1000000000.0);
	}

	// MB's laptop on 2U1 data (145 seqs) queryBNMarginal "N0" for different number of threads:
	// 0: 15 secs; 1: 16 secs; 2: 12 secs; 3: 11 secs; 4: 11 secs; 5: 10 secs; 6: 11 secs; 7: 12 secs; 8: 11 secs

	/**
	 * Infer gap/base character of each partial order alignment graph structure at each internal node of the phylogenetic tree using marginal inference.
	 *
	 * @param phyloNode		ancestral node to perform marginal reconstruction of
	 */
	private void queryBNMarginal(String phyloNode) throws InterruptedException {
		marginalDistributions = new EnumDistrib[pogAlignment.getNumNodes()];
//		long startTime = System.nanoTime();

		// infer base/gap of each aligned node

		List<Integer> nodeIDs = new ArrayList<>();
		nodeIDs.add(pogAlignment.getInitialNodeID());
		nodeIDs.addAll(pogAlignment.getNodeIDs());

		// First, create a queue with which all nodes will be inferred (by NodeID)
		Queue<Integer> nodeQueue = new PriorityQueue<>();
		for (Integer nodeId : nodeIDs)
			nodeQueue.add(nodeId);

		// Create a thread coordinator, set the number of threads it should be using
		MarginalInferenceExecutor batch = null;
		if (this.threads > 1)
			batch = new MarginalInferenceExecutor(this.threads);

		// Here's where all results will end up
		Map<Integer, EnumDistrib> result = new HashMap<>(nodeIDs.size());

		// Populate the coordinator with jobs, and run them successively;
		// results are kept for later (associated with IDs) but the actual internals of the inference are not
		while (!nodeQueue.isEmpty()) {
			Integer nodeId = nodeQueue.poll(); // next job
			pogAlignment.setCurrent(nodeId);
			// perform character inference if not the dummy initial node
			if (!nodeId.equals(pogAlignment.getInitialNodeID())) {
				VarElim ve = new VarElim();
				PhyloBNet phyloNet = createCharacterNetwork();
				ve.instantiate(phyloNet.getBN());
				// get character variable representing current phylogenetic tree node using marginal probability
				EnumVariable charNode = null;
				for (EnumVariable internalNode : phyloNet.getInternal())
					if (internalNode.getName().equalsIgnoreCase(marginalNode))
						charNode = internalNode;
				if (charNode != null) {
					if (this.threads <= 1)
						result.put(nodeId, getMarginalDistrib(ve, charNode));
					else
						batch.addMarginalInference(nodeId, ve, charNode);
				}
				if (batch != null && nodeQueue.isEmpty()) { // job list complete OR this was the last, so need to run batch
					Map<Integer, EnumDistrib> rets = batch.run();
					result.putAll(rets);
				}
			}
		}
		for (Map.Entry<Integer, EnumDistrib> entry : result.entrySet()) {
			marginalDistributions[entry.getKey()] = entry.getValue();
		}

		nodeIDs.add(pogAlignment.getFinalNodeID()); // add dummy final node for identifying backwards transitions
		for (Integer nodeId : nodeIDs) {
			pogAlignment.setCurrent(nodeId);
			Character base = '-';

			// get gap or character inference at phylogenetic nodes
			Map<String, Integer[]> phyloTransition = getPhyloTransitions();

			// store inferred transitions
			List<Integer> transitionIds = new ArrayList<>();
			Integer[] transitions = phyloTransition.get(phyloNode);
			for (int i = 0; i < transitions.length; i++)
				transitionIds.add(transitions[i]);

			if (!ancestralInferences.containsKey(phyloNode))
				ancestralInferences.put(phyloNode, new ArrayList<>());
			if (!nodeId.equals(pogAlignment.getFinalNodeID()) && nodeId >= 0 && marginalDistributions[nodeId] != null)
				base = (char) marginalDistributions[nodeId].getMax();
			ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), base, transitionIds));
		}

		// save sequence information in internal nodes of the phylogenetic tree
		ancestralSeqLabels = new ArrayList<>();
		ancestralSeqLabels.add(marginalNode);
//		long elapsedTimeNs = System.nanoTime() - startTime;x
//		System.out.printf("Elapsed time in secs: %5.3f\n", elapsedTimeNs/1000000000.0);
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
			d_marg = null;//EnumDistrib.uniform(Enumerable.aacid);
        }
        return d_marg;
    }

	/**
	 * Get the pairs of parent / child nodes that differ in insertion / deletion content
	 *
	 * @return				marginal distribution of query node using ve
	 */
	public HashMap<String, List> getIndelDifferences(int length){

		HashMap<String, List> pairs = new HashMap<>();
		Map<String, String> ancestralDict = this.getAncestralDict();

		System.out.println(ancestralDict);



		for (String ancestor: this.getAncestralSeqLabels()){
			String inference = ancestralDict.get(ancestor);

			Pattern pattern = Pattern.compile("\\-{" + length + ",}");
			Matcher matcher = pattern.matcher((inference));
			while (matcher.find()) {
				String parent = this.getParent(ancestor).getLabel().toString();

				System.out.println(matcher.start() + " " + matcher.end());



				String parent_inference = ancestralDict.get(parent);

				String parent_substring = parent_inference.substring(matcher.start(), matcher.end());

				Map<Integer, Integer> position_mapping = new HashMap<Integer, Integer>();

				int count = 0;
				for (int actual = matcher.start(); actual <= matcher.end(); actual ++){

					position_mapping.put(count, actual);
					count ++;


				}

				System.out.println(position_mapping);
				System.out.println("Found inference string is");
				System.out.println(inference);

				System.out.println("Parent inference string is");
				System.out.println(parent_inference);

				System.out.println("Exact locations are");
				System.out.println(inference.substring(matcher.start(), matcher.end()));
				System.out.println(parent_substring);



				Pattern parent_pattern = Pattern.compile("\\w{" + length + ",}");
				Matcher parent_matcher = parent_pattern.matcher(parent_substring);
				while (parent_matcher.find()) {


					System.out.println(parent_matcher.start() + " " + parent_matcher.end());

					System.out.println(parent_substring.substring(parent_matcher.start(), parent_matcher.end()));

					String pair_key = parent + ":" + ancestor;

					if (pairs.containsKey(pair_key)){
						pairs.get(pair_key).add(position_mapping.get(parent_matcher.start()) + ":" + position_mapping.get(parent_matcher.end()));
					}
					else {
						ArrayList<String> pair_list = new ArrayList<String>();
						pair_list.add(position_mapping.get(parent_matcher.start()) + ":" + position_mapping.get(parent_matcher.end()));
						pairs.put(pair_key, pair_list);


					}


				}



			}

		}
		return pairs;
	}

	public JSONArray getIndelDifferencesJSON(int length){

		HashMap<String, List> indelDifferences = this.getIndelDifferences(length);

		JSONArray pairs = new JSONArray();

		for (Map.Entry<String, List> list : indelDifferences.entrySet()) {

			JSONObject pair = new JSONObject();
			pair.put("label", list.getKey().toString());
			pair.put("value", list.getValue());
			pairs.put(pair);
		}
		return pairs;

	}


    /**
     * Helper class to store changes to an ancestral graph node
     * 
     * Information:
     * 		- POG structure index 
     * 		- Inferred base character: base character that is inferred or '-' to represent a gap (i.e. that the node needs to be deleted when updating the structure)
     */
    public class Inference {
    	Integer pogId;
    	char base;
    	List<Integer> transitions;

    	public Inference(Integer id, char ch, List<Integer> tr){
    		pogId = id;
    		base = ch;
    		transitions = tr;
    	}
    	public String toString(){
    		return pogId + "->" + base;
    	}
    }
}