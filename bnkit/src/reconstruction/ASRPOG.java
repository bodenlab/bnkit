package reconstruction;


import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.PhyloBNet;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.gap;
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
	private double[] rates = null; 										// Rates at positions in alignment

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile		filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 */
	public ASRPOG(String alignmentFile, String treeFile, boolean jointInference) {
		performASR("", treeFile, alignmentFile, jointInference, false);
	}

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile		filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param maxParsimonyGaps	flag for indicating gap inference using maximum parsimony (true: maximum parsimony or false: maximum likelihood)
	 */
	public ASRPOG(String alignmentFile, String treeFile, boolean jointInference, boolean maxParsimonyGaps) {
		performASR("", treeFile, alignmentFile, jointInference, maxParsimonyGaps);
	}

	/**
	 * Infer ancestral sequences given an alignment file (fasta or aln).
	 *
	 * @param alignmentFile	filepath to the sequence alignment (expected extension .fa, .fasta or .aln)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param marginalNode	node label for maginal inference
	 */
	public ASRPOG(String alignmentFile, String treeFile, String marginalNode) {
		this.marginalNode = marginalNode;
		performASR("", treeFile, alignmentFile, false, false);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, boolean jointInference) {
		performASR(pog, treeFile, sequenceFile, jointInference, false);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog				POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile			filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile		filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param jointInference	flag for indicating joint inference (true: 'joint' or false: 'marginal')
	 * @param maxParsimonyGaps	flag for indicating gap inference using maximum parsimony (true: maximum parsimony or false: maximum likelihood)
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean maxParsimonyGaps) {
		performASR(pog, treeFile, sequenceFile, jointInference, maxParsimonyGaps);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog			POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile	filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param marginalNode	node label for maginal inference
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, String marginalNode) {
		this.marginalNode = marginalNode;
		performASR(pog, treeFile, sequenceFile, false, false);
	}

	/**
	 * Construct phylogenetic tree structure, load partial order alignment graph, and infer ancestral sequences based on inference type
	 *
	 * @param pog			POG dot string or filepath to the partial order alignment graph (expected extension .dot)
	 * @param treeFile		filepath to the phylogenetic tree (expected extension .nwk)
	 * @param sequenceFile	filepath to the sequences (expected extension .fa, .fasta or .aln)
	 * @param marginalNode	node label for maginal inference
	 * @param maxParsimonyGaps	flag for indicating gap inference using maximum parsimony (true: maximum parsimony or false: maximum likelihood)
	 */
	public ASRPOG(String pog, String treeFile, String sequenceFile, String marginalNode, boolean maxParsimonyGaps) {
		this.marginalNode = marginalNode;
		performASR(pog, treeFile, sequenceFile, false, maxParsimonyGaps);
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
	public void saveSupportedAncestors(String filepath){
		Map<String, String> ancestralSeqs = new HashMap<>();
		if (marginalNode == null)
			for (String phyloNodeLabel : ancestralSeqLabels)
				ancestralSeqs.put(phyloNodeLabel, phyloTree.find(phyloNodeLabel).getSequence().toString());
		else
			ancestralSeqs.put(marginalNode, phyloTree.find(marginalNode).getSequence().toString());
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(filepath + "_recon.fa", false));
			for (String node : ancestralSeqs.keySet()) {
				bw.write(">" + node);
				bw.newLine();
				bw.write(ancestralSeqs.get(node));
				bw.newLine();
				bw.newLine();
			}
			bw.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
	}

	/**
	 * Save ASR information; ALN, tree, rates, marginal distribution (if applicable).
	 *
	 * @param filepath		filename to save ASR
	 * @param infSequence   sequences inferred (true) or marginal (false)
	 * @param format		format to save alignment to, "fasta" or "clustal"
	 */
	public void save(String filepath, boolean infSequence, String format) {
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
	public void saveALN(String filepath, String format) {
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
		try {
			if (format.equalsIgnoreCase("clustal")) {
				AlnWriter aw = new AlnWriter(filepath + ".aln");
				aw.save(allSeqs);
				aw.close();
			} else {
				FastaWriter fw = new FastaWriter(filepath + ".fa");
				fw.save(allSeqs);
				fw.close();
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	/**
	 * Save rates
	 *
	 * @param filename filepath to save rates
	 */
	public void saveRate(String filename) {
		if (rates == null)
			rates = new double[]{};
		try {
			Writer writer = new PrintWriter(filename, "UTF-8");
			for (Integer nodeId : pogAlignment.getNodeIDs()) {
				pogAlignment.setCurrent(nodeId);
				writer.write(pogAlignment.getCurrentId() + ":");
				for (Character base : pogAlignment.getCurrentBases())
					writer.write(base + ",");
				writer.write(" " + Double.toString(rates[nodeId]) + "\n");
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

	/**
	 * Save marginal distribution
	 *
	 * @param filename	filepath to save distribution
	 */
	public void saveDistrib(String filename){
		try {
			Writer writer = new PrintWriter(filename, "UTF-8");
			Object[] aacid = Enumerable.aacid_ext.getValues();
			Object[][] margMatrix = new Object[aacid.length][marginalDistributions.length + 1];
			writer.write("columns\t");
			for (int i = 1; i < marginalDistributions.length; i++) { //write header
				if (i == marginalDistributions.length - 1)
					writer.write(i + "\n");
				else
					writer.write(i + "\t");
			}
			for (int j = 0; j < aacid.length; j++) { //fill in row names
				margMatrix[j][0] = aacid[j];
			}
			for (int k = 1; k < marginalDistributions.length + 1; k++) {
				EnumDistrib distr = marginalDistributions[k-1];
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

	/**
	 * Save the phylogenetic tree with the ancestral nodes entered.
	 *
	 * @param filepath	filepath to save phylogenetic tree (.nwk)
	 */
	public void saveTree(String filepath) {
		try {
			Writer writer = new PrintWriter(filepath, "UTF-8");
			String newick = phyloTree.getRoot().toString();
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
	 * @param parsimony			flag to identify gaps in reconstruction using parsimony (true) or maximum likelihood (false)
	 */
	private void performASR(String pog, String treeFile, String sequenceFile, boolean jointInference, boolean parsimony) throws RuntimeException {
		loadData(treeFile, sequenceFile);
		if (pog.equals(""))	// load graph structure from alignment file
			pog = sequenceFile;
		pogAlignment = new POGraph(pog, sequenceFile);
		pog = null;

		// perform inference
		if (jointInference) {
			marginalDistributions = null;
			queryBNJoint(parsimony);
		} else if (marginalNode != null && phyloTree.find(marginalNode) != null) {
			queryBNMarginal(marginalNode, parsimony);
		} else {
			if (marginalNode == null)
				System.out.println("No node was specified for the marginal inference: inferring the root node");
			else
				throw new RuntimeException("Incorrect internal node label provided for marginal reconstruction: " + marginalNode + " tree: " + phyloTree.toString());
			marginalNode = phyloTree.getRoot().getLabel().toString();
			queryBNMarginal(phyloTree.getRoot().getLabel().toString(), parsimony);
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
		seq.set(chars);
		phyloTree.find(phyloNodeLabel).setSequence(seq);
		System.out.println(phyloNodeLabel + ": sequence " + seq.toString());
	}

	/**
	 * Generates the ancestor by applying the inference changes stored in the ancestralInferences log.
	 *
	 * @param label	Ancestral sequence label
	 * @return		POGStructure of ancestor
	 */
	private POGraph getAncestor(String label) {
		POGraph ancestor = new POGraph(pogAlignment);
		for (Inference inferredBase : ancestralInferences.get(label))
			if (ancestor.setCurrent(inferredBase.pogId))
				if (inferredBase.base == '-')
					ancestor.removeNode();
				else
					ancestor.setBase(inferredBase.base);

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
	private void loadData(String treeFile, String sequenceFile) {
		try {
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
			if (!eNodes.equals(seqNames))
				throw new RuntimeException("The sequence names in the provided alignment must all have a match" +
						" in the provided tree");

			// save sequence information in internal nodes of the phylogenetic tree
			for (EnumSeq.Gappy<Enumerable> extant : extantSequences)
				phyloTree.find(extant.getName()).setSequence(extant);

		} catch (Exception ex) {
			ex.printStackTrace();
		}
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


	/**
	 * Create gap/character Bayesian network for node in the partial order alignment graph. 
	 * 
	 * @return	Bayesian networks for node position in pogAlignment
	 */
	private PhyloBNet createGapNetwork(){
			PhyloBNet phyloBN = PhyloBNet.createGap(phyloTree, new gap());
			Map<Integer, Character> sequenceCharacterMapping = pogAlignment.getSequenceCharacterMapping();
			
			// for all extant sequences, if sequence is in this alignment, find where the location is in the phylogenetic tree and assign character,
			// otherwise assign gap marker
			for (int extantSeq = 0; extantSeq < extantSequences.size(); extantSeq++)
			    if (sequenceCharacterMapping.containsKey(extantSeq))
					// find where the sequence is in the BN and set the base character
					phyloBN.getBN().getNode(extantSequences.get(extantSeq).getName()).setInstance('C');
				else
					// sequence is not part of character inference, set to gap
					phyloBN.getBN().getNode(extantSequences.get(extantSeq).getName()).setInstance('G');
		
			return phyloBN;
		}
	
	/**
	 * Infer gap/base character of each partial order alignment graph structure at each internal node of the phylogenetic tree using joint inference.
	 *
	 * @param gapParsimony	flag to identify gaps using parsimony (true) or maximum likelihood (false)
	 */
	private void queryBNJoint(boolean gapParsimony){
		rates = new double[pogAlignment.getNumNodes()]; //Rate matrix

		// infer base/gap of each aligned node 
		List<Integer> nodeIDs = pogAlignment.getNodeIDs();
		for (Integer nodeId : nodeIDs) {
			pogAlignment.setCurrent(nodeId);

			if (pogAlignment.getCurrentBase() != null && pogAlignment.getSequenceCharacterMapping().size() == extantSequences.size())
				continue;

			Variable.Assignment[] charAssignments = null;
			VarElim ve = new VarElim();

			// get gap or character inference at phylogenetic nodes
			Map<String, Character> phyloGaps = new HashMap<>();
			if (gapParsimony) {
				// construct string[] of extant names and object[] of C/G at the appropriate leaf nodes
				// i.e. all sequences in the current node will have a character, every other extant will be a gap
				String[] extants = new String[extantSequences.size()];
				Character[] bases = new Character[extantSequences.size()];
				Map<Integer, Character> nodeSeqs = pogAlignment.getSequenceCharacterMapping();
				for (int seqId = 0; seqId < extantSequences.size(); seqId++) {
					extants[seqId] = extantSequences.get(seqId).getName();
					if (nodeSeqs.containsKey(seqId))
						bases[seqId] = 'C';
					else
						bases[seqId] = 'G';
				}
				phyloTree.setContentByParsimony(extants, bases);
				// get assigned internal bases
				for (String phyloNode : ancestralSeqLabels)
					phyloGaps.put(phyloNode, (char) phyloTree.find(phyloNode).getValue());
			} else {
				ve.instantiate(createGapNetwork().getBN());
				Variable.Assignment[] gapAssignments = getJointAssignment(ve);
				for (Variable.Assignment varassign : gapAssignments)
					phyloGaps.put(varassign.var.getName(), (char) varassign.val);
			}
			if (phyloGaps.values().contains('C')) {
				PhyloBNet charNet = createCharacterNetwork();
				rates[nodeId] = charNet.getRate();
				ve.instantiate(charNet.getBN());
				charAssignments = getJointAssignment(ve);
			}


			// for each node in the phylogenetic tree, if character is inferred at position, set inferred base, otherwise gap is inferred, remove from alignment
			for (String phyloNode : phyloGaps.keySet()) {
				Character base = null;
				if (phyloGaps.get(phyloNode) == 'C') {
					// infer base character
					for (Variable.Assignment varassign : charAssignments)
						if (phyloNode.equals(varassign.var.getName())) {
							base = (char) varassign.val;
							break;
						}
					// base will be null (but not set as a 'gap') if an ancestor needed to be pruned, and this node was
					// subsequently pruned with it
					if (base != null) {
						if (!ancestralInferences.containsKey(phyloNode))
							ancestralInferences.put(phyloNode, new ArrayList<>());
						ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), base));

					}
				}

				if (base == null){
					// gap, remove from alignment
					if (!ancestralInferences.containsKey(phyloNode))
						ancestralInferences.put(phyloNode, new ArrayList<>());
					ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), '-'));
				}
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
	 * @param gapParsimony	flag to identify gaps using parsimony (true) or maximum likelihood (false)
	 */
	private void queryBNMarginal(String phyloNode, boolean gapParsimony) {
		marginalDistributions = new EnumDistrib[pogAlignment.getNumNodes()];

		// infer base/gap of each aligned node
		List<Integer> nodeIDs = pogAlignment.getNodeIDs();
		for (Integer nodeId : nodeIDs) {
			pogAlignment.setCurrent(nodeId);

			VarElim ve = new VarElim();

			// check if gap or character inference at phylogenetic node
			Character base = null;
			if (gapParsimony) {
				// construct string[] of extant names and object[] of C/G at the appropriate leaf nodes
				// i.e. all sequences in the current node will have a character, every other extant will be a gap
				String[] extants = new String[extantSequences.size()];
				Character[] bases = new Character[extantSequences.size()];
				Map<Integer, Character> nodeSeqs = pogAlignment.getSequenceCharacterMapping();
				for (int seqId = 0; seqId < extantSequences.size(); seqId++) {
					extants[seqId] = extantSequences.get(seqId).getName();
					if (nodeSeqs.containsKey(seqId))
						bases[seqId] = 'C';
					else
						bases[seqId] = 'G';
				}
				phyloTree.setContentByParsimony(extants, bases);
				// get assigned internal bases
				base = (Character) phyloTree.find(phyloNode).getValue();
			} else {
				PhyloBNet phyloNet = createGapNetwork();
				ve.instantiate(phyloNet.getBN());
				for (EnumVariable pNode : phyloNet.getInternal())
					if (pNode.getName().equalsIgnoreCase(phyloNode)) {
						EnumDistrib gapAssignments = getMarginalDistrib(ve, pNode);
						base = (Character) gapAssignments.getMax();
						break;
					}
			}

			// if character is inferred at position, set inferred base, otherwise gap is inferred, remove from alignment
			if (base == 'C') {
				// infer base character
				base = null;

				PhyloBNet phyloNet = createCharacterNetwork();
				ve.instantiate(phyloNet.getBN());
				// get character variable representing current phylogenetic tree node using marginal probability
				EnumVariable charNode = phyloNet.getInternal().get(0);
				for (int nodeIndex = 1; !charNode.getName().equals(phyloNode); nodeIndex++)
					charNode = phyloNet.getInternal().get(nodeIndex);
				marginalDistributions[nodeId] = getMarginalDistrib(ve, charNode);
				base = (char)marginalDistributions[nodeId].getMax();

				if (!ancestralInferences.containsKey(phyloNode))
					ancestralInferences.put(phyloNode, new ArrayList<>());
				ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), base));
			} else {
				// gap, remove from alignment
				double[] empty = new double[Enumerable.aacid.size()];
				for (int d = 0; d < Enumerable.aacid.size(); d++)
					empty[d] = 0.0;
				marginalDistributions[nodeId] = new EnumDistrib(Enumerable.aacid, empty);
				if (!ancestralInferences.containsKey(phyloNode))
					ancestralInferences.put(phyloNode, new ArrayList<>());
				ancestralInferences.get(phyloNode).add(new Inference(pogAlignment.getCurrentId(), '-'));
			}
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
     * Helper class to store changes to an ancestral graph node
     * 
     * Information:
     * 		- POG structure index 
     * 		- Inferred base character: base character that is inferred or '-' to represent a gap (i.e. that the node needs to be deleted when updating the structure)
     */
    private class Inference {
    	int pogId;
    	char base;

    	public Inference(int id, char ch){
    		pogId = id;
    		base = ch;
    	}
    	public String toString(){
    		return pogId + "->" + base;
    	}
    }
}