package reconstruction;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 *  Reconstruct ancestral sequences using information stored in a partial order alignment graph. Output is given as a partial order
 *  alignment structure that illustrates possible indel positions and inferred mutations in the ancestral sequence.
 *
 *  @author Marnie
 *
 */
public class RunASRPOG {

	/**
	 *
	 *  @param args
	 *
	 *  1.  string representation (i.e. from POGraph.toString()) or filepath of partial order alignment graph structure (filetype .dot).
	 *  	If this is not specified, then a sequence fasta file (-s flag) must be specified. If the specified sequence file is in
	 *  	.aln format, the reconstruction will use this as an alignment, otherwise the sequences will be aligned using a partial
	 *  	order alignment graph for the reconstruction.
	 *
	 *  2.	filepath of phylogenetic tree (filetype: .nwk)
	 *
	 *  Optional flags:
	 *  			-s		filepath of sequences (filetype: .aln, .fa or .fasta)
	 *  			-o 		output filepath to save reconstructed partial order graphs of internal node of the phylogenetic tree
	 * 				-inf	inference type, 'marginal' or 'joint'. Default: joint. If 'marginal', specify the node name after 'marginal', default: root node.
	 * 			    -p      number of threads
	 * 				-mp		use maximum parsimony to infer gaps positions. By default, maximum likelihood is used within a Bayesian network framework.
	 * 				-msa	generate dot file in output directory representing multiple sequence alignment of input sequences or partial order alignment graph. Default: no msa dot file is generated
	 * 				-dot	generate dot file in output directory ancestral node sequence
	 * 				-align	perform sequence alignment prior to reconstruction, otherwise assumes sequences are aligned
	 *
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException, InterruptedException {

		if (args.length > -1) {
			ASRPOG asr = null;

			String inference = "joint";
			String marginalNode = null;
			String outputPath = "";
			String sequencePath = "";
			String treePath = "";
			String poagRepresentation = "";
			String alignmentType = "";
			String perturbNode = "";
			String model = null;

			System.out.println("RUNNING");
			int numThreads = 5;
			int indelLength = 0;
			if (!args[0].contains("-"))
				if (!args[0].contains(".nwk"))
					poagRepresentation = args[0];
				else
					treePath = args[0];
			boolean dotFile = false;
			boolean msaFile = false;
			boolean performAlignment = false;
			boolean checkBranchIsolation = false;
			boolean perturbAncestors = false;
			boolean indelDifferences = false;

			// parse parameters
			for (int arg = 0; arg < args.length; arg++) {
				if (args[arg].endsWith(".fa") || args[arg].endsWith(".fasta") || args[arg].endsWith(".aln"))
					sequencePath = args[arg];
				else if (args[arg].endsWith(".nwk"))
					treePath = args[arg];
				else if (args[arg].endsWith(".dot"))
					poagRepresentation = args[arg];
				else if (args[arg].equalsIgnoreCase("-o"))
					outputPath = args[arg + 1];
				else if (args[arg].equalsIgnoreCase("-help")) {
					usage("GRASP ancestral sequence predictions");
					return;
				} else if (args[arg].equalsIgnoreCase("-p"))
					numThreads = Integer.parseInt(args[arg + 1]);
				else if (args[arg].equalsIgnoreCase("-inf")) {
					inference = args[arg + 1];
					if (inference.equalsIgnoreCase("marginal") && arg + 2 < args.length && !args[arg + 2].startsWith("-"))
						marginalNode = args[arg + 2];
				} else if (args[arg].equalsIgnoreCase("-msa"))
					msaFile = true;
				else if (args[arg].equalsIgnoreCase("-dot"))
					dotFile = true;
				else if (args[arg].equalsIgnoreCase("-model"))
					model = args[arg + 1];
				else if (args[arg].equalsIgnoreCase("-align")) {
					performAlignment = true;
					alignmentType = args[arg + 1];
				}
				else if (args[arg].equalsIgnoreCase("-perturb")) {
					perturbAncestors = true;
					perturbNode = args[arg + 1];
				}

				else if (args[arg].equalsIgnoreCase("-indeldiff")) {
					indelDifferences = true;
					indelLength = Integer.parseInt(args[arg + 1]);

				}


			}

			// exit if the phylogenetic tree has not been specified, or the partial order alignment structure and sequence filepath both have not been specified (need one or the other)
			if (treePath.isEmpty())
				usage("Filepath to the phylogenetic tree must be provided as an input parameter.");
			if (poagRepresentation.isEmpty() && sequencePath.isEmpty())
				usage("A partial order alignment graph structure or filepath must be input as a parameter, or a sequence FASTA filepath must be specified using the [-s] parameter.");

			if (poagRepresentation.isEmpty()) {
				if (performAlignment) { // generate a partial order alignment graph if the alignment has not been specified
					asr = new ASRPOG(sequencePath, treePath, inference.equalsIgnoreCase("joint"), true, model, numThreads);
				} else if (marginalNode != null)
					asr = new ASRPOG(sequencePath, treePath, sequencePath, marginalNode, false, model, numThreads);
				else
					asr = new ASRPOG(sequencePath, treePath, inference.equalsIgnoreCase("joint"), false, model, numThreads);
			} else if (marginalNode != null)
				asr = new ASRPOG(poagRepresentation, treePath, sequencePath, marginalNode, performAlignment, model, numThreads);
			else
				asr = new ASRPOG(poagRepresentation, treePath, inference.equalsIgnoreCase("joint"), performAlignment, model, numThreads);

			if (!outputPath.isEmpty()) {
//				outputPath += "/";
				if (dotFile)
					asr.saveGraph(outputPath);
				if (msaFile)
					asr.saveMSAGraph(outputPath);
				//asr.saveSupportedAncestors(outputPath);
				//asr.saveGraph(outputPath);
				// FIXME: next call not saving right distribs
				asr.saveDistrib(outputPath  + marginalNode);
				if (inference.equalsIgnoreCase("joint"))
					asr.save(outputPath, true, "fasta");
				else
					asr.save(outputPath, false, "fasta");
			}

			if (indelDifferences){
				System.out.println("Get the indel differences");
				HashMap<String, List> pairs = (asr.getIndelDifferences(indelLength));
				System.out.println("Here are the pairs");
				System.out.println(pairs);
			}

//			if (perturbAncestors){
//
//				ArrayList<EnumSeq.Gappy<Enumerable>>  allSeqs = asr.getInferenceFromPerturbation(outputPath, perturbNode, sequencePath, treePath, inference, performAlignment, model, numThreads);
//				asr.outputALN(allSeqs, outputPath + "_all_ancestors", "fasta");
//
//
//			}


			if (checkBranchIsolation){


				ArrayList<String> ancestralNodes = new ArrayList<String>() {{
					add("N1");
					add("N3");
					add("N4");
					add("N5");

				}};



				Map<String, String> ancestralLabels = asr.getAncestralDict();


				// Code for using a specific list of Nodes
				for (String node : ancestralNodes){
					BranchIsolation branchIsolation = new BranchIsolation(asr, ancestralLabels, treePath, sequencePath, node, performAlignment, model, numThreads);


				}

				// Code for using one specific node
//				BranchIsolation branchIsolation = new BranchIsolation(asr, ancestralLabels, treePath, sequencePath, "N5", mp);


				// Code for using all nodes
//				for (String node : ancestralLabels.keySet()){
//					BranchIsolation branchIsolation = new BranchIsolation(asr, ancestralLabels, treePath, sequencePath, node, mp);
//				}


			}

		} else {
			usage("");
		}



	}

	private static void exit(String message){
		System.err.println(message + "\n");
		System.exit(1);
	}

	private static void usage(String message){
		if (!message.isEmpty())
			System.out.println(message + "\n");
		System.out.println("Usage: <poag_structure/poag_file.dot> <tree_file.nwk>");
		System.out.println("Usage: <tree_file.nwk> <sequence_file.fasta/.aln/.fa> [-o/-p/-inf/-dot/-msa/-align/-model]\n");
		System.out.println("Optional flags:");
		System.out.println("\t-o 		output filepath to save reconstruction");
		System.out.println("\t-p		number of threads to use. Default: 1");
		System.out.println("\t-inf		inference type, 'marginal' or 'joint'. Default: Joint");
		System.out.println("\t-msa		generate dot file in output directory representing multiple sequence alignment of input sequences or partial order alignment graph. Default: no msa dot file is generated");
		System.out.println("\t-dot		generate dot file in output directory representing ancestral node sequence. Default: no dot file is generated");
		System.out.println("\t-align	\tperform sequence alignment prior to reconstruction. Assumes sequences are aligned if this flag is not specified");
		System.out.println("\t-model	\tevolutionary model to use for reconstruction (JTT, Dayhoff, LG or WAG). Default: JTT");
		System.exit(1);
	}
}
