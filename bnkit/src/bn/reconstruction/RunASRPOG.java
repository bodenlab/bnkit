package bn.reconstruction;

import java.io.IOException;


/**
 *  Reconstruct ancestral sequences using information stored in a partial order alignment graph. Output is given as a partial order
 *  alignment structure that illustrates possible indel positions and inferred mutations in the ancestral sequence.
 *
 *  @author marnie
 *
 */
public class RunASRPOG {

	/**
	 *
	 *  @param args
	 *
	 *  1.  string representation (i.e. from POGraph.toString()) or filepath of partial order alignment graph structure (filetype .txt).
	 *  	If this is not specified, then a sequence fasta file (-s flag) must be specified. If the specified sequence file is in
	 *  	.aln format, the reconstruction will use this as an alignment, otherwise the sequences will be aligned using a partial
	 *  	order alignment graph for the reconstruction.
	 *
	 *  Required flags:
	 *  			-t		filepath of phylogenetic tree (filetype: .nwk)
	 *  Optional flags:
	 *  			-s		filepath of sequences (filetype: .aln, .fa or .fasta)
	 *  			-o 		output filepath to save reconstructed partial order graphs of internal node of the phylogenetic tree
	 * 				-p		inference type, 'marginal' or 'joint'. Default: joint
	 * 				-msa	generate dot file in output directory representing multiple sequence alignment of input sequences or partial order alignment graph. Default: no msa dot file is generated
	 * 				-dot	generate dot file in output directory ancestral node sequence
	 *
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		if (args.length > 1) {
			ASRPOG asr = null;
			
			String inference = "joint";
			String outputPath = "";
			String sequencePath = "";
			String treePath = "";
			String poagRepresentation = "";
			if (!args[0].contains("-"))
				poagRepresentation = args[0];
			boolean dotFile = false;
			boolean msaFile = false;
			
			// parse parameters
			for (int arg = 0; arg < args.length; arg++) {
				if (args[arg].equals("-s"))
					sequencePath = args[arg + 1];
				else if (args[arg].equals("-o"))
					outputPath = args[arg + 1];
				else if (args[arg].equals("-t"))
					treePath = args[arg + 1];
				else if (args[arg].equals("-p"))
					inference = args[arg + 1];
				msaFile = !msaFile && args[arg].equals("-msa");
				dotFile = !dotFile && args[arg].equals("-dot");
			}
			
			// exit if the phylogenetic tree has not been specified, or the partial order alignment structure and sequence filepath both have not been specified (need one or the other)
			if (treePath.isEmpty())
				exit("Filepath to the phylogenetic tree must be specified using the [-t] parameter.");
			if (poagRepresentation.isEmpty() && sequencePath.isEmpty())
				exit("A partial order alignment graph structure or filepath must be input as a parameter, or a sequence fasta filepath must be specified using the [-s] parameter.");
			
			// generate a partial order alignment graph if the alignment has not been specified
			if (poagRepresentation.isEmpty()) {
				asr = new ASRPOG(sequencePath, treePath, inference);
			} else
				asr = new ASRPOG(poagRepresentation, treePath, inference);
			
			if (!outputPath.isEmpty()) {
				if (dotFile)
					asr.saveGraph(outputPath);
				if (msaFile)
					asr.saveMSAGraph(outputPath);
				asr.saveSupportedAncestors(outputPath);
				asr.saveGraph(outputPath);
				if (inference.equalsIgnoreCase("joint"))
					asr.save(outputPath, true);
				else
					asr.save(outputPath, false);
			}

		} else if (args.length > 10) {
			exit("Too many arguments specified.");
		} else {
			exit("");
		}

	}

	private static void exit(String message){
		if (!message.isEmpty())
			System.out.println(message);
		System.out.println("Usage: <poag_structure/poag_file.txt> -t <tree_file.nwk>");
		System.out.println("Usage: -t <tree_file.nwk> -s <sequence_file.fasta/.aln/.fa> [-o/-p/-dot/-msa]");
		System.exit(1);
	}
}
