package bn.reconstruction;

/**
 * Created by aesseb on 08-Feb-16.
 */
public class runASR_package {

    /**
     * @param args
     * tree_file - Newick string representation of
     * aln_file - alignment file - clustal or fasta
     * ID - identifier for Reconstruction
     */
    public static void main(String[] args) {
        if(args.length < 4) {
            System.out.println("Usage: <tree_file> <aln_file> <Inference> <ID> ");
            System.out.println("Usage: <tree_file> <aln_file> <Inference> <nodeLabel> <ID> ");
            System.exit(1);
        } else if (args.length == 4) {
            ASR_package asr = new ASR_package(args[0], args[1], args[2]);
            if (args[2].equals("Joint"))
                asr.save(args[3], true);
            else
                asr.save(args[3], false);
//            Analysis test = new Analysis(asr); //the constructor currently handles all steps required
        } else if (args.length == 5) {
            ASR_package asr = new ASR_package(args[0], args[1], args[2], args[3]);
            if (args[2].equals("Joint"))
                asr.save(args[4], true);
            else
                asr.save(args[4], false);
        } else {
            System.out.println("Usage: <tree_file> <aln_file> <Inference> <ID> ");
            System.out.println("Usage: <tree_file> <aln_file> <Inference> <nodeLabel> <ID> ");
            System.exit(1);
        }
    }

}
