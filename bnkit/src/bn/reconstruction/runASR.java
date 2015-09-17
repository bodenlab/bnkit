/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bn.reconstruction;


/**
 *
 * @author Alex
 */
public class runASR {
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        if(args.length < 2) {
            System.out.println("Usage: <tree_file> <aln_file>");
        } else if (args.length == 2) {
            ASR asr = new ASR(args[0], args[1]);
            asr.save("JSONoutput.txt");
        } else {
            ASR asr = new ASR(args[0], args[1]);
            asr.save(args[2]);
        }
    }
}

//"../test_tree.txt" "../test_aln.aln"
//"../test_tree.txt" "../test_aln_jitter.txt"
//"../CYP2F.nwk" "../CYP2F.aln"
//"../soxb_tree.nwk" "../soxb_clustal.aln"
//"../voordeckers_tree.txt" "../voordeckers_aln.txt"