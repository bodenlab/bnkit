/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bn.reconstruction;

import bn.prob.EnumDistrib;
import java.util.Arrays;

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
        }
        
        ASR asr = new ASR(args[0], args[1]);
        String result = asr.getAsrSeq();
        EnumDistrib[] margDistr = asr.getMarginDistribs();
        double[] rate = asr.getRates();
        System.out.println("Ancestral sequence");
        System.out.println(result);
        System.out.println("Raw Marginal Distributions per position");
        System.out.println(Arrays.asList(margDistr));
        System.out.println("Single MargDistrib Column");
        System.out.println(margDistr[0]);
        System.out.println("Full Rate list");
        System.out.println(Arrays.asList(rate));
        System.out.println("Single rate column");
        System.out.println(rate[0]);
    }
}
