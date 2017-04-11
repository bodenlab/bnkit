/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package poag;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import api.PartialOrderGraph;

/**
 *
 * @author ariane
 */
public class POAG {
    
    PartialOrderGraph poag; 
    PathGen pathGen;

    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String dotPath = "/home/ariane/Documents/bodenlab/bnkit/bnkit/src/test/resources/testPOGraphMSAMed.dot"; //"/home/ariane/NetBeansProjects/POAG/src/poag/testPOGraphLarge.dot"; //new PartialOrderGraph("/home/ariane/Documents/stemformatics/biojs_alignment/data/defaultMSA.dot");
        POAG pg = new POAG(dotPath);
    }
    
    
    public POAG (String dotPath) {
        poag = new PartialOrderGraph(dotPath);
        // Call poagJSON on a poag 
        POAGJson poagJson = new POAGJson(poag);
        System.err.println(poagJson.toJSON());
    }
}
    
    
    
   
    
    