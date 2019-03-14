package api;

import dat.EnumSeq;
import dat.Enumerable;
import dat.POGraph;
import dat.POGraph.*;
import dat.PhyloTree;
import dat.EnumSeq.*;
import dat.file.AlnReader;
import dat.file.FastaReader;

import java.io.IOException;
import java.util.List;

/**
 * Created by mikael on 24/5/18.
 */
public class BasicExample {

    static String alnfile = "/Users/mikael/simhome/ASR/basic_example.aln";
    static String nwkfile = "/Users/mikael/simhome/ASR/basic_example.nwk";

    public static void main(String[] args) {
        List<Gappy<Enumerable>> aln = null;
        PhyloTree nwk = null;
        try {
            aln = EnumSeq.Gappy.loadClustal(alnfile, Enumerable.aacid_ext);
            System.out.println("Successfully loaded " + aln.size() + " aligned sequences");
            POGraph pog = new POGraph(aln);
            System.out.println("Successfully converted the alignment into a partial-order graph with " + pog.getNumNodes() + " positions (including dummy start and end states)");
            nwk = new PhyloTree().loadNewick(nwkfile);
            System.out.println("Successfully loaded " + nwk.toNodesBreadthFirst().length + " nodes arranged as a tree");
        } catch (IOException e) {
            System.err.println("Failed loading something");
            System.exit(1);
        }
    }
}
