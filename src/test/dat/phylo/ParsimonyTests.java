package dat.phylo;

import asr.Parsimony;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Newick;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * Created by mikael on 17/8/17.
 * Updated for Tree.java 10/6/20.
 */
public class ParsimonyTests {

    static Tree[] trees;
    static EnumSeq.Alignment[] alns;


    Tree tree1 = Newick.parse("((((x01,x02)X01_02,x03)X01_03,(x04,((x05,x06)X05_06,(x07,x08)X07_08)X05_08)X04_08)X01_08,(x09,(x10,x11)X10_11)X09_11)X01_11;");

    @BeforeAll
    public static void setUp() throws Exception {
        try {
            trees = new Tree[] { // MB: also successfully tried a 150-seq tree with alignment (not in test/resources)
                    Tree.load("bnkit/src/test/resources/large.nwk", "newick"),
                    Tree.load("bnkit/src/test/resources/default.nwk", "newick"),
            };
            alns = new EnumSeq.Alignment[] {
                    new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("bnkit/src/test/resources/large.aln", Enumerable.aacid)),
                    new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("bnkit/src/test/resources/default.aln", Enumerable.aacid)),
            };
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    @Test
    public void testParsimony1() throws Exception {
        for (int i = 0; i < trees.length; i ++) {
            for (int col = 0; col < alns[i].getWidth(); col++) {
                Map<Object, Object> assign = new HashMap<>();
                for (int n = 0; n < alns[i].getNames().length; n ++)
                    assign.put(alns[i].getNames()[n], alns[i].getColumn(col)[n]);
                TreeInstance ti = trees[i].getInstance(assign);
                Parsimony tip = new Parsimony(trees[i]);
                tip.SET_ONE_TARGET_PARSIMONY = true;
                tip.SET_RANDOM_PARSIMONY = true;
                double[] scores = tip.forward(ti);
                double best = Double.POSITIVE_INFINITY;
                int bestidx = 0;
                if (scores != null) {
                    for (int s = 0; s < scores.length; s++) {
//                            System.out.printf("\t%2.0f", scores[s]);
                        if (scores[s] < best) {
                            best = scores[s];
                            bestidx = s;
                        }
                    }
                }
/*
                double score = tip.getScore(bestidx);
                    System.out.printf("\t%3.0f\n", score);
                if (i == 1)
                    System.out.println("\t" + trees[i].printValues());
                assertEquals(score, best); // make sure the independently calculated score is the same as that determined with parsimony
 */
            }
        }
    }

    private int getBPIndex(Tree t, BranchPoint bp) {
        int idx = -1;
        for (int i = 0; i < t.getSize(); i++) {
            if (bp == t.getBranchPoint(i)) {
                idx = i;
                break;
            }
        }
        return idx;
    }

    @Test
    public void testParsimony2() throws Exception {
        String[] names = new String[]{"x01", "x02", "x03", "x04", "x05", "x06", "x07", "x08", "x09", "x10", "x11"};
        Object[] s = new Object[]{"A", "B", "C", "D", "E", "F", "G"};
        /* MB: manually checked the three inits below */
        Object[] init1 = new Object[]{s[0], s[0], s[2], s[3], s[1], s[1], s[1], s[1], s[4], s[5], s[6]};
        Object[] init2 = new Object[]{s[0], s[0], s[2], s[3], s[1], s[1], s[1], s[1], s[4], s[1], s[0]};
        Object[] init3 = new Object[]{s[0], s[0], s[2], s[3], s[1], s[1], s[1], s[1], s[4], s[1], s[5]};
        TreeInstance ti = tree1.getInstance(names, init1);
        Parsimony tip = new Parsimony(ti.getTree());
        tip.SET_ONE_TARGET_PARSIMONY = false;
        tip.SET_RANDOM_PARSIMONY = false;
        tip.infer(ti, false);
        BranchPoint n = tree1.find("X01_03");
        int bpidx = getBPIndex(tree1, n);
        if (bpidx >= 0 && n != null) {
            for (int i = 0; i < s.length; i++)
                assertEquals(tip.getOptimal(bpidx).contains(s[i]), true);
        }
        n = tree1.find("X04_08");
        bpidx = getBPIndex(tree1, n);
        if (bpidx >= 0 && n != null) {
            for (int i = 0; i < s.length; i++)
                assertEquals(tip.getOptimal(bpidx).contains(s[i]), true);
        }
        TreeInstance.LABEL_INCLUDES_INDEX = true;
        System.out.println("init1: " + tip);
        try { Newick.parse(tip.toString()).save("bnkit/src/test/resources/init1.nwk", "nwk"); } catch (IOException e) {}
        ti = tree1.getInstance(names, init2);
        tip = new Parsimony(ti.getTree());
        tip.SET_ONE_TARGET_PARSIMONY = false;
        tip.SET_RANDOM_PARSIMONY = false;
        tip.infer(ti, false);
        n = tree1.find("X01_03");
        bpidx = getBPIndex(tree1, n);
        if (bpidx >= 0 && n != null) {
            assertEquals(tip.getOptimal(bpidx).contains(s[0]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[1]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[2]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[3]), false);
        }
        n = tree1.find("X04_08");
        bpidx = getBPIndex(tree1, n);
        if (bpidx >= 0 && n != null) {
            assertEquals(tip.getOptimal(bpidx).contains(s[0]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[1]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[2]), false);
            assertEquals(tip.getOptimal(bpidx).contains(s[3]), true);
        }
        System.out.println("init2: " + tip);
        try { Newick.parse(tip.toString()).save("bnkit/src/test/resources/init2.nwk", "nwk"); } catch (IOException e) {}
        ti = tree1.getInstance(names, init3);
        tip = new Parsimony(ti.getTree());
        tip.SET_ONE_TARGET_PARSIMONY = false;
        tip.SET_RANDOM_PARSIMONY = false;
        tip.infer(ti, false);
        n = tree1.find("X01_03");
        bpidx = getBPIndex(tree1, n);
        if (bpidx >= 0 && n != null) {
            assertEquals(tip.getOptimal(bpidx).contains(s[0]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[1]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[2]), true);
        }
        n = tree1.find("X04_08");
        bpidx = getBPIndex(tree1, n);
        if (bpidx >= 0 && n != null) {
            assertEquals(tip.getOptimal(bpidx).contains(s[0]), false);
            assertEquals(tip.getOptimal(bpidx).contains(s[1]), true);
            assertEquals(tip.getOptimal(bpidx).contains(s[2]), false);
            assertEquals(tip.getOptimal(bpidx).contains(s[3]), false);
        }
        System.out.println("init3: " + tip);
        try { Newick.parse(tip.toString()).save("bnkit/src/test/resources/init3.nwk", "nwk"); } catch (IOException e) {}

    }

}
