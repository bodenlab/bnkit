import static org.junit.jupiter.api.Assertions.assertEquals;

import dat.EnumSeq;
import dat.Enumerable;
import dat.PhyloTree;
import java.io.IOException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;


/**
 * Created by mikael on 17/8/17.
 * @deprecated
 */
public class ParsimonyTests {

    static PhyloTree[] trees;
    static EnumSeq.Alignment[] alns;


    PhyloTree tree1 = new PhyloTree().parseNewick("((((x01,x02)X01_02,x03)X01_03,(x04,((x05,x06)X05_06,(x07,x08)X07_08)X05_08)X04_08)X01_08,(x09,(x10,x11)X10_11)X09_11)X01_11;");

    @BeforeAll
    public static void setUp() throws Exception {
        try {
            trees = new PhyloTree[] { // MB: also successfully tried a 150-seq tree with alignment (not in test/resources)
                    new PhyloTree().loadNewick("bnkit/src/test/resources/large.nwk"),
                    new PhyloTree().loadNewick("bnkit/src/test/resources/default.nwk"),
//                    PhyloTree.loadNewick("src/test/resources/edge1.nwk"),
            };
            alns = new EnumSeq.Alignment[] {
                    new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("bnkit/src/test/resources/large.aln", Enumerable.aacid)),
                    new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("bnkit/src/test/resources/default.aln", Enumerable.aacid)),
//                    new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("src/test/resources/edge2.aln", Enumerable.aacid)),
            };
            for (int i = 0; i < trees.length; i ++)
                trees[i].setAlignment(alns[i]);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    @Test
    public void testParsimony1() throws Exception {
        PhyloTree.Node.SET_ONE_TARGET_PARSIMONY = true;
        PhyloTree.Node.SET_RANDOM_PARSIMONY = true;
        for (int i = 0; i < trees.length; i ++) {
//            System.out.println("Tree " + i);
            for (int col = 0; col < alns[i].getWidth(); col++) {
//                System.out.println("Col " + col);
                double same = -1;
                for (int r = 0; r < 20; r++) {
                    trees[i].setContentByParsimony(alns[i].getNames(), alns[i].getColumn(col));
                    double[] scores = trees[i].getRoot().getScores();
                    double best = Double.POSITIVE_INFINITY;
                    if (scores != null) {
                        for (int s = 0; s < scores.length; s++) {
//                            System.out.printf("\t%2.0f", scores[s]);
                            if (scores[s] < best)
                                best = scores[s];
                        }
                    }
                    double score = trees[i].getRoot().getParsimonyScore();
/*                    System.out.printf("\t%3.0f\n", score);
                    if (i == 1)
                        System.out.println("\t" + trees[i].printValues()); */
                    if (same == -1)
                        same = best;
                    assertEquals(score, best); // make sure the independently calculated score is the same as that determined with parsimony
                    assertEquals(same, best); // make sure the score is the same each time for the same col
                }
            }
        }
    }

    @Test
    public void testParsimony2() throws Exception {
        PhyloTree.Node.SET_ONE_TARGET_PARSIMONY = false;
        PhyloTree.Node.SET_RANDOM_PARSIMONY = false;
        String[] names = new String[] {"x01","x02","x03","x04","x05","x06","x07","x08","x09","x10","x11"};
        Object[] s = new Object[] {"A","B","C","D","E","F","G"};
        /* MB: manually checked the three inits below */
        Object[] init1 = new Object[] {s[0],s[0],s[2],s[3],s[1],s[1],s[1],s[1],s[4],s[5],s[6]};
        Object[] init2 = new Object[] {s[0],s[0],s[2],s[3],s[1],s[1],s[1],s[1],s[4],s[1],s[0]};
        Object[] init3 = new Object[] {s[0],s[0],s[2],s[3],s[1],s[1],s[1],s[1],s[4],s[1],s[5]};
        tree1.setContentByParsimony(names, init1);
        PhyloTree.Node n = tree1.find("N2_X01_03");
        if (n != null) {
            for (int i = 0; i < s.length; i ++)
                assertEquals(n.getValues().contains(s[i]), true);
        }
        n = tree1.find("N4_X04_08");
        if (n != null) {
            for (int i = 0; i < s.length; i ++)
                assertEquals(n.getValues().contains(s[i]), true);
        }
//        System.out.println(tree1.printValues());
        tree1.setContentByParsimony(names, init2);
        n = tree1.find("N2_X01_03");
        if (n != null) {
            assertEquals(n.getValues().contains(s[0]), true);
            assertEquals(n.getValues().contains(s[1]), true);
            assertEquals(n.getValues().contains(s[2]), true);
            assertEquals(n.getValues().contains(s[3]), false);
        }
        n = tree1.find("N4_X04_08");
        if (n != null) {
            assertEquals(n.getValues().contains(s[0]), true);
            assertEquals(n.getValues().contains(s[1]), true);
            assertEquals(n.getValues().contains(s[2]), false);
            assertEquals(n.getValues().contains(s[3]), true);
        }
//        System.out.println(tree1.printValues());
        tree1.setContentByParsimony(names, init3);
        n = tree1.find("N2_X01_03");
        if (n != null) {
            assertEquals(n.getValues().contains(s[0]), true);
            assertEquals(n.getValues().contains(s[1]), true);
            assertEquals(n.getValues().contains(s[2]), true);
        }
        n = tree1.find("N4_X04_08");
        if (n != null) {
            assertEquals(n.getValues().contains(s[0]), false);
            assertEquals(n.getValues().contains(s[1]), true);
            assertEquals(n.getValues().contains(s[2]), false);
            assertEquals(n.getValues().contains(s[3]), false);
        }
//        System.out.println(tree1.printValues());
    }

}
