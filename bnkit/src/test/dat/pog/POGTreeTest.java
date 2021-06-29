package dat.pog;

import asr.Parsimony;
import asr.Prediction;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

class POGTreeTest {

    static EnumSeq.Alignment aln1 = null;
    static EnumSeq.Alignment aln2 = null;
    static POGTree pogt1 = null;
    static POGTree pogt2 = null;
    static Tree tree1 = null;
    static Tree tree2 = null;

    @BeforeAll
    static void setPogt1() {
        try {
            aln1 = new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("bnkit/src/test/resources/default.aln", Enumerable.aacid));
            aln2 = new EnumSeq.Alignment(EnumSeq.Gappy.loadFasta("bnkit/src/test/resources/simmons1.aln", Enumerable.aacid, '-'));
            tree1 = Tree.load("bnkit/src/test/resources/default.nwk", "newick");
            tree2 = Tree.load("bnkit/src/test/resources/simmons1.nwk", "newick");
            pogt1 = new POGTree(aln1, tree1);
            pogt2 = new POGTree(aln2, tree2);
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    @Test
    void getNodeInstance() {
        for (int i = 0; i < aln1.getWidth(); i ++) {
            TreeInstance ti = pogt1.getNodeInstance(i);
            int count = 0;
            for (int j = 0; j < ti.getTree().getSize(); j ++)
                count += ti.getInstance(j) == null ? 0 : 1;
            assertEquals(aln1.getOccupancy(i), count);
        }
    }

    @Test
    void parsimonyAncestors2() {
        Prediction ap = Prediction.PredictBySICP(pogt2);
        for (int idx : pogt2.getTree()) {
            if (!pogt2.getTree().isLeaf(idx)) { // ancestor
                Object ancID = pogt2.getTree().getLabel(idx);
                POGraph ancestor = ap.getAncestor(ancID);
                try {
                    ancestor.saveToDOT("/Users/mikael/Downloads/pog" + ancID + ".dot");
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    @Test
    void parsimonyAncestors1() {
        EnumSeq.Gappy[] aln = getParsimonyGappedAlignment(pogt1);
        Prediction ap = Prediction.PredictByParsimony(pogt1);
        for (int idx : pogt1.getTree()) {
            if (pogt1.getTree().getChildren(idx).length > 0) { // ancestor
                Object ancID = pogt1.getTree().getLabel(idx);
                POGraph ancestor = ap.getAncestor(ancID);
                Random r = new Random(System.currentTimeMillis());
                for (EnumSeq.Gappy seq : aln) {
                    if (seq.getName().equals(ancID.toString())) {
                        int ptr = -1;
                        for (int i = 0; i < seq.length(); i++) {
                            Object sym = seq.get(i);
                            if (sym.equals('C') || (sym.equals('?') && r.nextBoolean())) {
                                int[] next = ancestor.getForward(ptr);
                                boolean found = false;
                                for (int m : next) {
                                    if (m == i)
                                        found = true;
                                }
                                if (!found)
                                    System.out.println(seq.getName() + ": " + seq + "\t@\t" + ptr + "\tto\t" + i);
                                assertTrue(found);
                                ptr = i;
                            }
                        }
                        System.out.println(seq.getName() + ": " + seq + "\tis fine");
                    }
                }
            }
        }
    }

    public EnumSeq.Gappy[] getParsimonyGappedAlignment(POGTree pogTree) {
        int nPos = pogTree.getPositions(); //
        IdxTree tree = pogTree.getTree();
        EnumSeq.Gappy[] aln = new EnumSeq.Gappy[tree.getSize()];
        TreeInstance[] ti = new TreeInstance[nPos];
        Parsimony[] pi = new Parsimony[nPos];
        for (int i = 0; i < nPos; i ++)
            ti[i] = pogTree.getNodeInstance(i, true);
        for (int i = 0; i < nPos; i ++) {
            Parsimony p = new Parsimony(ti[i].getTree());
            //p.SET_RANDOM_PARSIMONY = true; // default is false
            //p.SET_ONE_TARGET_PARSIMONY = true; // default is false
            p.infer(ti[i], false);
            pi[i] = p;
        }
        for (int j = 0; j < tree.getSize(); j ++) {
            EnumSeq.Gappy seq = new EnumSeq.Gappy(Enumerable.gap_ext);
            Object[] arr = new Object[nPos];
            for (int i = 0; i < nPos; i++) {
                if (pi[i].getOptimal(j).contains(Boolean.FALSE) && pi[i].getOptimal(j).contains(Boolean.TRUE))  // admissible
                    arr[i] = '?';
                else if (pi[i].getOptimal(j).contains(Boolean.FALSE))
                    arr[i] = 'C';
                else
                    arr[i] = 'G';
                System.err.print(Enumerable.gap_ext.isValid(arr[i]) ? "" : arr[i]);
            }
            seq.set(arr);
            seq.setName(tree.getBranchPoint(j).getID().toString());
            aln[j] = seq;
        }
        return aln;
    }

}