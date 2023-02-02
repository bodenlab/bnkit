package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.prob.GaussianDistrib;
import dat.Enumerable;
import dat.file.Newick;
import dat.phylo.*;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

class MaxLhoodJointTest {

    IdxTree createTree() {
        BranchPoint root = new BranchPoint("Root");
        BranchPoint anc_left = new BranchPoint("Anc_left", root, 0.21);
        root.addChild(anc_left);
        BranchPoint anc_right = new BranchPoint("Anc_right", root, 0.17);
        root.addChild(anc_right);
        BranchPoint leaf1 = new BranchPoint("Leaf1", anc_left, 0.22);
        BranchPoint leaf2 = new BranchPoint("Leaf2", anc_left, 0.12);
        anc_left.addChild(leaf1);
        anc_left.addChild(leaf2);
        BranchPoint leaf3 = new BranchPoint("Leaf3", anc_right, 0.09);
        BranchPoint leaf4 = new BranchPoint("Leaf4", anc_right, 0.15);
        anc_right.addChild(leaf3);
        anc_right.addChild(leaf4);
        return new IdxTree(new BranchPoint[] {root, anc_left, leaf1, leaf2, anc_right, leaf3, leaf4});
    }

    Tree tree = Newick.parse("((Leaf1:0.02,Leaf2:0.12)Anc_left:0.19,(Leaf3:0.18,Leaf4:0.15)Anc_right:0.17)Root;");
            //createTree();
    Tree mini1 = Newick.parse("((Leaf1:0.09,Leaf2:0.11)N1:0.07,Leaf3:0.12)N0");
    Tree mini2 = Newick.parse("((Leaf1:0.03,Leaf2:0.05)N1:0.13,Leaf3:0.12)N0");
    Object A = 'A';
    Object C = 'C';
    Object G = 'G';
    Object T = 'T';

    @Test
    void infer1() {
        PhyloBN pbn = PhyloBN.create(tree, SubstModel.createModel("JC"));
        MaxLhoodJoint mlj = new MaxLhoodJoint(pbn);
        TreeInstance ti = new TreeInstance(tree, new Object[] {null, null, A, null, null, C, null});
        MaxLhoodJoint.Inference inf = mlj.infer(ti);
        System.out.println("Input: \t" + ti);
        System.out.println("Output:\t" + inf);
        assertEquals(A, inf.getTreeInstance().getInstance(tree.getIndex("Leaf1")));
        assertEquals(A, inf.getTreeInstance().getInstance(tree.getIndex(1)));
        assertEquals(C, inf.getTreeInstance().getInstance(tree.getIndex(2)));
        assertEquals(A, inf.getTreeInstance().getInstance(tree.getIndex(0)));
    }

    @Test
    void infer1b() {
        SubstModel jc2 = new JC(1, new Object[]{A,C});
        PhyloBN pbn1 = PhyloBN.create(mini1, jc2);
        MaxLhoodJoint mlj1 = new MaxLhoodJoint(pbn1);
        MaxLhoodMarginal mlm1 = new MaxLhoodMarginal(0, pbn1);
        MaxLhoodMarginal mlm1b = new MaxLhoodMarginal(1, pbn1);
        TreeInstance ti1 = new TreeInstance(mini1, new Object[] {null, null, A, null, C});
        MaxLhoodJoint.Inference inf1 = mlj1.infer(ti1);
        mlm1.decorate(ti1);
        mlm1b.decorate(ti1);
        System.out.println(mlm1.getDecoration(0));
        System.out.println(mlm1b.getDecoration(1));
        System.out.println("Input: \t" + ti1);
        System.out.println("Output:\t" + inf1);
        PhyloBN pbn2 = PhyloBN.create(mini2, jc2);
        MaxLhoodJoint mlj2 = new MaxLhoodJoint(pbn2);
        MaxLhoodMarginal mlm2 = new MaxLhoodMarginal(0, pbn2);
        TreeInstance ti2 = new TreeInstance(mini2, new Object[] {null, null, A, null, C});
        MaxLhoodJoint.Inference inf2 = mlj2.infer(ti2);
        mlm2.decorate(ti2);
        System.out.println(mlm2.getDecoration(0));
        System.out.println("Input: \t" + ti2);
        System.out.println("Output:\t" + inf2);
        assertEquals(A, inf1.getTreeInstance().getInstance(tree.getIndex(0)));
        assertEquals(C, inf2.getTreeInstance().getInstance(tree.getIndex(0)));
        assertEquals(A, inf1.getTreeInstance().getInstance(tree.getIndex(1)));
        assertEquals(A, inf2.getTreeInstance().getInstance(tree.getIndex(1)));
    }

    @Test
    void infer2b() {
        PhyloBN pbn1 = PhyloBN.withGDTs(mini1, new JC(1, new Object[]{A,C}), 1, true, 1L);
        pbn1.setMasterGDT(new Object[]{A,C}, new GaussianDistrib[]{new GaussianDistrib(1.45, 0.1), new GaussianDistrib(2.55, 0.1)});
        MaxLhoodJoint mlj1 = new MaxLhoodJoint(pbn1);
        MaxLhoodMarginal mlm1 = new MaxLhoodMarginal(0, pbn1);
        TreeInstance ti1 = new TreeInstance(mini1, new Object[] {null, null, 1.45, null, 2.55});
        MaxLhoodJoint.Inference inf1 = mlj1.infer(ti1);
        mlm1.decorate(ti1);
        System.out.println(mlm1.getDecoration(0));
        PhyloBN pbn2 = PhyloBN.withGDTs(mini2, new JC(1, new Object[]{A,C}), 1, true, 1L);
        pbn2.setMasterGDT(new Object[]{A,C}, new GaussianDistrib[]{new GaussianDistrib(1.45, 0.1), new GaussianDistrib(2.55, 0.1)});
        MaxLhoodJoint mlj2 = new MaxLhoodJoint(pbn2);
        TreeInstance ti2 = new TreeInstance(mini2, new Object[] {null, null, 1.45, null, 2.55});
        MaxLhoodJoint.Inference inf2 = mlj2.infer(ti2);
        System.out.println("Input: \t" + ti1);
        System.out.println("Output:\t" + inf1);
        System.out.println("Input: \t" + ti2);
        System.out.println("Output:\t" + inf2);
        assertEquals(A, inf1.getTreeInstance().getInstance(tree.getIndex(0)));
        assertEquals(C, inf2.getTreeInstance().getInstance(tree.getIndex(0)));
        assertEquals(A, inf1.getTreeInstance().getInstance(tree.getIndex(1)));
        assertEquals(A, inf2.getTreeInstance().getInstance(tree.getIndex(1)));
    }
    @Test
    void infer3a() {
        PhyloBN pbn1 = PhyloBN.withGDTs(tree, new JC(1, new Object[]{A,C}), 1, true, 1L);
        pbn1.setMasterGDT(new Object[]{A,C}, new GaussianDistrib[]{new GaussianDistrib(1.45, 0.1), new GaussianDistrib(2.55, 0.1)});
        MaxLhoodJoint mlj1 = new MaxLhoodJoint(pbn1);
        MaxLhoodMarginal mlm1 = new MaxLhoodMarginal(0, pbn1);
        TreeInstance ti = new TreeInstance(tree, new Object[] {null, null, 1.5, 1.4, null, 2.5, 2.6});
        MaxLhoodJoint.Inference inf = mlj1.infer(ti);
        mlm1.decorate(ti);
        System.out.println(mlm1.getDecoration(0));
        System.out.println("Input: \t" + ti);
        System.out.println("Output:\t" + inf);
        assertEquals(inf.getTreeInstance().getInstance(tree.getIndex(1)), A);
        assertEquals(inf.getTreeInstance().getInstance(tree.getIndex(2)), C);
        assertEquals(inf.getTreeInstance().getInstance(tree.getIndex(0)), C);
    }
    @Test
    void infer3b() {
        PhyloBN pbn1 = PhyloBN.withGDTs(tree, new JC(1, new Object[]{A,C}), 1, true, 1L);
        pbn1.setMasterGDT(new Object[]{A,C}, new GaussianDistrib[]{new GaussianDistrib(1.45, 0.1), new GaussianDistrib(2.55, 0.1)});
        MaxLhoodJoint mlj1 = new MaxLhoodJoint(pbn1);
        MaxLhoodMarginal mlm1 = new MaxLhoodMarginal(0, pbn1);
        TreeInstance ti = new TreeInstance(tree, new Object[] {null, null, 1.5, null, null, null, 2.6});
        MaxLhoodJoint.Inference inf = mlj1.infer(ti);
        mlm1.decorate(ti);
        System.out.println(mlm1.getDecoration(0));
        System.out.println("Input: \t" + ti);
        System.out.println("Output:\t" + inf);
        assertEquals(inf.getTreeInstance().getInstance(tree.getIndex(1)), A);
        assertEquals(inf.getTreeInstance().getInstance(tree.getIndex(2)), C);
        assertEquals(inf.getTreeInstance().getInstance(tree.getIndex(0)), C);
    }

    @Test
    void joint1() {
        SubstModel model = SubstModel.createModel("JC");
        Enumerable alpha = new Enumerable(model.getDomain().getValues());
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        double SCALEDIST = 1.0;
        int NLEAVES = 8; //
        for (int SEED = 0; SEED < 100; SEED ++) {
            // make a random tree with N - 1 ancestors for a binary tree, 2 x N - 1 variables in total
            Tree tree = Tree.Random(NLEAVES, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, 2, 2);
            tree.adjustDistances(SCALEDIST);
            int[] ancidxs = tree.getAncestors();
            int[] leafidxs = tree.getLeaves();
            // generate a random set of leaf states
            Random rand = new Random(SEED);
            Object[] allstates = new Object[tree.getSize()];
            Object[] leafstates = new Object[leafidxs.length];
            for (int i = 0; i < leafstates.length; i ++) {
                leafstates[i] = alpha.get(rand.nextInt(alpha.size()));
                allstates[leafidxs[i]] = leafstates[i];
            }
            // use bnkit to do the inference...
            MaxLhoodJoint mlj = new MaxLhoodJoint(tree, model);
            MaxLhoodJoint.Inference inf = mlj.infer(new TreeInstance(tree, allstates));
            TreeInstance ti_bnkit = inf.getTreeInstance();
            // to be benchmarked against the naive product of all conditional probs...
            // nominate all possible ancestors
            Object[][] ancstates = new Object[(int)Math.pow(alpha.size(), ancidxs.length)][];
            double[] jointp = new double[(int)Math.pow(alpha.size(), ancidxs.length)];
            for (int j  = 0; j < ancstates.length; j ++)
                ancstates[j] = alpha.getWord4Key(j, ancidxs.length);
            // now, main loop: for each ancestor state...
            int maxidx = 0;
            for (int j  = 0; j < ancstates.length; j ++) {
                for (int i = 0; i < ancstates[j].length; i++)
                    allstates[ancidxs[i]] = ancstates[j][i];
                double joint = 1; // joint prob of the states
                for (int idx : tree) { // loop through all nodes in the tree
                    int child = idx;
                    int parent = tree.getParent(child);
                    double p = parent < 0 ? model.getProb(allstates[child]) : model.getProb(allstates[child], allstates[parent], tree.getDistance(child));
                    joint *= p;
                }
                jointp[j] = joint;
                if (joint > jointp[maxidx])
                    maxidx = j;
            }
            for (int i = 0; i < ancstates[maxidx].length; i++)
                allstates[ancidxs[i]] = ancstates[maxidx][i];
            TreeInstance ti_naive = new TreeInstance(tree, allstates);
            Object[] ancs_bnkit = new Object[ancidxs.length];
            for (int i = 0; i < ancidxs.length; i ++)
                ancs_bnkit[i] = ti_bnkit.getInstance(ancidxs[i]);
            int bnkit_idx = alpha.getKey4Word(ancs_bnkit);
            if (bnkit_idx != maxidx) {
                System.out.println(ti_bnkit);
                System.out.println(" [" + bnkit_idx + "] -LogL = " + -Math.log(jointp[bnkit_idx]));
                System.out.println(ti_naive);
                System.out.println(" [" + maxidx + "] -LogL = " + -Math.log(jointp[maxidx])); // + "\t < " + " [" + 42 + "] -LogL = " + -Math.log(jointp[42]));
            }

            for (int idx : ancidxs) {
//                assertEquals(ti_naive.getInstance(idx), ti_bnkit.getInstance(idx));
            }
        }
    }
}