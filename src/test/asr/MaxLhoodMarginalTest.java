package asr;

import bn.CountTable;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;
import bn.prob.MixtureDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import json.JSONObject;
import org.junit.jupiter.api.Test;
import util.MilliTimer;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class MaxLhoodMarginalTest {

    PhyloBN pbn_gdt, pbn_cpt, pbn_gdt_ext, pbn_cpt_ext;
    // This is a tree Sanjana generated from the DHAD9112 tree, modified a little to have shorter labels and distances with fewer decimals. Here's it's Newick string if required:
    // (((S001:0.14,(S002:0.16,S003:0.1)N3:0.08)N2:0.03,S004:0.12)N1:0.14,((S005:0.28,(S006:0.12,S007:0.14)N6:0.13)N5:0.06,(S008:0.2,(S009:0.12,S010:0.19)N8:0.07)N7:0.11)N4:0.06)N0:0.0;
    IdxTree tree = IdxTree.fromJSON(new JSONObject("{\"Parents\":[-1,0,1,2,2,4,4,1,0,8,9,9,11,11,8,14,14,16,16],\"Labels\":[\"0\",\"1\",\"2\",\"S001\",\"3\",\"S002\",\"S003\",\"S004\",\"4\",\"5\",\"S005\",\"6\",\"S006\",\"S007\",\"7\",\"S008\",\"8\",\"S009\",\"S010\"],\"Distances\":[0,0.14,0.03,0.14,0.08,0.16,0.10,0.12,0.06,0.06,0.28,0.13,0.12,0.14,0.11,0.20,0.07,0.12,0.19],\"Branchpoints\":19}\n"));
    String[] headers = {   "S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"};
    // The values assigned to the tree above, tabulated as per headers; they were intended to group nicely per clade
    Double[][] rows1 = {{   3.63,  3.81,  2.89,  3.81,  2.54,  2.76,  3.79,  3.70,  1.94,  3.97}};
    Double[][] rows1tst = {{3.6,   3.8,   2.9,   3.8,   2.5,   2.8,   3.8,   3.7,   1.9,   4.0}};
    String[][] rows2 = {{   "b",   "c",   "a",   "c",   "a",   "a",   "b",   "b",   "a",   "c"}};
    String[][] rows2tst = {{"b",   "c",   "a",   "c",   "a",   "a",   "b",   "b",   "a",   "c"}};
    String[] states = {"A", "B"};
    String[] states3 = {"A", "B", "C"};
    double GAMMA = 1.0;
    SubstModel model2 = new JC(GAMMA, states);
    SubstModel model3 = new JC(GAMMA, states3);

    void setup_gdt() {
        pbn_gdt = PhyloBN.withGDTs(tree, model2, 1);
        pbn_gdt.trainEM(headers, rows1, 1L);
        // BNBuf.save(pbn_gdt.getBN(), "gdt10.xml");
    }

    void setup_gdt_ext() {
        pbn_gdt_ext = PhyloBN.withGDTs(tree, model2, 1, false, 1L);
        pbn_gdt_ext.trainEM(headers, rows1, 1L);
        // BNBuf.save(pbn_gdt.getBN(), "gdtext10.xml");
    }

    void setup_cpt() {
        Set<String> alphas = new HashSet<>();
        for (int i = 0; i < rows2[0].length; i ++)
            alphas.add(rows2[0][i]);
        String[] alphax = new String[alphas.size()];
        alphas.toArray(alphax);
        pbn_cpt = PhyloBN.withCPTs(tree, model2, alphax, 1);
        pbn_cpt.trainEM(headers, rows2, 1L);
        // BNBuf.save(pbn_cpt.getBN(), "cpt10.xml");
    }
    void setup_cpt_ext() {
        Set<String> alphas = new HashSet<>();
        for (int i = 0; i < rows2[0].length; i ++)
            alphas.add(rows2[0][i]);
        String[] alphax = new String[alphas.size()];
        alphas.toArray(alphax);
        pbn_cpt_ext = PhyloBN.withCPTs(tree, model3, alphax, 1, false, 1L);
        pbn_cpt_ext.trainEM(headers, rows2, 1L);
        // BNBuf.save(pbn_cpt.getBN(), "cptext10.xml");
    }

    @Test
    void getDecoration1a() {
        // create a BN with a topology from tree, adding GDT leaves to all extants
        setup_gdt();
        // Go through all branch points in the tree
        // Expect that the root represents the majority class (bp-idx >= 8)
        // and that the minority ancestors bp-idx 1-4 belong to opposites of majority ancestors bp-idx 8-16
        int root_max_index = 0;
        for (int idx : tree) {
            // if the branch point is a parent, it is an ancestor
            if (tree.isParent(idx)) {
                Object label = tree.getLabel(idx);
                // set-up inference for that ancestor
                MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(idx, pbn_gdt);
                // look at the data and match the given values to particular branch points in a tree-instance
                TreeInstance ti = tree.getInstance(headers, rows1tst[0]);
                // put the values on the tree and infer the nominated ancestor
                inf.decorate(ti);
                /* can save the values to a file if desired:
                try {
                    Newick.save(ti, "c10.nwk");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                 */
                // retrieve the value, which happens to be a distribution of enumerable
                EnumDistrib d = inf.getDecoration(idx);
                System.out.println(idx + "\t" + label + "\t" + d);
                if (idx == 0)
                    root_max_index = d.getMaxIndex();
                else {
                    assertTrue(idx < 8 ? root_max_index != d.getMaxIndex() : root_max_index == d.getMaxIndex());
                }
            }
        }
    }

    @Test
    void getDecoration1b() {
        // create a BN with a topology from tree, adding GDT leaves to all extants and ancestors
        setup_gdt_ext();
        // Go through all branch points in the tree
        // Expect that the root represents the middle prediction
        double root_avg = 0;
        double minority_avg = 0;
        int minority_cnt = 0;
        double majority_avg = 0;
        int majority_cnt = 0;
        for (int idx : tree) {
            // if the branch point is a parent, it is an ancestor
            if (tree.isParent(idx)) {
                Object label = tree.getLabel(idx);
                // set-up inference for that ancestor
                MaxLhoodMarginal<MixtureDistrib> inf = new MaxLhoodMarginal<>(idx, pbn_gdt_ext);
                // look at the data and match the given values to particular branch points in a tree-instance
                TreeInstance ti = tree.getInstance(headers, rows1tst[0]);
                // put the values on the tree and infer the nominated ancestor
                inf.decorate(ti);
                // retrieve the value, which happens to be a distribution of enumerable
                MixtureDistrib d = inf.getDecoration(idx);
                System.out.println(idx + "\t" + label + "\t" + d);
                if (idx == 0) {
                    for (int i = 0; i < 30; i ++)
                        root_avg += (((Double) d.sample()) / 30);
                } else {
                    if (idx < 8) {
                        minority_cnt += 1;
                        for (int i = 0; i < 30; i ++)
                            minority_avg += (((Double) d.sample()) / 30);
                    } else {
                        majority_cnt += 1;
                        for (int i = 0; i < 30; i ++)
                            majority_avg += (((Double) d.sample()) / 30);
                    }
                }
            }
        }
        minority_avg /= minority_cnt;
        majority_avg /= majority_cnt;
        assertTrue(minority_avg < root_avg && root_avg < majority_avg);
    }

    @Test
    void getDecoration2a() {
        setup_cpt();
        int root_max_index = 0;
        for (int idx : tree) {
            if (tree.isParent(idx)) {
                Object label = tree.getLabel(idx);
                MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(idx, pbn_cpt);
                TreeInstance ti = tree.getInstance(headers, rows2tst[0]);
                inf.decorate(ti);
                EnumDistrib d = inf.getDecoration(idx);
                System.out.println(idx + "\t" + label + "\t" + d);
                if (idx == 0)
                    root_max_index = d.getMaxIndex();
                else {
                    assertTrue(idx < 8 ? root_max_index != d.getMaxIndex() : root_max_index == d.getMaxIndex());
                }
            }
        }
    }
    @Test
    void getDecoration2b() {
        // Extended nodes on all ancestor nodes
        // three different labels lining up with clades and three latent states, so
        // expecting EM to produce roughly the same number of ancestors for each latent state
        setup_cpt_ext();
        int[] cnt = new int[3];
        for (int idx : tree) {
            if (tree.isParent(idx)) {
                Object label = tree.getLabel(idx);
                MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(idx, pbn_cpt_ext);
                TreeInstance ti = tree.getInstance(headers, rows2tst[0]);
                inf.decorate(ti);
                EnumDistrib d = inf.getDecoration(idx);
                System.out.println(idx + "\t" + label + "\t" + d);
                if (idx == 0)
                    ;
                else {
                    cnt[d.getMaxIndex()] += 1;
                }
            }
        }
        assertTrue(Math.abs(cnt[0] - cnt[1]) <= 1 && Math.abs(cnt[1] - cnt[2]) <= 1);
    }

    /**
     * Test of exact marginal inference in phylogenetic trees, comparing against a naive product, exhaustively determining the prob of all combinations.
     * Testing accuracy and reporting time efficiency.
     * Naive_ML_Leaves=8	13758
     * Naive_ML_Leaves=7	 2996
     * Naive_ML_Leaves=6	  707
     * Naive_ML_Leaves=5	  145
     * Naive_ML_Leaves=4	   33
     * Naive_ML_Leaves=3	    7
     * BNKit_ML_Leaves=8	   14
     * BNKit_ML_Leaves=7	   12
     * BNKit_ML_Leaves=6	   12
     * BNKit_ML_Leaves=5	    7
     * BNKit_ML_Leaves=4	   23
     * BNKit_ML_Leaves=3	   30
     */
    @Test
    void marginal1() {
        SubstModel model = SubstModel.createModel("JC"); // create a Jukes-Cantor substitution model
        Enumerable alpha = new Enumerable(model.getDomain().getValues()); // the character states used by the above model (ACGT)
        // we will create random trees, so the following parameters will control what they look like
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        double SCALEDIST = 1.0; // calibrate the leaf to root distance to be around 1
        MilliTimer timer = new MilliTimer();
        // the naive product is copied from MaxLhoodJointTest
        for (int NLEAVES = 3; NLEAVES < 9; NLEAVES ++) {
            timer.start("Leaves=" + NLEAVES);
            for (int SEED = 0; SEED < 100; SEED++) { // generate a tree for each random seed
                // make a random tree with N - 1 ancestors for a binary tree, 2 x N - 1 variables in total
                Tree tree = Tree.Random(NLEAVES, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, 2, 2);
                tree.adjustDistances(SCALEDIST);
                int[] ancidxs = tree.getAncestors();    // indices that recover all internal nodes; always start with "0"
                int[] leafidxs = tree.getLeaves();      // indices for all leaves
                // generate a random set of leaf states
                Random rand = new Random(SEED);
                Object[] allstates = new Object[tree.getSize()];            // holder for states of all nodes (indexed by tree)
                Object[] leafstates = new Object[leafidxs.length];          // holder for the subset of "leaf" states
                for (int i = 0; i < leafstates.length; i++) {              // for each leaf...
                    leafstates[i] = alpha.get(rand.nextInt(alpha.size()));  // randomly select one of the character states (defined in the JC model)
                    allstates[leafidxs[i]] = leafstates[i];                 // patch the holder of all states
                }
                // use bnkit to do the inference...
                timer.start("BNKit_ML_Leaves=" + NLEAVES);
                MaxLhoodMarginal mlj = new MaxLhoodMarginal(0, tree, model);         // create inference object
                mlj.decorate(new TreeInstance(tree, allstates));
                EnumDistrib ed_bnkit = (EnumDistrib) mlj.getDecoration(0);
                timer.stop("BNKit_ML_Leaves=" + NLEAVES);
                // to be benchmarked against marginalising the naive product of all conditional probs...
                // nominate all possible ancestors
                timer.start("Naive_ML_Leaves=" + NLEAVES);
                Object[][] ancstates = new Object[(int) Math.pow(alpha.size(), ancidxs.length)][];   // holds all combinations of ancestor states
                for (int j = 0; j < ancstates.length; j++)                                        // they are enumerated here...
                    ancstates[j] = alpha.getWord4Key(j, ancidxs.length);                            // assigned the state...
                double[] jointp = new double[(int) Math.pow(alpha.size(), ancidxs.length)];          // and each will be assigned a probability based on tree and leaf states
                // now, main loop: for each ancestor state...
                for (int j = 0; j < ancstates.length; j++) {  //
                    for (int i = 0; i < ancstates[j].length; i++)   // assign the ancestor state to the holder of all states
                        allstates[ancidxs[i]] = ancstates[j][i];    // only ancestor states are over-written (leaf states have been assigned already)
                    double joint = 1;                               // start calc joint prob of the states
                    for (int idx : tree) {  // loop through all nodes in the tree
                        int child = idx;    // the perspective is that the current node is a child (of a possible parent node)
                        int parent = tree.getParent(child); // this is the parent
                        double p;       // determine the probability of the implied substitution
                        if (parent < 0) {   // no parent, so use a-priori prob of child state...
                            p = model.getProb(allstates[child]);    // taken from the substitution model
                        } else {            // parent indeed...
                            double[][] probs = model.getProbs(tree.getDistance(child));     // the model has them all...
                            p = model.getProb(allstates[child], allstates[parent], probs);  // so extract the appropriate substitution
                        }
                        joint *= p; // the prob of the combination of states is the product of their (conditional) probabilities
                    }
                    jointp[j] = joint;              // put the prob in place for the list of all ancestor combinations
                } // all ancestor combinations have been calculated...
                // marginalise; sum-out all variables except N0
                double[] counts = new double[alpha.size()];
                for (int j = 0; j < ancstates.length; j++)
                    counts[alpha.getIndex(ancstates[j][0])] += jointp[j]; // count the character seen in the first position (which is N0)
                EnumDistrib ed_naive = new EnumDistrib(alpha, counts); // create a distribution for N0
                // done marginalising...
                timer.stop("Naive_ML_Leaves=" + NLEAVES);
                assertArrayEquals(ed_naive.get(), ed_bnkit.get(), 0.01);
            }
            timer.stop("Leaves=" + NLEAVES);
        }
        timer.report(true);
    }

}