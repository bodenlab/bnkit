package asr;

import bn.BNode;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.factor.AbstractFactor;
import bn.factor.Factor;
import bn.factor.Factorize;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;
import dat.EnumTable;
import dat.Enumerable;
import dat.Variable;
import dat.file.Newick;
import dat.phylo.*;
import org.junit.jupiter.api.Test;
import util.MilliTimer;

import java.util.*;

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

    Tree tree = Newick.parse("((Leaf1:0.02,Leaf2:0.12)Anc_left:0.19,(Leaf3:0.18,Leaf4:0.15)Anc_right:0.191)Root;");
            //createTree();
    Tree mini1 = Newick.parse("((Leaf1:0.09,Leaf2:0.11)N1:0.07,Leaf3:0.12)N0");
    Tree mini2 = Newick.parse("((Leaf1:0.03,Leaf2:0.05)N1:0.13,Leaf3:0.12)N0");
    Object A = 'A';
    Object C = 'C';
    Object G = 'G';
    Object T = 'T';

    int NSEEDs = 100;
    int NPICKS = 3; // how many variables to pick values for

    @Test
    void infer0() {
        Tree t = tree;
        SubstModel model = SubstModel.createModel("WAG");
        PhyloBN pbn = PhyloBN.create(t, model);
        MaxLhoodJoint mlj = new MaxLhoodJoint(pbn);
        for (int SEED = 0; SEED < NSEEDs; SEED ++) {
            Random rand = new Random(SEED + 1);
            BNode[] bnodes = new BNode[t.getSize()];
            List<Integer> mylist = new ArrayList<>();
            for (int bpidx : t) {
                bnodes[bpidx] = pbn.getBNode(bpidx);
                mylist.add(bpidx); // here ordered
            }
            Collections.shuffle(mylist, rand);
            int[] rank = new int[bnodes.length]; // the "rank" of a node
            int cnt = 0;
            for (int idx : mylist)
                rank[idx] = cnt ++;

            // set instances for the variables that have been "picked" (indexed by their node index)
            Map<Variable, Object> inst4map = new HashMap<>();
            Object[] inst4arr = new Object[mylist.size()];
            for (int i = 0; i < NPICKS; i++) {
                int pick = mylist.get(i);
                Object inst = model.getDomain().getValues()[rand.nextInt(model.getDomain().size())];
                inst4map.put(bnodes[pick].getVariable(), inst);
                inst4arr[pick] = inst;
            }

            int[] free = new int[bnodes.length - NPICKS];
            int idx = 0;
            for (int i = NPICKS; i < bnodes.length; i++) {
                int not_picked = mylist.get(i);
                free[idx ++] = not_picked;
            }
            TreeInstance ti = new TreeInstance(t, inst4arr);
            // System.out.println("Input: \t" + ti);
            MaxLhoodJoint.Inference inf = mlj.infer(ti);
            // System.out.println("Output:\t" + inf);

            Object[][] states = new Object[(int) Math.pow(model.getDomain().size(), free.length)][]; // holds all combinations of "free" states
            for (int j = 0; j < states.length; j++)                                                  // they are enumerated here...
                states[j] = model.getDomain().getWord4Key(j, free.length);                           // assigned the state...
            double[] jointp = new double[(int) Math.pow(model.getDomain().size(), free.length)];     // and each will be assigned a probability based on tree and leaf states

            int maxp_idx = -1;
            for (int j = 0; j < jointp.length; j ++) {          // for every possible "free" state...
                // calculate the joint probability of each permutation of assignments, which is to...
                jointp[j] = 1;
                for (int i = 0; i < bnodes.length; i++) { // calculate the product of all conditional probabilities
                    BNode node = bnodes[i];
                    int paridx = t.getParent(i);
                    // Four possibilities
                    if (inst4arr[i] == null) {
                        // 1: variable is not set, parent variable is not set (or the node is not conditioned, i.e. root)
                        Object varinstance = states[j][rank[i] - NPICKS];
                        if (paridx == -1) { // root (1a)
                            jointp[j] *= node.get(varinstance);
                        } else if (inst4arr[paridx] == null) { // not root (1b), and the parent is not set
                            jointp[j] *= node.get(varinstance, states[j][rank[paridx] - NPICKS]);
                        } else { //
                            // 2: variable is not set, but parent is/condition is set
                            jointp[j] *= node.get(varinstance, inst4arr[paridx]);
                        }
                    } else {
                        // 3: variable is set
                        Object varinstance = inst4arr[i];
                        if (paridx == -1) { // root (3a)
                            jointp[j] *= node.get(varinstance);
                        } else if (inst4arr[paridx] == null) { // not root (3b), has parent that is not set
                            jointp[j] *= node.get(varinstance, states[j][rank[paridx] - NPICKS]);
                        } else {
                            // 4: variable is set, parent is set
                            jointp[j] *= node.get(varinstance, inst4arr[paridx]);
                        }
                    }
                }
                if (maxp_idx == -1)
                    maxp_idx = j;
                else if (jointp[j] > jointp[maxp_idx])
                    maxp_idx = j;
            }
            for (int i = 0; i < bnodes.length; i++) { // show values that give the maximum product
//                System.out.print("Node " + bnodes[i].getName() + " = ");
                if (inst4arr[i] == null) {
//                    System.out.println(states[maxp_idx][rank[i] - NPICKS]);
                    assertEquals(inf.getTreeInstance().getInstance(i), states[maxp_idx][rank[i] - NPICKS]);
                } else {
//                    System.out.println(inst4arr[i] + " *");
                }
            }
        }
    }


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

    /**
     * Test of exact joint inference in phylogenetic trees, comparing against a naive product, exhaustively determining the prob of all combinations.
     * Testing accuracy and reporting time efficiency.
     * Data collected from MacBook Pro M1 2022:
     * Naive_ML_Leaves=8	    27335
     * Naive_ML_Leaves=7	     5949
     * Naive_ML_Leaves=6	     1312
     * Naive_ML_Leaves=5	      300
     * Naive_ML_Leaves=4	       88
     * Naive_ML_Leaves=3	       39
     * BNKit_ML_Leaves=8	       10
     * BNKit_ML_Leaves=7	        4
     * BNKit_ML_Leaves=6	        4
     * BNKit_ML_Leaves=5	        5
     * BNKit_ML_Leaves=4	       11
     * BNKit_ML_Leaves=3	       22
     * BNKit_ML_CPT_Leaves=8	   13
     * BNKit_ML_CPT_Leaves=7	    9
     * BNKit_ML_CPT_Leaves=6	   11
     * BNKit_ML_CPT_Leaves=5	   10
     * BNKit_ML_CPT_Leaves=4	   12
     * BNKit_ML_CPT_Leaves=3	   26
     * BNKit_ML_GDT_Leaves=8	    7
     * BNKit_ML_GDT_Leaves=7	    5
     * BNKit_ML_GDT_Leaves=6	    6
     * BNKit_ML_GDT_Leaves=5	    6
     * BNKit_ML_GDT_Leaves=4	   10
     * BNKit_ML_GDT_Leaves=3	   15
     */
    @Test
    void joint1() {
        SubstModel model = SubstModel.createModel("JC"); // create a Jukes-Cantor substitution model
        Enumerable alpha = new Enumerable(model.getDomain().getValues()); // the character states used by the above model (ACGT)
        // we will create random trees, so the following parameters will control what they look like
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        double SCALEDIST = 1.0; // calibrate the leaf to root distance to be around 1
        MilliTimer timer = new MilliTimer();

        for (int NLEAVES = 3; NLEAVES < 8; NLEAVES ++) {
            // int NLEAVES = 7; // this is how many leaf nodes each tree will have
            timer.start("Leaves=" + NLEAVES);
            for (int SEED = 0; SEED < 100; SEED++) { // generate a tree for each random seed
                // make a random tree with N - 1 ancestors for a binary tree, 2 x N - 1 variables in total
                Tree tree = Tree.Random(NLEAVES, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, 2, 2);
                tree.adjustDistances(SCALEDIST);
                int[] ancidxs = tree.getAncestors();    // indices that recover all internal nodes
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
                MaxLhoodJoint mlj = new MaxLhoodJoint(tree, model);         // create inference object
                MaxLhoodJoint.Inference inf = mlj.infer(new TreeInstance(tree, allstates)); // infer
                timer.stop("BNKit_ML_Leaves=" + NLEAVES);
                TreeInstance ti_bnkit = inf.getTreeInstance();              // recover the states from inference
                // to be benchmarked against the naive product of all conditional probs...
                // nominate all possible ancestors
                timer.start("Naive_ML_Leaves=" + NLEAVES);
                Object[][] ancstates = new Object[(int) Math.pow(alpha.size(), ancidxs.length)][];   // holds all combinations of ancestor states
                for (int j = 0; j < ancstates.length; j++)                                        // they are enumerated here...
                    ancstates[j] = alpha.getWord4Key(j, ancidxs.length);                            // assigned the state...
                double[] jointp = new double[(int) Math.pow(alpha.size(), ancidxs.length)];          // and each will be assigned a probability based on tree and leaf states
                // now, main loop: for each ancestor state...
                int maxidx = 0; // remember which state (by index) has the max prob
                String[] trace = new String[ancstates.length];  // document the calcs done for each ancestor state combination
                for (int j = 0; j < ancstates.length; j++) {  //
                    for (int i = 0; i < ancstates[j].length; i++)   // assign the ancestor state to the holder of all states
                        allstates[ancidxs[i]] = ancstates[j][i];    // only ancestor states are over-written (leaf states have been assigned already)
                    double joint = 1;                               // start calc joint prob of the states
                    StringBuilder sb = new StringBuilder();         // document the calcs
                    for (int idx : tree) {  // loop through all nodes in the tree
                        int child = idx;    // the perspective is that the current node is a child (of a possible parent node)
                        int parent = tree.getParent(child); // this is the parent
                        double p;       // determine the probability of the implied substitution
                        if (parent < 0) {   // no parent, so use a-priori prob of child state...
                            p = model.getProb(allstates[child]);    // taken from the substitution model
                            sb.append(String.format("P(N%d=%s)=%5.3f", child, allstates[child], p));
                        } else {            // parent indeed...
                            //double[][] probs = model.getProbs(tree.getDistance(child));     // the model has them all...
                            p = model.getProb(allstates[child], allstates[parent], tree.getDistance(child));  // so extract the appropriate substitution
                            sb.append(String.format("P(N%d=%s|N%d=%s)=%5.3f", child, allstates[child], parent, allstates[parent], p));
                        }
                        joint *= p; // the prob of the combination of states is the product of their (conditional) probabilities
                        if (idx < tree.getSize() - 1)
                            sb.append(" x ");
                    }
                    trace[j] = sb.toString();
                    jointp[j] = joint;              // put the prob in place for the list of all ancestor combinations
                    if (joint > jointp[maxidx])     // check if max
                        maxidx = j;
                } // all ancestor combinations have been calculated...
                timer.stop("Naive_ML_Leaves=" + NLEAVES);
                for (int i = 0; i < ancstates[maxidx].length; i++) // place the best combination in "allstates" (which already has the leaf states)
                    allstates[ancidxs[i]] = ancstates[maxidx][i];
                TreeInstance ti_naive = new TreeInstance(tree, allstates);  // put them in a tree, so we can easily look at them
                Object[] ancs_bnkit = new Object[ancidxs.length];           // holder for the best combination with bnkit's ML inference (above)
                for (int i = 0; i < ancidxs.length; i++)
                    ancs_bnkit[i] = ti_bnkit.getInstance(ancidxs[i]);       // take them from the tree instance from before
                int bnkit_idx = alpha.getKey4Word(ancs_bnkit);              // bnkit's inference would correspond to a state enumeration index
                if (bnkit_idx != maxidx) {                                  // if they are different (which they should not be)...
                    System.out.println(ti_bnkit);                           // print out some useful info for debugging
                    System.out.println(" [" + bnkit_idx + "] -LogL = " + -Math.log(jointp[bnkit_idx]));
                    System.out.println(ti_naive);
                    System.out.println(" [" + maxidx + "] -LogL = " + -Math.log(jointp[maxidx])); // + "\t < " + " [" + 42 + "] -LogL = " + -Math.log(jointp[42]));
                }
                for (int idx : ancidxs) { // just look at each ancestor, to test the same outcome with naive as with bnkit ML
                    assertEquals(ti_naive.getInstance(idx), ti_bnkit.getInstance(idx));
                }

                // next test accessory CPT nodes to those in the tree
                String[] lowercase_alpha = new String[]{"a", "c", "g", "t"};   // use lowercase to distinguish the accessory leaf states from states in the actual tree
                Enumerable lowercase_enum = new Enumerable(lowercase_alpha);    // create an alphabet from them
                // next create a BN from tree, with each leaf node having an extra, accessory node hanging of it
                PhyloBN pbn = PhyloBN.withCPTs(tree, model, lowercase_alpha, 1, true, SEED + 1);
                Enumerable leaf_alpha = model.getDomain();                      // the latent states
                pbn.setMasterCPT(leaf_alpha.getValues(), new EnumDistrib[]{    // set the (shared/master) conditional prob table for the accessory nodes
                        new EnumDistrib(lowercase_enum, 0.997, 0.001, 0.001, 0.001), // A is parent
                        new EnumDistrib(lowercase_enum, 0.001, 0.997, 0.001, 0.001), // C is parent
                        new EnumDistrib(lowercase_enum, 0.001, 0.001, 0.997, 0.001), // G is parent
                        new EnumDistrib(lowercase_enum, 0.001, 0.001, 0.001, 0.997)  // T is parent
                });
                Object[] lowercase_states = new Object[tree.getSize()];         // holder for accessory states
                for (int i = 0; i < leafstates.length; i++)                    // we're pinching the leaf states from the original, already tested trees/BNs
                    lowercase_states[leafidxs[i]] = leafstates[i].toString().toLowerCase();
                TreeInstance ti_inp = new TreeInstance(tree, lowercase_states); // from them, create a tree instance that can be used as input to bnkit inference
                timer.start("BNKit_ML_CPT_Leaves=" + NLEAVES);
                MaxLhoodJoint mlj_cpt = new MaxLhoodJoint(pbn);                 // setup of ML inference
                MaxLhoodJoint.Inference inf_cpt = mlj_cpt.infer(ti_inp);        // bnkit inference...
                timer.stop("BNKit_ML_CPT_Leaves=" + NLEAVES);
                TreeInstance ti_cpt = inf_cpt.getTreeInstance();                // extract output from inference
                // inferred states should be identical to the max (naive) product of all conditional probs that was determined before
                Object[] ancs_cpt = new Object[ancidxs.length];                 // collect the ancestor states from bnkit's output
                for (int i = 0; i < ancidxs.length; i++)
                    ancs_cpt[i] = ti_cpt.getInstance(ancidxs[i]);
                int cpt_idx = alpha.getKey4Word(ancs_cpt);                      // determine the index of bnkit's favourite state combination
                if (cpt_idx != maxidx) {                                        // debug: were they not the same? they should be
                    System.out.println(ti_cpt);
                    System.out.println(" [" + cpt_idx + "] -LogL = " + -Math.log(jointp[cpt_idx]));
                    System.out.println(ti_naive);
                    System.out.println(" [" + maxidx + "] -LogL = " + -Math.log(jointp[maxidx])); // + "\t < " + " [" + 42 + "] -LogL = " + -Math.log(jointp[42]));
                }
                for (int idx : ancidxs) {                                       // junit testing
                    assertEquals(ti_naive.getInstance(idx), ti_cpt.getInstance(idx));
                }

                // next test accessory GDT nodes to those in the tree
                PhyloBN pbn_gdt = PhyloBN.withGDTs(tree, model, 1, true, SEED + 1);
                GaussianDistrib[] gds = new GaussianDistrib[]{        // spike the table so that parent state is fixed (in principle)
                        new GaussianDistrib(0, 0.01),   // A is parent
                        new GaussianDistrib(1, 0.01),   // C is parent
                        new GaussianDistrib(2, 0.01),   // G is parent
                        new GaussianDistrib(3, 0.01)    // T is parent
                };
                pbn_gdt.setMasterGDT(leaf_alpha.getValues(), gds);
                Object[] double_states = new Object[tree.getSize()];    // holder of accessory values
                for (int i = 0; i < leafstates.length; i++)            // set them so each parent state equals that of the original tree (above)
                    double_states[leafidxs[i]] = gds[leaf_alpha.getIndex(leafstates[i])].sample();
                TreeInstance ti_inp2 = new TreeInstance(tree, double_states);
                timer.start("BNKit_ML_GDT_Leaves=" + NLEAVES);
                MaxLhoodJoint mlj_gdt = new MaxLhoodJoint(pbn_gdt);
                MaxLhoodJoint.Inference inf_gdt = mlj_gdt.infer(ti_inp2);
                timer.stop("BNKit_ML_GDT_Leaves=" + NLEAVES);
                TreeInstance ti_gdt = inf_gdt.getTreeInstance();
                // to be benchmarked against the naive product of all conditional probs...
                Object[] ancs_gdt = new Object[ancidxs.length];
                for (int i = 0; i < ancidxs.length; i++)
                    ancs_gdt[i] = ti_gdt.getInstance(ancidxs[i]);
                int gdt_idx = alpha.getKey4Word(ancs_gdt);
                if (gdt_idx != maxidx) {
                    System.out.println(ti_gdt);
                    System.out.println(" [" + gdt_idx + "] -LogL = " + -Math.log(jointp[gdt_idx]));
                    System.out.println(ti_naive);
                    System.out.println(" [" + maxidx + "] -LogL = " + -Math.log(jointp[maxidx])); // + "\t < " + " [" + 42 + "] -LogL = " + -Math.log(jointp[42]));
                }
                for (int idx : ancidxs) {
                    assertEquals(ti_naive.getInstance(idx), ti_gdt.getInstance(idx));
                }
            }
            timer.stop("Leaves=" + NLEAVES);
        }
        timer.report(true);
    }

    int NTHREADS = 10;
    @Test // reconstruction when multithreaded
    void jointThreaded() {
        SubstModel model = SubstModel.createModel("LG"); // create a Jukes-Cantor substitution model
        Enumerable alpha = new Enumerable(model.getDomain().getValues()); // the character states used by the above model (ACGT)
        // we will create random trees, so the following parameters will control what they look like
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        double SCALEDIST = 1.0; // calibrate the leaf to root distance to be around 1

        for (int NLEAVES = 3; NLEAVES < 5; NLEAVES ++) { // this is how many leaf nodes each tree will have
            TreeDecor[] inf = new TreeDecor[NSEEDs];            // number of seeds we try is also how many inferences we will carry out
            TreeInstance[] treeinstances = new TreeInstance[NSEEDs];
            Object[][] allstates = new Object[NSEEDs][];
            for (int SEED = 0; SEED < NSEEDs; SEED++) { // generate a tree for each random seed
                // make a random tree with N - 1 ancestors for a binary tree, 2 x N - 1 variables in total
                Tree tree = Tree.Random(NLEAVES, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, 2, 2);
                tree.adjustDistances(SCALEDIST);
                allstates[SEED] = new Object[tree.getSize()];            // holder for states of all nodes (indexed by tree)
                Object[] instantiated = new Object[tree.getSize()];            // holder for states of instantiated nodes (indexed by tree)
                int[] leafidxs = tree.getLeaves();      // indices for all leaves
                // generate a random set of leaf states
                Random rand = new Random(SEED);
                Object[] leafstates = new Object[leafidxs.length];          // holder for the subset of "leaf" states
                for (int i = 0; i < leafstates.length; i++) {              // for each leaf...
                    leafstates[i] = alpha.get(rand.nextInt(alpha.size()));  // randomly select one of the character states (defined in the JC model)
                    allstates[SEED][leafidxs[i]] = leafstates[i];                 // patch the holder of all states
                    instantiated[leafidxs[i]] = leafstates[i];                 // patch the holder of all states
                }
                // use bnkit to do the inference...
                inf[SEED] = new MaxLhoodJoint(tree, model);        //   configure inference
                treeinstances[SEED] = new TreeInstance(tree, instantiated);
                // infer the joint prob naively ...
                int[] ancidxs = treeinstances[SEED].getTree().getAncestors();
                Object[][] ancstates = new Object[(int) Math.pow(alpha.size(), ancidxs.length)][];   // holds all combinations of ancestor states
                for (int j = 0; j < ancstates.length; j++)                                        // they are enumerated here...
                    ancstates[j] = alpha.getWord4Key(j, ancidxs.length);                            // assigned the state...
                double[] jointp = new double[(int) Math.pow(alpha.size(), ancidxs.length)];          // and each will be assigned a probability based on tree and leaf states
                // now, main loop: for each ancestor state...
                int maxidx = 0; // remember which state (by index) has the max prob
                String[] trace = new String[ancstates.length];  // document the calcs done for each ancestor state combination
                for (int j = 0; j < ancstates.length; j++) {  //
                    for (int i = 0; i < ancstates[j].length; i++)   // assign the ancestor state to the holder of all states
                        allstates[SEED][ancidxs[i]] = ancstates[j][i];    // only ancestor states are over-written (leaf states have been assigned already)
                    double joint = 1;                               // start calc joint prob of the states
                    StringBuilder sb = new StringBuilder();         // document the calcs
                    for (int idx : tree) {  // loop through all nodes in the tree
                        int child = idx;    // the perspective is that the current node is a child (of a possible parent node)
                        int parent = tree.getParent(child); // this is the parent
                        double p;       // determine the probability of the implied substitution
                        if (parent < 0) {   // no parent, so use a-priori prob of child state...
                            p = model.getProb(allstates[SEED][child]);    // taken from the substitution model
                            sb.append(String.format("P(N%d=%s)=%5.3f", child, allstates[SEED][child], p));
                        } else {            // parent indeed...
                            // double[][] probs = model.getProbs(tree.getDistance(child));     // the model has them all...
                            p = model.getProb(allstates[SEED][child], allstates[SEED][parent], tree.getDistance(child));  // so extract the appropriate substitution
                            sb.append(String.format("P(N%d=%s|N%d=%s)=%5.3f", child, allstates[SEED][child], parent, allstates[SEED][parent], p));
                        }
                        joint *= p; // the prob of the combination of states is the product of their (conditional) probabilities
                        if (idx < tree.getSize() - 1)
                            sb.append(" x ");
                    }
                    trace[j] = sb.toString();
                    jointp[j] = joint;              // put the prob in place for the list of all ancestor combinations
                    if (joint > jointp[maxidx])     // check if max
                        maxidx = j;
                } // all ancestor combinations have been calculated...
                for (int i = 0; i < ancstates[maxidx].length; i++) // place the best combination in "allstates" (which already has the leaf states)
                    allstates[SEED][ancidxs[i]] = ancstates[maxidx][i];
            }
            ThreadedDecorators threadpool = new ThreadedDecorators(inf, treeinstances, NTHREADS);
            try {
                Map<Integer, TreeDecor> ret = threadpool.runBatch();
                for (int SEED = 0; SEED < NSEEDs; SEED++) {           // for each seed...
                    int[] ancidxs = treeinstances[SEED].getTree().getAncestors();
                    for (int idx : ancidxs) {                      // for each ancestor...
                        assertEquals(ret.get(SEED).getDecoration(idx), allstates[SEED][idx]);  //     extract state
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
}