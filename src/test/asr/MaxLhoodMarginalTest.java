package asr;

import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import bn.prob.MixtureDistrib;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.TreeInstance;
import json.JSONObject;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
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
}