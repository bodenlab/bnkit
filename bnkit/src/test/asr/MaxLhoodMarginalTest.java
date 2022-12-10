package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.file.BNBuf;
import bn.prob.EnumDistrib;
import dat.file.Newick;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import json.JSONObject;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class MaxLhoodMarginalTest {

    PhyloBN pbn_gdt, pbn_cpt;
    IdxTree tree = IdxTree.fromJSON(new JSONObject("{\"Parents\":[-1,0,1,2,2,4,4,1,0,8,9,9,11,11,8,14,14,16,16],\"Labels\":[\"0\",\"1\",\"2\",\"S001\",\"3\",\"S002\",\"S003\",\"S004\",\"4\",\"5\",\"S005\",\"6\",\"S006\",\"S007\",\"7\",\"S008\",\"8\",\"S009\",\"S010\"],\"Distances\":[0,0.14,0.03,0.14,0.08,0.16,0.10,0.12,0.06,0.06,0.28,0.13,0.12,0.14,0.11,0.20,0.07,0.12,0.19],\"Branchpoints\":19}\n"));
    String[] headers = {   "S009","S005","S002","S006","S003","S001","S008","S010","S004","S007"};
    Double[][] rows1 = {{   3.63,  3.81,  2.89,  3.81,  2.54,  2.76,  3.79,  3.70,  1.94,  3.97}};
    Double[][] rows1tst = {{3.6,   3.8,   2.9,   3.8,   2.5,   2.8,   3.8,   3.7,   1.9,   4.0}};
    String[][] rows2 = {{   "b",   "c",   "a",   "c",   "a",   "a",   "b",   "b",   "a",   "c"}};
    String[][] rows2tst = {{"b",   "c",   "a",   "c",   "a",   "a",   "b",   "b",   "a",   "c"}};
    String[] states = {"A", "B"};
    double GAMMA = 1.0;
    SubstModel model = new JC(GAMMA, states);

    void setup_gdt() {
        pbn_gdt = PhyloBN.withGDTs(tree, model, 1);
        pbn_gdt.trainEM(headers, rows1, 1L);
        // BNBuf.save(pbn_gdt.getBN(), "gdt10.xml");
    }

    void setup_cpt() {
        Set<String> alphas = new HashSet<>();
        for (int i = 0; i < rows2[0].length; i ++)
            alphas.add(rows2[0][i]);
        String[] alphax = new String[alphas.size()];
        alphas.toArray(alphax);
        pbn_cpt = PhyloBN.withCPTs(tree, model, alphax, 1);
        pbn_cpt.trainEM(headers, rows2, 1L);
        // BNBuf.save(pbn_cpt.getBN(), "cpt10.xml");
    }

    @Test
    void getDecoration1() {
        setup_gdt();
        for (int idx : tree) {
            if (tree.isParent(idx)) {
                Object label = tree.getLabel(idx);
                MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(idx, pbn_gdt);
                TreeInstance ti = tree.getInstance(headers, rows1tst[0]);
                inf.decorate(ti);
                try {
                    Newick.save(ti, "c10.nwk");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                EnumDistrib d = inf.getDecoration(idx);
                System.out.println(idx + "\t" + label + "\t" + d);
            }
        }
    }

    @Test
    void getDecoration2() {
        setup_cpt();
        for (int idx : tree) {
            if (tree.isParent(idx)) {
                Object label = tree.getLabel(idx);
                MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(idx, pbn_cpt);
                TreeInstance ti = tree.getInstance(headers, rows2tst[0]);
                try {
                    Newick.save(ti, "d10.nwk");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                inf.decorate(ti);
                EnumDistrib d = inf.getDecoration(idx);
                System.out.println(idx + "\t" + label + "\t" + d);
            }
        }
    }
}