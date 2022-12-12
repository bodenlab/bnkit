/*
 *  bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2015  M. Boden et al.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package bn.example;

import bn.BNet;
import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.PhyloBNet;
import bn.ctmc.matrix.JTT;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.PhyloTree;
import dat.PhyloTree.Node;
import dat.Variable.Assignment;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author mikael
 */
public class ASRwGamma {
    static String file_tree = "/Users/mikael/simhome/ASR/CYP2/CYP2F.nwk";
    static String file_aln = "/Users/mikael/simhome/ASR/CYP2/CYP2F.aln";
    
    static PhyloTree tree = new PhyloTree();
    static List<EnumSeq.Gappy<Enumerable>> seqs;
    static EnumSeq.Alignment<Enumerable> aln;

    public static double check(EnumDistrib[] distribs) {
        Object[] shouldBe = new Object[] {'S','G','D','M','D','H','A'};
        int[] columns = new int[] {33,119,160,210,212,278,421};
        double cnt = 0;
        for (int i = 0; i < columns.length; i ++) {
            cnt += distribs[columns[i]].get(shouldBe[i]);
//            if (distribs[columns[i]].getMax().equals(shouldBe[i]))
//                cnt ++;
        }
        return cnt;
    }
    
    public static void main(String[] args) {
        
        int N_Gamma = 100; // number of samples for estimating joint with gamma
        double[] alphas = new double[] {1.67, 1.69, 1.71, 1.73};
        try {
            
//            for (int a = 0; a < alphas.length; a ++) {
                //System.out.println("Alpha = " + alphas[a]);

                
            tree = tree.loadNewick(file_tree);
            Node[] nodes = tree.toNodesBreadthFirst();
            List<String> indexForNodes = new ArrayList<>(); // Newick string for subtree
            Map<String, String> mapForNodes = new HashMap<>(); // Shortname --> Newick string for subtree
            for (Node n : nodes) {
                indexForNodes.add(replacePunct(n.toString()));
                mapForNodes.put(n.getLabel().toString(), replacePunct(n.toString()));
            }
            
            seqs = EnumSeq.Gappy.loadClustal(file_aln, Enumerable.aacid);
            aln = new EnumSeq.Alignment<>(seqs);

            String[] names = aln.getNames();
            List<String> labels = new ArrayList<>();
            
            for (int i = 0; i < names.length; i ++) {
                int index = names[i].indexOf("/"); // in this aln file, names have been amended
                if (index > 0)
                    labels.add(names[i].substring(0, index));
                else
                    labels.add(names[i]);
            }
            
            PhyloBNet[] pbnets = new PhyloBNet[aln.getWidth()];
            
            Object[][] asr_matrix = new Object[indexForNodes.size()][aln.getWidth()]; // single joint reconstruction for tree
            
            double[][] counts = new double[aln.getWidth()][Enumerable.aacid.size()];   // counts from which distributions of joint reconstructions for root are estimated
            EnumDistrib[] distribs = new EnumDistrib[aln.getWidth()];

            // 1. prepare BNs with constant rate--to gather data for estimating a gamma distribution
            
            for (int col = 0; col < aln.getWidth(); col ++) {
                Object[] gaps = aln.getGapColumn(col); // array with true for gap, false for symbol
                Object[] column = aln.getColumn(col);  // array for symbols, null for gaps
                tree.setContentByParsimony(names, gaps);
                PhyloBNet pbn;
                pbn = PhyloBNet.create(tree, new JTT());
                pbnets[col] = pbn;

                // set variables according to alignment
                for (int i = 0; i < labels.size(); i ++) {
                    String shortname = labels.get(i);
                    String longname = mapForNodes.get(shortname);
                    if (longname != null) {
                        BNode bnode = pbn.getBN().getNode(longname);
                        bnode.setInstance(column[i]);
                    }
                }
            }

            BNode root = null;
            String asr_root = null;
            
            double[] R = new double[aln.getWidth()];
            for (int col = 0; col < aln.getWidth(); col ++) {
                PhyloBNet pbn = pbnets[col];
                BNet bn = pbn.getBN();
                root = pbn.getRoot();
                VarElim ve = new VarElim();
                ve.instantiate(bn);

                List<EnumVariable> intern = pbn.getInternal();

                int purged_leaves = pbn.purgeGaps();
                int collapsed_nodes = pbn.collapseSingles();
                //System.out.println("Col " + col + "\tPurged: " + purged_leaves + " + " + collapsed_nodes);
                Query q_joint = ve.makeMPE();
                CGTable r_joint = (CGTable)ve.infer(q_joint);
                Assignment[] ass = r_joint.getMPE();
                for (Assignment a0 : ass) {
                    EnumVariable asr_var = (EnumVariable)a0.var;
                    Object asr_val = a0.val;
                    int index = indexForNodes.indexOf(replacePunct(asr_var.getName()));
                    if (index >= 0) 
                        asr_matrix[index][col] = asr_val;
                    BNode node = bn.getNode(asr_var);
                    node.setInstance(asr_val);
                }
                R[col] = pbn.getRate();
            }
            
            List<EnumSeq.Gappy<Enumerable>> asrs = new ArrayList<>();
            for (int row = 0; row < asr_matrix.length; row ++) {
                Object[] asr_obj = asr_matrix[row];
                EnumSeq.Gappy<Enumerable> myasr = new EnumSeq.Gappy<>(Enumerable.aacid_alt);
                myasr.set(asr_obj);
                myasr.setName(indexForNodes.get(row));
                asrs.add(myasr);
            }
            
            String rootname = root.getVariable().getName();
            
            EnumSeq.Alignment aln_asr = new EnumSeq.Alignment(asrs);
            for (int i = 0; i < aln_asr.getHeight(); i ++) {
                EnumSeq.Gappy<Enumerable> asr_seq = aln_asr.getEnumSeq(i);
                String nodename = asr_seq.getName();
                if (rootname.equals(nodename))
                    asr_root = asr_seq.toString();
                //System.out.println(asr_seq.getName() + "\t" + asr_seq.toString());
            }
            
            /*
            System.out.println("Joint reconstruction (with rate 1.0): " + asr_root);
            System.out.println("Actual rates from joint reconstruction used to estimate gamma distribution");
            for (int col = 0; col < aln.getWidth(); col ++) {
                System.out.println(col + "\t" + asr_root.charAt(col) + "\t" + R[col]);
            } 
            */
            double alpha = GammaDistrib.getAlpha(R);
            
            // ******
            //alpha = alphas[a];
            // ******
            alpha = 1.73;
            double beta = 1 / alpha;
            System.out.println("Gamma alpha = " + alpha + " beta = " + beta);
            GammaDistrib gd = new GammaDistrib(alpha, 1.0/beta);
            double mean = 0.0;
//           System.out.println("Sample (showing only first 10)");
            double[] rates = new double[N_Gamma];
            for (int i = 0; i < N_Gamma; i ++) {
                rates[i] = gd.sample();
                mean += rates[i];
//                if (i < 10)
//                   System.out.println(i + "\t" + rates[i]);
            }
            System.out.println("Mean\t" + mean / N_Gamma + " in the limit it should be 1.0");

            
            // 2. Repeat for samples drawn from gamma
            
            for (int g = 0; g < N_Gamma; g ++) {
                
                for (int col = 0; col < aln.getWidth(); col ++) {
                    Object[] gaps = aln.getGapColumn(col); // array with true for gap, false for symbol
                    Object[] column = aln.getColumn(col);  // array for symbols, null for gaps
                    tree.setContentByParsimony(names, gaps);
                    PhyloBNet pbn;
                    pbn = PhyloBNet.create(tree, new JTT(), rates[g]); // rate sampled from reconstructions when rate was 1.0
                    pbnets[col] = pbn;

                    // set variables according to alignment
                    for (int i = 0; i < labels.size(); i ++) {
                        String shortname = labels.get(i);
                        String longname = mapForNodes.get(shortname);
                        if (longname != null) {
                            BNode bnode = pbn.getBN().getNode(longname);
                            bnode.setInstance(column[i]);
                        }
                    }
                }

                for (int col = 0; col < aln.getWidth(); col ++) {
                    PhyloBNet pbn = pbnets[col];
                    BNet bn = pbn.getBN();
                    root = pbn.getRoot();
                    VarElim ve = new VarElim();
                    ve.instantiate(bn);

                    List<EnumVariable> intern = pbn.getInternal();

                    int purged_leaves = pbn.purgeGaps();
                    int collapsed_nodes = pbn.collapseSingles();
                    //System.out.println("Col " + col + "\tPurged: " + purged_leaves + " + " + collapsed_nodes);
                    Query q_joint = ve.makeMPE();
                    CGTable r_joint = (CGTable)ve.infer(q_joint);
                    Assignment[] ass = r_joint.getMPE();
                    for (Assignment a0 : ass) {
                        EnumVariable asr_var = (EnumVariable)a0.var;
                        Object asr_val = a0.val;
                        int index = indexForNodes.indexOf(replacePunct(asr_var.getName()));
                        if (index >= 0) 
                            asr_matrix[index][col] = asr_val;
                        BNode node = bn.getNode(asr_var);
                        node.setInstance(asr_val);
                    }
                }

                asrs = new ArrayList<>();
                for (int row = 0; row < asr_matrix.length; row ++) {
                    Object[] asr_obj = asr_matrix[row];
                    EnumSeq.Gappy<Enumerable> myasr = new EnumSeq.Gappy<>(Enumerable.aacid_alt);
                    myasr.set(asr_obj);
                    myasr.setName(indexForNodes.get(row));
                    asrs.add(myasr);
                }

                rootname = root.getVariable().getName();

                EnumSeq.Gappy<Enumerable> sampled_root = null;
                aln_asr = new EnumSeq.Alignment(asrs);
                for (int i = 0; i < aln_asr.getHeight(); i ++) {
                    EnumSeq.Gappy<Enumerable> asr_seq = aln_asr.getEnumSeq(i);
                    String nodename = asr_seq.getName();
                    if (rootname.equals(nodename))
                        sampled_root = asr_seq;
                    // System.out.println(asr_seq.getName() + "\t" + asr_seq.toString());
                }
                
                for (int col = 0; col < aln.getWidth(); col ++) {
                    counts[col][Enumerable.aacid.getIndex(sampled_root.get()[col])] += 1.0 / N_Gamma;
                }
                
            }
            
            System.out.print("\t  ");
            for (int j = 0; j < Enumerable.aacid.size(); j ++) {
                System.out.print(Enumerable.aacid.get(j) + "    ");
            }
            System.out.println();
            for (int col = 0; col < aln.getWidth(); col ++) {
                distribs[col] = new EnumDistrib(Enumerable.aacid, counts[col]);
                if (distribs[col].getMax().equals(asr_root.charAt(col)))
                    System.out.println(col + "\t" + distribs[col] + "\t" + distribs[col].getMax());
                else
                    System.out.println(col + "\t" + distribs[col] + "\t" + distribs[col].getMax()+ "\t" + asr_root.charAt(col));
            }
            double cnt = check(distribs);
            System.out.println("Correct = " + cnt);
            // ******
            //}
        } catch (IOException ex) {
            ex.printStackTrace();

        }
    }
    
    private static String replacePunct(String str) {
        return str.replace('.', '_');
    }
}
