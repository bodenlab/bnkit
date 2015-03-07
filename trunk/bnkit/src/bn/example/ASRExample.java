/*
 *  bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2014  M. Boden et al.
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
import bn.factor.Factorize;
import bn.file.BNBuf;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.PhyloTree;
import dat.PhyloTree.Node;
import dat.Variable.Assignment;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author mikael
 */
public class ASRExample {
    static String file_tree = "/Users/mikael/simhome/ASR/gap_example.nwk";
    static String file_aln = "/Users/mikael/simhome/ASR/gap_example.aln";
    
    static PhyloTree tree;
    static List<EnumSeq.Gappy<Enumerable>> seqs;
    static EnumSeq.Alignment<Enumerable> aln;
    
    public static void main(String[] args) {
        try {
            tree = PhyloTree.loadNewick(file_tree);
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
            
            Object[][] asr_matrix = new Object[indexForNodes.size()][aln.getWidth()]; // joint reconstruction for tree

            for (int col = 0; col < aln.getWidth(); col ++) {
                Object[] gaps = aln.getGapColumn(col); // array with true for gap, false for symbol
                Object[] column = aln.getColumn(col);  // array for symbols, null for gaps
                tree.setContentByParsimony(names, gaps);
                PhyloBNet pbn = PhyloBNet.create(tree, new JTT());
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

            double[][] margin_probs = new double[Enumerable.aacid.size()][aln.getWidth()]; // marginal reconstruction for root
            EnumDistrib[] margin_distribs = new EnumDistrib[aln.getWidth()];
            
            BNode root = null;
            String asr_root = null;
            
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
                Query q_marg = ve.makeQuery(root.getVariable());
                CGTable r_marg = (CGTable)ve.infer(q_marg);
                EnumDistrib d_marg = (EnumDistrib)r_marg.query(root.getVariable());
                margin_distribs[col] = d_marg;
                Query q_joint = ve.makeMPE();
                CGTable r_joint = (CGTable)ve.infer(q_joint);
                Assignment[] a = r_joint.getMPE();
                for (Assignment a0 : a) {
                    EnumVariable asr_var = (EnumVariable)a0.var;
                    Object asr_val = a0.val;
                    int index = indexForNodes.indexOf(replacePunct(asr_var.getName()));
                    if (index >= 0) 
                        asr_matrix[index][col] = asr_val;
                    BNode node = bn.getNode(asr_var);
                    node.setInstance(asr_val);
                }
//                BNBuf.save(pbn.getBN(), "/Users/mikael/test.xml");
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
                System.out.println(asr_seq.getName() + "\t" + asr_seq.toString());
            }
            
            System.out.println("Joint reconstruction: " + asr_root);
            System.out.print("  ");
            for (int j = 0; j < Enumerable.aacid.size(); j ++) {
                System.out.print(Enumerable.aacid.get(j) + "    ");
            }
            System.out.println();
            for (int col = 0; col < aln.getWidth(); col ++) {
                System.out.println(margin_distribs[col] + "\t" + asr_root.charAt(col));
            }            
            
        } catch (IOException ex) {
            ex.printStackTrace();

        }
    }
    
    private static String replacePunct(String str) {
        return str.replace('.', '_');
    }
}
