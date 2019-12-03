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
import bn.Distrib;
import bn.Predef;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.QueryResult;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.*;
import bn.node.CPT;
import dat.EnumVariable;
import dat.Variable;
import dat.Variable.Assignment;

/**
 *
 * @author mikael
 */
public class SubstExample {
    static SubstModel model = new JTT();
    // The value n in the PAM-n matrix represents the number of mutations per 100 amino acids
    // r is the rate of accepted mutation accumulation in mutations per amino acid site per million years
    // static double r_base = 2.2e-9 * 1e6; // per base pair per million years, in mammalian genomes according to Kumar,  803–808, doi: 10.1073/pnas.022629899
    // static double r = r_base * 3; // per amino acid
    static double T = 0.1;
    // static double K = T * 2 * r; // T = K/2r where K is the number of mutations per amino acid, and 

    // For probability of insertion and deletions, see
    // Empirical and Structural Models for Insertions and Deletions in the Divergent Evolution of Proteins
    // Steven A. Benner, Mark A. Cohen, Gaston H. Gonnet, doi:10.1006/jmbi.1993.1105

    static EnumVariable aa_a1, aa_b1, aa_c1, aa_x1, aa_x2;
    static SubstNode node_a1, node_b1, node_c1, node_x1, node_x2;
    static BNet bn;

    public SubstExample() {
        System.out.println("Evolutionary time: " + T);
        aa_a1 = Predef.AminoAcid("A1");
        aa_b1 = Predef.AminoAcid("B1");
        aa_c1 = Predef.AminoAcid("C1");
        aa_x1 = Predef.AminoAcid("X1");
        aa_x2 = Predef.AminoAcid("X2");
        node_a1 = new SubstNode(aa_a1, aa_x1, model, T);
        node_b1 = new SubstNode(aa_b1, aa_x1, model, T);
        node_x1 = new SubstNode(aa_x1, aa_x2, model, T);
        node_c1 = new SubstNode(aa_c1, aa_x2, model, T);
        node_x2 = new SubstNode(aa_x2, model);
        bn = new BNet();
        bn.add(node_x1, node_a1, node_b1, node_c1, node_x2);
        node_a1.setInstance('F');
        node_b1.setInstance('V');
        node_c1.setInstance('I');
    }

    public void runJoint() {
        VarElim ve = new VarElim();
        ve.instantiate(bn);
        Query q = ve.makeMPE();
        CGTable qr = (CGTable)ve.infer(q);
        System.out.println(qr);
        double joint = qr.getFactor();
        System.out.println("\t" + joint);
        Assignment[] as = qr.getMPE();
        for (Assignment a : as)
            System.out.println("\t" + a);
        double loglikelihood = ve.logLikelihood();
        System.out.println("LL = " +loglikelihood + "\tL = " + Math.exp(loglikelihood));
    }

    public void runMarginal(Variable aa) {
        VarElim ve = new VarElim();
        ve.instantiate(bn);
        Query q = ve.makeQuery(aa);
        CGTable qr = (CGTable)ve.infer(q);
        qr.display();
        Distrib d = qr.query(aa);
        System.out.println(d);
    }
    public static void main(String[] args) {
        SubstExample example = new SubstExample();
        example.runJoint();
        example.runMarginal(aa_x1);
        example.runMarginal(aa_x2);
    }
        
}
