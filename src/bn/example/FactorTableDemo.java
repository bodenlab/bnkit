/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package bn.example;

import bn.BNet;
import bn.BNode;
import bn.node.CPT;
import bn.node.GDT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import bn.prob.GaussianDistrib;
import bn.Predef;
import dat.Variable;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.factor.AbstractFactor;
import bn.factor.Factorize;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author mikael
 */
public class FactorTableDemo {

    /**
     * Demo of how FactorTable can be used in variable elimination inference.
     * Example amended from Russell and Norvig (2003; p. 493-494).
     * @param args
     */
    public static void main(String[] args) {

        // Define variables
        EnumVariable B = Predef.Boolean("Burglary");
        EnumVariable E = Predef.Boolean("Earthquake");
        Variable S     = Predef.Real("Seismic signal");
        EnumVariable A = Predef.Boolean("Alarm");
        EnumVariable J = Predef.Boolean("John calls");
        EnumVariable M = Predef.Boolean("Mary calls");

        // Define nodes (connecting the variables into an acyclic graph, i.e. the structure)
        CPT b = new CPT(B);
        CPT e = new CPT(E);
        CPT a = new CPT(A,    B, E);
        CPT j = new CPT(J,    A);
        CPT m = new CPT(M,    A);
        GDT s = new GDT(S,    E);
        
        // Parameterise the nodes using our "expertise"
        b.put(new EnumDistrib(Enumerable.bool, 0.001, 0.999));
        e.put(new EnumDistrib(Enumerable.bool, 0.002, 0.998));
        a.put(new EnumDistrib(Enumerable.bool, 0.95, 0.05), true, true);
        a.put(new EnumDistrib(Enumerable.bool, 0.94, 0.06), true, false);
        a.put(new EnumDistrib(Enumerable.bool, 0.29, 0.71), false, true);
        a.put(new EnumDistrib(Enumerable.bool, 0.001, 0.999), false, false);
        j.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), true);
        j.put(new EnumDistrib(Enumerable.bool, 0.05, 0.95), false);
        m.put(new EnumDistrib(Enumerable.bool, 0.70, 0.30), true);
        m.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), false);
        s.put(new GaussianDistrib(6.0, 3.0), true);
        s.put(new GaussianDistrib(2.0, 3.0), false);
        
        b.print();
        e.print();
        a.print();
        j.print();
        m.print();
        s.print();
        
        // Put all the nodes into the Bayesian network data structure
        // The BNet class manages efficient access to the nodes, based on the structure.
        BNet bn = new BNet();
	bn.add(b, e, a, j, m, s);

        // Once in the BNet, variables (through the nodes) can be instantiated to values from their respective domains.
        j.setInstance(true);
        m.setInstance(true);
        //s.setInstance(5.5);
        
        // Variable elimination works by factorising CPTs, and then by performing products and variable sum-outs in
        // an order that heuristically is computationally efficient.
        
        Map<Variable, Object> evidence = new HashMap<>(); // evidence
        for (BNode node : bn.getOrdered()) {
            Variable var = node.getVariable();
            Object val = node.getInstance(); // will be null if not instantiated
            evidence.put(var, val);
        }

        // First we make each CPT into a FactorTable, considering variables that are instantiated.
        // We assume that all nodes are involved in the inference, though that is not always going to be true.
        // In fact, BNet has a method for creating a new BNet instance that does not contain nodes that are
        // irrelevant to a particular query.
        AbstractFactor ft_b = b.makeDenseFactor(evidence);
        System.out.println("Factor B");
        ft_b.display();
        
        AbstractFactor ft_e = e.makeDenseFactor(evidence);
        System.out.println("Factor E");
        ft_e.display();
        
        AbstractFactor ft_a = a.makeDenseFactor(evidence);
        System.out.println("Factor A");
        ft_a.display();
        
        AbstractFactor ft_j = j.makeDenseFactor(evidence);
        System.out.println("Factor J");
        ft_j.display();
        
        AbstractFactor ft_m = m.makeDenseFactor(evidence);
        System.out.println("Factor M");
        ft_m.display();
        
        AbstractFactor ft_s = s.makeDenseFactor(evidence);
        System.out.println("Factor S");
        ft_s.display();
        
        

        // To produce a JPT, all FactorTables relevant to the query need to enter into the product
        // at some point. From a theoretical point of view the order has no impact.
        // But the order in which these products are done can impact on the number
        // of operations that are required. Smaller, overlapping FTs give smaller 
        // products. We thus combine those that have variables in common. The topological order
        // may also be very helpful to use. See Dechter's paper.
        AbstractFactor ft = Factorize.getProduct(ft_b, ft_e);
        System.out.println("Factor B * E");
        ft.display();

        ft = Factorize.getProduct(ft, ft_s);
        System.out.println("Factor (B * E) * S");
        ft.display();

        ft = Factorize.getProduct(ft, ft_a);
        System.out.println("Factor (B * E) * A");
        ft.display();

        ft = Factorize.getProduct(ft, ft_j);
        System.out.println("Factor ((B * E) * A) * J");
        ft.display();

        ft = Factorize.getProduct(ft, ft_m);
        System.out.println("Factor (((B * E) * A) * J) * M)");
        ft.display();

        // Next we sum-out variables that have not been specified/instantiated or are part of the query.
        // Note that we could have done this earlier, and this would have resulted in smaller products.
        // Variable elimination also optimises when variables are summed-out.
        ft = Factorize.getMargin(ft, new EnumVariable[] {A, E});
        ft.display();

        // Normalise the FT
        ft = Factorize.getNormal(ft);
        
        System.out.println("Factor (((B * E) * A) * J) * M) - (A, E)");
        ft.display();

        VarElim ve = new VarElim();
        ve.instantiate(bn);
        Query q = ve.makeQuery(B);
        CGTable qr = (CGTable) ve.infer(q);
        qr.display();

    }
}

