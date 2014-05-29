/*
 * Copyright (C) 2014 mikael
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
import bn.CPT;
import bn.Distrib;
import bn.EnumDistrib;
import bn.EnumVariable;
import bn.Enumerable;
import bn.GDT;
import bn.GaussianDistrib;
import bn.JPT;
import bn.Predef;
import bn.Variable;
import bn.Variable.Assignment;
import bn.alg.CGTable;
import bn.alg.CGVarElim;
import bn.alg.Query;
import bn.alg.QueryResult;

/**
 *
 * @author mikael
 */
public class GDTExample2 {
        /**
     * Demo of how FactorTable can be used in variable elimination inference.
     * Example amended from Russell and Norvig (2003; p. 493-494).
     * @param args
     */
    public static void main(String[] args) {

        // Define variables
        EnumVariable B = Predef.Boolean("Burglary");
        EnumVariable E = Predef.Boolean("Earthquake");
        Variable S     = Predef.Real("Seismic");
        Variable R     = Predef.Real("Rumble");
        EnumVariable A = Predef.Boolean("Alarm");
        Variable N     = Predef.Real("Noise");
//        EnumVariable N = Predef.Boolean("Noise");
        EnumVariable C = Predef.Boolean("Concert");
        EnumVariable J = Predef.Boolean("John calls");
        EnumVariable M = Predef.Boolean("Mary calls");

        // Define nodes (connecting the variables into an acyclic graph, i.e. the structure)
        CPT b = new CPT(B);
        CPT e = new CPT(E);
        CPT c = new CPT(C);
        CPT a = new CPT(A,    B, E);
        CPT j = new CPT(J,    A);
        CPT m = new CPT(M,    A, C);
        GDT s = new GDT(S,    E);
        GDT r = new GDT(R,    B, E);
        GDT n = new GDT(N,    A, C);
//        CPT n = new CPT(N,    A, C);
        
        // Parameterise the nodes using our "expertise"
        b.put(new EnumDistrib(Enumerable.bool, 0.001, 0.999));
        e.put(new EnumDistrib(Enumerable.bool, 0.002, 0.998));
        c.put(new EnumDistrib(Enumerable.bool, 0.03, 0.97));
        a.put(new EnumDistrib(Enumerable.bool, 0.95, 0.05), true, true);
        a.put(new EnumDistrib(Enumerable.bool, 0.94, 0.06), true, false);
        a.put(new EnumDistrib(Enumerable.bool, 0.29, 0.71), false, true);
        a.put(new EnumDistrib(Enumerable.bool, 0.001, 0.999), false, false);
        j.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), true);
        j.put(new EnumDistrib(Enumerable.bool, 0.05, 0.95), false);
        m.put(new EnumDistrib(Enumerable.bool, 0.50, 0.50), true, true);
        m.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), false, true);
        m.put(new EnumDistrib(Enumerable.bool, 0.70, 0.30), true, false);
        m.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), false, false);
        s.put(new GaussianDistrib(6.0, 2.0), true);
        s.put(new GaussianDistrib(2.0, 2.0), false);
        r.put(new GaussianDistrib(0.65, 0.25), true, true);
        r.put(new GaussianDistrib(0.30, 0.25), true, false);
        r.put(new GaussianDistrib(0.50, 0.25), false, true);
        r.put(new GaussianDistrib(0.20, 0.25), false, false);
        n.put(new GaussianDistrib(150, 40), true, true);
        n.put(new GaussianDistrib(100, 40), true, false);
        n.put(new GaussianDistrib(100, 40), false, true);
        n.put(new GaussianDistrib(50, 40), false, false);
//        n.put(new EnumDistrib(Enumerable.bool, 0.8, 0.2), true, true);
//        n.put(new EnumDistrib(Enumerable.bool, 0.5, 0.5), true, false);
//        n.put(new EnumDistrib(Enumerable.bool, 0.5, 0.5), false, true);
//        n.put(new EnumDistrib(Enumerable.bool, 0.1, 0.9), false, false);
                
        b.print();
        e.print();
        c.print();
        a.print();
        j.print();
        m.print();
        s.print();
        r.print();
        n.print();
        
        // Put all the nodes into the Bayesian network data structure
        // The BNet class manages efficient access to the nodes, based on the structure.
        BNet bn = new BNet();
	bn.add(b, e, c, a, j, m, s, r, n);

        // Once in the BNet, variables (through the nodes) can be instantiated to values from their respective domains.
        b.setInstance(false);
        j.setInstance(true);
        //m.setInstance(true);
        c.setInstance(false);
        //n.setInstance(100.0);
        //s.setInstance(6.5);
        //r.setInstance(0.55);
        
        // Variable elimination works by factorising CPTs, and then by performing products and variable sum-outs in
        // an order that heuristically is computationally efficient.
        
        CGVarElim ve = new CGVarElim();
        ve.instantiate(bn);
        Query q = ve.makeQuery(new Variable[] {M,E,N,R,S});
        //Query q = ve.makeQuery(new Variable[] {S});
        CGTable qr = (CGTable) ve.infer(q);
        qr.display();
        qr.displaySampled();
        Distrib d = qr.query(R);
        System.out.println(d);
        d = qr.query(N, new Assignment[] {
            Variable.assign(R, 0.01),
            Variable.assign(M, true),
            });
        System.out.println(d);
        d = qr.query(S, new Assignment[] {
            Variable.assign(M, true),
            });
        System.out.println(d);
        
    }

}