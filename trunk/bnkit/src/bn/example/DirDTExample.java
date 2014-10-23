/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bn.example;
import bn.node.CPT;
import bn.node.DirDT;
import dat.EnumVariable;
import dat.Variable;
import dat.Enumerable;
import bn.prob.DirichletDistrib;
import bn. *;
import bn.alg.*;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author mikael
 */
public class DirDTExample {
    
    public static void main(String[] args) {

        // Define variables
        EnumVariable G = Predef.Nominal(new String[] {"Male", "Female"}, "Gender");
        Enumerable colours = new Enumerable(new String[] {"Pink", "Green", "Blue"});
        Enumerable sports = new Enumerable(new String[] {"Netball", "Soccer", "Rugby"});
        Variable C     = Predef.Distrib(colours, "Colours");
        Variable S     = Predef.Distrib(sports,  "Sports");

        // Define nodes (connecting the variables into an acyclic graph, i.e. the structure)
        CPT g = new CPT(G);
        DirDT c = new DirDT(C,    G);
        DirDT s = new DirDT(S,    G);

        // Parameterise the nodes using our "expertise"
        g.put(new EnumDistrib(new Enumerable(new String[] {"Male", "Female"}), 0.49, 0.51));
        c.put(new DirichletDistrib(colours, new double[] {3.0, 5.0, 7.0}), "Male");
        c.put(new DirichletDistrib(colours, new double[] {7.0, 2.0, 3.0}), "Female");
        s.put(new DirichletDistrib(sports,  new double[] {2.0, 5.0, 5.0}), "Male");
        s.put(new DirichletDistrib(sports,  new double[] {9.0, 4.0, 3.0}), "Female");
        
        BNet bn = new BNet();
        bn.add(g,c,s);
        c.setInstance(new EnumDistrib(colours, new double[] {0.3, 0.3, 0.4})); // primarily blue, so probably male...
        VarElim inf = new VarElim();
        inf.instantiate(bn);
        Query q = inf.makeQuery(G,S);
        CGTable r = (CGTable) inf.infer(q);
        r.display();
        Distrib d1 = r.query(S);
        System.out.println("Prob of sports: " + d1);
        double[] means = new double[sports.size()];
        int NSAMPLE = 20;
        for (int i = 0; i < NSAMPLE; i ++) {
            EnumDistrib d1_sample = (EnumDistrib)d1.sample();
            for (int j = 0; j < sports.size(); j ++) 
                means[j] += d1_sample.get(j) / NSAMPLE;
            System.out.println("\t" + (i+1) + "\t" + d1_sample);
        }
        System.out.print("\tMean\t");
        for (int j = 0; j < sports.size(); j ++)
            System.out.print(String.format("%5.2f", means[j]));
        System.out.println();
        Distrib d2 = r.query(G);
        System.out.println("Prob of gender: " + d2);
        
        // re-training
        DirichletDistrib colours_male = new DirichletDistrib(colours, new double[] {3.0, 5.0, 7.0});
        DirichletDistrib colours_female = new DirichletDistrib(colours, new double[] {7.0, 2.0, 3.0});
        DirichletDistrib sports_male = new DirichletDistrib(sports,  new double[] {2.0, 5.0, 5.0});
        DirichletDistrib sports_female = new DirichletDistrib(sports,  new double[] {9.0, 4.0, 3.0});

        Object[][] data = new Object[200][3]; // data set
        for (int i = 0; i < 200; i ++) {
            if (i % 2 == 0) { // even so male
                if (i % 20 == 0) data[i][0] = "Male";
                EnumDistrib e1 = (EnumDistrib)colours_male.sample();
                data[i][1] = e1;
                EnumDistrib e2 = (EnumDistrib)sports_male.sample();
                data[i][2] = e2;
            } else { // female
                if (i % 20 == 1) data[i][0] = "Female";
                EnumDistrib e1 = (EnumDistrib)colours_female.sample();
                data[i][1] = e1;
                EnumDistrib e2 = (EnumDistrib)sports_female.sample();
                data[i][2] = e2;
            }
        }

        g.put(new EnumDistrib(new Enumerable(new String[] {"Male", "Female"}), 0.25, 0.75));
        c.put(new DirichletDistrib(colours, 1), "Male");
        c.put(new DirichletDistrib(colours, 1), "Female");
        s.put(new DirichletDistrib(sports,  1), "Male");
        s.put(new DirichletDistrib(sports,  1), "Female");

        EM em = new EM(bn);
        em.EM_MAX_ROUNDS = 1;
        em.EM_PRINT_STATUS = false;
        for (int round = 0; round < 20; round ++) {
            em.train(data, new Variable[] {G, C, S}, 0);
            c.setInstance(new EnumDistrib(colours, new double[] {0.3, 0.3, 0.4})); // primarily blue, but the gender node is latent...
            inf = new VarElim();
            inf.instantiate(bn);
            q = inf.makeQuery(G,S);
            r = (CGTable) inf.infer(q);
            d2 = r.query(G);
            System.out.println("Prob of gender: " + d2);
        }
        r.display();
        d1 = r.query(S);
        System.out.println("Prob of sports: " + d1);
    }
}
