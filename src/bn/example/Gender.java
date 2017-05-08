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
import bn.prob.EnumDistrib;
import bn.node.CPT;
import bn.node.GDT;
import dat.EnumVariable;
import dat.Variable;
import dat.Enumerable;
import bn.prob.GaussianDistrib;
import bn. *;
import bn.alg.*;

/**
 *
 * @author mikael
 */
public class Gender {
    
    public static void main(String[] args) {

        // Define variables
        EnumVariable G = Predef.Nominal(new String[] {"Male", "Female"}, "Gender");
        Variable H     = Predef.Real("Height");
        Variable W     = Predef.Real("Weight");

        // Define nodes (connecting the variables into an acyclic graph, i.e. the structure)
        CPT g = new CPT(G);
        GDT h = new GDT(H,    G);
        GDT w = new GDT(W,    G);

        // Parameterise the nodes using our "expertise"
        g.put(new EnumDistrib(new Enumerable(new String[] {"Male", "Female"}), 0.45, 0.55));
        h.put(new GaussianDistrib(180.0, 15.0), "Male");
        h.put(new GaussianDistrib(165.0, 15.0), "Female");
        w.put(new GaussianDistrib(80.0, 15.0), "Male");
        w.put(new GaussianDistrib(60.0, 15.0), "Female");
        
        BNet bn = new BNet();
        bn.add(g,h,w);
        h.setInstance(175.0);
        CGVarElim inf = new CGVarElim();
        inf.instantiate(bn);
        //ApproxInference inf = new ApproxInference();
        //inf.setIterations(2000);
        //inf.instantiate(bn);
        Query q = inf.makeQuery(G,W);
        CGTable r = (CGTable) inf.infer(q);
        r.display();
        Distrib d = r.query(W);
        double sum = 0;
        int nElem = 50;
        int minw = 50;
        int maxw = 100;
        double[] hist = new double[nElem];
        for (int i = 0; i < 1000; i ++) {
            double s = (Double)d.sample();
            if (s > minw)
                hist[(int)((s - minw) * ((maxw - minw) / nElem))] += 1;
            sum += s;
        }
        System.out.println("Avg weight = " + sum / 1000.0);
        for (int i = 0; i < hist.length; i ++) {
            System.out.println((i * ((maxw - minw) / nElem)) + minw + "kg\t" + hist[i]);
        }
    }
}
