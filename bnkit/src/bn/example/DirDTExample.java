/*
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
import bn. *;
import bn.alg.*;

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
    }
}
