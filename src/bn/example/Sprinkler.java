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
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import bn.JPT;
import bn.Predef;
import bn.alg.ApproxInference;
import bn.alg.CGVarElim;
import bn.alg.Query;
import bn.alg.QueryResult;
import bn.file.BNBuf;

/**
 *
 * @author mikael
 */
public class Sprinkler {
        
        /**
         * Example from Russell and Norvig 2002 p 510.
         * @param args 
         */
	public static void main(String[] args) {
		EnumVariable CLOUDY = Predef.Boolean("Cloudy");
		EnumVariable SPRINKLER = Predef.Boolean("Sprinkler");
		EnumVariable RAIN = Predef.Boolean("Rain");
		EnumVariable WETGRASS = Predef.Boolean("WetGrass");
		
		CPT cloudy = new CPT(CLOUDY);
		cloudy.put(new EnumDistrib(Enumerable.bool, 0.5, 0.5));
		CPT sprinkler = new CPT(SPRINKLER, CLOUDY); // we think cloud cover may influence sprinkler activity
		sprinkler.put(new EnumDistrib(Enumerable.bool, 0.10, 0.90), true); 		
                sprinkler.put(new EnumDistrib(Enumerable.bool, 0.50, 0.50), false); 		
		CPT rain = new CPT(RAIN, CLOUDY);	   // clouds cause rain 	
		rain.put(new EnumDistrib(Enumerable.bool, 0.80, 0.20), true); 		
                rain.put(new EnumDistrib(Enumerable.bool, 0.20, 0.80), false); 		
                CPT wetGrass = new CPT(WETGRASS, SPRINKLER, RAIN); // wetness of grass is caused by sprinkler and rain
		wetGrass.put(new EnumDistrib(Enumerable.bool, 0.99, 0.01), true, true); 		
                wetGrass.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), true, false); 		
                wetGrass.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), false, true); 		
                wetGrass.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), false, false); 		
		// construct the network
		BNet bn=new BNet();
		bn.add(cloudy,sprinkler,rain,wetGrass);
		
		cloudy.print();
		rain.print();
		wetGrass.print(); // print the CPT

		sprinkler.setInstance(true);
		
//		CGVarElim inf = new CGVarElim();
                ApproxInference inf = new ApproxInference();
		inf.instantiate(bn);
                Query q = inf.makeQuery(RAIN);
                QueryResult qr = inf.infer(q);
		JPT jpt=qr.getJPT();
		jpt.display();
		
		BNBuf.save(bn, "data/sprinkler0.xml");
	}
	

}
