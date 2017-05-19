package bn.example;

import bn.BNet;
import bn.BNode;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import bn.node.GDT;
import bn.prob.GaussianDistrib;
import bn.JPT;
import bn.Predef;
import dat.Variable;
import bn.alg.ApproxInference;
import bn.alg.CGTable;
import bn.alg.CGVarElim;
import bn.alg.Query;
import bn.alg.QueryResult;
import bn.file.BNBuf;

public class GDTExample {

	public static void main(String[] args) {
		EnumVariable CLOUDY = Predef.Boolean("Cloudy");
		EnumVariable SPRINKLER = Predef.Boolean("Sprinkler");
		Variable HUMIDITY = Predef.Real("Humidity");
		EnumVariable RAIN = Predef.Boolean("Rain");
		EnumVariable SUN = Predef.Nominal(new String[] {"Sunny", "Partly sunny", "No sun"}, "Sun");
		EnumVariable WETGRASS = Predef.Boolean("WetGrass");
		
		CPT cloudy = new CPT(CLOUDY);
		cloudy.put(new EnumDistrib(Enumerable.bool, 0.5, 0.5));
		CPT sprinkler = new CPT(SPRINKLER, CLOUDY); // we think cloud cover may influence sprinkler activity
		sprinkler.put(new EnumDistrib(Enumerable.bool, 0.10, 0.90), true); 		
                sprinkler.put(new EnumDistrib(Enumerable.bool, 0.50, 0.50), false); 		
                GDT humidity = new GDT(HUMIDITY, RAIN); // we expect humidity to be caused by rain (not sprinklers) 
		humidity.put(new GaussianDistrib(70.0, 90.0), true);
		humidity.put(new GaussianDistrib(30.0, 90.0), false);
		CPT rain = new CPT(RAIN, CLOUDY);	   // clouds cause rain 	
		rain.put(new EnumDistrib(Enumerable.bool, 0.25, 0.75), true); 		
                rain.put(new EnumDistrib(Enumerable.bool, 0.00, 1.00), false); 		
                CPT sun = new CPT(SUN, CLOUDY);	   // clouds cause rain 	
		sun.put(new EnumDistrib(new Enumerable(new String[] {"Sunny", "Partly sunny", "No sun"}), 0.0, 0.25, 0.75), true);
                sun.put(new EnumDistrib(new Enumerable(new String[] {"Sunny", "Partly sunny", "No sun"}), 1.0, 0.0, 0.0), false);
                CPT wetGrass = new CPT(WETGRASS, SPRINKLER, RAIN); // wetness of grass is caused by sprinkler and rain
		wetGrass.put(new EnumDistrib(Enumerable.bool, 0.99, 0.01), true, true); 		
                wetGrass.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), true, false); 		
                wetGrass.put(new EnumDistrib(Enumerable.bool, 0.90, 0.10), false, true); 		
                wetGrass.put(new EnumDistrib(Enumerable.bool, 0.01, 0.99), false, false); 		
		// construct the network
		BNet bn=new BNet();
		bn.add(cloudy);
		bn.add(sprinkler);
		bn.add(rain);
		bn.add(sun);
		bn.add(humidity);
		bn.add(wetGrass);
		
		cloudy.print();
		rain.print();
		humidity.print();
		wetGrass.print(); // print the CPT

		humidity.setInstance(55.0);
		sun.setInstance("Partly sunny");
		sprinkler.setInstance(false);
		
		CGVarElim inf1 = new CGVarElim();
		inf1.instantiate(bn);
                Query q1 = inf1.makeQuery(WETGRASS);
                CGTable qr1 = (CGTable)inf1.infer(q1);
		qr1.display();
                
                ApproxInference inf2 = new ApproxInference();
                inf2.setIterations(10000);
		inf2.instantiate(bn);
                Query q2 = inf2.makeQuery(WETGRASS);
                CGTable qr2 = (CGTable)inf2.infer(q2);
		qr2.display();
		
		BNBuf.save(bn, "data/sprinkler1.xml");
	}
	
	public static void main2(String[] args) {
		BNet bn = BNBuf.load("data/sprinkler1.xml");
		
		BNode cloudy = bn.getNode("Cloudy");
		BNode sprinkler = bn.getNode("Sprinkler");
		BNode humidity = bn.getNode("Humidity");
		BNode wetgrass = bn.getNode("WetGrass");
		BNode sun = bn.getNode("Sun");
		
		humidity.setInstance(55.0);
		sun.setInstance("Partly sunny");
		//cloudy.setInstance(true);
		sprinkler.setInstance(false);
		
		CGVarElim ve = new CGVarElim();
		ve.instantiate(bn);
		JPT jpt=ve.infer(wetgrass).getJPT();
		jpt.display();
		
	}

	

}
