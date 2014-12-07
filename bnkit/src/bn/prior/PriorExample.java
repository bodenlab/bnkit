package bn.prior;

import bn.BNet;
import bn.Predef;
import bn.alg.CGTable;
import bn.alg.EM;
import bn.alg.MAP;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;

public class PriorExample {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		complexExample();
	}
	
	/**
	 * This example is for a single node
	 */
	public static void rootExample() {
		EnumVariable sun = Predef.Boolean("Sunrise");
		CPTPrior cpt = new CPTPrior(sun);
		DirichletDistribPrior betaDistrib = new DirichletDistribPrior(sun.getDomain(), new double[] {0.7,0.3}, 2);
		cpt.setPrior(betaDistrib);
		
		System.out.println(betaDistrib.toString());
		cpt.countInstance(null, true);
		cpt.countInstance(null, true);
		cpt.countInstance(null, false);
		cpt.countInstance(null, false);
		cpt.maximizeInstance();
		System.out.println(cpt.getDistrib().toString());
	}
	
	/**
	 * This example has two nodes, sun rise node and rain node.
	 * each of them has prior.
	 */
	public static void complexExample() {
		EnumVariable sun = Predef.Boolean("Sunrise");
		CPTPrior cpt = new CPTPrior(sun);
		DirichletDistribPrior betaDistrib = new DirichletDistribPrior(sun.getDomain(), new double[] {0.7,0.3}, 2);
		cpt.setPrior(betaDistrib);
		
		EnumVariable rain = Predef.Boolean("rain");
		CPTPrior cpt2 = new CPTPrior(rain, sun);
		// prior for sun rise. 90% it won't rain
		DirichletDistribPrior sunRiseBetaDistrib = new DirichletDistribPrior(sun.getDomain(), new double[] {0.1,0.9}, 1);
		// prior for sun not rise. 60% it will rain
		DirichletDistribPrior sunNotRiseBetaDistrib = new DirichletDistribPrior(sun.getDomain(), new double[] {0.6,0.4}, 1);
		
		/**
		 * Node that when set prior, the order of parent value should 
		 * be the same as the one used in constructor of CPTPrior
		 */
		cpt2.setPrior(new Object[] {true}, sunRiseBetaDistrib);
		cpt2.setPrior(new Object[] {false}, sunNotRiseBetaDistrib);
		BNet bn = new BNet();
		bn.add(cpt, cpt2);
		
		// learn data through EM
		Boolean[][] data = new Boolean[10][];
		data[0] = new Boolean[] {true, false};
		data[1] = new Boolean[] {true, false};
		data[2] = new Boolean[] {true, true};
		data[3] = new Boolean[] {true, true};
		data[4] = new Boolean[] {true, true};
		data[5] = new Boolean[] {false, false};
		data[6] = new Boolean[] {false, false};
		data[7] = new Boolean[] {false, true};
		data[8] = new Boolean[] {false, true};
		data[9] = new Boolean[] {false, true};
		
		/**
		 * MAP can be applied when data is fully observed (No missing label)
		 * if no prior is given, MAP would treat the prior as uniformPrior
		 * This is when ML would be applied. 
		 */
		MAP map = new MAP(bn);
		map.train(data, new Variable[] {sun, rain}, 1);
		
		VarElim ve = new VarElim();
        ve.instantiate(bn);
        
        cpt2.setInstance(true);
        Query q1 = ve.makeQuery(sun);
        CGTable r1 = (CGTable)ve.infer(q1);
        r1.display();
        
        cpt2.setInstance(false);
        Query q2 = ve.makeQuery(sun);
        CGTable r2 = (CGTable)ve.infer(q2);
        r2.display();
	}

}
