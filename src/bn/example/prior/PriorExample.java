package bn.example.prior;

import java.util.Random;

import bn.BNet;
import bn.Predef;
import bn.alg.CGTable;
import bn.alg.EM;
import bn.alg.MAP;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.node.CPT;
import bn.prior.CPTPrior;
import bn.prior.DirichletDistribPrior;
import bn.prior.PriorBNode;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;

public class PriorExample {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Double a = Double.POSITIVE_INFINITY;
		System.out.println( 1 / a * 1);
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
	public static void complexExample(Boolean[][] data) {
		EnumVariable sun = Predef.Boolean("Sunrise");
		CPTPrior cpt = new CPTPrior(sun);
		DirichletDistribPrior betaDistrib = new DirichletDistribPrior(sun.getDomain(), new double[] {0.7,0.3}, 2);
		cpt.setPrior(betaDistrib);
		
		EnumVariable rain = Predef.Boolean("rain");
		CPTPrior cpt2 = new CPTPrior(rain, sun);
		// prior for sun rise. 90% it won't rain
		DirichletDistribPrior sunRiseBetaDistrib = new DirichletDistribPrior(rain.getDomain(), new double[] {0.1,0.9}, 1);
		// prior for sun not rise. 60% it will rain
		DirichletDistribPrior sunNotRiseBetaDistrib = new DirichletDistribPrior(rain.getDomain(), new double[] {0.6,0.4}, 1);
		
		/**
		 * Node that when set prior, the order of parent value should 
		 * be the same as the one used in constructor of CPTPrior
		 */
		cpt2.setPrior(new Object[] {true}, sunRiseBetaDistrib);
		//cpt2.setPrior(new Object[] {false}, sunNotRiseBetaDistrib);
		BNet bn = new BNet();
		bn.add(cpt, cpt2);
		
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
	
	public static void newSample(Boolean[][] data) {
		EnumVariable sun = Predef.Boolean("Sunrise");
		CPT cpt = new CPT(sun);
		PriorBNode priorNode = new PriorBNode(cpt);
		DirichletDistribPrior betaDistrib = new DirichletDistribPrior(sun.getDomain(), new double[] {0.7,0.3}, 2);
		priorNode.setPrior(betaDistrib);
		
		EnumVariable rain = Predef.Boolean("rain");
		CPT cpt2 = new CPT(rain, sun);
		PriorBNode priorNode2 = new PriorBNode(cpt2);
		// prior for sun rise. 90% it won't rain
		DirichletDistribPrior sunRiseBetaDistrib = new DirichletDistribPrior(rain.getDomain(), new double[] {0.1,0.9}, 1);
		// prior for sun not rise. 60% it will rain
		DirichletDistribPrior sunNotRiseBetaDistrib = new DirichletDistribPrior(rain.getDomain(), new double[] {0.6,0.4}, 1);
		DirichletDistribPrior uni = new DirichletDistribPrior(rain.getDomain(), new double[] {1, 1}, 1);
		//set the uniform distribution, used when there is no prior privoded for some condition
		priorNode2.setUniformPrior(DirichletDistribPrior.getUniformDistrib(rain.getDomain()));
		/**
		 * Node that when set prior, the order of parent value should 
		 * be the same as the one used in constructor of CPTPrior
		 */
		priorNode2.setPrior(new Object[] {true}, sunRiseBetaDistrib);
		//priorNode2.setPrior(new Object[] {false}, sunNotRiseBetaDistrib);
		BNet bn = new BNet();
		bn.add(priorNode, priorNode2);
		
		// learn data through EM
		
		
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
	
	public static Boolean[][] getDataSet() {
		Boolean[][] data = new Boolean[100][];
		Random rn = new Random();
		for(int i = 0; i < data.length; i++) {
			Boolean r1 = (rn.nextInt() % 2 == 0? true : false);
			Boolean r2 = (rn.nextInt() % 2 == 0? true : false);
			data[i] = new Boolean[] {r1 , r2};
		}
		return data;
	}
	
	public static void sampleSize() {
		EnumVariable s1 = Predef.Boolean("s1");
		CPTPrior cpt = new CPTPrior(s1);
		DirichletDistribPrior betaDistrib = new DirichletDistribPrior(s1.getDomain(), 1000);
		cpt.setPrior(betaDistrib);
		
		EnumVariable s2 = Predef.Boolean("s2");
		CPTPrior cpt2 = new CPTPrior(s2, s1);
		DirichletDistribPrior BetaDistrib1 = new DirichletDistribPrior(s2.getDomain(), 1000);
		DirichletDistribPrior BetaDistrib2 = new DirichletDistribPrior(s2.getDomain(), 1000);
		cpt2.setPrior(new Object[] {true}, BetaDistrib1);
		cpt2.setPrior(new Object[] {false}, BetaDistrib2);
		int n = 5000;
		Boolean[][] data = new Boolean[8 * n][];
		
		for(int i = 0; i < n; i++){
			data[0 + i * 8] = new Boolean[] {true, false};
			data[1 + i * 8] = new Boolean[] {true, true};
			data[2 + i * 8] = new Boolean[] {false, true};
			data[3 + i * 8] = new Boolean[] {false, false};
			data[4 + i * 8] = new Boolean[] {false, true};
			data[5 + i * 8] = new Boolean[] {false, true};
			data[6 + i * 8] = new Boolean[] {true, false};
			data[7 + i * 8] = new Boolean[] {false, false};
		}
		/*
		data[0] = new Boolean[] {false, true};
		data[1] = new Boolean[] {false, false};
		data[2] = new Boolean[] {false, true};
		data[3] = new Boolean[] {true, true};
		*/
		BNet bn = new BNet();
		bn.add(cpt, cpt2);
		EM em = new EM(bn);
		em.train(data, new Variable[] {s1, s2}, 1);
		VarElim ve = new VarElim();
        ve.instantiate(bn);
		Query q1 = ve.makeQuery(s2);
        CGTable r1 = (CGTable)ve.infer(q1);
        r1.display();
	}

}
