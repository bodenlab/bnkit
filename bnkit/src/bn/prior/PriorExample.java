package bn.prior;

import bn.Predef;
import dat.EnumVariable;

public class PriorExample {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
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

}
