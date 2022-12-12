package bn.example.prior;

import dat.EnumVariable;
import bn.BNet;
import bn.Predef;
import bn.prior.CPTPrior;
import bn.prior.DirichletDistribPrior;

public class JulianExample {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		EnumVariable parentPosition = Predef.AminoAcid("parent");
		CPTPrior cpt = new CPTPrior(parentPosition);
		DirichletDistribPrior prior = 
				new DirichletDistribPrior(parentPosition.getDomain(), 
						new double[] {
					                  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
									  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
									  0.1,0.1,0.1,0.1
									  }, 
						20);
		EnumVariable childPosition = Predef.AminoAcid("child");
		CPTPrior cpt2 = new CPTPrior(childPosition);
		cpt.setPrior(prior);
		for(int i = 0; i < parentPosition.getDomain().size(); i++) {
			Object parentValue =  parentPosition.getDomain().get(i);
			cpt2.setPrior(new Object[] {parentValue}, prior);
		}
		
		BNet bn = new BNet();
		bn.add(cpt, cpt2);
		
		// training

	}

}
