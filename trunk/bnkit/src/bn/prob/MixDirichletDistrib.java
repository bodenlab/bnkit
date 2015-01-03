package bn.prob;

import java.io.Serializable;

import dat.Enumerable;

import bn.Distrib;
import bn.prior.DirichletDistribPrior;

/**
 * Mixture Dirichlet distribution is weighted sum of different Dirichlet distribution
 * MDir = w1 * Dir1 + w2 * Dir2 + ... 
 * This class includes learning parameters from data
 * @author wangyufei
 *
 */

public class MixDirichletDistrib  extends MixtureDistrib implements Serializable{

	/**
	 * construct a mixture Dirichlet model from a single component
	 * @param d1
	 * @param weight1
	 */
	public MixDirichletDistrib(DirichletDistrib d1, double weight1) {
		super(d1, weight1);
	}
	
	/**
	 * given the domain of Dirichlet distribution
	 * build an empty Mixture model
	 * @param domain
	 * @param ComponentNum
	 */
	public MixDirichletDistrib(Enumerable domain, int ComponentNum) {
		super(new DirichletDistrib(domain), 1.0);
		for(int i = 0; i < ComponentNum - 1; i++) {
			super.addDistrib(new DirichletDistrib(domain), 1.0);
		}
	}
	
	/**
	 * add either a Dirichlet distribution or mixDirichlet distribution
	 */
	public double addDistrib(Distrib d2, double weight2) {
		if(d2 instanceof DirichletDistrib || d2 instanceof MixDirichletDistrib) {
			return super.addDistrib(d2, weight2);
		}
		throw new RuntimeException("only accept DirichletDistrib or MixDirichletDistrib");
		
	}
	
	public void learnParameters(Object[] data, Double[] weight) {
		// TODO
	}

	

}
