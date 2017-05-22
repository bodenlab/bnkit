package bn.prior;

import java.io.Serializable;

import bn.Distrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;

/**
 * Gamma distribution is used as prior of Gaussian distribution
 * here, Gamma distribution is described as parameters "shape" k and 
 * parameters scale Lambda.
 * 
 * Gamma prior is useful for likelihood of Guassian distribution with 
 * unknown variance. All Prior uses the same interface so please @see bn.example.prior for
 * more information
 * @author wangyufei
 *
 */
public class GammaDistribPrior extends GammaDistrib implements Prior {

	// currently, GammaDistrib is only used as conjudge prior for GaussianDistrib
	private GaussianDistrib likelihoodDistrib;
	private double oldK;
	private double oldLambda;
	
	/**
	 * 
	 * @param k, describes shape
	 * @param lambda, describes scale
	 */
	public GammaDistribPrior(double k, double lambda) {
		super(k, lambda);
		likelihoodDistrib = null;
		oldK = k;
		oldLambda = lambda;
	}
	
	/**
     * MAP algorithm for Gaussian Distribution with known mean
     */
	@Override
	public void learn(Object[] data, double[] prob) {
		if(likelihoodDistrib == null) {
			System.err.println("likelihood distribution should be specificed");
			return;
		}
		Double[] training = (Double[]) data;
		
		// learn K
		setK(getK() + ((double)data.length) / 2);
		
		// learn lambda
		double dataVariance = 0.0;
		double mean = likelihoodDistrib.getMean();
		for(Double point: training) {
			dataVariance += Math.pow(point - mean, 2);
		}
		setLambda(getLambda() + dataVariance / 2);
	}

	@Override
	public void setEstimatedDistrib(Distrib distrib) {
		try {
			likelihoodDistrib = (GaussianDistrib) distrib;
		} catch(ClassCastException e) {
			System.out.println("the likelihood for Gamma prior should be Guassian distribution");
		}
		
	}

	@Override
	public Distrib getEstimatedDistrib() {
		double precision = (getK() - 1) * (1 / getLambda());
		likelihoodDistrib.sigmaSquared = 1 / precision;
		likelihoodDistrib.sigma = Math.sqrt(likelihoodDistrib.sigmaSquared);
		return likelihoodDistrib;
	}

	@Override
	public void resetParameters() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void learnPrior(Object[] data, double[] prob) {
		// TODO Auto-generated method stub	
	}
	
	public static GammaDistribPrior getUniformDistrib() {
		return new GammaDistribPrior(0.0, 0.0);
	}



}
