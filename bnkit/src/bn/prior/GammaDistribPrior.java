package bn.prior;

import java.io.Serializable;

import bn.Distrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;


public class GammaDistribPrior extends GammaDistrib implements Prior, Serializable{

	// currently, GammaDistrib is only used as conjudge prior for GaussianDistrib
	private GaussianDistrib likelihoodDistrib;
	private double oldK;
	private double oldLambda;
	
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
		double oldK = getK();
		setK(oldK + ((double)data.length) / 2);
		
		// learn lambda
		double dataVariance = 0.0;
		double mean = likelihoodDistrib.getMean();
		for(Double point: training) {
			dataVariance += Math.pow(point - mean, 2);
		}
		setLambda(getLambda() + dataVariance / 2);
	}

	@Override
	public void setLikelihoodDistrib(Distrib distrib) {
		try {
			likelihoodDistrib = (GaussianDistrib) distrib;
		} catch(ClassCastException e) {
			System.out.println("the likelihood for Gamma prior should be Guassian distribution");
		}
		
	}

	@Override
	public Distrib getMAPDistrib() {
		double precision = (getK() - 1) * getLambda();
		likelihoodDistrib.sigmaSquared = 1 / precision;
		likelihoodDistrib.sigma = Math.sqrt(likelihoodDistrib.sigmaSquared);
		return likelihoodDistrib;
	}

	@Override
	public void resetParameters() {
		// TODO Auto-generated method stub
		
	}



}
