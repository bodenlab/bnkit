package bn.prior;

import java.io.Serializable;
import java.util.Random;

import bn.Distrib;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;

public class GaussianDistribPrior extends GaussianDistrib implements Prior, Serializable{
    
	// currently, Gaussian is only used as conjudge prior for GaussianDistrib
	private GaussianDistrib likelihoodDistrib;
	
    public GaussianDistribPrior(double mean, double variance) {
    	super(mean, variance);
    	likelihoodDistrib = null;
    }

	@Override
	public void learn(Object[] data, double[] prob) {
		Double[] trainingData = (Double[]) data;
		double sum = 0.0;
		double posteriorMean;
		double posteriorVariance;
		
		if(likelihoodDistrib == null) {
			System.err.println("likelihood distribution should be specificed");
			return;
		}
		
		
		for(Double point: trainingData) {
			sum += point.doubleValue();
		}
		posteriorMean = (mu / sigmaSquared + sum / likelihoodDistrib.getVariance()) 
				/ ( 1 / sigmaSquared + trainingData.length / likelihoodDistrib.getVariance());
		posteriorVariance = 1 / ( 1 / sigmaSquared + trainingData.length / likelihoodDistrib.getVariance());
		setMean(posteriorMean);
		setVariance(posteriorVariance);
	}


	@Override
	public void setLikelihoodDistrib(Distrib distrib) {
		try {
			likelihoodDistrib = (GaussianDistrib) distrib;
		} catch(ClassCastException e) {
			System.out.println("the likelihood for Gaussian prior should be Guassian distribution");
		}
		
	}

	@Override
	public Distrib getMAPDistrib() {
		// the mode for Gaussian is mu
		likelihoodDistrib.mu = getMean();
		return likelihoodDistrib;
	}


}
