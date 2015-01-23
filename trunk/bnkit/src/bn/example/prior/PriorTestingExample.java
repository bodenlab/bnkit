package bn.example.prior;

import bn.prior.GaussianDistribPrior;
import bn.prob.GaussianDistrib;

/**
 * class used to testing if our prior & uniform prior is right
 * @author wangyufei
 *
 */
public class PriorTestingExample {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int nSample = 2000;
		GaussianDistrib guassian1 = new GaussianDistrib(100,10);
		Double[] sample = new Double[nSample];
		double[] prob = new double[nSample];
		for(int i = 0; i < nSample; i++) {
			sample[i] = guassian1.sample();
			prob[i] = 1.0;
		}
		GaussianDistribPrior prior = GaussianDistribPrior.getUniformDistrib();
		prior.setLikelihoodDistrib(new GaussianDistrib(0.1, 10));
		prior.learn(sample, prob);
		GaussianDistrib g2 = (GaussianDistrib) prior.getBayesDistrib();
		System.out.println("mean: " + g2.getMean() + "  variance: " + g2.getVariance());
	}

}
