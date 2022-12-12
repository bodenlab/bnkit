package bn.example.prior;

import bn.prior.GammaDistribPrior;
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
		GaussianDistrib guassian1 = new GaussianDistrib(1000,100);
		Double[] sample = new Double[nSample];
		double[] prob = new double[nSample];
		for(int i = 0; i < nSample; i++) {
			sample[i] = guassian1.sample();
			prob[i] = 1.0;
		}
		GuassianPriorTesting(sample, prob);
		GammaPriorTesting(sample, prob);
	}
	
	// testing with GuassianPrior, trying to get mean
	public static void GuassianPriorTesting(Double[] sample, double[] prob) {
		GaussianDistribPrior prior = GaussianDistribPrior.getUniformDistrib();
		prior.setEstimatedDistrib(new GaussianDistrib(0.1, 100));
		prior.learn(sample, prob);
		GaussianDistrib g2 = (GaussianDistrib) prior.getEstimatedDistrib();
		System.out.println(g2.toString());
	}
	
	// testing with GammaPrior, tryig to get Variance
	public static void GammaPriorTesting(Double[] sample, double[] prob) {
		GammaDistribPrior prior = GammaDistribPrior.getUniformDistrib();
		prior.setEstimatedDistrib(new GaussianDistrib(1000, 0.001));
		prior.learn(sample, prob);
		GaussianDistrib g2 = (GaussianDistrib) prior.getEstimatedDistrib();
		System.out.println(g2.toString());
		
	}

}
