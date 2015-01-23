package bn.prior;

import java.io.Serializable;

import bn.Distrib;

/**
 * interface for conjugate prior
 * Note: it is still on design stage, never use it before everything works all right
 * WORKING STATUS
 * @author wangyufei
 *
 */
public interface Prior extends Serializable {
	/**
	 * distribution can learn from the data
	 * and change its own parameters
	 * this is a process from prior to posterior
	 * @param data
	 */
	public void learn(Object[] data, double[] prob);
	
	/**
	 * set likelihood distribution
	 * @param distrib
	 */
	public void setLikelihoodDistrib(Distrib distrib);
	
	/**
	 * get the MAP result distribution
	 * @return
	 */
	public Distrib getBayesDistrib();
	
	/**
	 * reset the parameters to the initial value
	 * used in EM
	 */
	public void resetParameters();
	/**
	 * This is the interface to learn prior parameters
	 * from raw dataset
	 * @param data raw data
	 * @param prob count for that data
	 */
	public void learnPrior(Object[] data, double[] prob);
	
}
