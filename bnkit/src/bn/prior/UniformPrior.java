package bn.prior;

import java.io.Serializable;
import java.util.Arrays;

import dat.Enumerable;

import bn.Distrib;
import bn.prob.EnumDistrib;

/**
 * used to calculate the MLE when the prior is uniform for all x
 * That is no prior
 * @author wangyufei
 *
 */
public class UniformPrior implements Prior, Serializable {

	private Distrib likelihoodDistrib; 
	
	@Override
	public void learn(Object[] data, double[] prob) {
		// FIXME bad design pattern
		if(likelihoodDistrib instanceof EnumDistrib) {
			EnumDistrib enumDistrib= (EnumDistrib)likelihoodDistrib;
			training(enumDistrib, data, prob);
		}

	}
	
	private void training(EnumDistrib enumDistrib, Object[] data, double[] prob) {
		Enumerable domain = enumDistrib.getDomain();
		double[] probs = new double[domain.size()];
		Arrays.fill(probs, 0.0);
		
		for(int i = 0; i < data.length; i++) {
			Object point = data[i];
			probs[domain.getIndex(point)] += prob[i];
		}
		
		enumDistrib.set(probs);
		
	}

	@Override
	public void setLikelihoodDistrib(Distrib distrib) {
		likelihoodDistrib = distrib;
	}

	@Override
	public Distrib getMAPDistrib() {
		return likelihoodDistrib;
	}

}
