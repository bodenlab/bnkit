package bn.prior;

import java.io.Serializable;
import java.util.Arrays;

import bn.Distrib;
import bn.prob.DirichletDistrib;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import dat.Enumerable;

public class DirichletDistribPrior extends DirichletDistrib implements Prior, Serializable {
	
    // currently, Dirichlet is only used as conjudge prior for EnumDistrib
	private EnumDistrib likelihoodDistrib;
	private double[] originalAlpha;
    
	
	public DirichletDistribPrior(Enumerable domain, double[] p, double m) {
        super(domain, p, m);
        likelihoodDistrib = null;
        originalAlpha = new double[p.length];
        for(int i = 0; i < p.length; i++) {
        	originalAlpha[i] = p[i] * m;
        }
    }
    
    
    private void setPosterior(double[] alpha) {
    	setPrior(alpha);
    }
    
    
    /**
     * learn from data
     * just add count and old alpha
     * the each item in prob is 
     * corresponding to each data point 
     */
	@Override
	public void learn(Object[] data, double[] prob) {
		Enumerable domain = (Enumerable) getDomain();
		double[] alpha = getAlpha();
		double[] countVector = new double[domain.size()];
		double[] newAlpha = new double[alpha.length];
		
		if(likelihoodDistrib == null) {
			System.err.println("likelihood distribution should be specificed");
			return;
		}
		Arrays.fill(countVector, 0.0);
		/**
		 * get the count vector for each variable
		 */
		for(int i = 0; i < data.length; i++) {
			Object point = data[i];
			countVector[domain.getIndex(point)] += prob[i];
		}
		/**
		 * get new alpha value
		 */
		for(int i = 0; i < alpha.length; i++) {
			newAlpha[i] = alpha[i] + countVector[i];
		}
		setPosterior(newAlpha);
		
	}


	@Override
	public void setLikelihoodDistrib(Distrib distrib) {
		try {
			likelihoodDistrib = (EnumDistrib) distrib;
		} catch(ClassCastException e) {
			System.out.println("the likelihood for Dirichlet prior should be enum distribution");
		}
	}

	/**
	 * MAP for Dirichlet Distribution
	 * Pi = (alphai - 1) / (sum(alpha) - K)
	 */
	@Override
	public Distrib getMAPDistrib() {
		Enumerable domain = (Enumerable) getDomain();
		double[] probs = new double[domain.size()];
		Arrays.fill(probs, 0.0);
		for(int i = 0; i < domain.size(); i++) {
			probs[i] = (getAlpha()[i] - 1) / (getSum() - domain.size());
		}
		likelihoodDistrib.set(probs);
		likelihoodDistrib.normalise();
		return likelihoodDistrib;
	}

	@Override
	public void resetParameters() {
		setPrior(originalAlpha);
	}



}
