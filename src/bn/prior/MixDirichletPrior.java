package bn.prior;

import java.util.Arrays;

import dat.Enumerable;
import bn.Distrib;
import bn.prob.DirichletDistrib;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import bn.prob.MixDirichletDistrib;

/**
 * This algorithm comes from 
 * Sj�lander, K., Karplus, K., e.l.(1996). 
 * Computer applications in the biosciences: CABIOS, 12(4), 327-345.
 * 
 * @Note this prior is only used for the estimation of enum distribution in BN
 *       To learn the parameters of this prior, learnPrior can be useful
 *       You don't have to use the method in MixDirichlet.java!
 *       If you use Mixture Dirichlet distribution for other purpose, you should
 *       use MixDirichletDistrib instead
 * @author wangyufei
 *
 */
public class MixDirichletPrior extends MixDirichletDistrib implements Prior {
	
	double[][] alpha;
	double[] m;
	private EnumDistrib likelihoodDistrib;
	private double[] countVector;
	private double scale;
	
	/**
	 * constructure a mixture Dirichlet Prior given domain and component size
	 * @param domain enumerable variable for each Dirichlet component
	 * @param component the number of components in the "mixture"
	 */
	public MixDirichletPrior(Enumerable domain, int component) {
		super(domain, component);
		scale = 1.0;
		alpha = new double[component][domain.size()];
		m = new double[component];
		countVector = new double[domain.size()];
		Arrays.fill(countVector, 0.0);
		for(int i = 0; i < component; i++) {
            DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i); 
            System.arraycopy(dirichlet.getAlpha(), 0, alpha[i], 0, domain.size());
            m[i] = this.getWeights(i);
		}	
	}
	
	/**
	 * construct a mixture Dirichlet Prior given one single Dirichlet and component number
	 * each component will be identical to the given Dirichlet distribution
	 * @param d1
	 * @param component
	 */
	public MixDirichletPrior(DirichletDistrib d1, int component) {
		super(d1, 1 / (float)component);
		Enumerable domain = (Enumerable)d1.getDomain();
		scale = 1.0;
		m = new double[component];
		alpha = new double[component][domain.size()];
		countVector = new double[domain.size()];
		Arrays.fill(countVector, 0.0);
		for(int i = 0; i < component; i++) {
			System.arraycopy(d1.getAlpha(), 0, alpha[i], 0, domain.size());
			m[i] = 1 / (float)component;
		}
		for(int i = 1; i < component; i++) {
			DirichletDistrib dirichlet = new DirichletDistrib(domain, d1.getAlpha().clone(), 1);
			this.addDistrib(dirichlet, m[i]);
		}
	}
	
	/**
	 * calculate the probability of the count vector of the training data
	 * given the different Dirichlet component represented
	 * by array of alpha value
	 * @param alpha
	 * @return the probability of the count vector
	 */
	private double probCountVector(double[] alpha) {
		double result = 0.0;
		double sumCountVector = 0.0;
		double sumAlpha = 0;
		if(countVector.length != alpha.length) {
			throw new RuntimeException("the length of count vector and alpha should be same");
		}
		for(int i = 0; i < countVector.length; i++) {
			sumCountVector += countVector[i];
			sumAlpha += alpha[i] * scale;
		}
		result = GammaDistrib.lgamma(sumCountVector + 1) + GammaDistrib.lgamma(sumAlpha) - GammaDistrib.lgamma(sumCountVector + sumAlpha);
		for(int i = 0; i < countVector.length; i++) {
			result += GammaDistrib.lgamma(countVector[i] + alpha[i] * scale) - GammaDistrib.lgamma(countVector[i] + 1) - GammaDistrib.lgamma(alpha[i] * scale);
		}
		return result;
	}

	/**
	 * learn data
	 */
	@Override
	public void learn(Object[] data, double[] prob) {
		Enumerable domain = getDomain();
		
		if(likelihoodDistrib == null) {
			System.err.println("likelihood distribution should be specificed");
			return;
		}
		/**
		 * get the count vector for each variable
		 */
		for(int i = 0; i < data.length; i++) {
			Object point = data[i];
			countVector[domain.getIndex(point)] += prob[i];
		}
	}

	@Override
	public void setEstimatedDistrib(Distrib distrib) {
		try {
			likelihoodDistrib = (EnumDistrib) distrib;
		} catch(ClassCastException e) {
			System.out.println("the likelihood for Mixture Dirichlet prior should be enum distribution");
		}

	}

	@Override
	public Distrib getEstimatedDistrib() {
		Enumerable domain = getDomain();
		double[] prob = new double[this.getMixtureSize()];
		double[] alphaSums = new double[this.getMixtureSize()];
		double[] dist = new double[domain.size()];
		double probSum = 0.0;
		double countSum = 0.0;
		
		for(int i = 0; i < domain.size(); i++) {
			countSum += countVector[i];
		}
		Arrays.fill(alphaSums, 0);
		for(int i = 0; i < this.getMixtureSize(); i++) {
			DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
			double[] alpha = dirichlet.getAlpha();
			prob[i] = this.getWeights(i) * Math.exp(probCountVector(alpha));
			probSum += prob[i];
			for(int j = 0; j < domain.size(); j++) {
				alphaSums[i] += alpha[j] * scale; 
			}
		}
		Arrays.fill(dist, 0);
		double componentProb = 0;
		for(int i = 0; i < domain.size(); i++) {
			for(int j = 0; j < this.getMixtureSize();j++) {
				DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(j);
				componentProb = prob[j] / probSum;
				componentProb *= ((this.countVector[i] + dirichlet.getAlpha()[i] * scale) / (countSum + alphaSums[j]));
				dist[i] += componentProb;
			}
		}
		likelihoodDistrib.set(dist);
		return likelihoodDistrib;
	}
	
	public void setAlphaScale(double NewScale) {
		scale = NewScale;
	}

	@Override
	public void resetParameters() {
		this.setWeights(m);
		for(int i = 0; i < this.getMixtureSize(); i++) {
			DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i); 
			dirichlet.setPrior(alpha[i]);
		}
		Arrays.fill(countVector, 0.0);
	}
	/**
	 * learn parameters mixture Dirichlet distribution from the data
	 * here, the probability of data are treated as 1  
	 */
	@Override
	public void learnPrior(Object[] data, double[] prob) {
		int[][] learningData = new int[data.length][];
		for(int i = 0; i < data.length; i++) {
			try {
				learningData[i] = (int[])data[i];
			} catch (ClassCastException e) {
				throw new ClassCastException("only accept int[][] data type");
			}
		}
		this.learnParameters(learningData);
		for(int i = 0; i < this.getMixtureSize(); i++) {
            DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i); 
            System.arraycopy(dirichlet.getAlpha(), 0, alpha[i], 0, this.getDomain().size());
            m[i] = this.getWeights(i);
		}
	}
	
	/**
	 * the same interface for learnPrior
	 * only learning data is necessary. 
	 * @param data
	 */
    public void learnPrior(Object[] data){
        this.learnPrior(data, new double[] {1});
    }
    
    public static MixDirichletPrior getUniformDistrib(Enumerable domain, int component) {
    	double[] alpha = new double[domain.size()];
		Arrays.fill(alpha, 1.0);
    	DirichletDistrib dirichlet = new DirichletDistrib(domain, alpha, 1);
    	MixDirichletPrior uniformPrior = new MixDirichletPrior(dirichlet, component);
    	return uniformPrior;
    }

}
