package bn.prior;

import java.util.Arrays;

import bn.Distrib;
import bn.prob.DirichletDistrib;
import bn.prob.EnumDistrib;
import dat.Enumerable;

/**
 * Dirichlet distribution prior is useful for likelihood of Enum distribution.
 * All priors use the same interface so please @see bn.example.prior for more
 * information.
 *
 * @author wangyufei
 *
 */
public class DirichletDistribPrior extends DirichletDistrib implements Prior {

    // currently, Dirichlet is only used as conjudge prior for EnumDistrib
    private EnumDistrib likelihoodDistrib;
    private double[] originalAlpha;
    private double scale;

    /**
     * give the enumerable variable for enum distribution and the alpha
     *
     * @param domain enumerable variable
     * @param same_alpha this value will be assigned to all alpha of Dirichlet
     * distribution
     */
    public DirichletDistribPrior(Enumerable domain, double same_alpha) {
        super(domain, same_alpha);
        likelihoodDistrib = null;
        scale = 1.0;
        originalAlpha = new double[domain.size()];
        for (int i = 0; i < domain.size(); i++) {
            originalAlpha[i] = same_alpha;
        }
    }

    /**
     * do have to give the parameters Used when training for prior is necessary
     *
     * @param domain
     */
    public DirichletDistribPrior(Enumerable domain) {
        super(domain, 0.0);
        scale = 1.0;
        likelihoodDistrib = null;
        originalAlpha = new double[domain.size()];
        for (int i = 0; i < domain.size(); i++) {
            originalAlpha[i] = 0.0;
        }
    }

    /**
     *
     * @param domain the enumerable for enum distribution
     * @param p array of alpha values of Dirichlet distribution
     * @param m the scaling variable for the array above
     */
    public DirichletDistribPrior(Enumerable domain, double[] p, double m) {
        super(domain, p, m);
        likelihoodDistrib = null;
        scale = 1.0;
        originalAlpha = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            originalAlpha[i] = p[i] * m;
        }
    }

    private void setPosterior(double[] alpha) {
        setPrior(alpha);
    }

    /**
     * learn from data just add count and old alpha the each item in prob is
     * corresponding to each data point
     */
    @Override
    public void learn(Object[] data, double[] prob) {
        Enumerable domain = (Enumerable) getDomain();
        double[] alpha = getAlpha();
        double[] countVector = new double[domain.size()];
        double[] newAlpha = new double[alpha.length];

        if (likelihoodDistrib == null) {
            System.err.println("likelihood distribution should be specificed");
            return;
        }
        Arrays.fill(countVector, 0.0);
        /**
         * get the count vector for each variable
         */
        for (int i = 0; i < data.length; i++) {
            Object point = data[i];
            countVector[domain.getIndex(point)] += prob[i];
        }
        /**
         * get new alpha value
         */
        for (int i = 0; i < alpha.length; i++) {
            newAlpha[i] = alpha[i] + countVector[i];
        }
        setPosterior(newAlpha);

    }

    @Override
    public void setEstimatedDistrib(Distrib distrib) {
        try {
            likelihoodDistrib = (EnumDistrib) distrib;
        } catch (ClassCastException e) {
            System.out.println("the likelihood for Dirichlet prior should be enum distribution");
        }
    }

    /**
     * MAP for Dirichlet Distribution Pi = (alphai - 1) / (sum(alpha) - K)
     */
    @Override
    public Distrib getEstimatedDistrib() {
        Enumerable domain = (Enumerable) getDomain();
        double[] probs = new double[domain.size()];
        Arrays.fill(probs, 0.0);
        for (int i = 0; i < domain.size(); i++) {
            probs[i] = getAlpha()[i] / getSum();
        }
        likelihoodDistrib.set(probs);
        likelihoodDistrib.normalise();
        return likelihoodDistrib;
    }

    @Override
    public void resetParameters() {
        setPrior(originalAlpha);
    }

    /**
     * learn parameters of Dirichlet distribution from the dataset
     */
    @Override
    public void learnPrior(Object[] data, double[] prob) {
        int[][] rawData = new int[data.length][];
        for (int i = 0; i < data.length; i++) {
            Integer[] vector = (Integer[]) data[i];
            for (int j = 0; j < vector.length; j++) {
                rawData[i][j] = vector[j].intValue();
            }
        }
        double[] alpha = DirichletDistrib.getAlpha(rawData, prob);
        setPrior(alpha);
    }

    /**
     *
     * @param domain the enum distribution
     * @return the uniform Distribution of Dirichlet distribution
     */
    public static DirichletDistribPrior getUniformDistrib(Enumerable domain) {
        double[] alpha = new double[domain.size()];
        Arrays.fill(alpha, 1.0);
        return new DirichletDistribPrior(domain, alpha, 1.0);
    }

}
