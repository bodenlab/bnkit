package bn.ctmc;

import dat.Enumerable;

public class PIPSubstModel extends SubstModel{

    double mu;
    double lambda;

    /**
     * Create time reversible evolutionary model.
     *
     * @param F        stationary base frequencies
     * @param S        Symmetric, un-scaled version of Q matrix Q_ij = s_ij*pi_j as defined by PAML
     * @param alphabet the values that substitutable variables can take, listed strictly in the order of the array and matrix
     */
    public PIPSubstModel(double[] F, double[][] S, Enumerable alphabet, double mu, double lambda) {
        super(F, S, alphabet);
        this.mu = mu;
        this.lambda = lambda;
    }

    public double getInsertionProb(double treeLength) {

        double factor = 1.0 / (treeLength + (1.0 / mu));
        return factor * (1.0 / mu);

    }

    public double getInsertionProb(double t, double treeLength) {

        double factor = 1.0 / (treeLength + (1.0 / mu));
        return factor * t;
    }

    public double survivalProb(double t, boolean isRoot) {

        if (isRoot) {
            return 1.0;
        } else {
            return (1.0 - Math.exp(-mu * t)) / (mu * t);
        }
    }

    public double pureSurvivalProb(double t) {
        return Math.exp(-mu * t);
    }

    @Override
    public double getProb(Object X, Object Y, double t) {
        if (X.equals('-')) {
            return mu;
        } else if (Y.equals('-')) {
            return 0.0;
        } else {
            return super.getProb(X, Y, t);
        }
    }

    @Override
    public String getName() {
        return "PIP";
    }
}
