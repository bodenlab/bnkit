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

    public double getMu() {
        return mu;
    }

    public double getLambda() {
        return lambda;
    }

    public PIPSubstModel(double[] F, double[][] S, Enumerable alphabet, double mu, double lambda, boolean symmetric, boolean normalise) {
        super(F, S, alphabet, symmetric, normalise);
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
        if (Y.equals('-')) {
            // parent is gap — absorbing state
            return X.equals('-') ? 1.0 : 0.0;
        } else if (X.equals('-')) {
            // child is gap, parent is real — deletion occurred
            return (1.0 - Math.exp(-mu * t));

        } else {
            // both real characters — standard substitution
            return super.getProb(X, Y, t);
        }
    }

    /**
     * Get probability P(X=x)
     * @param X
     * @return
     */
    @Override
    public double getProb(Object X) {
        if (X.equals("-")) {
            return 0.0;
        } else {
            return super.getProb(X);
        }
    }

    @Override
    public String getName() {
        return "PIP";
    }
}
