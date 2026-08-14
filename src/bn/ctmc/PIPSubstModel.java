package bn.ctmc;

import bn.math.Matrix;
import dat.Enumerable;

public class PIPSubstModel extends SubstModel{

    double mu;
    double lambda;
    double[] origF;
    private final Enumerable originalAlpha;

    /**
     * Create time reversible evolutionary model.
     *
     * @param F        stationary base frequencies
     * @param S        Symmetric, un-scaled version of Q matrix Q_ij = s_ij*pi_j as defined by PAML
     * @param alphabet the values that substitutable variables can take, listed strictly in the order of the array and matrix
     */
//    public PIPSubstModel(double[] F, double[][] S, Enumerable alphabet, double mu, double lambda) {
//        super(F, S, alphabet);
//        this.mu = mu;
//        this.lambda = lambda;
//    }

    public PIPSubstModel(double[] F, double[][] IRM, Enumerable alphabet, double mu, double lambda,
                         boolean symmetric, boolean normalise, boolean copy) {

        super(F, IRM, alphabet, symmetric, normalise);
        this.mu = mu;
        this.lambda = lambda;
        this.origF = F;
        this.originalAlpha = alphabet;

        if (!copy) {
            Character[] gap_alphabet = addGapToAlphabet();
            this.alpha = new Enumerable(gap_alphabet);
            this.F = addGapToStationaryFreqs(F);

            int numChars = F.length;
            double[][] R_EPS = new double[numChars + 1][numChars + 1];

            // copy normalised substitution rates from this.R (set by super)
            for (int i = 0; i < numChars; i++) {
                for (int j = 0; j < numChars; j++) {
                    R_EPS[i][j] = this.R[i][j];
                }
            }

            for (int i = 0; i < numChars; i++) {
                R_EPS[i][numChars] = mu;
                R_EPS[i][i] -= mu;
            }

            for (int j = 0; j < numChars + 1; j++) {
                R_EPS[numChars][j] = 0.0;
            }

            this.R = R_EPS;
            this.Rexp = new Matrix.Exp(R);
        }
    }

    public Enumerable getOriginalAlpha() {
        return originalAlpha;
    }

    private Character[] addGapToAlphabet() {
        Character[] gapAlphabet = new Character[alpha.size() + 1];
        for (int i = 0; i < alpha.size(); i++) {
            gapAlphabet[i] = (Character) alpha.get(i);
        }
        gapAlphabet[alpha.size()] = '-';

        return gapAlphabet;
    }

    /**
     * Adjusts stationary frequencies to account for gaps
     *
     *  @return array of modified stationary frequencies with gap stationary frequency added
     */
    private double[] addGapToStationaryFreqs(double[] F) {

        double[] fGap = new double[alpha.size()];

        if (lambda < 0 || mu < 0) {
            throw new IllegalArgumentException("mu + lambda must be >= 0");
        }
        // no indels - zero prob of gaps
        for (int i = 0; i < alpha.size() - 1; i++) {
            fGap[i] = F[i];
        }
        fGap[alpha.size() - 1] = 0.0;

        return fGap;
    }


    public PIPSubstModel(double[] F, double[][] IRM, Enumerable alphabet, double mu, double lambda) throws IllegalArgumentException{
        this(F, IRM, alphabet, mu, lambda, true, true, false);
    }



    public double getMu() {
        return mu;
    }

    public double getLambda() {
        return lambda;
    }

    public PIPSubstModel(double[] F, double[][] S, Enumerable alphabet,
                         double mu, double lambda, boolean symmetric, boolean normalise) {
        this(F, S, alphabet, mu, lambda, symmetric, normalise, false);
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
//        } else if (X.equals('-')) {
//            // child is gap, parent is real — deletion occurred
//            return (1.0 - Math.exp(-mu * t));
        } else {
            // both real characters — standard substitution
            return super.getProb(X, Y, t);
        }
    }


    @Override
    public String getName() {
        return "PIP";
    }
}
