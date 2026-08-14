package bn.ctmc;

import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.JTTGap;
import dat.Enumerable;
import bn.math.Matrix.Exp;

/**
 * This is a gap augmented conditional probability table for CTMC
 * based on discrete alphabets. This is based directly on the
 * paper by <a href="https://doi.org/10.1371/journal.pcbi.1000172">Eddy and Rivas, 2008</a>,
 * which describes a non-reversible generative (birth-death) evolutionary
 * model for insertions and deletions. Makes the key assumption that
 * indel events occur one residue at a time. When the insertion (λ) and
 * deletion (μ) rates are zero, the model will return the same values as the
 * standard CTMC table.
 *
 * @author Sebastian
 */
public class GapSubstModel extends SubstModel {

    final double mu; // deletion rate
    final double lambda; // insertion rate
    final double[] origF;

    public GapSubstModel(double[] F, double[][] IRM, Enumerable alphabet, double mu, double lambda,
                         boolean symmetric, boolean normalise, boolean copy) {

        super(F, IRM, alphabet, symmetric, normalise);
        this.mu = mu;
        this.lambda = lambda;
        this.origF = F;

        if (!copy) {
            Character[] gap_alphabet = addGapToAlphabet();
            this.alpha = new Enumerable(gap_alphabet);
            this.F = addGapToStationaryFreqs(F);

            double[][] R_EPS = new double[F.length + 1][F.length + 1];
            for (int i = 0; i < IRM.length; i ++) {
                if (IRM[i].length != F.length)
                    throw new IllegalArgumentException("IRM must be a square matrix");

                // first make the matrix symmetric
                if (symmetric) {
                    for (int j = i + 1; j < IRM[i].length; j++) {
                        double s = IRM[i][j];
                        R_EPS[i][j] = s * F[j];
                        R_EPS[j][i] = s * F[i];
                    }
                } else {
                    for (int j = 0; j < IRM[i].length; j ++) {
                        if (i == j)
                            continue;
                        R_EPS[i][j] = IRM[i][j];
                    }
                }
            }

            int numChars = F.length;
            // Applies R - μI (substitutions with deletions)
            for (int i = 0; i < numChars; i++) {
                // last column all μ (deletions to gap)
                R_EPS[i][numChars] = this.mu;

//                for (int j = 0; j < numChars; j++) {
//                    if (i == j) {
//                        R_EPS[i][j] -= this.mu;
//                    }
//                }
            }

            for (int j = 0; j < numChars; ++j) {
                R_EPS[numChars][j] = 0.0;
            }
            // sums bottom row to zero
            R_EPS[numChars][numChars] = 0.0;

            //double[][] R_EPS = constructIndelR(F, IRM);

            // make sum of rows == 0
            int dim = R_EPS.length;
            for (int i = 0; i < dim; i++) {
                double sum = 0.0;
                for (int j = 0; j < dim; j++) {
                    if (i != j)
                        sum += R_EPS[i][j];
                }
                R_EPS[i][i] = -sum;
            }

            // normalise
            if (normalise) {
                double sum = 0.0;
                for (int i = 0; i < dim; i++) {
                    sum += -R_EPS[i][i]*this.F[i];
                }
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++)
                        R_EPS[i][j] = R_EPS[i][j]/sum;
                }
            }

            this.R = R_EPS;
            this.Rexp = new Exp(R);
        }
    }


    public GapSubstModel(double[] F, double[][] IRM, Enumerable alphabet, double mu, double lambda) throws IllegalArgumentException{
        this(F, IRM, alphabet, mu, lambda, true, true, false);
    }


    public GapSubstModel deepCopy() {
        return new GapSubstModel(getF().clone(), getR().clone(), getDomain(), mu, lambda, false,
                false, true);
    }


    /**
     * Constructs the indel augmented rate matrix R_EPS
     * as per Eddy & Rivas (2008).
     * @return a double array of shape (numChars + 1) x (numChars + 1) to account for gaps
     */
    private double[][] constructIndelR(double[] F, double[][] IRM) {

        double[][] R_EPS = new double[F.length + 1][F.length + 1];
        for (int i = 0; i < IRM.length; i ++) {
            if (IRM[i].length != F.length)
                throw new IllegalArgumentException("IRM must be a square matrix");

            // first make the matrix symmetric
            for (int j = i + 1; j < IRM[i].length; j++) {
                double s = IRM[i][j];
                R_EPS[i][j] = s * F[j];
                R_EPS[j][i] = s * F[i];
            }
        }

        int numChars = F.length;
        // Applies R - μI (substitutions with deletions)
        for (int i = 0; i < numChars; i++) {
            // last column all μ (deletions to gap)
            R_EPS[i][numChars] = this.mu;

            for (int j = 0; j < numChars; j++) {
                if (i == j) {
                    R_EPS[i][j] -= this.mu;
                }
            }
        }

        //  Bottom row: λ * π_j (insertions from gap)
        for (int j = 0; j < numChars; ++j) {
            R_EPS[numChars][j] = F[j] * this.lambda;
        }
        // sums bottom row to zero
        R_EPS[numChars][numChars] = -this.lambda;

        return R_EPS;
    }

    /**
     * Adds gap character to the alphabet
     * @return Character array with gap added
     */
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


    public double getMu() {
        return this.mu;
    }

    public double getLambda() {
        return this.lambda;
    }

    @Override
    public String getName() {
        return "GapSubstModel";
    }

    @Override
    public double getProb(Object X) {

        // gaps are encoded by nulls in EnumSeq.Gappy<Enumerable>
        if (X == null) {
            X = '-';
        }

        int index_X = alpha.getIndex(X);
        return F[index_X];
    }


    /**
     * Function for the rate of deletion. For the case of
     * no insertions or deletions (mu=lambda=0) the function
     * is defined as zero.
     * (Eddy & Rivas, 2008) - Equation 6
     * @param t branch length
     * @return probability of deletion
     */
    public double gammaT(double time) {
        if (this.mu == 0 && this.lambda == 0) {
            return 0.0;
        }

        double deletionProp = this.mu / (this.mu + this.lambda);
        double indelProp = 1 - Math.exp(-((this.mu + this.lambda) * time));

        return deletionProp * indelProp;
    }

    /**
     * Function for the rate of insertion. For the case of
     * no insertions or deletions (mu=lambda=0) the function
     * is defined as zero.
     * (Eddy & Rivas, 2008) - Equation 6
     * @param time branch length
     * @return probability of insertion
     */
    public double ksiT(double time) {
        if (this.mu == 0 && this.lambda == 0) {
            return 0.0;
        }
        double insertionRate = this.lambda / (this.mu + this.lambda);
        double indelProp = 1 - Math.exp(-((this.mu + this.lambda) * time));

        return insertionRate * indelProp;
    }

    /**
     * Probability of observing some child residue, x, given the parent, y,
     * for a particular branch length, t. Re-weights the probability based on
     * observing this residue given it is not an insertion.
     * (Eddy & Rivas, 2008) - Equation 22: P(x|y,t) = (1 - ksi) * P(x|y)
     * @param child child state
     * @param parent parent state
     * @param time branch length
     * @return P(x|y,t)
     */
    public double getProbGapAugmented(Character child, Character parent, double time) {
        double probNoInsertion = 1 - ksiT(time);
        double conditionalProb = getProb(child, parent, time);

        return conditionalProb * probNoInsertion;
    }

    /**
     * Probability of a gap occurring given the parent
     * contained content, i. Is effectively the probability of
     * no insertion (1 - ksi) occurring jointly with a deletion event (gamma).
     * (Eddy & Rivas, 2008) - Equation 23
     * @param t branch length
     * @return P(-| i, t)
     */
    public double getProbOfGap(double t) {
        double probNoInsertion = 1 - ksiT(t);
        double probDeletion = gammaT(t);

        return probNoInsertion * probDeletion;
    }

    /**
     * Probability of an insertion occurring for a particular
     * character given the parent had no content. Weights the
     * stationary probability of that character by the insertion
     * probability function, ksi.
     * (Eddy & Rivas, 2008) - Equation 24. = ksi * pi
     * @param t branch length
     * @param state residue character
     * @return  P(j| -, t)
     */
    public double getProbOfInsertion(double time, Object state) {

        double insertionProb = ksiT(time);
        int index_X = alpha.getIndex(state);
        double stationaryFreqResidue = origF[index_X];
        return insertionProb * stationaryFreqResidue;
    }

    public static void main(String[] argv) {
        SubstModel JTT = new JTT();
        GapSubstModel JTT_Gap = new JTTGap(0.05, 0.05);
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);

        // Stationary probabilities
        System.out.println(JTT.getProb('A'));
        System.out.println(JTT_no_indel.getProb('A'));
        System.out.println(JTT_Gap.getProb('A'));

        // Conditional probabilities
        System.out.println(JTT.getProb('A', 'Y', 0.5));
        System.out.println(JTT_no_indel.getProb('A', 'Y', 0.5));
        System.out.println(JTT_Gap.getProb('A', 'Y', 0.5));

        // Gap probabilities
        System.out.println(JTT_Gap.getProb('-'));
        System.out.println(JTT_no_indel.getProb('-'));

        // Gap specific functions
        System.out.println(JTT_Gap.ksiT(0.05));
        System.out.println(JTT_Gap.gammaT(0.05));
        System.out.println(JTT_Gap.getProbOfGap(0.05));
        System.out.println(JTT_Gap.getProbOfInsertion(0.05, 'A'));
        System.out.println(JTT_Gap.getProbOfGap(0.05));
    }
}

