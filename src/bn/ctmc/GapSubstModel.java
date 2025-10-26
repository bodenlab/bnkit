package bn.ctmc;

import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.JTTGap;
import dat.Enumerable;
import bn.math.Matrix.Exp;

/**
 * This is a gap augmented conditional probability table for CTMC
 * based on discrete alphabets. This is based directly on the
 * paper by Eddy and Rivas (https://doi.org/10.1371/journal.pcbi.1000172)
 * which describes a non-reversible generative (birth-dirth) evolutionary
 * model for insertions and deletions. Makes the key assumption that
 * indel events occur one residue at a time. When the insertion (lambda) and
 * deletion (mu) rates are zero, the model will return the same values as the
 * standard CTMC table.
 *
 * @author Sebastian
 */
public class GapSubstModel extends SubstModel {

    final double mu; // deletion rate
    final double lambda; // insertion rate

    public GapSubstModel(double[] F, double[][] S, Enumerable alphabet, double mu, double lambda)
            throws IllegalArgumentException{
        // perform normal setup
        super(F, S, alphabet);
        this.mu = mu;
        this.lambda = lambda;

        double[][] R_EPS = new double[F.length + 1][F.length + 1];
        int K = alphabet.size();
        // Top-left block: R - μI (substitutions with deletions)
        for (int i = 0; i < K; i++) {
            // Top-right column: μ (deletions to gap)
            R_EPS[i][K] = mu;

            for (int j = 0; j < K; j++) {
                R_EPS[i][j] = this.R[i][j];
                if (i == j) {
                    R_EPS[i][j] -= mu;
                }
            }
        }

        //  Bottom-left row: λ * π_j (insertions from gap)
        for (int j = 0; j < K; ++j) {
            R_EPS[K][j] = F[j] * lambda;
        }
        R_EPS[K][K] = -lambda;

        this.R = R_EPS;
        this.Rexp = new Exp(R_EPS);

        // update the alphabet to have an extra character
        Character[] gap_alphabet = new Character[alphabet.size() + 1];
        for (int i = 0; i < alphabet.size(); i++) {
            gap_alphabet[i] = (Character) alphabet.get(i);
        }
        gap_alphabet[alphabet.size()] = '-';
        this.alpha = new Enumerable(gap_alphabet);
        double[] F_gap = new double[alpha.size()];
        if (mu + lambda < 0) {
            throw new IllegalArgumentException("mu + lambda must be >= 0");
        } else if (mu + lambda > 0) {
            // need to adjust stationary freqs by deletion ratio
            // first just adjust real letters
            double gap_freq = mu / (mu + lambda);
            for (int i = 0; i < alpha.size() - 1; i++) {
                F_gap[i] = (1 - gap_freq) * F[i];
            }
            // stationary prob for an indel
            F_gap[alpha.size() - 1] = gap_freq;
        } else if (mu + lambda == 0) {
            // no indels - zero prob of gaps
            for (int i = 0; i < alpha.size() - 1; i++) {
                F_gap[i] = F[i];
            }
            F_gap[alpha.size() - 1] = 0.0;
        } else {
            throw new IllegalArgumentException("mu and lambda combination not supported");
        }
        this.F = F_gap;

    }

    @Override
    public String getName() {
        return "GapSubstModel";
    }

    /**
     * Function for the rate of deletion. For the case of
     * no insertions or deletions (mu=lambda=0) the function
     * is defined as zero.
     * (Eddy & Rivas, 2008) - Equation 6
     * @param t branch length
     * @return probability of deletion
     */
    public double gamma_t(double t) {
        if (this.mu == 0 && this.lambda == 0) {
            return 0.0;
        }

        double deletion_prop = this.mu / (this.mu + this.lambda);
        double indel_prop = 1 - Math.exp(-((this.mu + this.lambda) * t));

        return deletion_prop * indel_prop;
    }

    /**
     * Function for the rate of insertion. For the case of
     * no insertions or deletions (mu=lambda=0) the function
     * is defined as zero.
     * (Eddy & Rivas, 2008) - Equation 6
     * @param t branch length
     * @return probability of insertion
     */
    public double ksi_t(double t) {
        if (this.mu == 0 && this.lambda == 0) {
            return 0.0;
        }
        double insertion_rate = this.lambda / (this.mu + this.lambda);
        double indel_prop = 1 - Math.exp(-((this.mu + this.lambda) * t));

        return insertion_rate * indel_prop;
    }

    /**
     * Probability of observing some child residue, x, given the parent, y,
     * for a particular branch length, t. Re-weights the probability based on
     * observing this residue given it is not an insertion.
     * (Eddy & Rivas, 2008) - Equation 22: P(x|y,t) = (1 - ksi) * P(x|y)
     * @param child_x child state
     * @param parent_y parent state
     * @param t branch length
     * @return P(x|y,t)
     */
    public double prob_j_given_i_t(Character child_x, Character parent_y, double t) {
        double prob_no_insertion = 1 - ksi_t(t);
        double prob_j_given_i = getProb(child_x, parent_y, t);

        return prob_j_given_i * prob_no_insertion;
    }

    /**
     * Probability of a gap occurring given the parent
     * contained content, i. Is effectively the probability of
     * no insertion (1 - ksi) occurring jointly with a deletion event (gamma).
     * (Eddy & Rivas, 2008) - Equation 23
     * @param t branch length
     * @return P(-| i, t)
     */
    public double prob_gap_given_i_t(double t) {
        double prob_no_insertion = 1 - ksi_t(t);
        double prob_deletion = gamma_t(t);

        return prob_no_insertion * prob_deletion;
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
    public double prob_j_given_gap_t(double t, Character state) {
        double insertion_prob = ksi_t(t);
        double stationary_freq_residue = getProb(state);

        return insertion_prob * stationary_freq_residue;
    }

    public static void main(String[] argv) {
        SubstModel JTT = new JTT();
        GapSubstModel JTT_Gap = new JTTGap(0.05, 0.05);
        GapSubstModel JTT_no_indel = new JTTGap(0, 0);

        // Stationary probabilities
        System.out.println(JTT.getProb('A'));
        System.out.println(JTT_no_indel.getProb('A'));
        System.out.println(JTT_Gap.getProb('A') + '\n');

        // Conditional probabilities
        System.out.println(JTT.getProb('A', 'Y', 0.5));
        System.out.println(JTT_no_indel.getProb('A', 'Y', 0.5));
        System.out.println(JTT_Gap.getProb('A', 'Y', 0.5) + '\n');

        // Gap probabilities
        System.out.println(JTT_Gap.getProb('-'));
        System.out.println(JTT_no_indel.getProb('-') + '\n');

        // Gap specific functions
        System.out.println(JTT_Gap.ksi_t(0.05));
        System.out.println(JTT_Gap.gamma_t(0.05));
        System.out.println(JTT_Gap.prob_gap_given_i_t(0.05));
        System.out.println(JTT_Gap.prob_j_given_gap_t(0.05, 'A'));
        System.out.println(JTT_Gap.prob_gap_given_i_t(0.05));
    }
}

