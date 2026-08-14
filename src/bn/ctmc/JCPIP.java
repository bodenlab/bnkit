package bn.ctmc;

import bn.ctmc.matrix.JC;
import dat.Enumerable;

public class JCPIP extends PIPSubstModel {

    public static double[] F = {0.25, 0.25, 0.25, 0.25};

    public static double[][] Q = {
            //A     T      G     C
            { 0,    0.25,  0.25, 0.25},
            { 0.25, 0,     0.25, 0.25},
            { 0.25, 0.25,  0,    0.25},
            { 0.25, 0.25,  0.25, 0   }
    };

    public static Character[] S = {'A','C','G','T'};

    /**
     * Create time reversible evolutionary model.
     *
     * @param F        stationary base frequencies
     * @param S        Symmetric, un-scaled version of Q matrix Q_ij = s_ij*pi_j as defined by PAML
     * @param alphabet the values that substitutable variables can take, listed strictly in the order of the array and matrix
     * @param mu
     * @param lambda
     */
    public JCPIP(double[] F, double[][] S, Enumerable alphabet, double mu, double lambda) {
        super(F, S, alphabet, mu, lambda);
    }

    public JCPIP(double alpha, Enumerable domain, double mu, double lambda) {
        super(JC.F(domain.size()), JC.Q(alpha, domain.size()), domain, mu, lambda,false, false);
    }

    public JCPIP(double mu, double lambda) {
        super(F, Q, new Enumerable(S), mu, lambda, false, false);
    }
}
