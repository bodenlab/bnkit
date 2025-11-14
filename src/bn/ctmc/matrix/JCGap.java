package bn.ctmc.matrix;

import bn.ctmc.GapSubstModel;
import dat.Enumerable;

public class JCGap extends GapSubstModel {

    public static double[] F = {0.25, 0.25, 0.25, 0.25};

    public static double[][] Q = {
            //A     T      G     C
            { 0,    0.25,  0.25, 0.25},
            { 0.25, 0,     0.25, 0.25},
            { 0.25, 0.25,  0,    0.25},
            { 0.25, 0.25,  0.25, 0   }
    };

    public static Character[] S = {'A','C','G','T'};

    public JCGap(double mu, double lambda) {
        super(F, Q, new Enumerable(S), mu, lambda, false,
                false, false);
    }


    public JCGap(double[][] IRM, double mu, double lambda, boolean symmetric, boolean normalise) {
        super(F, IRM, new Enumerable(S), mu, lambda, symmetric, normalise, false);
    }


    @Override
    public String getName() {
        return "JCGap";
    }
}
