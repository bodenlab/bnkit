package bn.ctmc.matrix;

import bn.ctmc.GapSubstModel;

public class YangGap extends GapSubstModel {

    public static Character[] S = {'A','C','G','T'};

    public static double[] F = {
            // A      C   	 G      T
            0.308, 0.185, 0.199, 0.308};
    public static double[][] Q = {
            //A      C      G      T
            {0.000, 0.139, 0.566, 0.115}, //
            {0.232, 0.000, 0.223, 0.882},
            {0.877, 0.207, 0.000, 0.215},
            {0.115, 0.530, 0.139, 0.000}};

    public YangGap(double mu, double lambda) {
        super(F, Q, new dat.Enumerable(S), mu, lambda);
    }

    public YangGap(double[] empiricalFreqs, double mu, double lambda) {
        super(empiricalFreqs, Q, new dat.Enumerable(S), mu, lambda);
    }

    @Override
    public String getName() {
        return "YangGap";
    }
}
