package bn.prob;

import java.util.*;

import smile.stat.distribution.*;

public class GammaMixtureEM {

    public static void main(String[] args) {
        int k = 2; // number of components
        double[] data = {0.13, 0.21, 0.09, 0.3, 0.14, 0.45, 0.03, 0.06, 0.01, 0.05, 0.04, 0.07, 0.35, 0.08, 0.18, 0.08, 0.38, 0.25, 0.33};
        Arrays.sort(data);
        double[] subset = new double[data.length/k];
        GammaDistribution[] gamma = new GammaDistribution[k];
        for (int i = 0; i < k; i++) {
            System.arraycopy(data, i*subset.length, subset, 0, subset.length);
            gamma[i] = GammaDistribution.fit(subset);
            for (int j = 0; j < subset.length; j++) {
                System.out.print(subset[j] + " ");
            }
            System.out.println(" rendered: " + gamma[i]);
            GammaDistrib gd = GammaDistrib.fitMLE(subset);
            System.out.println("\tMLE: " + gd);
        }
        Mixture mixture = ExponentialFamilyMixture.fit(data,
            new Mixture.Component(0.5, gamma[0]),
            new Mixture.Component(0.5, gamma[1])
        );
        for (int i = 0; i < k; i++) {
            ;
            double a = ((GammaDistribution)(mixture.components[i].distribution)).k;
            double b = ((GammaDistribution)(mixture.components[i].distribution)).theta;
            System.out.println("Gamma " + i +"\tShape (k, a in Matlab)= " + a + "\tScale (theta, b in Matlab)= " + b);
            System.out.println("y" + (i+1) + " = gampdf(x, " + a + ", " + b + ") * " + mixture.components[i].priori + ";");

        }
        System.out.println(mixture);
    }
}