/*
 * Copyright (c) 2010-2021 Haifeng Li. All rights reserved.
 *
 * Smile is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Smile is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Smile.  If not, see <https://www.gnu.org/licenses/>.
 */

package smile.stat.distribution;

import smile.math.MathEx;
import java.util.Arrays;
import java.io.Serial;

/**
 * The finite mixture of distributions from exponential family.
 * The EM algorithm can be used to learn the mixture model from data.
 * EM is particularly useful when the likelihood is an exponential family: the
 * E-step becomes the sum of expectations of sufficient statistics, and the
 * M-step involves maximizing a linear function. In such a case, it is usually
 * possible to derive closed form updates for each step.
 *
 * @author Haifeng Li
 */
public class ExponentialFamilyMixture extends Mixture {
    @Serial
    private static final long serialVersionUID = 2L;

    /** The log-likelihood when the distribution is fit on a sample data. */
    public final double L;
    /** The BIC score when the distribution is fit on a sample data. */
    public final double bic;

    /**
     * Constructor.
     * @param components a list of exponential family distributions.
     */
    public ExponentialFamilyMixture(Component... components) {
        this(0.0, 1, components);
    }

    /**
     * Constructor.
     * @param components a list of discrete exponential family distributions.
     * @param L the log-likelihood.
     * @param n the number of samples to fit the distribution.
     */
    ExponentialFamilyMixture(double L, int n, Component... components) {
        super(components);

        for (Component component : components) {
            if (!(component.distribution instanceof ExponentialFamily)) {
                throw new IllegalArgumentException("Component " + component + " is not of exponential family.");
            }
        }

        this.L = L;
        this.bic = L - 0.5 * length() * Math.log(n);
    }

    /**
     * Fits the mixture model with the EM algorithm.
     * @param components the initial configuration of mixture. Components may have
     *                   different distribution form.
     * @param x the training data.
     * @return the distribution.
     */
    public static ExponentialFamilyMixture fit(double[] x, Component... components) {
        return fit(x, components, 0.0, 500, 1E-4);
    }

    /**
     * Fits the mixture model with the EM algorithm.
     *
     * @param components the initial configuration.
     * @param x the training data.
     * @param gamma the regularization parameter. Although regularization works
     *              well for high dimensional data, it often reduces the model
     *              to too few components. For one-dimensional data, gamma should
     *              be 0 in general.
     * @param maxIter the maximum number of iterations.
     * @param tol the tolerance of convergence test.
     * @return the distribution.
     */
    public static ExponentialFamilyMixture fit(double[] x, Component[] components, double gamma, int maxIter, double tol) {
        if (x.length < components.length / 2) {
            throw new IllegalArgumentException("Too many components compared to data size.");
        }

        if (gamma < 0.0 || gamma > 0.2) {
            throw new IllegalArgumentException("Invalid regularization factor gamma: " + gamma);
        }

        int n = x.length;
        int k = components.length;

        double[][] posteriori = new double[k][n];
        double L = Double.NEGATIVE_INFINITY;
        double diff = Double.MAX_VALUE;

        for (int iter = 1; iter <= maxIter && diff > tol; iter++) {
            // ----------- E-step -----------
            for (int i = 0; i < k; i++) {
                Component c = components[i];
                for (int j = 0; j < n; j++) {
                    posteriori[i][j] = c.priori * c.distribution.p(x[j]);
                }
            }

            // Normalize posterior probabilities
            for (int j = 0; j < n; j++) {
                double total = 0.0;
                for (int i = 0; i < k; i++) {
                    total += posteriori[i][j];
                }

                if (total <= 1e-300) {
                    for (int i = 0; i < k; i++) {
                        posteriori[i][j] = 1.0 / k;
                    }
                } else {
                    for (int i = 0; i < k; i++) {
                        posteriori[i][j] /= total;
                    }
                }

                if (gamma > 0) {
                    for (int i = 0; i < k; i++) {
                        posteriori[i][j] *= (1 + gamma * MathEx.log2(posteriori[i][j]));
                        if (Double.isNaN(posteriori[i][j]) || posteriori[i][j] < 0.0) {
                            posteriori[i][j] = 0.0;
                        }
                    }
                }
            }

            // ----------- M-step -----------
            double[] newPrioris = new double[k];
            for (int i = 0; i < k; i++) {
                Component updated = ((ExponentialFamily) components[i].distribution).M(x, posteriori[i]);

                // If it's a GammaDistribution, extract shape and scale from toString
                if (updated.distribution instanceof GammaDistribution) {
                    GammaDistribution g = (GammaDistribution) updated.distribution;

                    // Parse shape and scale from toString()
                    String[] parts = g.toString().replace("Gamma Distribution(", "").replace(")", "").split(",");
                    double parsedScale = Double.parseDouble(parts[0].trim());
                    double parsedShape = Double.parseDouble(parts[1].trim());

                    // Apply lower bound to scale
                    double minScale = 1e-6;
                    if (parsedScale < minScale) {
                        parsedScale = minScale;
                    }

                    // Recreate GammaDistribution with safe parameters
                    g = new GammaDistribution(parsedShape, parsedScale);
                    components[i] = new Component(updated.priori, g);
                } else {
                    components[i] = updated;
                }

                newPrioris[i] = components[i].priori;
            }

            // Normalize priors
            double Z = Arrays.stream(newPrioris).sum();
            if (Z <= 1e-300) {
                throw new IllegalStateException("Total weight Z is zero or very small. EM collapsed.");
            }

            for (int i = 0; i < k; i++) {
                components[i] = new Component(newPrioris[i] / Z, components[i].distribution);
            }

            // ----------- Log-likelihood -----------
            double loglikelihood = 0.0;
            for (double xi : x) {
                double px = 0.0;
                for (Component c : components) {
                    px += c.priori * c.distribution.p(xi);
                }
                if (px > 1e-300) {
                    loglikelihood += Math.log(px);
                } else {
                    loglikelihood += -1e6;
                }
            }

            diff = loglikelihood - L;
            L = loglikelihood;
        }

        return new ExponentialFamilyMixture(L, x.length, components);
    }

}
