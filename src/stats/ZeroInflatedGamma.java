package stats;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import bn.Distrib;
import bn.prob.GammaDistrib;
import smile.stat.distribution.GammaDistribution;

/**
 * A zero-inflated gamma distribution is a mixture of:
 * a point mass at zero (for excess zeros), and
 * a gamma distribution (for the positive continuous part).
 *
 * @author Chongting Zhao, Mikael Boden
 */
public class ZeroInflatedGamma implements RateModel, Distrib {
    private GammaDistrib gamma;
    private double pi;     // Zero-inflation probability
    private double shape;  // Gamma shape parameter (aka alpha, or k)
    private double scale;  // Gamma scale parameter; note that rate (aka beta or lambda) is the inverse of scale (aka theta)
    private Random rand;

    /**
     * Constructor for Zero-Inflated Gamma distribution.
     *
     * @param pi    Zero-inflation probability
     * @param shape Gamma shape parameter (aka alpha, or k)
     * @param scale Gamma scale parameter; note that rate (aka beta or lambda) is the inverse of scale (aka theta)
     * @param seed  Random seed for sampling
     */
    public ZeroInflatedGamma(double pi, double shape, double scale, long seed) {
        if (pi < 0 || pi > 1) {
            throw new IllegalArgumentException("Invalid zero-inflation probability p: " + pi);
        }
        if (shape <= 0 || scale <= 0) {
            throw new IllegalArgumentException("Shape and scale must be positive.");
        }
        this.gamma = new GammaDistrib(shape, 1.0/scale, seed);
        this.pi = pi;
        this.shape = shape;
        this.scale = scale; // here...
        this.rand = new Random(seed);
    }

    /**
     * Constructor for Zero-Inflated Gamma distribution.
     *
     * @param pi    Zero-inflation probability
     * @param shape Gamma shape parameter
     * @param scale Gamma scale parameter; note that rate is the inverse of scale
     */
    public ZeroInflatedGamma(double pi, double shape, double scale) {
        this(pi, shape, scale, System.currentTimeMillis());
    }

    public void setSeed(long seed) {
        rand = new Random(seed);
    }

    /**
     * Get the zero-inflation probability p.
     * @return p (probability of zero inflation)
     */
    public double getZIP() {
        return pi;
    }

    /**
     * Get the shape parameter (k).
     * @return shape parameter
     */
    public double getShape() {
        return shape;
    }

    /**
     * Get the scale parameter (θ).
     * @return scale parameter
     */
    public double getScale() {
        return scale;
    }

    /**
     * Probability mass function (PMF)
     *
     * @param x the value
     * @return the probability of x
     */
    public double p(double x) {
        if (x < 0) {
            return 0.0;
        } else if (x == 0) {
            return pi;
        } else {
            return (1 - pi) * gamma.get(x);
        }
    }

    /**
     * Estimates the parameters (p, k, θ) from data using MLE / MOM.
     *
     * @param data sample data
     * @return a fitted ZeroInflatedGamma distribution
     */
    public static ZeroInflatedGamma fitMLE(double[] data) {
        return fitMLE(data, System.currentTimeMillis());
    }

    /**
     * Estimates the parameters (p, k, θ) from data using MLE / MOM.
     *
     * @param data sample data
     * @param seed random seed
     * @return a fitted ZeroInflatedGamma distribution
     */
    public static ZeroInflatedGamma fitMLE(double[] data, long seed) {
        if (data.length == 0) {
            throw new IllegalArgumentException("Data cannot be empty.");
        }

        int zeroCount = (int) Arrays.stream(data).filter(x -> x == 0.0).count();
        double p = (double) zeroCount / data.length;

        double[] nonZeroData = Arrays.stream(data).filter(x -> x > 0).toArray();
        if (nonZeroData.length == 0) {
            return new ZeroInflatedGamma(1.0, 1.0, 1.0, seed);
        }

        double mean = Arrays.stream(nonZeroData).average().orElse(0.0);
        double variance = Arrays.stream(nonZeroData)
                .map(x -> Math.pow(x - mean, 2))
                .sum() / nonZeroData.length;

        double shape = mean * mean / variance;
        double scale = variance / mean;

        return new ZeroInflatedGamma(p, shape, scale, seed);
    }

    /**
     * Retrieve the probability of the value
     *
     * @param value
     * @return the probability
     */
    @Override
    public double get(Object value) {
        double x = (Double) value;
        return p(x);
    }

    /**
     * Sample from the Zero-Inflated Gamma distribution.
     *
     * @return the sampled value
     */
    public Double sample() {
        if (rand.nextDouble() < pi) {
            return 0.0;
        }
        return gamma.sample();
    }

    @Override
    public double cdf(double rate) {
        throw new UnsupportedOperationException("CDF for ZeroInflatedGamma is not supported yet.");
        /* // CDF for regular Gamma:
        if (rate < 0) {
            return 0.0;
        } else {
            return Gamma.regularizedIncompleteGamma(getShape(), rate / getScale());
        }
         */
    }

    /**
     * Return a string representation.
     *
     * @return a formatted string
     */
    public String getTrAVIS() {
        return String.format("ZeroInflatedGamma:%.3f,%.3f,%.3f", pi, shape, scale);
    }

    /**
     * Calculate the log likelihood of the data given the model/distribution
     *
     * @param data dataset
     * @return the log-likelihood of the data given the model
     */
    @Override
    public double getLogLikelihood(double[] data) {
        double logL = 0.0;
        for (double x : data) {
            if (x == 0.0) {
                logL += Math.log(pi);
            } else if (x > 0.0) {
                logL += Math.log(1.0 - pi) + Math.log(gamma.get(x));
            } // ignore negative values (not supported by the model)
        }
        return logL;
    }

    /**
     * Compound mixture of zero point mass and mixture of Gamma distributions.
     */
    public static class Mixture implements RateModel, Distrib {
        public double zeromass;
        public GammaDistrib[] distribs;
        public double[] priors;
        private final Random rand;

        /**
         * Constructs a mixture of Zero-Inflated Gamma distributions.
         * @param zeromass the probability of zero
         * @param gammadistribs array of Gamma components
         * @param priors mixture weights (must sum to 1 - zeromass)
         * @param seed random seed for sampling
         */
        public Mixture(double zeromass, GammaDistrib[] gammadistribs, double[] priors, long seed) {
            this.zeromass = zeromass;
            if (gammadistribs.length != priors.length) {
                throw new IllegalArgumentException("Gamma components and priors must have same length.");
            }
            this.distribs = gammadistribs;
            this.priors = priors;
            this.rand = new Random(seed);
        }

        @Override
        public void setSeed(long seed) {
            rand.setSeed(seed);
            for (GammaDistrib comp : distribs) {
                comp.setSeed(seed);
            }
        }

        /**
         * Probability mass function for the mixture.
         *
         * @param x value
         * @return probability
         */
        @Override
        public double p(double x) {
            double prob = 0.0;
            if (x == 0.0)
                return zeromass;
            for (int i = 0; i < distribs.length; i++) {
                prob += priors[i] * distribs[i].p(x);
            }
            return prob;
        }

        @Override
        public double cdf(double x) {
            double result = x >= 0 ? zeromass : 0.0;
            for (int i = 0; i < distribs.length; i++) {
                result += priors[i] * distribs[i].cdf(x);
            }
            return result;
        }

        @Override
        public double get(Object value) {
            double x = (Double) value;
            return p(x);
        }

        /**
         * Sample from the mixture.
         *
         * @return sampled value
         */
        @Override
        public Double sample() {
            double r = rand.nextDouble();
            double cum = zeromass;
            if (r < cum)
                return 0.0;
            for (int i = 0; i < priors.length; i++) {
                cum += priors[i];
                if (r < cum) {
                    return distribs[i].sample();
                }
            }
            return distribs[distribs.length - 1].sample();
        }

        /**
         * Log-likelihood of data under the mixture.
         *
         * @param data dataset
         * @return log-likelihood
         */
        @Override
        public double getLogLikelihood(double[] data) {
            double logL = 0.0;
            for (double x : data) {
                double px = p(x);
                if (px > 0) {
                    logL += Math.log(px);
                }
            }
            return logL;
        }

        /**
         * Get the number of mixture components.
         * @return number of components
         */
        public int getComponentCount() {
            return distribs.length;
        }

        /**
         * Get the i-th mixture component.
         * @param i index
         * @return ZeroInflatedGamma component
         */
        public GammaDistrib getComponent(int i) {
            return distribs[i];
        }

        /**
         * Get the mixture weights.
         * @return array of priors
         */
        public double[] getPriors() {
            return Arrays.copyOf(priors, priors.length);
        }

        /**
         * Fits a compound mixture of zero point mass and mixture of Gamma distributions to the specified data set
         * using Maximum Likelihood Estimation (MLE) via the EM algorithm.
         *
         * @param data the data array to fit
         * @param components the number of mixture Gamma components
         * @param seed the random seed for initialization
         * @return a fitted Mixture model
         */
        public static Mixture fitMLE(double[] data, int components, long seed) {
            if (data.length == 0) {
                throw new IllegalArgumentException("Data cannot be empty.");
            }
            GammaDistrib[] distribs = new GammaDistrib[components];
            double[] priors = new double[components];
            int zeroCount = (int) Arrays.stream(data).filter(x -> x == 0.0).count();
            double zeromass = (double) zeroCount / data.length;
            double[] nonZeroData = Arrays.stream(data).filter(x -> x > 0).toArray();
            if (nonZeroData.length == 0) {
                for (int i = 0; i < components; i++) {
                    distribs[i] = new GammaDistrib(1.0, 1.0, seed + i);
                    priors[i] = 0.0;
                }
                return new Mixture(1.0, distribs, priors, seed);
            }

            int n = nonZeroData.length;
            double[][] resp = new double[n][components];
            // Initialize priors and components
            Arrays.fill(priors, 1.0 / components);
            for (int k = 0; k < components; k++) {
                // Randomly assign data to clusters for initial fit
                int start = k * n / components;
                int end = (k + 1) * n / components;
                double[] subset = Arrays.copyOfRange(nonZeroData, start, end);
                distribs[k] = GammaDistrib.fitMLE(subset, seed + k);
            }

            // EM algorithm
            int maxIter = 100;
            for (int iter = 0; iter < maxIter; iter++) {
                // E-step: compute responsibilities
                for (int i = 0; i < n; i++) {
                    double sum = 0;
                    for (int k = 0; k < components; k++) {
                        resp[i][k] = priors[k] * distribs[k].p(nonZeroData[i]);
                        sum += resp[i][k];
                    }
                    for (int k = 0; k < components; k++) {
                        resp[i][k] /= sum;
                    }
                }

                // M-step: update priors and fit each component
                for (int k = 0; k < components; k++) {
                    // Weighted data for component k
                    List<Double> weightedData = new ArrayList<>();
                    for (int i = 0; i < n; i++) {
                        for (int r = 0; r < (int) (resp[i][k] * 100); r++) {
                            weightedData.add(nonZeroData[i]);
                        }
                    }
                    if (weightedData.size() > 0) {
                        distribs[k] = GammaDistrib.fitMLE(weightedData, seed + k);
                    }
                    // Update priors
                    double sumResp = 0;
                    for (int i = 0; i < n; i++) {
                        sumResp += resp[i][k];
                    }
                    priors[k] = sumResp / n;
                }
            }
            for (int k = 0; k < components; k++) {
                priors[k] -= zeromass / components;
            }
            return new Mixture(zeromass, distribs, priors, seed);
        }

        public String getTrAVIS() {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("MixtureZIG:%.3f,", zeromass));
            for (int i = 0; i < distribs.length; i++) {
                sb.append(String.format("%.3f,%.3f,%.3f", distribs[i].getShape(), distribs[i].getScale(), priors[i]));
                if (i < distribs.length - 1)
                    sb.append(",");
            }
            return sb.toString();
        }
    }

    public static void main(String[] args) {
        double[] X = {0.000001, 11.2, 8.3, 13.1, 15.9, 11.5, 11.4, 12.3, 11.9, 5.5, 0.001, 1.2, 2.3, 3.1, 0, 0, 5.9, 1.5, 1.4, 2.3, 1.9, 3.5, 2.3, 2.1, 2.9, 0.5, 0.4, 0.3, 0.9, 0.15, 0.01, 0.22, 0.23, 0.123, 0, 0};

        ZeroInflatedGamma zig = ZeroInflatedGamma.fitMLE(X, 1);
        System.out.println(zig.getTrAVIS());
        ZeroInflatedGamma.Mixture zigm2 = ZeroInflatedGamma.Mixture.fitMLE(X, 2, 1);
        System.out.println(zigm2.getTrAVIS());

        System.out.println("Sampled values:");
        for (int i = 0; i < 10; i++) {
            System.out.printf("%.3f\t%.3f\n", zig.sample(), zigm2.sample());
        }
    }


}

