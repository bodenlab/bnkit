package stats;

import java.util.Arrays;
import java.util.Random;

import bn.Distrib;
import bn.prob.GammaDistrib;

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
        if (data.length == 0) {
            throw new IllegalArgumentException("Data cannot be empty.");
        }

        int zeroCount = (int) Arrays.stream(data).filter(x -> x == 0.0).count();
        double p = (double) zeroCount / data.length;

        double[] nonZeroData = Arrays.stream(data).filter(x -> x > 0).toArray();
        if (nonZeroData.length == 0) {
            return new ZeroInflatedGamma(1.0, 1.0, 1.0, 42);
        }

        double mean = Arrays.stream(nonZeroData).average().orElse(0.0);
        double variance = Arrays.stream(nonZeroData)
                .map(x -> Math.pow(x - mean, 2))
                .sum() / nonZeroData.length;

        double shape = mean * mean / variance;
        double scale = variance / mean;

        return new ZeroInflatedGamma(p, shape, scale, 42);
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

    public static void main(String[] args) {
        double[] data = {0.0, 1.2, 3.4, 0.0, 2.1, 0.5, 0.0, 4.2, 3.3, 0.0, 1.7};

        // 训练 ZIG 模型
        ZeroInflatedGamma zig = ZeroInflatedGamma.fitMLE(data);
        System.out.println("Fitted Model: " + zig.pi);
        System.out.println("Fitted Model: " + zig.shape);
        System.out.println("Fitted Model: " + zig.scale);
        

        // 采样测试
        System.out.println("Sampled values:");
        for (int i = 0; i < 10; i++) {
            System.out.printf("%.3f ", zig.sample());
        }
    }
}

