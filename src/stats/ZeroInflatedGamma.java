package stats;

import java.util.Arrays;
import java.util.Random;
import bn.prob.GammaDistrib;

/**
 * Zero-Inflated Gamma (ZIG) Distribution
 * <p>
 * This distribution combines a probability mass at zero with a Gamma-distributed component.
 *
 * @author Chongting Zhao
 */
public class ZeroInflatedGamma {
    private double p;      // Zero-inflation probability
    private double shape;  // Gamma shape parameter
    private double scale;  // Gamma scale parameter
    private Random rand;
    private long seed;

    /**
     * Constructor for Zero-Inflated Gamma distribution.
     *
     * @param p     Zero-inflation probability
     * @param shape Gamma shape parameter
     * @param scale Gamma scale parameter
     * @param seed  Random seed for sampling
     */
    public ZeroInflatedGamma(double p, double shape, double scale, long seed) {
        if (p < 0 || p > 1) {
            throw new IllegalArgumentException("Invalid zero-inflation probability p: " + p);
        }
        if (shape <= 0 || scale <= 0) {
            throw new IllegalArgumentException("Shape and scale must be positive.");
        }

        this.p = p;
        this.shape = shape;
        this.scale = 1/scale;
        this.seed = seed;
        this.rand = new Random(seed);
    }
    /**
     * Get the zero-inflation probability p.
     * @return p (probability of zero inflation)
     */
    public double getP() {
        return p;
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
            return p;
        } else {
            GammaDistrib gamma = new GammaDistrib(shape, scale);
            return (1 - p) * gamma.get(x);
        }
    }


    /**
     * Estimates the parameters (p, k, θ) from data using MLE / MOM.
     *
     * @param data sample data
     * @return a fitted ZeroInflatedGamma distribution
     */
    public static ZeroInflatedGamma fit(double[] data) {
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
     * Sample from the Zero-Inflated Gamma distribution.
     *
     * @return the sampled value
     */
    public double sample() {
        if (rand.nextDouble() < p) {
            return 0.0;
        }
        GammaDistrib gamma = new GammaDistrib(shape, scale);
        return gamma.sample();
    }

    /**
     * Return a string representation.
     *
     * @return a formatted string
     */
    public String getTrAVIS() {
        return String.format("ZeroInflatedGamma p=%.4f, shape=%.4f, scale=%.4f", p, shape, scale);
    }

    public static void main(String[] args) {
        double[] data = {0.0, 1.2, 3.4, 0.0, 2.1, 0.5, 0.0, 4.2, 3.3, 0.0, 1.7};

        // 训练 ZIG 模型
        ZeroInflatedGamma zig = ZeroInflatedGamma.fit(data);
        System.out.println("Fitted Model: " + zig.p);
        System.out.println("Fitted Model: " + zig.shape);
        System.out.println("Fitted Model: " + zig.scale);
        

        // 采样测试
        System.out.println("Sampled values:");
        for (int i = 0; i < 10; i++) {
            System.out.printf("%.3f ", zig.sample());
        }
    }
}

