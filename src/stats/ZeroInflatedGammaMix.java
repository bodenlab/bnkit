package stats;
import smile.stat.distribution.ExponentialFamilyMixture;
import smile.stat.distribution.GammaDistribution;
import smile.stat.distribution.Mixture;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class ZeroInflatedGammaMix {    private double p;  // Zero-inflation probability
    private ExponentialFamilyMixture gammaMixture;
    private Random rand;
    private long seed;

    /** Constructor */
    public ZeroInflatedGammaMix(double p, ExponentialFamilyMixture gammaMixture, long seed) {
        if (p < 0 || p > 1) {
            throw new IllegalArgumentException("Invalid zero-inflation probability p: " + p);
        }
        this.p = p;
        this.gammaMixture = gammaMixture;
        this.seed = seed;
        this.rand = new Random(seed);
    }

    public double getP() {
        return p;
    }

    public ExponentialFamilyMixture getGammaMixture() {
        return gammaMixture;
    }

    /** Probability mass function (PMF) */
    public double p(double x) {
        if (x < 0) {
            return 0.0;
        } else if (x == 0) {
            return p;
        } else {
            return (1 - p) * gammaMixture.p(x);
        }
    }

    /** Sample a value */
    public double sample() {
        if (rand.nextDouble() < p) {
            return 0.0;
        }

        double r = rand.nextDouble();
        double cumulative = 0.0;
        for (Mixture.Component c : gammaMixture.components) {
            cumulative += c.priori;
            if (r <= cumulative) {
                return sampleGamma((GammaDistribution) c.distribution);
            }
        }

        return sampleGamma((GammaDistribution) gammaMixture.components[gammaMixture.length() - 1].distribution);
    }
    private double getShape(GammaDistribution gamma) {
        return Double.parseDouble(gamma.toString().split(",")[1].replace(")", "").trim());
    }


    private double getScale(GammaDistribution gamma) {
        return Double.parseDouble(gamma.toString().split(",")[0]
                .replace("Gamma Distribution(", "")
                .trim());
    }

    private double sampleGamma(GammaDistribution gamma) {
        double shape = getShape(gamma);
        double scale = getScale(gamma);
        if (shape <= 0.0 || scale <= 0.0) {
            throw new IllegalArgumentException("Shape and scale must be positive.");
        }

        if (shape < 1.0) {
            // Boosting method: Gamma(alpha) = Gamma(alpha + 1) * U^(1/alpha)
            double u = rand.nextDouble();
            GammaDistribution boosted = new GammaDistribution(shape + 1.0, scale);
            return sampleGamma(boosted) * Math.pow(u, 1.0 / shape);
        }

        // Marsaglia and Tsang method for shape >= 1
        double d = shape - 1.0 / 3.0;
        double c = 1.0 / Math.sqrt(9.0 * d);

        while (true) {
            double x = rand.nextGaussian(); // Normal(0,1)
            double v = 1.0 + c * x;
            if (v <= 0.0) continue;
            v = v * v * v;
            double u = rand.nextDouble();

            if (u < 1.0 - 0.0331 * x * x * x * x) {
                return scale * d * v;
            }
            if (Math.log(u) < 0.5 * x * x + d * (1.0 - v + Math.log(v))) {
                return scale * d * v;
            }
        }
    }

    /** Fit Zero-Inflated Gamma Mixture model from data */
    public static ZeroInflatedGammaMix fit(double[] data) {
        if (data.length == 0) {
            throw new IllegalArgumentException("Data cannot be empty.");
        }

        int zeroCount = (int) Arrays.stream(data).filter(x -> x == 0.0).count();
        double p = (double) zeroCount / data.length;

        double[] nonZeroData = Arrays.stream(data).filter(x -> x > 0).toArray();
        if (nonZeroData.length == 0) {
            return new ZeroInflatedGammaMix(1.0, null, 42);
        }

        // Initialize 2 gamma components
        Mixture.Component[] components = new Mixture.Component[2];
        components[0] = new Mixture.Component(0.5, new GammaDistribution(2.0, 2.0));
        components[1] = new Mixture.Component(0.5, new GammaDistribution(5.0, 1.0));

        ExponentialFamilyMixture gammaMixture = ExponentialFamilyMixture.fit(nonZeroData, components);

        return new ZeroInflatedGammaMix(p, gammaMixture, 42);
    }
    /**
     * Fit a Zero-Inflated Gamma Mixture model with specified number of components.
     *
     * @param data training data
     * @param numComponents number of Gamma components in the mixture
     * @return fitted ZeroInflatedGammaMix
     */
    public static ZeroInflatedGammaMix fit(double[] data, int numComponents) {
        if (data.length == 0) {
            throw new IllegalArgumentException("Data cannot be empty.");
        }
        if (numComponents < 1) {
            throw new IllegalArgumentException("Number of components must be >= 1.");
        }


        int zeroCount = (int) Arrays.stream(data).filter(x -> x == 0.0).count();
        double p = (double) zeroCount / data.length;


        double[] nonZeroData = Arrays.stream(data).filter(x -> x > 0).toArray();
        if (nonZeroData.length == 0) {
            return new ZeroInflatedGammaMix(1.0, null, 42);
        }


        Arrays.sort(nonZeroData);
        Mixture.Component[] components = new Mixture.Component[numComponents];
        for (int i = 0; i < numComponents; i++) {
            int qIndex = (int) ((i + 1.0) / (numComponents + 1) * nonZeroData.length);
            double mean = Math.max(nonZeroData[qIndex], 1e-3);  // 避免为0
            double shape = 2.0;
            double scale = mean / shape;
            shape = Math.max(shape, 0.01);
            scale = Math.max(scale, 1e-4);
            components[i] = new Mixture.Component(1.0 / numComponents, new GammaDistribution(shape, scale));
        }


        ExponentialFamilyMixture gammaMixture;
        try {
            gammaMixture = ExponentialFamilyMixture.fit(nonZeroData, components);
        } catch (Exception e) {
            System.err.println("Gamma mixture fit failed: " + e.getMessage());
            return new ZeroInflatedGammaMix(p, null, 42);
        }


        List<Mixture.Component> valid = new ArrayList<>();
        for (Mixture.Component c : gammaMixture.components) {
            GammaDistribution g = (GammaDistribution) c.distribution;
            double shape = g.k;
            double scale = g.theta;
            double weight = c.priori;

            if (Double.isNaN(shape) || Double.isNaN(scale) || shape > 100 || scale < 1e-5 || weight < 1e-3) {
                continue;
            }


            shape = Math.min(Math.max(shape, 0.01), 100);
            scale = Math.min(Math.max(scale, 1e-4), 10);
            valid.add(new Mixture.Component(weight, new GammaDistribution(shape, scale)));
        }


        if (valid.isEmpty()) {
            return new ZeroInflatedGammaMix(p, null, 42);
        } else {

            double totalWeight = valid.stream().mapToDouble(c -> c.priori).sum();
            Mixture.Component[] normalized = new Mixture.Component[valid.size()];
            for (int i = 0; i < valid.size(); i++) {
                Mixture.Component c = valid.get(i);
                double newWeight = c.priori / totalWeight;
                GammaDistribution g = (GammaDistribution) c.distribution;
                normalized[i] = new Mixture.Component(newWeight,
                        new GammaDistribution(g.k, g.theta));
            }

            return new ZeroInflatedGammaMix(p, new ExponentialFamilyMixture(normalized), 42);
        }
    }

    /** Summary output */
    public String getTrAVIS() {
        return String.format("ZeroInflatedGammaMix p=%.4f, mixture components=%d", p, gammaMixture.length());
    }

    /** 测试用主程序 */
    public static void main(String[] args) {
        double[] data = {0.0, 1.2, 3.4, 0.0, 2.1, 0.5, 0.0, 4.2, 3.3, 0.0, 1.7};

        ZeroInflatedGammaMix zig = ZeroInflatedGammaMix.fit(data);
        System.out.println("Fitted Model p: " + zig.getP());
        System.out.println("Fitted Model #Components: " + zig.getGammaMixture().length());

        System.out.println("Components Details:");
        int idx = 1;
        for (Mixture.Component c : zig.getGammaMixture().components) {
            GammaDistribution gamma = (GammaDistribution) c.distribution;
            String[] parts = gamma.toString().split(",");
            double scale = Double.parseDouble(parts[0].replace("Gamma Distribution(", "").trim());
            double shape = Double.parseDouble(parts[1].replace(")", "").trim());
            double weight = c.priori;

            System.out.printf("  Component %d: weight = %.4f, shape = %.4f, scale = %.4f%n", idx, weight, shape, scale);
            idx++;
        }

        System.out.println("Sampled values:");
        for (int i = 0; i < 10; i++) {
            System.out.printf("%.3f ", zig.sample());
        }
    }
}
