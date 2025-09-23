package stats;

import bn.prob.GammaDistrib;

import java.util.Arrays;
import java.util.Random;

public class Zipf implements IndelModel {

    private double s; // Exponent parameter for the Zipf distribution
    private Random rand;
    private int maxK; // Maximum value for observations in the sampling range

    /**
     * Defines a Zipf distribution
     *
     * @param s The exponent parameter of the distribution
     * @param seed The random seed for sampling
     * @param maxK The maximum value for the sampling range
     */
    public Zipf(double s, long seed, int maxK) {
        this.s = s;
        this.maxK = maxK;
        this.rand = new Random(seed);
    }

    /**
     * Defines a Zipf distribution
     *
     * @param s The exponent parameter of the distribution
     * @param maxK The maximum value for the sampling range
     */
    public Zipf(double s, int maxK) {
        this(s, System.currentTimeMillis(), maxK);
    }

    /** Default max value to be observed */
    public static int DEFAULT_MAXK = 1000; //
    /**
     * Defines a Zipf distribution
     *
     * @param s The exponent parameter of the distribution
     */
    public Zipf(double s) {
        this(s, System.currentTimeMillis(), DEFAULT_MAXK);
    }

    public String toString() {
        return "Zipf(lambda=" + s + ", maxK=" + maxK + ")";
    }

    @Override

    public String getTrAVIS() {
        if (this.maxK == DEFAULT_MAXK)
            return String.format("Zipf:%.3f", s);
        else
            return String.format("Zipf:%.3f,%d", s, maxK);
    }

    @Override
    public void setSeed(long seed) {
        this.rand = new Random(seed);
    }

    /**
     * Probability mass function (PMF)
     *
     * @param k The rank
     * @return The probability of rank k
     */
    public double p(int k) {
        return 1.0 / Math.pow(k, s) / zeta(s);
    }

    /**
     * Cumulative distribution function (CDF)
     *
     * @param k The rank up to which the cumulative probability is calculated
     * @return The cumulative probability of ranks up to and including k
     */
    public double cdf(int k) {
        double sum = 0;
        for (int i = 1; i <= k; i++) {
            sum += p(i);
        }
        return sum;
    }

    /**
     * Samples from the Zipf distribution
     *
     * @return A rank sampled from the distribution
     */
    public int sample() {
        double toss = rand.nextDouble();
        double sum = 0;
        for (int i = 1; i <= maxK; i++) {
            sum += p(i);
            if (sum > toss) {
                return i;
            }
        }
        return maxK;
    }

    /**
     * Computes an approximation of the Riemann zeta function value
     *
     * @param s The exponent parameter
     * @return The approximate value of the zeta function
     */
    private double zeta(double s) {
        double sum = 0;
        for (int i = 1; i <= maxK; i++) {
            sum += 1.0 / Math.pow(i, s);
        }
        return sum;
    }


    // Generalized harmonic number H_{N,s}
    public static double harmonic(int N, double s) {
        double sum = 0.0;
        for (int k = 1; k <= N; k++) {
            sum += Math.pow(k, -s);
        }
        return sum;
    }

    // Compute log-likelihood of Zipf(s, N)
    public static double logLikelihood(int[] data, double s, int N) {
        double sumLogs = 0.0;
        for (int x : data) {
            sumLogs += Math.log(x);
        }

        double harmonic = harmonic(N, s);

        return -s * sumLogs - data.length * Math.log(harmonic);
    }

    public double getLogLikelihood(int[] data) {
        return logLikelihood(data, this.s, this.maxK);
    }


    // Derivative of H_{N,s} with respect to s
    public static double harmonicDerivative(int N, double s) {
        double sum = 0.0;
        for (int k = 1; k <= N; k++) {
            sum -= Math.log(k) * Math.pow(k, -s);
        }
        return sum;
    }

    // Second derivative of H_{N,s} with respect to s
    public static double harmonicSecondDerivative(int N, double s) {
        double sum = 0.0;
        for (int k = 1; k <= N; k++) {
            sum += Math.pow(Math.log(k), 2) * Math.pow(k, -s);
        }
        return sum;
    }

    // Log-likelihood derivative wrt s
    public static double scoreFunction(int[] data, int N, double s) {
        double avgLogX = Arrays.stream(data).mapToDouble(Math::log).average().orElse(0.0);
        double H = harmonic(N, s);
        double Hprime = harmonicDerivative(N, s);
        return - (Hprime / H) - avgLogX;
    }

    /**
     * Estimate the parameter value of Zipf by MLE and return Zipf instance.
     * Does this via Newton-Raphson.
     * @param data dataset
     * @return an instance of Zipf with parameter value
     */
    public static Zipf fitMLE(int[] data) {
        return fitMLE(data, Arrays.stream(data).max().getAsInt());
    }

    /**
     * Estimate the parameter value s of the Zipf distribution by MLE
     * The implementation uses Newton–Raphson
     * @param data data set
     * @param N max value of observations
      */
    public static Zipf fitMLE(int[] data, int N) {
        double s = 2.0; // initial guess
        for (int iter = 0; iter < 100; iter++) {
            double H = harmonic(N, s);
            double Hprime = harmonicDerivative(N, s);
            double Hsecond = harmonicSecondDerivative(N, s);

            double avgLogX = Arrays.stream(data).mapToDouble(Math::log).average().orElse(0.0);

            double score = - (Hprime / H) - avgLogX;
            double fisherInfo = - ((Hsecond * H - Hprime * Hprime) / (H * H)); // d/ds of (H'/H)

            double step = score / fisherInfo;
            s -= step;

            if (Math.abs(step) < 1e-8) break; // convergence
        }
        return new Zipf(s, N);
    }


    // compute log posterior up to additive constant
    public static double logPosterior(double s, int[] data, int N, double alpha, double beta) {
        if (s <= 0) return Double.NEGATIVE_INFINITY;

        // compute sufficient statistics
        double sumLogX = 0.0;
        for (int x : data) {
            sumLogX += Math.log(x);
        }
        int n = data.length;

        // likelihood part
        double H = harmonic(N, s);
        double logLikelihood = -s * sumLogX - n * Math.log(H);

        // prior part
        double logPrior = GammaDistrib.logPrior(s, alpha, beta);

        return logLikelihood + logPrior;
    }

    /**
     * Estimate the parameter value of Zipf by MAP and return Zipf instance.
     * Does this via Golden section search.
     * @param data dataset
     * @return an instance of Zipf with parameter value
     */
    public static Zipf fitMAP(int[] data, double alpha, double beta, double lower, double upper) {
        int N = Arrays.stream(data).max().getAsInt();
        double s = goldenSectionSearch(data, N, alpha, beta, lower, upper,1e-6, 1000);
        return new Zipf(s, N);
    }


    // Golden-section search for maximum
    public static double goldenSectionSearch(int[] data, int N, double alpha, double beta,
                                             double lower, double upper, double tol, int maxIter) {
        double gr = (Math.sqrt(5) + 1) / 2; // golden ratio

        double c = upper - (upper - lower) / gr;
        double d = lower + (upper - lower) / gr;

        int iter = 0;
        while (Math.abs(c - d) > tol && iter < maxIter) {
            double fc = logPosterior(c, data, N, alpha, beta);
            double fd = logPosterior(d, data, N, alpha, beta);

            if (fc > fd) {
                upper = d;
                d = c;
                c = upper - (upper - lower) / gr;
            } else {
                lower = c;
                c = d;
                d = lower + (upper - lower) / gr;
            }
            iter++;
        }
        return (lower + upper) / 2;
    }

    /**
     * Find the distribution with maximum data likelihood
     * @param indel_data
     * @return the best model
     */
    public static IndelModel bestfit(int[] indel_data) {
        return fitMLE(indel_data);
    }

    public static void main(String[] args) {
        // Example dataset
        int[] data = {1, 2, 2, 3, 5, 1, 10};
        int N = Arrays.stream(data).max().getAsInt();

        // Gamma prior hyperparameters
        double alpha = 2.0;
        double beta = 0.1;

        // Search interval for s
        double lower = 0.1;
        double upper = 5.0;

        // Run golden-section search
        System.out.printf("MAP gives " + fitMAP(data, alpha, beta, lower, upper));
    }

    public static void main1(String[] args) {
        // Example dataset (Zipf-like samples)
        int[] data = {1, 2, 2, 3, 5, 1, 10};

        // N = maximum observed value
        int N = Arrays.stream(data).max().getAsInt();

        Zipf zipf = fitMLE(data, N);
        System.out.println("MLE gives " + zipf);
    }

    public static void main0(String[] args) {
        // Example: Using different s and maxK values
        for (double s = 1.0; s < 5.5; s += 0.5) {
            Zipf zipf = new Zipf(s, 1, 1000); // Set maxK to 5000
            System.out.println("s = " + s);
            int[] cnt = new int[zipf.maxK + 1];
            for (int i = 0; i < 1000; i++) {
                cnt[zipf.sample()] += 1;
            }
            for (int j = 1; j <= 20; j++) {
                System.out.print(cnt[j] + "\t");
            }
            System.out.println();
        }
    }
}
