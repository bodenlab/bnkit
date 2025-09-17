package stats;

import bn.prob.GammaDistrib;

import java.util.Arrays;
import java.util.Random;
public class ZeroTruncatedPoisson implements IndelModel {

    private double lambda;
    private Random rand;


    /**
     * Define a zero-truncated Poisson distribution
     *
     * @param lambda the average number of events per interval
     * @param seed the random seed in case it is used for sampling
     */
    public ZeroTruncatedPoisson(double lambda, long seed) {
        this.lambda = lambda;
        this.rand = new Random(seed);
    }

    /**
     * Define a zero-truncated Poisson distribution
     *
     * @param lambda the average number of events per interval
     */
    public ZeroTruncatedPoisson(double lambda) {
        this(lambda, System.currentTimeMillis());
    }

    public String toString() {
        return "ZeroTruncatedPoisson(lambda=" + lambda + ")";
    }

    @Override
    public String getTrAVIS() {
        return String.format("ZeroTruncatedPoisson:%.3f", lambda);
    }

    @Override
    public void setSeed(long seed) {
        this.rand = new Random(seed);
    }

    /**
     * The probability mass function
     * @param k the number of events in the interval
     * @return the probability of k events
     */
    public double p(double k) {
        return Math.exp(k * Math.log(lambda) - lambda - lgamma(k + 1)) / (1 - Math.exp(-lambda));
    }

    /**
     * The probability mass function.
     * The implementation is numerically stable and uses the log gamma
     * (https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution)
     * @param k the number of events in the interval
     * @return the probability of k events
     */
    public double p(int k) {
        return Math.exp(k * Math.log(lambda) - lambda - lgamma(k + 1)) / (1 - Math.exp(-lambda));
    }

    /**
     * The cumulative probability function.
     * The implementation calls the PMF for all values i from 0 to floor(k)
     * @param k the number of events in the interval
     * @return the cumulative probability of k events
     * https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution
     */
    public double cdf(int k) {
        double sum = 0;
        for (int i = 0; i <= k; i ++) {
            sum += p(i);
            if (sum >= 1.0) // happens only due to poor numerical precision
                return 1.0;
        }
        return sum;
    }

    /**
     * The cumulative probability function.
     * The implementation calls the PMF for all values i from 0 to floor(k)
     * @param k the number of events in the interval
     * @return the cumulative probability of k events
     * https://en.wikipedia.org/wiki/Poisson_distribution
     */
    public double cdf(double k) {
        int my_k = (int)Math.floor(k);
        double sum = 0;
        for (int i = 0; i <= my_k; i ++)
            sum += p(i);
        return sum;
    }

    int DEFAULT_MAXK = 1000;

    /**
     *
     * @return
     */
    public int sample() {
        double toss = rand.nextDouble();
        double sum = 0;
        for (int i = 1; i <= DEFAULT_MAXK; i ++) {
            sum += p(i);
            if (sum > toss)
                return i;
        }
        return DEFAULT_MAXK;
    }

    /**
     * Estimate lambda for Zero-Truncated Poisson via Newton-Raphson
     * @param data dataset
     * @return instance of ZeroTruncatedPoisson with parameter value set by MLE
     */
    public static ZeroTruncatedPoisson fitMLE(int[] data) {
        double mean = Arrays.stream(data).average().orElse(0.0);

        // Initial guess: start with mean (approx for small lambda)
        double lambda = mean > 1 ? mean - 0.5 : 0.5;

        for (int iter = 0; iter < 100; iter++) {
            double expNegLambda = Math.exp(-lambda);
            double f = lambda / (1.0 - expNegLambda) - mean; // equation f(λ)=0

            // derivative f'(λ)
            double denom = (1.0 - expNegLambda);
            double fprime = (denom - lambda * expNegLambda) / (denom * denom);

            double step = f / fprime;
            lambda -= step;

            if (lambda <= 0) lambda = 1e-6; // keep positive
            if (Math.abs(step) < 1e-10) break;
        }
        return new ZeroTruncatedPoisson(lambda);
    }

    // Log-likelihood for ZTP
    public static double logLikelihood(int[] data, double lambda) {
        double logL = 0.0;
        for (int x : data) {
            logL += x * Math.log(lambda) - lambda - logFactorial(x);
        }
        logL -= data.length * Math.log(1.0 - Math.exp(-lambda));
        return logL;
    }

    public double getLogLikelihood(int[] data) {
        return logLikelihood(data, this.lambda);
    }

    // log(x!) using logGamma
    public static double logFactorial(int x) {
        return lgamma(x + 1.0);
    }

    /**
     * Returns an approximation of the log of the Gamma function of x. Laczos
     * Approximation Reference: Numerical Recipes in C
     * http://www.library.cornell.edu/nr/cbookcpdf.html
     */
    public static double lgamma(double x) {
        return Poisson.lgamma(x);
    }

    // Log posterior (up to constant)
    public static double logPosterior(double lambda, int[] data, double alpha, double beta) {
        return logLikelihood(data, lambda) + GammaDistrib.logPrior(lambda, alpha, beta);
    }

    /**
     * Estimate the parameter value of ZTP by MAP and return ZeroTruncatedPoisson instance.
     * Does this via Golden section search.
     * @param data dataset
     * @return an instance of ZeroTruncatedPoisson with parameter value
     */
    public static ZeroTruncatedPoisson fitMAP(int[] data, double alpha, double beta, double lower, double upper) {
        double lambda = goldenSectionSearch(data, alpha, beta, lower, upper,1e-6, 1000);
        return new ZeroTruncatedPoisson(lambda);
    }


    // Golden-section search for maximum
    public static double goldenSectionSearch(int[] data, double alpha, double beta,
                                             double lower, double upper, double tol, int maxIter) {
        double gr = (Math.sqrt(5) + 1) / 2; // golden ratio

        double c = upper - (upper - lower) / gr;
        double d = lower + (upper - lower) / gr;

        int iter = 0;
        while (Math.abs(c - d) > tol && iter < maxIter) {
            double fc = logPosterior(c, data, alpha, beta);
            double fd = logPosterior(d, data, alpha, beta);

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
        // Example dataset (ZTP observations)
        int[] data = {1,2,1,3,2,4,1,2,3,1,5};

        // Gamma prior hyperparameters
        double alpha = 2.0;
        double beta = 1.0;

        // Search interval for lambda
        double lower = 0.01;
        double upper = 10.0;

        // Run golden-section search
        ZeroTruncatedPoisson ztp = fitMAP(data, alpha, beta, lower, upper);
        System.out.printf("MAP gives " + ztp + "\n");
        for (int i = 0; i < data.length; i ++) {
            System.out.printf("%d ", ztp.sample());
        }
    }

    public static void main1(String[] args) {
        // Example dataset with no zeros
        int[] data = {1,2,1,3,2,4,1,2,3,1,5};

        // Run golden-section search
        ZeroTruncatedPoisson ztp = fitMLE(data);
        System.out.printf("MLE gives " + ztp);
        for (int i = 0; i < data.length; i ++) {
            System.out.printf("%d ", ztp.sample());
        }
        double logL = logLikelihood(data, ztp.lambda);
        System.out.println("Log-likelihood at λ̂ = " + logL);
    }

    public static void main0(String[] args) {
        for (int lambda = 1; lambda < 10; lambda ++) {
            ZeroTruncatedPoisson ZeroTruncatedPoisson = new ZeroTruncatedPoisson(lambda, 0);
            System.out.println("Lambda = " + lambda);
            int[] cnt = new int[ZeroTruncatedPoisson.DEFAULT_MAXK + 1];
            for (int i = 0; i < 100; i ++)
                cnt[ZeroTruncatedPoisson.sample()] += 1;
            for (int j = 0; j < 20; j ++)
                System.out.print(cnt[j] + "\t");
            System.out.println();
        }
    }

}
