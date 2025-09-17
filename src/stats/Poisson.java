package stats;

import java.util.Arrays;
import java.util.Random;

/**
 * The discrete Poisson distribution
 * Created by mikael on 31/3/17.
 */
public class Poisson implements IndelModel {

    private double lambda;
    private Random rand = null;


    /**
     * Define a Poisson distribution
     *
     * @param lambda the average number of events per interval
     * @param seed the random seed in case it is used for sampling
     */
    public Poisson(double lambda, long seed) {
        this.lambda = lambda;
        this.rand = new Random(seed);
    }

    /**
     * Define a Poisson distribution
     *
     * @param lambda the average number of events per interval
     */
    public Poisson(double lambda) {
        this(lambda, System.currentTimeMillis());
    }

    public String toString() {
        return "Poisson(lambda=" + lambda + ")";
    }

    @Override
    public String getTrAVIS() {
        return String.format("Poisson:%.3f", lambda);
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
        return Math.exp(k * Math.log(lambda) - lambda - lgamma(k + 1));
    }

    /**
     * The probability mass function.
     * The implementation is numerically stable and uses the log gamma
     * (https://en.wikipedia.org/wiki/Poisson_distribution)
     * @param k the number of events in the interval
     * @return the probability of k events
     */
    public double p(int k) {
        return Math.exp(k * Math.log(lambda) - lambda - lgamma(k + 1));
    }

    /**
     * The cumulative probability function.
     * The implementation calls the PMF for all values i from 0 to floor(k)
     * @param k the number of events in the interval
     * @return the cumulative probability of k events
     * https://en.wikipedia.org/wiki/Poisson_distribution
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
        for (int i = 0; i <= DEFAULT_MAXK; i ++) {
            sum += p(i);
            if (sum > toss)
                return i;
        }
        return DEFAULT_MAXK;
    }

    /**
     * Estimate lambda via MLE
     * @param data dataset
     * @return an instance of Poisson with parameter value from MLE
     */
    // Estimate Poisson parameter λ via MLE (sample mean)
    public static Poisson fitMLE(int[] data) {
        double lambda = Arrays.stream(data).average().orElse(0.0);
        return new Poisson(lambda);
    }

    /** Container class for the posterior; when prior is Gamma, the posterior is Gamma
     *
     */
    static class Posterior {
        double alpha;   // shape
        double beta;    // rate
        double mean;
        double variance;
        double map;

        Posterior(double alpha, double beta, int sum, int n) {
            this.alpha = alpha + sum;    // posterior shape
            this.beta = beta + n;        // posterior rate
            this.mean = this.alpha / this.beta;
            this.variance = this.alpha / (this.beta * this.beta);
            if (this.alpha > 1) {
                this.map = (this.alpha - 1) / this.beta;
            } else {
                this.map = 0.0; // mode at 0 if shape ≤ 1
            }
        }
    }

    // Compute posterior given prior and data
    public static Posterior posterior(double alphaPrior, double betaPrior, int[] data) {
        int sum = 0;
        for (int x : data) sum += x;
        int n = data.length;
        return new Posterior(alphaPrior, betaPrior, sum, n);
    }

    /**
     * Estimate the parameter value of Poisson by MAP and return Poisson instance.
     * @param data dataset
     * @return an instance of Poisson with parameter value
     */
    public static Poisson fitMAP(int[] data, double alpha, double beta) {
        return new Poisson(posterior(alpha, beta, data).map);
    }

    /**
     * Find the distribution with maximum data likelihood
     * @param indel_data
     * @return the best model
     */
    public static IndelModel bestfit(int[] indel_data) {
        return fitMLE(indel_data);
    }

    // Example usage
    public static void main(String[] args) {
        // Prior Gamma(alpha=2, beta=1)
        double alphaPrior = 2.0;
        double betaPrior = 2.0;

        // Observed Poisson data
        // Example dataset (Poisson-like counts)
        int[] data = {2, 3, 4, 2, 1, 0, 3, 2, 5, 4, 3, 2};

        Poisson poisson = fitMAP(data, alphaPrior, betaPrior);
        System.out.println("MAP gives " + poisson);
        for (int i = 0; i < data.length; i ++) {
            System.out.printf("%d ", poisson.sample());
        }

    }

    public static void main1(String[] args) {
        // Example dataset (Poisson-like counts)
        int[] data = {2, 3, 4, 2, 1, 0, 3, 2, 5, 4, 3, 2};

        Poisson poisson = fitMLE(data);
        System.out.println("MLE gives " + poisson);
        for (int i = 0; i < data.length; i ++) {
            System.out.printf("%d ", poisson.sample());
        }
    }

    public static void main0(String[] args) {
        for (int lambda = 1; lambda < 10; lambda ++) {
            Poisson poisson = new Poisson(lambda, 0);
            System.out.println("Lambda = " + lambda);
            int[] cnt = new int[poisson.DEFAULT_MAXK + 1];
            for (int i = 0; i < 100; i ++)
                cnt[poisson.sample()] += 1;
            for (int j = 0; j < 20; j ++)
                System.out.print(cnt[j] + "\t");
            System.out.println();
        }
    }
    /**
     * Returns an approximation of the log of the Gamma function of x. Laczos
     * Approximation Reference: Numerical Recipes in C
     * http://www.library.cornell.edu/nr/cbookcpdf.html
     */
    public static double lgamma(double x) {
        double[] cof = { 76.18009172947146, -86.50532032941677, 24.01409824083091,
                -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
        double y, z, ser, tmp;
        y = x;
        tmp = x + 5.5;
        tmp -= ((x + 0.5) * Math.log(tmp));
        ser = 1.000000000190015;
        for (int j = 0; j < 6; j += 1) {
            y += 1;
            ser += (cof[j] / y);
        }
        return (-tmp + Math.log(2.5066282746310005 * ser / x));
    }

    // Log-likelihood of Poisson for given lambda
    public static double logLikelihood(int[] data, double lambda) {
        double logL = 0.0;
        for (int x : data) {
            logL += x * Math.log(lambda) - lambda - logFactorial(x);
        }
        return logL;
    }


    // Compute log(x!) using logGamma(x+1)
    // Lanczos approximation
    public static double logFactorial(int x) {
        return lgamma(x + 1.0);
    }

    public double getLogLikelihood(int[] data) {
        return logLikelihood(data, this.lambda);
    }




}
