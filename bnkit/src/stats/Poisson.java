package stats;

import java.util.Random;

/**
 * The discrete Poisson distribution
 * Created by mikael on 31/3/17.
 */
public class Poisson {

    private double lambda;
    private Random rand = null;
    private long seed = System.currentTimeMillis();

    /**
     * Define a Poisson distribution
     *
     * @param lambda the average number of events per interval
     */
    public Poisson(double lambda) {
        this.lambda = lambda;
    }


    /**
     * Define a Poisson distribution
     *
     * @param lambda the average number of events per interval
     * @param seed the random seed in case it is used for sampling
     */
    public Poisson(double lambda, long seed) {
        this.lambda = lambda;
        this.seed = seed;
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

    int MAX_K = 1000;

    /**
     *
     * @return
     */
    public int sample() {
        if (rand == null)
            rand = new Random(seed);
        double toss = rand.nextDouble();
        double sum = 0;
        for (int i = 0; i <= MAX_K; i ++) {
            sum += p(i);
            if (sum > toss)
                return i;
        }
        return MAX_K;
    }

    public static void main(String[] args) {
        for (int lambda = 1; lambda < 10; lambda ++) {
            Poisson poisson = new Poisson(lambda, 0);
            System.out.println("Lambda = " + lambda);
            int[] cnt = new int[poisson.MAX_K + 1];
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


}
