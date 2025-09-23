package stats;
import bn.prob.GammaDistrib;

import java.util.Arrays;
import java.util.Random;

public class Lavalette implements IndelModel {

    private double a;  // Exponent parameter of the Lavalette distribution
    private int maxK;  // Maximum rank value (N in the formula)
    private Random rand;

    public Lavalette(double a, long seed, int maxK) {
        this.a = a;
        this.maxK = maxK;
        this.rand = new Random(seed);
    }

    public Lavalette(double a, int maxK) {
        this(a, System.currentTimeMillis(), maxK);
    }

    public static int DEFAULT_MAXK = 1000; //

    public Lavalette(double a) {
        this(a, System.currentTimeMillis(), DEFAULT_MAXK);
    }

    public String toString() {
        return "Lavalette(lambda=" + a + ", maxK=" + maxK + ")";
    }

    // Probability mass function (PMF) for Lavalette distribution
    public double p(int r) {
        return Math.pow(maxK + 1 - r, a) * Math.pow(r, -a) / lavaletteZeta();
    }

    // Calculate the normalization constant for Lavalette distribution
    private double lavaletteZeta() {
        double sum = 0.0;
        for (int r = 1; r <= maxK; r++) {
            sum += Math.pow(maxK + 1 - r, a) * Math.pow(r, -a);
        }
        return sum;
    }

    // Sample from Lavalette distribution
    public int sample() {
        double toss = rand.nextDouble();
        double sum = 0.0;
        for (int r = 1; r <= maxK; r++) {
            sum += p(r);
            if (sum > toss) {
                return r;
            }
        }
        return maxK;
    }

    // Cumulative distribution function (CDF) for Lavalette distribution
    public double cdf(int r) {
        if (r < 1) return 0.0;
        if (r >= maxK) return 1.0;

        double sum = 0.0;
        for (int i = 1; i <= r; i++) {
            sum += p(i);
        }
        return sum;
    }

    @Override
    public String getTrAVIS() {
        if (this.maxK == DEFAULT_MAXK)
            return String.format("Lavalette:%.3f", a);
        else
            return String.format("Lavalette:%.3f,%d", a, maxK);
    }

    public double getLogLikelihood(int[] data) {
        return Lavalette.logLikelihood(data, this.a, this.maxK);
    }

    @Override
    public void setSeed(long seed) {
        this.rand = new Random(seed);
    }

    // Compute Z(a)
    // Normalization constant: sum_{k=1}^N ((N+1-k)/k)^a
    public static double Z(int N, double a) {
        double sum = 0.0;
        for (int k = 1; k <= N; k++) {
            double ratio = (double)(N + 1 - k) / k;
            sum += Math.pow(ratio, a);
        }
        return sum;
    }

    // Compute log-likelihood of Lavalette(a, N)
    public static double logLikelihood(int[] data, double a, int N) {
        double sumLogs = 0.0;
        for (int x : data) {
            sumLogs += Math.log((double)(N + 1 - x) / x);
        }

        double denom = Z(N, a);

        return a * sumLogs - data.length * Math.log(denom);
    }

    // Derivative Z'(a)
    public static double Zprime(int N, double a) {
        double sum = 0.0;
        for (int k = 1; k <= N; k++) {
            double ratio = (double)(N + 1 - k) / k;
            sum += Math.log(ratio) * Math.pow(ratio, a);
        }
        return sum;
    }

    // Second derivative Z''(a)
    public static double Zsecond(int N, double a) {
        double sum = 0.0;
        for (int k = 1; k <= N; k++) {
            double ratio = (double)(N + 1 - k) / k;
            sum += Math.pow(Math.log(ratio), 2) * Math.pow(ratio, a);
        }
        return sum;
    }

    /**
     * Estimate the parameter value of Lavalette by MLE.
     * Does this via Newton-Raphson.
     * @param data dataset
     * @return an instance of Lavalette with parameter value
     */
    public static Lavalette fitMLE(int[] data) {
        return fitMLE(data, Arrays.stream(data).max().getAsInt());
    }

    /**
     * Estimate the parameter value of Lavalette by MLE.
     * Does this via Newton-Raphson.
     * @param data dataset
     * @param N max rank
     * @return an instance of Lavalette with parameter value
     */
    // Estimate parameter a
    public static Lavalette fitMLE(int[] data, int N){
        int n = data.length;
        // Precompute average log term
        double avgLogTerm = Arrays.stream(data)
                .mapToDouble(x -> Math.log((double)(N + 1 - x) / x))
                .average()
                .orElse(0.0);

        double a = 1.0; // initial guess
        for (int iter = 0; iter < 100; iter++) {
            double Z = Z(N, a);
            double Zp = Zprime(N, a);
            double Zpp = Zsecond(N, a);

            double score = n * (avgLogTerm - Zp / Z);
            double info = -n * ((Zpp * Z - Zp * Zp) / (Z * Z)); // second derivative of log-likelihood

            double step = score / info;
            a -= step;

            if (Math.abs(step) < 1e-8) break; // convergence
        }
        return new Lavalette(a, N);
    }

    // Log posterior (up to constant)
    public static double logPosterior(double s, int[] data, int N, double alpha, double beta) {
        if (s <= 0) return Double.NEGATIVE_INFINITY;

        double sumLogTerms = 0.0;
        for (int x : data) {
            sumLogTerms += Math.log(N + 1 - x);
        }
        int n = data.length;

        double norm = Z(N, s);
        if (norm <= 0) return Double.NEGATIVE_INFINITY;

        double logLikelihood = s * sumLogTerms - n * Math.log(norm);
        double logPrior = GammaDistrib.logPrior(s, alpha, beta);

        return logLikelihood + logPrior;
    }

    /**
     * Estimate the parameter value of Lavalette by MAP and return Lavalette instance.
     * Does this via Golden section search.
     * @param data dataset
     * @return an instance of Lavalette with parameter value
     */
    public static Lavalette fitMAP(int[] data, double alpha, double beta, double lower, double upper) {
        int N = Arrays.stream(data).max().getAsInt();
        double s = goldenSectionSearch(data, N, alpha, beta, lower, upper,1e-6, 1000);
        return new Lavalette(s, N);
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
        int[] data = {1,1,2,3,2,5,1,2,4,2,1,3,6,2,1};
        int N = Arrays.stream(data).max().getAsInt();

        // Gamma prior hyperparameters
        double alpha = 2.0;
        double beta = 1.0;

        // Search interval for s
        double lower = 0.1;
        double upper = 5.0;


        // Run golden-section search
        Lavalette lava = fitMAP(data, alpha, beta, lower, upper);
        System.out.printf("MAP gives " + lava);
        for (int i = 0; i < data.length; i ++) {
            System.out.printf("%d ", lava.sample());
        }
    }

    public static void main1(String[] args) {
        // Example dataset (ranks from Lavalette-like distribution)
        int[] data = {1,1,2,3,2,5,1,2,4,2,1,3,6,2,1};

        // Maximum observed rank
        int N = Arrays.stream(data).max().getAsInt();

        Lavalette lava = fitMLE(data, N);
        System.out.println("MLE gives " + lava);
        for (int i = 0; i < data.length; i ++) {
            System.out.printf("%d ", lava.sample());
        }
    }

    public static void main0(String[] args) {
        for (double a = 1.0; a < 5.5; a += 0.5) {
            Lavalette lavalette = new Lavalette(a, 0, 100);
            System.out.println("a = " + a);
            int[] cnt = new int[lavalette.maxK + 1];
            for (int i = 0; i < 1000; i++) {
                cnt[lavalette.sample()] += 1;
            }
            for (int j = 1; j <= 20; j++) {
                System.out.print(cnt[j] + "\t");
            }
            System.out.println();
        }
    }
}

