package stats;
import java.util.Random;

public class Lavalette implements IndelModel {

    private double a;  // Exponent parameter of the Lavalette distribution
    private int maxK;  // Maximum rank value (N in the formula)
    private Random rand;
    private long seed;

    public Lavalette(double a,long seed, int maxK) {
        this.a = a;
        this.maxK = maxK;
        this.seed = seed;
        this.rand = new Random(seed);
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

    public static void main(String[] args) {
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

