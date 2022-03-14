package stats;

import java.util.Collection;

/**
 * Compute the z-score.
 * @author m.boden
 */
public class ZScore {
    private final double mean;
    private final double var;
    private final double sd;

    /**
     * Prepare the Z-score by supplying the population from which the sample
     * mean and deviation is determined
     *
     * @param x the population
     */
    public ZScore(double[] x) {
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i];
        }
        mean = sum / x.length;
        double dev = 0;
        for (int i = 0; i < x.length; i++) {
            dev += (x[i] - mean) * (x[i] - mean);
        }
        if (dev == 0) {
            throw new RuntimeException("Invalid population: zero variance");
        }
        var = dev / x.length;
        sd = Math.sqrt(var);
    }
    
    public ZScore(double mean, double var) {
        this.mean = mean;
        this.var = var;
        this.sd = Math.sqrt(var);
    }

    private static double[] toArr(Collection<Double> x) {
        double[] arr = new double[x.size()];
        int i = 0;
        for (double value : x) {
            arr[i] = value;
            i ++;
        }
        return arr;
    }
    
    public ZScore(Collection<Double> x) {
        this(toArr(x));
    }
    
    public double getMean() {
        return mean;
    }

    public double getSD() {
        return sd;
    }

    public double getVariance() {
        return var;
    }

    /**
     * Compute the z-score from the mean and standard deviation.
     *
     * @param x the value for which a z-score is requested
     * @param mean the mean
     * @param sd the standard deviation
     * @return the z-score for x
     */
    public static double getZ(double x, double mean, double sd) {
        return (x - mean) / sd;
    }

    /**
     * Compute the z-score on basis of the population provided through the
     * constructor
     *
     * @param x the value for which a z-score is requested
     * @return the z-score for x
     */
    public double getZ(double x) {
        return (x - mean) / sd;
    }
    
    /**
     * Compute P-value on basis of Z.
     * @param x observed value
     * @param aboveMean true if test is cumulative probability above mean, false if below mean
     * @return the P-value
     */
    public double getP(double x, boolean aboveMean) {
        double z = getZ(x);
        if (aboveMean)
            return 1.0 - NormalDistribution.f(z);
        else
            return NormalDistribution.f(z);
    }

    public static void main(String[] args) {
        double[] X = {12.1, 11.2, 12.3, 11.8, 11.2, 12.3, 11.1, 13.2, 12.3, 11.6, 10.8};
        ZScore zscore = new ZScore(X);
        double x = 13.6;
        double z = zscore.getZ(x);
        System.out.println(z + "\t" + NormalDistribution.f(z));
    }
	
}
