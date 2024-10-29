package stats;

public interface IndelModel {
    /**
     * Samples a value from the distribution
     * @return the sampled value
     */
    int sample();

    /**
     * Computes the probability mass function (PMF) for a given value
     * @param k the value
     * @return the probability of k
     */
    double p(int k);

    /**
     * Computes the cumulative distribution function (CDF) for a given value
     * @param k the value
     * @return the cumulative probability up to k
     */
    double cdf(int k);

    /**
     * Returns a string representation of the distribution as it should be specified on the TrAVIS command line
     * @return the text string
     */
    String getTrAVIS();

}
