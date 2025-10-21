package stats;

import bn.prob.GammaDistrib;

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

    public void setSeed(long seed);

    /**
     * Calculate the log likelihood of the data given the model/distribution
     * @param data dataset
     * @return the log-likelihood of the data given the model
     */
    public double getLogLikelihood(int[] data);


    public static double[] parseParams(String str) throws RuntimeException {
        try {
            String[] parts = str.split(",");
            double[] result = new double[parts.length];
            for (int i = 0; i < parts.length; i++) {
                result[i] = Double.parseDouble(parts[i].trim());
            }
            return result;
        } catch (NumberFormatException e) {
            throw new RuntimeException("Failed to parse doubles: " + str, e);
        }
    }

    static IndelModel create(String distrib_name, String params) {
        double[] params_arr = parseParams(params);
        switch (distrib_name) {
            case "ZTP":
            case "ZeroTruncatedPoisson":
            case "zerotruncatedpoisson":
            case "ztp":
                if (params_arr.length == 1)
                    return new ZeroTruncatedPoisson(params_arr[0]);
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            case "Poisson":
            case "poisson":
                if (params_arr.length == 1)
                    return new Poisson(params_arr[0]);
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            case "Lavalette":
            case "lavalette":
                if (params_arr.length == 1)
                    return new Lavalette(params_arr[0]);
                else if (params_arr.length == 2)
                    return new Lavalette(params_arr[0], (int) params_arr[1]);
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            case "Zipf":
            case "zipf":
                if (params_arr.length == 1)
                    return new Zipf(params_arr[0]);
                else if (params_arr.length == 2)
                    return new Zipf(params_arr[0], (int) params_arr[1]);
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            default:
                throw new RuntimeException("Invalid distribution " + distrib_name);
        }
    }

    /**
     * Find the distribution with maximum data likelihood
     * @param indel_data
     * @return the best model
     */
    static IndelModel bestfit(int[] indel_data, long seed) {
        IndelModel[] models = new IndelModel[] {
                Lavalette.fitMLE(indel_data, seed),
                Zipf.fitMLE(indel_data, seed),
                ZeroTruncatedPoisson.fitMLE(indel_data, seed),
                Poisson.fitMLE(indel_data, seed),
        };
        double best_ll = Double.MIN_VALUE;
        IndelModel best_model = null;
        int best_idx = 0;
        for (int i = 0; i < models.length; i++) {
            IndelModel model = models[i];
            double ll = model.getLogLikelihood(indel_data);
            if (ll > best_ll) {
                best_ll = ll;
                best_idx = i;
            }
        }
        return models[best_idx];
    }

    static IndelModel bestfit(String distrib_name, int[] indel_data, long seed) {
        switch (distrib_name) {
            case "ZTP":
            case "ZeroTruncatedPoisson":
            case "zerotruncatedpoisson":
            case "ztp":
                return ZeroTruncatedPoisson.bestfit(indel_data, seed);
            case "Poisson":
            case "poisson":
                return Poisson.bestfit(indel_data, seed);
            case "Lavalette":
            case "lavalette":
                return Lavalette.bestfit(indel_data, seed);
            case "Zipf":
            case "zipf":
                return Zipf.bestfit(indel_data, seed);
            default:
                throw new RuntimeException("Invalid distribution " + distrib_name);
        }
    }

}
