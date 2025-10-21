package stats;

import bn.prob.GammaDistrib;

/**
 * Interface for rate models, e.g. for substitution, or insertions/deletions.
 * Used by TrAVIS.
 */
public interface RateModel {

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

    public static RateModel create(String distrib_name, String params) {
        double[] params_arr = parseParams(params);
        switch (distrib_name) {
            case "ZIG":
            case "ZeroInflatedGamma":
            case "zeroinflatedgamma":
            case "zig":
                if (params_arr.length == 3)
                    return new ZeroInflatedGamma(params_arr[0], params_arr[1], params_arr[2]);
                else if (params_arr.length == 2)
                    return new ZeroInflatedGamma(1.0, params_arr[0], params_arr[1]);
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            case "Gamma":
            case "gamma":
                if (params_arr.length == 2)
                    return new GammaDistrib(params_arr[0], 1 / params_arr[1]);
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            case "MixtureGamma":
            case "mixturegamma":
                if (params_arr.length % 3 == 0) { // must be evenly divisible by 3
                    GammaDistrib[] gammas = new GammaDistrib[params_arr.length / 3];
                    double[] priors = new double[params_arr.length / 3];
                    for (int i = 0; i < params_arr.length; i += 3) { // three params for each gamma distrib
                        gammas[i/3] = new GammaDistrib(params_arr[i], 1 / params_arr[i + 1]);
                        priors[i/3] = params_arr[i + 2];
                    }
                    return new GammaDistrib.Mixture(gammas, priors);
                }
                throw new RuntimeException("Failed to parse parameters \"" + params + "\" for nominated distribution " + distrib_name);
            default:
                throw new RuntimeException("Invalid distribution " + distrib_name);
        }
    }

    /**
     * Find the distribution with maximum data likelihood
     * @param rate_data
     * @return the best model
     */
    static RateModel bestfit(double[] rate_data, long seed) {
        RateModel[] models = new RateModel[] {
                ZeroInflatedGamma.fitMLE(rate_data, seed),
                ZeroInflatedGamma.Mixture.fitMLE(rate_data, 2, seed),
                ZeroInflatedGamma.Mixture.fitMLE(rate_data, 3, seed),
        };
        double best_ll = Double.MIN_VALUE;
        RateModel best_model = null;
        int best_idx = 0;
        for (int i = 0; i < models.length; i++) {
            RateModel model = models[i];
            double ll = model.getLogLikelihood(rate_data);
            if (ll > best_ll) {
                best_ll = ll;
                best_idx = i;
            }
        }
        return models[best_idx];
    }

    static RateModel bestfit(String distrib_name, double[] rate_data, long seed) {
        switch (distrib_name) {
            case "ZIG":
            case "ZeroInflatedGamma":
            case "zeroinflatedgamma":
            case "zig":
                return ZeroInflatedGamma.fitMLE(rate_data, seed);
            case "ZIG2":
            case "zig2":
            case "MixtureZIG":
            case "MixtureZeroInflatedGamma":
            case "mixturezeroinflatedgamma":
            case "mixturezig":
            case "MixtureZIG2":
            case "MixtureZeroInflatedGamma2":
            case "mixturezeroinflatedgamma2":
            case "mixturezig2":
                return ZeroInflatedGamma.Mixture.fitMLE(rate_data, 2, seed);
            case "ZIG3":
            case "zig3":
            case "MixtureZIG3":
            case "MixtureZeroInflatedGamma3":
            case "mixturezeroinflatedgamma3":
            case "mixturezig3":
                return ZeroInflatedGamma.Mixture.fitMLE(rate_data, 3, seed);
            case "Gamma":
            case "gamma":
                return GammaDistrib.fitMLE(rate_data, seed);
            case "MixtureGamma":
            case "mixturegamma":
            case "MixtureGamma2":
            case "mixturegamma2":
                return GammaDistrib.Mixture.fitMLE(rate_data, 2, seed);
            case "MixtureGamma3":
            case "mixturegamma3":
                return GammaDistrib.Mixture.fitMLE(rate_data, 3, seed);
            default:
                throw new RuntimeException("Invalid distribution " + distrib_name);
        }
    }

    /**
     * Samples a rate value from the distribution
     * @return the sampled value
     */
    Double sample();

    /**
     * Computes the probability mass function (PMF) for a given value
     * @param rate the rate value
     * @return the probability of k
     */
    double p(double rate);

    /**
     * Computes the cumulative distribution function (CDF) for a given value
     * @param rate the rate value
     * @return the cumulative probability up to rate
     */
    double cdf(double rate);

    void setSeed(long seed);

    /**
     * Returns a string representation of the distribution as it should be specified on the TrAVIS command line
     * @return the text string
     */
    String getTrAVIS();

    /**
     * Calculate the log likelihood of the data given the model/distribution
     * @param data dataset
     * @return the log-likelihood of the data given the model
     */
    public double getLogLikelihood(double[] data);


}
