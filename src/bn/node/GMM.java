package bn.node;

import java.util.Arrays;
import java.util.Random;

public class GMM {
    private int numComponents;
    private double[] weights;
    private double[][] means;
    private double[][] variances;

    public GMM(int numComponents) {
        this.numComponents = numComponents;
        this.weights = new double[numComponents];
        this.means = new double[numComponents][];
        this.variances = new double[numComponents][];
    }

    public void train(double[][] data, int maxIterations, double tolerance) {
        int numSamples = data.length;
        int numFeatures = data[0].length;

        // Initialize model parameters
        Random random = new Random();
        for (int i = 0; i < numComponents; i++) {
            weights[i] = 1.0 / numComponents;
            means[i] = new double[numFeatures];
            variances[i] = new double[numFeatures];
            for (int j = 0; j < numFeatures; j++) {
                means[i][j] = random.nextDouble(); // Random initialization
                variances[i][j] = random.nextDouble(); // Random initialization
            }
        }

        double prevLogLikelihood = Double.NEGATIVE_INFINITY;
        for (int iteration = 0; iteration < maxIterations; iteration++) {
            // E-step
            double[][] responsibilities = new double[numSamples][numComponents];
            for (int i = 0; i < numSamples; i++) {
                double sum = 0.0;
                for (int j = 0; j < numComponents; j++) {
                    responsibilities[i][j] = weights[j] * calculateGaussian(data[i], means[j], variances[j]);
                    sum += responsibilities[i][j];
                }
                for (int j = 0; j < numComponents; j++) {
                    responsibilities[i][j] /= sum; // normalising the responsibilities across components
                }
            }

            // M-step
            double[] Nk = new double[numComponents]; // a total of "responsibilities" (one for each component)
            for (int j = 0; j < numComponents; j++) {
                for (int i = 0; i < numSamples; i++) {
                    Nk[j] += responsibilities[i][j];
                }
                for (int f = 0; f < numFeatures; f++) {
                    double meanSum = 0.0;
                    double varianceSum = 0.0;
                    for (int i = 0; i < numSamples; i++) {
                        meanSum += responsibilities[i][j] * data[i][f];
                        varianceSum += responsibilities[i][j] * Math.pow(data[i][f] - means[j][f], 2);
                    }
                    means[j][f] = meanSum / Nk[j];
                    variances[j][f] = varianceSum / Nk[j];
                }
                weights[j] = Nk[j] / numSamples;
            }

            // Compute log-likelihood
            double logLikelihood = 0.0;
            for (int i = 0; i < numSamples; i++) {
                double sampleLogLikelihood = 0.0;
                for (int j = 0; j < numComponents; j++) {
                    sampleLogLikelihood += weights[j] * calculateGaussian(data[i], means[j], variances[j]);
                }
                logLikelihood += Math.log(sampleLogLikelihood);
            }

            // Convergence check
            if (Math.abs(logLikelihood - prevLogLikelihood) < tolerance) {
                break;
            }
            prevLogLikelihood = logLikelihood;
        }
    }

    private double calculateGaussian(double[] x, double[] mean, double[] variance) {
        double result = 1.0;
        for (int i = 0; i < x.length; i++) {
            result *= Math.exp(-0.5 * Math.pow(x[i] - mean[i], 2) / variance[i]) / (Math.sqrt(2 * Math.PI * variance[i]));
        }
        return result;
    }

    public static void main(String[] args) {
        int numComponents = 3;
        GMM gmm = new GMM(numComponents);

        // Sample data (you should replace this with your actual data)
        double[][] data = {
                {2.0},
                {3.5},
                {1.5},
                {7.0},
                {8.0},
                {2.2},
                {7.7}
        };

        double[][] data2 = {
                {8.547969513147832}, {11.055365021610807},{9.515919623038354},{7.309429864631631},{5.155113714709337},{5.9934818128526},{9.737441853198932},{5.866697379742977},{8.730212309448154},{5.107871539173301},{7.091657649202554},{6.005496193011316},{5.784554136706749},{7.24665515021407},{7.588302408954381},{4.654811488291314}};

        int maxIterations = 100;
        double tolerance = 1e-4;
        gmm.train(data2, maxIterations, tolerance);

        // Print the learned parameters
        System.out.println("Weights: " + Arrays.toString(gmm.weights));
        System.out.println("Means: " + Arrays.deepToString(gmm.means));
        System.out.println("Variances: " + Arrays.deepToString(gmm.variances));
    }
}
