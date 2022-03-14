package bn.prob;

import bn.Distrib;
import dat.Enumerable;
import java.io.Serializable;
import java.util.Random;

/**
 * The univariate Gaussian density function. Implementation is based on
 * Bayesian Logic (BLOG) inference engine version 0.3 (see copyright message
 * below).
 */
public class GaussianDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;

    public Double mu;
    public Double sigma; // standard deviation
    public Double sigmaSquared; // variance

    private double normConst;

    @Override
    public double get(Object input) {
        double y = this.getDensity((Double) input);
        return y;
    }

    @Override
    public String toString() {
        return String.format("%4.2f;%4.2f", mu, sigmaSquared);
    }

    public double getMean() {
        return mu;
    }

    public double getVariance() {
        return sigmaSquared;
    }

    private double logNormConst;

    private final double ROOT_2PI = Math.sqrt(2 * Math.PI);
    private final double LOG_ROOT_2PI = 0.5 * (Math.log(2) + Math.log(Math.PI));

    private Random rand = new Random();

    /**
     * Creates a univariate Gaussian distribution with the given fixed mean and
     * variance.
     * @param mean
     * @param variance
     */
    public GaussianDistrib(double mean, double variance) {
        setMean(mean);
        setVariance(variance);
    }

    /**
     * Returns the density of this Gaussian distribution at the given value.
     * This method should only be called if the mean and variance were set in
     * the constructor (internal calls are ok if the private method
     * <code>initParams</code> is called first).
     * @param x input value
     * @return density of x
     */
    public double getDensity(double x) {
        return (Math.exp(-Math.pow((x - mu), 2) / (2 * sigmaSquared)) / normConst);
    }

    /**
     * Returns the natural log of the density of this Gaussian distribution at
     * the given value. This method should only be called if the mean and
     * variance were set in the constructor, or if <code>initParams</code> has
     * been called.
     * @param x input value
     * @return the log of the density of x
     */
    public double getLogDensity(double x) {
        return (-Math.pow((x - mu), 2) / (2 * sigmaSquared)) - logNormConst;
    }

    /**
     * Returns a value sampled from this Gaussian distribution. This method
     * should only be called if the mean and variance were set in the
     * constructor (internal calls are ok if the private method
     * <code>initParams</code> is called first).
     *
     * <p>
     * The implementation uses the Box-Muller transformation [G.E.P. Box and
     * M.E. Muller (1958) "A note on the generation of random normal deviates".
     * Ann. Math. Stat 29:610-611]. See also
     * http://www.cs.princeton.edu/introcs/26function/MyMath.java.html
     */
    @Override
    public Double sample() {
        double U = rand.nextDouble();
        double V = rand.nextDouble();
        return Double.valueOf(mu + (sigma * Math.sin(2 * Math.PI * V) * Math.sqrt((-2 * Math.log(U)))));
    }

    public void setMean(double mean) {
        mu = mean;
    }

    public void setVariance(double variance) {
        sigmaSquared = variance;

        if (sigmaSquared <= 0) {
            throw new IllegalArgumentException("Variance of UnivarGaussian distribution must be positive, "
                    + "not " + sigmaSquared);
        }
        sigma = Math.sqrt(sigmaSquared);
        normConst = sigma * ROOT_2PI;
        logNormConst = (0.5 * Math.log(sigmaSquared)) + LOG_ROOT_2PI;
    }

    public void randomize(int seed) {
        Random rand1 = new Random(seed);
        setMean(rand1.nextGaussian());
        setVariance(rand1.nextDouble());
    }

    public void setSeed(long seed) {
        rand = new Random(seed);
    }
    /**
     * Create a density resembling those in the specified samples.
     * @param samples samples
     * @return a new distribution which could serve as a starting point for representing a sub-group of those in the samples.
     */
    public static GaussianDistrib probe(double[] samples) {
        double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY, mean = 0;
        for (int i = 0; i < samples.length; i ++) {
            if (samples[i] < min)
                min = samples[i];
            if (samples[i] > max)
                max = samples[i];
            mean += samples[i] / samples.length;
        }
        Random myrand = new Random();
        int choose = myrand.nextInt(samples.length);
        if (samples[choose] > mean)
            return new GaussianDistrib(samples[choose] - mean, myrand.nextDouble() * (max - mean));
        else 
            return new GaussianDistrib(mean - samples[choose], myrand.nextDouble() * (mean - min));
    }
    
    /**
     * Create a density based on the specified samples.
     * @param samples samples
     * @return a new distribution for the samples.
     */
    public static GaussianDistrib estimate(double[] samples) {
        if (samples.length == 0)
            return null;
        double mean = 0;
        for (int i = 0; i < samples.length; i ++) {
            mean += samples[i] / samples.length;
        }
        double diff = 0;
        for (int i = 0; i < samples.length; i ++) {
            diff += (mean - samples[i]) * (mean - samples[i]);
        }
        if (diff == 0)
            return null;
        return new GaussianDistrib(mean, diff / samples.length);
    }
    
    /**
     * Example that finds a mixture of Gaussians using Gibbs sampling.
     * @param args 
     */
    public static void main(String[] args) {
        int clusters = 2;
        Enumerable label = new Enumerable(clusters);
        EnumDistrib C = EnumDistrib.random(label); // mixture distrib
        GaussianDistrib[] Q = new GaussianDistrib[clusters];
        double[] X = {2.9, 3.3, 4.3, 6.1, 6.4, 7.5, 8.1, 8.8};
        
        // initialisation
        Random r = new Random();

        int[] Z_sample = new int[X.length];
        EnumDistrib[] Z = new EnumDistrib[X.length];
        for (int k = 0; k < clusters; k ++) {
            Q[k] = GaussianDistrib.probe(X);
            Q[k].setVariance(1);
            System.out.println("P(Q|Z = " + k + ") = " + Q[k]);
        }
        double[] counts = new double[clusters];
        for (int i = 0; i < X.length; i ++) {
            Z[i] = new EnumDistrib(label);
            Z_sample[i] = r.nextInt(clusters);
            counts[Z_sample[i]] += 1.0;
            System.out.println("X_" + i + " = " + X[i] + ": Z = " + Z_sample[i]);
        }
        C.set(counts);
        for (int k = 0; k < clusters; k ++) {
//            System.out.println("P(Z = " + k + ") = " + C.get(k));
        }
            
        // training
        
        for (int round = 0; round < 10; round ++) {
            int i_z = r.nextInt(X.length);
            System.out.println("Hold-out X_" + i_z + " = " + X[i_z]);
            
            for (int k = 0; k < clusters; k ++) {
                System.out.println("\tSamples labelled Z = " + k + ": ");
                double sum = 0.0;
                int cnt = 0;
                for (int i = 0; i < X.length; i ++) {
                    if (i != i_z && k == Z_sample[i]) {
                        sum += X[i];
                        System.out.println("\t\tX_" + i + " = " + X[i]);
                        cnt ++;
                    }
                }
                double mean = sum / cnt;
                if (cnt > 0)
                    Q[k].setMean(mean);
                double diff = 0;
                for (int i = 0; i < X.length; i ++) {
                    if (i != i_z && k == Z_sample[i]) {
                        diff += (mean - X[i]) * (mean - X[i]);
                    }
                }
                double var = diff / cnt;
                //if (cnt > 1)
                //    Q[k].setVariance(var);
                System.out.println("\t\tNew P(Q|Z = " + k + ") = " + Q[k]);
            }
            
            int i = i_z;
            double[] p = new double[clusters];
            for (int k = 0; k < clusters; k ++) {
                p[k] = Q[k].getDensity(X[i]);
                //p[k] = Q[k].getDensity(X[i]) * C.get(k);
            }
            Z[i].set(p);
            Z_sample[i] = (int)Z[i].sample();
            System.out.println("\tX_" + i + " = " + X[i] + ": P(Z|Q = " + X[i] + ") = " + Z[i] + " =sample=> " + Z_sample[i]);
            counts = new double[clusters];
            for (i = 0; i < X.length; i ++) {
                counts[Z_sample[i]] += 1.0;
            }
            C.set(counts);
        }
        
        System.out.println("RESULT:");
        for (int i = 0; i < X.length; i ++) {
            System.out.println("X_" + i + " = " + X[i] + ": P(Z|Q = " + X[i] + ") = " + Z[i] + " =sample=> " + Z_sample[i]);
        }
        
    }
}
/*
 * Copyright (c) 2005, Regents of the University of California
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in
 *   the documentation and/or other materials provided with the
 *   distribution.  
 *
 * * Neither the name of the University of California, Berkeley nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior 
 *   written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
