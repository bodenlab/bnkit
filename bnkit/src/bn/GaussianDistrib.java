package bn;

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

    private final Random rand = new Random();

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

    private void setMean(double mean) {
        mu = mean;
    }

    private void setVariance(double variance) {
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
