package bn.prob;


import bn.Distrib;
import java.io.Serializable;
import java.util.*;

/**
 * A Gamma distribution with shape parameter k and scale parameter 1/lambda.
 * Defined as f(x) = (lambda*e^(-lambda*x)*(lambda*x)^(k - 1)) / Gamma(k) where
 * Gamma(k) = integral from 0 to infinity of t^(k-1) * e^(-t) dt
 * Implementation is partly based on Bayesian Logic (BLOG) inference engine version 0.7 
 * (see copyright message below). The digamma and trigamma implementations are 
 * under an Apache Commons License as displayed below.
 */

public class GammaDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;

    public static final double GAMMA = 0.577215664901532860606512090082;

    // limits for switching algorithm in digamma
    private static final double C_LIMIT = 49;
    private static final double S_LIMIT = 1e-5;

    private double lambda; // same as 1 / beta; see http://en.wikipedia.org/wiki/Gamma_distribution
    private double k; // alpha = k
    private Random rand = new Random();

    public void setSeed(long seed) {
        rand = new Random(seed);
    }

    /**
     * Creates a new Gamma distribution with parameters k and lambda.
     * @param k shape (same as a)
     * @param lambda rate, which is the inverse of scale (or 1/b)
     */
    public GammaDistrib(double k, double lambda) {
        this.k = k;
        this.lambda = lambda;
    }

    /**
     * Returns the probability of x under this distribution
     */
    public double get(Object xobj) {
        double x = (Double) xobj;
        return (lambda * Math.exp(-lambda * x) * Math.pow(lambda * x, k - 1) / gamma(k));
    }

    /**
     * Returns a double sampled according to this distribution. Uniformly fast for
     * all k > 0. (Reference: Non-Uniform Random Variate Generation, Devroye
     * http://cgm.cs.mcgill.ca/~luc/rnbookindex.html) Uses Cheng's rejection
     * algorithm (GB) for k>=1, rejection from Weibull distribution for 0 < k < 1.
     */
    public Double sample() {
        boolean accept = false;
        if (k >= 1) {
            // Cheng's algorithm
            double b = (k - Math.log(4));
            double c = (k + Math.sqrt(2 * k - 1));
            double lam = Math.sqrt(2 * k - 1);
            double cheng = (1 + Math.log(4.5));
            double u, v, x, y, z, r;
            do {
                u = rand.nextDouble();
                v = rand.nextDouble();
                y = ((1 / lam) * Math.log(v / (1 - v)));
                x = (k * Math.exp(y));
                z = (u * v * v);
                r = (b + (c * y) - x);
                if ((r >= ((4.5 * z) - cheng)) || (r >= Math.log(z))) {
                        accept = true;
                }
            } while (!accept);
            return (x / lambda);
        } else {
            // Weibull algorithm
            double c = (1 / k);
            double d = ((1 - k) * Math.pow(k, (k / (1 - k))));
            double u, v, z, e, x;
            do {
                u = rand.nextDouble();
                v = rand.nextDouble();
                z = -Math.log(u); // generating random exponential variates
                e = -Math.log(v);
                x = Math.pow(z, c);
                if ((z + e) >= (d + x)) {
                    accept = true;
                }
            } while (!accept);
            return (x / lambda);
        }
    }

    /*
     * Returns an approximation of the Gamma function of x r(x) = integral from 0
     * to infinity of (t^(x-1) * e^(-t) dt) with |error| < 2e-10. Laczos
     * Approximation Reference: Numerical Recipes in C
     * http://www.library.cornell.edu/nr/cbookcpdf.html
     */
    public static double gamma(double x) {
        return Math.exp(lgamma(x));
    }

    /*
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
    
    /*
     * used for parameter learning
     * */
    public double getK() {
    	return k;
    }
    
    public void setK(double _k) {
    	k = _k;
    }
    
    public double getLambda() {
    	return lambda;
    }
    
    public void setLambda(double _lambda) {
    	lambda = _lambda;
    }

    public double getAlpha() {
        return k;
    }
    
    public double getBeta() {
        return 1.0 / lambda;
    }
    
    public void setBeta(double beta) {
        lambda = 1.0 / beta;
    }
    
    public void setAlpha(double alpha) {
        k = alpha;
    }
    
    /**
     * The log likelihood of the data X given the gamma distribution.
     * Based on Thomas Minka "Estimating a Gamma distribution" 2002.
     * http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
     * @param X data
     * @return the log likelihood of the data; log p(X|alpha, beta)
     */
    public double logLikelihood(double[] X) {
        double a = getAlpha();
        double b = getBeta();
        int n = X.length;
        double x_mean = 0;
        double log_x_mean = 0;
        for (double xi : X) {
            x_mean += xi;
            log_x_mean += Math.log(xi);
        }
        x_mean /= n;
        log_x_mean /= n;
        return n * (a - 1.0) * log_x_mean - n * lgamma(a) - n * a * Math.log(b) - n * x_mean / b;
    }
    
    /**
     * Estimate the parameters of a gamma distribution from data.
     * Specifically the implementation estimates alpha, and beta is given by
     * beta = mean(x) / alpha.
     * It uses a fast approximation based on Thomas Minka "Estimating a Gamma distribution" 2002.
     * http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
     * @param X data
     * @return alpha parameter (same as lambda here)
     */
    public static double getAlpha(double[] X) {
        double delta = 1;
        int n = X.length;
        double x_mean = 0;
        double log_x_mean = 0;
        for (double xi : X) {
            x_mean += xi;
            log_x_mean += Math.log(xi);
        }
        x_mean /= n;
        log_x_mean /= n;
        double a = 0.5 / (Math.log(x_mean) - log_x_mean); // good starting point (see Minka 2002)
        double a_inv = 1.0 / a;
        for (int r = 0; r < 10; r ++) { // max 10 iterations, should converge in ~4
            double numerator = log_x_mean - Math.log(x_mean) + Math.log(a) - digamma(a);
            double denominator = a * a * (1.0 / a - trigamma(a));
            double next_inv = a_inv + numerator / denominator;
            double next_a = 1.0 / next_inv;
            delta = Math.abs(next_a - a);
            a = next_a;
            if (delta < .01) {
                //System.out.println("Converged after " + (r + 1) + " rounds");
                break;
            }
        }
        return a;
    }
    
    /**
     * Get beta that maximises likelihood as computed for a specified alpha.
     * @param X data
     * @param alpha
     * @return beta
     */
    public static double getBeta(double[] X, double alpha) {
        int n = X.length;
        double x_mean = 0;
        for (double xi : X) {
            x_mean += xi;
        }
        x_mean /= n;
        return x_mean / alpha;
    }
    
    /*
    The following two methods: digamma and trigamma are ...

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
    
    The code has been modified as a result of the incorporation into bnkit.
    */
    /**
     * <p>Computes the digamma function of x.</p>
     *
     * <p>This is an independently written implementation of the algorithm described in
     * Jose Bernardo, Algorithm AS 103: Psi (Digamma) Function, Applied Statistics, 1976.</p>
     *
     * <p>Some of the constants have been changed to increase accuracy at the moderate expense
     * of run-time.  The result should be accurate to within 10^-8 absolute tolerance for
     * x >= 10^-5 and within 10^-8 relative tolerance for x > 0.</p>
     *
     * <p>Performance for large negative values of x will be quite expensive (proportional to
     * |x|).  Accuracy for negative values of x should be about 10^-8 absolute for results
     * less than 10^5 and 10^-8 relative for results larger than that.</p>
     *
     * @param x Argument.
     * @return digamma(x) to within 10-8 relative or absolute error whichever is smaller.
     * @see <a href="http://en.wikipedia.org/wiki/Digamma_function">Digamma</a>
     * @see <a href="http://www.uv.es/~bernardo/1976AppStatist.pdf">Bernardo&apos;s original article </a>
     * @since 2.0
     */
    public static double digamma(double x) {
        if (x > 0 && x <= S_LIMIT) {
            // use method 5 from Bernardo AS103
            // accurate to O(x)
            return -GAMMA - 1 / x;
        }
        if (x >= C_LIMIT) {
            // use method 4 (accurate to O(1/x^8)
            double inv = 1 / (x * x);
            //            1       1        1         1
            // log(x) -  --- - ------ + ------- - -------
            //           2 x   12 x^2   120 x^4   252 x^6
            return Math.log(x) - 0.5 / x - inv * ((1.0 / 12) + inv * (1.0 / 120 - inv / 252));
        }
        return digamma(x + 1) - 1 / x;
    }

    /**
     * Computes the trigamma function of x.
     * This function is derived by taking the derivative of the implementation
     * of digamma.
     *
     * @param x Argument.
     * @return trigamma(x) to within 10-8 relative or absolute error whichever is smaller
     * @see <a href="http://en.wikipedia.org/wiki/Trigamma_function">Trigamma</a>
     * @see GammaDistrib#digamma(double)
     * @since 2.0
     */
    public static double trigamma(double x) {
        if (x > 0 && x <= S_LIMIT) {
            return 1 / (x * x);
        }
        if (x >= C_LIMIT) {
            double inv = 1 / (x * x);
            //  1    1      1       1       1
            //  - + ---- + ---- - ----- + -----
            //  x      2      3       5       7
            //      2 x    6 x    30 x    42 x
            return 1 / x + inv / 2 + inv / x * (1.0 / 6 - inv * (1.0 / 30 + inv / 42));
        }
        return trigamma(x + 1) + 1 / (x * x);
    }
    
    public static void main(String[] args) {
        double[] X = {0.000001, 11.2, 8.3, 13.1, 15.9, 11.5, 11.4, 12.3, 11.9, 5.5};
        double alpha = GammaDistrib.getAlpha(X);
        double beta = GammaDistrib.getBeta(X, alpha);
        //double beta = 1 / alpha; // force mean to be 1
        System.out.println("Setting Gamma distrib with alpha = " + alpha + " beta = " + beta);
        GammaDistrib gd = new GammaDistrib(alpha, 1/beta);
        double mean = 0.0;
        System.out.println("Sample");
        int N = 2000;
        for (int i = 0; i < N; i ++) {
            double y = gd.sample();
            mean += y;
            //System.out.println(i + "\t" + y);
        }
        System.out.println("Mean\t" + mean / N);
        
        
        
//        double[] x = {0, Double.MIN_VALUE, 1e-200, 1e-100, 1e-10, 0.00001, 0.1, 1.0, 100, 10000, 1e8, 1e10, 1e100, 1e200, Double.MAX_VALUE, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY};
//        for (double xx : x) {
//            System.out.println("Gamma(" + xx + ") = " + gamma(xx));
//            System.out.println("\tLog Gamma(" + xx + ") = " + lgamma(xx));
//            System.out.println("\tDi Gamma(" + xx + ") = " + digamma(xx));
//            System.out.println("\tTri Gamma(" + xx + ") = " + trigamma(xx));
//        }
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
