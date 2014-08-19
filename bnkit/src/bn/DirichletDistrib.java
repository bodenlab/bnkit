package bn;

import java.io.Serializable;
import java.util.Arrays;

/*
 * Implementation is based on Bayesian Logic (BLOG) inference engine version 0.7 
 * (see copyright message below).
 */

public class DirichletDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;
	
    private double[] alpha;
    private GammaDistrib[] gammas;
    private final Enumerable domain;
    
    /**
     * Constructs a Dirichlet distribution with the given domain and
     * parameter value for all dimensions.
     * 
     */
    public DirichletDistrib(Enumerable domain, double paramVal) {
        this.domain = domain;
        this.alpha = new double[domain.size()];
        Arrays.fill(this.alpha, paramVal);
        gammas = new GammaDistrib[this.alpha.length];
        for (int i = 0; i < gammas.length; i++) {
            gammas[i] = new GammaDistrib(alpha[i], 1);
        }
    }

    public DirichletDistrib(Enumerable domain, double... alpha) {
        if (domain.size() != alpha.length)
            throw new RuntimeException("Invalid distribution");
        this.domain = domain;
        this.alpha = alpha;
        gammas = new GammaDistrib[this.alpha.length];
        for (int i = 0; i < gammas.length; i++) {
            gammas[i] = new GammaDistrib(alpha[i], 1);
        }
    }
    
    /**
     * Returns the probability of a vector of values from this distribution.
     */
    public double get(Object dist) {
        EnumDistrib d = (EnumDistrib) dist;
        double[] pdist = d.get();
        double prob = 1.0;
        for (int i = 0; i < pdist.length; i++) {
            double count_term = alpha[i];
            double x = pdist[i];
            prob *= Math.pow(x, count_term - 1);
        }
        prob /= normalize(alpha);
        return prob;
    }
	
    /**
     * Returns an enumerable sampled from this distribution.
     */
    @Override
    public Object sample() {
        double sum = 0.0;
        double[] sample = new double[alpha.length];
        for (int i = 0; i < alpha.length; i++) {
            GammaDistrib component = gammas[i];
            double y = component.sample();
            sum += y;
            sample[i] = y;
        }

        for (int i = 0; i < alpha.length; i++) {
            sample[i] /= sum;
        }
        return new EnumDistrib(domain, sample);
    }

    /**
     * Computes the normalization constant for a Dirichlet distribution with
     * the given parameters.
     * 
     * @param a list of parameters of a Dirichlet distribution
     * @return the normalization constant for such a distribution
     */
    public static final double normalize(double[] params) {
        double denom = 0.0;
        double numer = 1.0;
        for (double param: params) {
            numer *= GammaDistrib.gamma(param);
            denom += param;
        }
        denom = GammaDistrib.gamma(denom);
        return numer / denom;
    }
    
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < alpha.length; i ++) {
            if (i < alpha.length - 1)
                sb.append(String.format("%4.2f,", alpha[i]));
            else
                sb.append(String.format("%4.2f;", alpha[i]));
        }
        return sb.toString();
    }
    
    public static void main(String[] args) {
        DirichletDistrib d = new DirichletDistrib(new Enumerable(3), 1, 2, 5);
        System.out.println(d);
        for (int i = 0; i < 10; i ++) {
            EnumDistrib d0 = (EnumDistrib) d.sample();
            double p = d.get(d0);
            System.out.println(i + ": " + d0 + " p = " + String.format("%4.2f;", p));
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
