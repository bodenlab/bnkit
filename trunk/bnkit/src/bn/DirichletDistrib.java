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
    public DirichletDistrib(Enumerable domain, double same_alpha) {
        this.domain = domain;
        this.alpha = new double[domain.size()];
        Arrays.fill(this.alpha, same_alpha);
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
     * Construct a Dirichlet distribution based on a magnitude and a multinomial
     * @param m The magnitude of the Dirichlet: sum_i alpha_i
     * @param p A probability distribution: p_i = alpha_i / m
     */
    public DirichletDistrib(Enumerable domain, double[] p, double m) {
        if (domain.size() != p.length)
            throw new RuntimeException("Invalid distribution");
        this.domain = domain;
        this.alpha = new double[p.length];
        for (int i = 0; i < p.length; i ++)
            this.alpha[i] = p[i] * m;
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
    
    
    /**
     * Retrieve the number of documents that contain exactly a specified number of 
     * occurrences of a specified word
     * @param counts the word count matrix [document][word]
     * @param word the word identifier [0..m]
     * @param nOccur the number of occurrences of the word we are looking for
     * @return the number of matching documents
     */
    private static int getObservations(int[][] counts, int word, int nOccur) {
        int nDoc = 0;
        for (int doc = 0; doc < counts.length; doc ++) {
            if (counts[doc][word] == nOccur)
                nDoc ++;
        }
        return nDoc;
    }
    
    /**
     * Retrieve the document count histogram for a specified word (indexed by absolute word count)
     * @param counts the word count matrix [document][word]
     * @param word the word identifier [0..m]
     * @return the numbers of matching documents as an array indexed by each word count
     */
    private static int[] getHistogram(int[][] counts, int word) {
        int maxWordCount = 0;
        for (int doc = 0; doc < counts.length; doc ++) {
            if (counts[doc][word] > maxWordCount)
                maxWordCount = counts[doc][word];
        }
        int[] ret = new int[maxWordCount + 1];
        for (int wcnt = 0; wcnt < ret.length; wcnt ++) {
            int nDoc = 0;
            for (int doc = 0; doc < counts.length; doc ++) {
                if (counts[doc][word] == wcnt)
                    nDoc ++;
            }
            ret[wcnt] = nDoc;
        }
        return ret;
    }
    
    
    /**
     * Retrieve the number of documents that are a specified number of (counted) words long.
     * Example, getObservationLengths(counts, 20) is the number of documents that are exactly 20 tokens long.	 
     * @return the number of matching documents
     */
    private static int getObservationLengths(int[][] counts, int nWords) {
        int nDoc = 0;
        for (int[] count : counts) {
            int wSum = 0;
            for (int word = 0; word < count.length; word ++) {
                wSum += count[word];
            }
            if (wSum == nWords)
                nDoc ++;
        }
        return nDoc;
    }
    
    /**
     * Retrieve the maximum number of words in a document.
     * Example, getObservationLengths(counts, 20) is the number of documents that are exactly 20 tokens long.	 
     * @return the number of matching documents
     */
    private static int getMaxObservations(int[][] counts) {
        int wMax = 0;
        for (int[] count : counts) {
            int wSum = 0;
            for (int word = 0; word < count.length; word ++) {
                wSum += count[word];
            }
            if (wSum > wMax)
                wMax = wSum;
        }
        return wMax;
    }
    
    
    /*
    The following method getMLE is based on the method 
    cc.mallet.types.Dirichlet#learnParameters, found in MALLET.
    This software is copyright (C) 2002 Univ. of Massachusetts Amherst, Computer Science Dept.
    "MALLET" (MAchine Learning for LanguagE Toolkit).
    http://www.cs.umass.edu/~mccallum/mallet
    This software is provided under the terms of the Common Public License,
    version 1.0, as published by http://www.opensource.org.	For further
    information, see the file `LICENSE' included with this distribution. 
    */

    
    /** 
     * Learn Dirichlet parameters using frequency histograms.
     * Currently this implementation has issues, including 
     * inefficiencies,
     * convergence criterion, and 
     * finding negative priors.
     * @param counts An array of count histograms. <code>observations[10][3]</code> could be the number of documents that contain exactly 3 tokens of word type 10.
     * @returns The sum of the learned parameters.
     */ 
    public double setMLE(int[][] counts) {
        return setMLE(counts, 1.00001, 1.0, 200);
    }

    /** 
     * Learn Dirichlet parameters using frequency histograms
     * 
     * @param counts An array of count histograms. <code>observations[10][3]</code> could be the number of documents that contain exactly 3 tokens of word type 10.
     * @param shape Gamma prior E(X) = shape * scale, var(X) = shape * scale<sup>2</sup>
     * @param scale 
     * @param numIterations 200 to 1000 generally insures convergence, but 1-5 is often enough to step in the right direction
     * @returns The sum of the learned parameters.
     */ 
    public double setMLE(int[][] counts, double shape, double scale, int numIterations) {
        double parametersSum = 0;
        //	Initialize the parameter sum
        for (int k=0; k < alpha.length; k++) {
            parametersSum += alpha[k];
        }
        double oldParametersK;
        double currentDigamma;
        double denominator;

        int nWords = counts[0].length;
        int nonZeroLimit;
        int[] nonZeroLimits = new int[nWords];
        Arrays.fill(nonZeroLimits, -1);

        // The histogram arrays go up to the size of the largest document,
        //	but the non-zero values will almost always cluster in the low end.
        //	We avoid looping over empty arrays by saving the index of the largest
        //	non-zero value.
        int[] histogram;
        for (int w = 0; w < nWords; w ++) {
            histogram = getHistogram(counts, w);
            int sum = 0;
            for (int cnt = 0; cnt < histogram.length; cnt++) {
                sum += histogram[cnt];
                if (histogram[cnt] > 0) {
                    nonZeroLimits[w] = cnt;
                }
            }
        }
        for (int iteration=0; iteration<numIterations; iteration++) {
            // Calculate the denominator
            denominator = 0;
            currentDigamma = 0;
            // Iterate over the histogram:
            int wMax = getMaxObservations(counts);
            for (int i = 1; i < wMax; i ++) {
                currentDigamma += 1 / (parametersSum + i - 1);
                denominator += getObservationLengths(counts, i) * currentDigamma;
            }
            // Bayesian estimation Part I
            denominator -= 1/scale;
            // Calculate the individual parameters
            parametersSum = 0;
            for (int k=0; k<alpha.length; k++) {
                // What's the largest non-zero element in the histogram?
                nonZeroLimit = nonZeroLimits[k];
                oldParametersK = alpha[k];
                alpha[k] = 0;
                currentDigamma = 0;
                histogram = getHistogram(counts, k); // the number of a counts for word k (index by word count)
                for (int i=1; i <= nonZeroLimit; i++) {
                    currentDigamma += 1 / (oldParametersK + i - 1);
                    alpha[k] += histogram[i] * currentDigamma;
                }
                // Bayesian estimation part II
                alpha[k] = (oldParametersK * alpha[k] + shape) / denominator;
                parametersSum += alpha[k];
            }
        }
        return parametersSum;
    }



// Below is Python code for MLE of Dir
// Source: https://github.com/maxsklar/research/tree/master/2014_05_Dirichlet/python
// Paper: http://arxiv.org/pdf/1405.0099.pdf
//
//    #!/usr/bin/python
//    #
//    # A library for finding the optimal dirichlet prior from counts
//    # By: Max Sklar
//    # @maxsklar
//    # https://github.com/maxsklar
//
//    # Copyright 2013 Max Sklar
//
//    # Find the "sufficient statistic" for a group of multinomials.
//    # Essentially, it's the average of the log probabilities
    public static double[] getSufficientStatistic(EnumDistrib[] multinomials, int k) {
        int N = multinomials.length;
        int K = k;
        double[] retVal = new double[K];
        for (EnumDistrib m : multinomials) {
            for (int cnt = 0; cnt < K; cnt ++) 
                retVal[cnt] += Math.log(m.get(cnt));
        }
        for (int cnt = 0; cnt < K; cnt ++)
            retVal[cnt] /= N;
        return retVal;
    }
//
//    # Find the log probability of the data for a given dirichlet
//    # This is equal to the log probabiliy of the data.. up to a linear transform
    private static double logProbForMultinomials(double[] alphaTrials, double[] ss) {
        double alpha_sum = 0;
        double lgamma_sum = 0;
        double pq = 0;
        for (int i = 0; i < alphaTrials.length; i ++) {
            alpha_sum += alphaTrials[i];
            double lgamma_alpha = GammaDistrib.lgamma(alphaTrials[i]);
            lgamma_sum += lgamma_alpha;
            pq += (alphaTrials[i] * ss[i]);
        }
        double retVal = GammaDistrib.lgamma(alpha_sum) - lgamma_sum + pq;
        return retVal;
    }
//
//    #Gives the derivative with respect to the log of prior.  This will be used to adjust the loss
    private static double[] getGradientForMultinomials(double[] alphaTrials, double[] ss) {
        int K = alphaTrials.length;
        double alpha_sum = 0;
        for (int i = 0; i < alphaTrials.length; i ++) {
            alpha_sum += alphaTrials[i];
        }
        double C = GammaDistrib.digamma(alpha_sum);
        double[] retVal = new double[K];
        Arrays.fill(retVal, C);
        for (int i = 0; i < alphaTrials.length; i ++) {
            retVal[i] += (C + ss[i] - GammaDistrib.digamma(alphaTrials[i]));
        }
        return retVal;
    }
//
//    #The hessian is actually the sum of two matrices: a diagonal matrix and a constant-value matrix.
//    #We'll write two functions to get both
    private static double priorHessianConst(double[] alpha) { 
        double alpha_sum = 0;
        for (int i = 0; i < alpha.length; i ++) {
            alpha_sum += alpha[i];
        }
        return -GammaDistrib.trigamma(alpha_sum);
    }
    private static double[] priorHessianDiag(double[] alpha) { 
        double[] retVal = new double[alpha.length];
        for (int i = 0; i < alpha.length; i ++) 
            retVal[i] = GammaDistrib.trigamma(alpha[i]);
        return retVal;
    }
//
//    # Compute the next value to try here
//    # http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf (eq 18)
    private static double[] getPredictedStep(double[] alpha, double hConst, double[] hDiag, double[] gradient) {
        int K = alpha.length;
        double[] retVal = new double[K];
        double numSum = 0.0;
        for (int doc = 0; doc < K; doc ++) 
            numSum += gradient[doc] / hDiag[doc];
        double denSum = 0.0;
        for (int doc = 0; doc < K; doc ++)  
            denSum += 1.0 / hDiag[doc];
        double b = numSum / ((1.0/hConst) + denSum);
        for (int doc = 0; doc < K; doc ++)
            retVal[doc] = (b - gradient[doc]) / hDiag[doc];
        return retVal;
    }
//
//    # Uses the diagonal hessian on the log-alpha values	
    private static double[] getPredictedStepAlt(double[] alpha, double hConst, double[] hDiag, double[] gradient) {
        int K = alpha.length;
        double[] retVal = new double[K];
        double Z = 0;
        for (int cnt = 0; cnt < K; cnt ++)
            Z += alpha[cnt] / (gradient[cnt] - alpha[cnt]*hDiag[cnt]);
        Z *= hConst;
        double[] Ss = new double[K];
        double sumSs = 0;
        for (int cnt = 0; cnt < K; cnt ++) {
            Ss[cnt] = 1.0 / (gradient[cnt] - alpha[cnt]*hDiag[cnt]) / (1 + Z);
            sumSs += Ss[cnt];
        }
        for (int doc = 0; doc < K; doc ++)
            retVal[doc] = gradient[doc] / (gradient[doc] - alpha[doc]*hDiag[doc]) * (1 - hConst * alpha[doc] * sumSs);
        return retVal;
    }
//
//    #The currentPriors and data are global, so we don't need to pass them in
    private static double getTotalLoss(double[] trialPriors, double[] ss) {
        return -1.0*logProbForMultinomials(trialPriors, ss);
    }
//
    private static double[] predictStepUsingHessian(double[] alpha, double[] gradient) {
        double totalHConst = priorHessianConst(alpha);
        double[] totalHDiag = priorHessianDiag(alpha);		
        return getPredictedStep(alpha, totalHConst, totalHDiag, gradient);
    }
//
    private static double[] predictStepLogSpace(double[] alpha, double[] gradient, double[] ss) {
        double totalHConst = priorHessianConst(alpha);
        double[] totalHDiag = priorHessianDiag(alpha);
        return getPredictedStepAlt(alpha, totalHConst, totalHDiag, gradient);
    }
//
//
//    # Returns whether it's a good step, and the loss	
    private static double testTrialPriors(double[] trialPriors, double[] ss) {
        for (double alpha : trialPriors) {
            if (alpha <= 0)
                return Double.POSITIVE_INFINITY;
        }
        return getTotalLoss(trialPriors, ss);
    }
//
    private static double sqVectorSize(double[] v) {
        double s = 0;
        for (int doc = 0; doc < v.length; doc ++) 
            s += Math.pow(v[doc], 2);
        return s;
    }
//
    public static double[] findDirichletPriors(double[] alphaStart, double[] ss) {
        double[] currentPriors = new double[alphaStart.length];
        System.arraycopy(alphaStart, 0, currentPriors, 0, currentPriors.length);
        double[] trialPriors = new double[currentPriors.length];
        // Only step in a positive direction, get the current best loss.
        double currentLoss = getTotalLoss(currentPriors, ss);
        double gradientToleranceSq = Math.pow(2, -20);
        double learnRateTolerance = Math.pow(2, -10);
        int count = 0;
        while (count < 1000) {
            count += 1;
            // Get the data for taking steps
            double[] gradient = getGradientForMultinomials(currentPriors, ss);
            double gradientSize = sqVectorSize(gradient);
            if (gradientSize < gradientToleranceSq) {
                // Converged with small gradient
                return currentPriors;
            }
            // First, try the second order method...
            double[] trialStep = predictStepUsingHessian(currentPriors, gradient);
            for (int a = 0; a < currentPriors.length; a ++)
                trialPriors[a] = currentPriors[a] + trialStep[a];
            double loss = testTrialPriors(trialPriors, ss);
            if (loss < currentLoss) {
                currentLoss = loss;
                // replace with new-found priors
                System.arraycopy(trialPriors, 0, currentPriors, 0, currentPriors.length);
                continue;
            }
            // Next, try first-order, iterative option...
            trialStep = predictStepLogSpace(currentPriors, gradient, ss);
            for (int a = 0; a < currentPriors.length; a ++)
                trialPriors[a] = currentPriors[a] * Math.exp(trialStep[a]);
            // un-used result ?
            // loss = testTrialPriors(trialPriors, ss);
            // Step in the direction of the gradient until there is a loss improvement
            loss = Double.POSITIVE_INFINITY;
            double learnRate = 1.0;
            while (loss > currentLoss) {
                for (int a = 0; a < currentPriors.length; a ++)
                    trialPriors[a] = currentPriors[a] + gradient[a] * learnRate;
                loss = testTrialPriors(trialPriors, ss);
                learnRate *= 0.9;
            }
            if (learnRate / 0.9 < learnRateTolerance)
                // Converged with small learn rate
                return currentPriors;

            currentLoss = loss;
            for (int a = 0; a < currentPriors.length; a ++)
                currentPriors[a] = trialPriors[a];
        }
        // Reached max iterations
        return currentPriors;
    }

//    def findDirichletPriorsFromMultinomials(multinomials, initAlphas):
//            ss = getSufficientStatistic(multinomials)
//            return findDirichletPriors(ss, initAlphas)

    
    public static void main(String[] args) {
        java.util.Random rand = new java.util.Random();
        int N = 30;
        int K = 3;
        Enumerable dom = new Enumerable(K);
        DirichletDistrib d0 = new DirichletDistrib(dom, 12, 2, 15);
        System.out.println(d0);
        int[][] counts = new int[N][K];
        EnumDistrib[] samples = new EnumDistrib[N];
        for (int i = 0; i < N; i ++) {
            EnumDistrib d = (EnumDistrib) d0.sample();
            samples[i] = d;
            double p = d0.get(d);
            System.out.println(i + ": " + d + " p = " + String.format("%4.2f;", p));
        }
        double[] ss = DirichletDistrib.getSufficientStatistic(samples, K);
        double[] alpha = DirichletDistrib.findDirichletPriors(new double[] {2,3,4}, ss);
        d0.alpha = alpha;
        System.out.println(d0);
        for (int i = 0; i < 30; i ++) {
            EnumDistrib d = (EnumDistrib) d0.sample();
            double p = d0.get(d);
            System.out.println(i + ": " + d + " p = " + String.format("%4.2f;", p));
        }
    }

}

// Below is R-code for MLE estimation of Dirichlet distribution params from data
// Source: http://web.hku.hk/~kaing/Dirichlet_app01.pdf
//
//    function(X, ind, N)
//    {
//        # Aa <- MLE.Dirichlet.EM.gradient(X, ind, N)
//        # --------------------------------------------------------
//        # Aim:    Finding the MLEs of the Dirichlet parameter
//        #         vector a = (a_1, ..., a_n)
//        # Method: Using the EM gradient algorithm (2.75)
//        # Input:  X:     an m x n observed data matrix
//        #         ind=1: using (2.71) as the initial values
//        #         ind=2: using (2.73) as the initial values
//        #         N:     the number of iterations required for
//    #
//    # Output: A ---- N by n matrix whose t-th row is aˆ(t)
//    #         a ---- the MLE of a = (a_1, ..., a_n)
//    # --------------------------------------------------------
//    m <- dim(X)[1]
//    n <- dim(X)[2]
//    one <- rep(1, n)
//    xmean <- apply(X, 2, mean)
//    G <- (apply(X, 2, prod))ˆ(1/m)
//    if(ind == 1) {
//          a <- min(X) * one
//    }
//    if(ind == 2) {
//          ga <-  - digamma(1)
//          b <- xmean
//          de0 <- sum(b * log(G))
//          aplus <- ((n - 1) * ga)/(sum(b * log(b)) - de0)
//          a <- b * aplus
//    }
//    A <- matrix(0, N, n)
//    for(tt in 1:N) {
//          aplus <- sum(a)
//          g <- digamma(aplus) - digamma(a) + log(G)
//          b <- trigamma(a)
//    the EM gradient algorithm
//    }
//        return(A, a)
//    }
//    a <- a + g/b A[tt, ]<-a

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
