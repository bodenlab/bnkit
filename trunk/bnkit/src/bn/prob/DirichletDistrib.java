package bn.prob;

import bn.Distrib;
import dat.Domain;
import bn.EnumDistrib;
import dat.Enumerable;
import java.io.Serializable;
import java.util.Arrays;

/**
 * Implementation of Dirichlet distribution.
 * The class implements functionality for constructing Dirichlet for 
 * generating enumerable distributions of a specified domain.
 * The class implements sampling and ML estimation from enumerable distributions.
 * Some functionality is adapted from other sources (see copyright messages below).
 */
public class DirichletDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;

    /** Dirichlet prior */
    private double[] alpha;
    private GammaDistrib[] gammas;
    /** Domain for enumerable distributions that can be generated from the Dirichlet and/or
     * learned from.     */
    private final Enumerable domain;
    
    /**
     * Construct a Dirichlet distribution with the given domain and
     * "alpha" parameter value for all dimensions.
     * @param domain the domain over which this distribution can generate distributions
     * @param same_alpha the alpha value that applies to all components
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

    /**
     * Construct a Dirichlet distribution with the given domain and
     * "alpha" parameter values.
     * @param domain the domain over which this distribution can generate distributions
     * @param alpha the alpha values that apply
     */
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
     * Set new prior parameters for distribution.
     * @param alpha new prior
     */
    public void setPrior(double[] alpha) {
        for (int i = 0; i < this.alpha.length; i ++) {
            this.alpha[i] = alpha[i];
            gammas[i] = new GammaDistrib(this.alpha[i], 1);
        }
    }
    
    /**
     * Find and set the optimal alpha parameters for a Dirichlet based on given instances of enumerable distribution.
     * @param dists array of enumerable distributions compatible with this Dirichlet
     */
    public void setPrior(EnumDistrib[] dists) {
        if (dists.length > 0) {
            double[] ss = DirichletDistrib.getSufficientStatistic(dists, alpha.length);
            setPrior(DirichletDistrib.findPrior(alpha, ss));
            return;
        }
        throw new RuntimeException("Invalid data for estimation of Dirichlet");
    }
    
    public Domain getDomain() {
        return domain;
    }
    
    /**
     * Returns the probability of a vector of values from this distribution.
     * @param dist an enumerable distribution
     * @return the probability
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
     * @return an enumerable distribution, chosen in proportion to its probability
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
    
    /**
     * String representation of the Dirichlet distribution, containing alpha values.
     * @return string representation of object
     */
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
    The following method findDirichletFromCounts is based on the method 
    cc.mallet.types.Dirichlet#learnParameters, found in MALLET.
    This software is copyright (C) 2002 Univ. of Massachusetts Amherst, Computer Science Dept.
    "MALLET" (MAchine Learning for LanguagE Toolkit).
    http://www.cs.umass.edu/~mccallum/mallet
    This software is provided under the terms of the Common Public License,
    version 1.0, as published by http://www.opensource.org.	For further
    information, see the file `LICENSE' available with the MALLET software distribution. 
    */
    
    /** 
     * Learn Dirichlet parameters using frequency histograms.
     * Currently this implementation has issues, including 
     * inefficiencies, convergence criterion, and finding negative priors.
     * @param counts An array of count histograms.
     * @param alphaStart initial alpha values
     * @return the new alpha values
     */ 
    public double[] findPriorFromCounts(int[][] counts, double[] alphaStart) {
        return findPriorFromCounts(counts, alphaStart, 1.00001, 1.0, 200);
    }

    /** 
     * Learn Dirichlet parameters using frequency histograms
     * 
     * @param counts An array of count histograms. <code>observations[10][3]</code> could be the number of documents that contain exactly 3 tokens of word type 10.
     * @param alphaStart initial alpha values
     * @param shape Gamma prior E(X) = shape * scale, var(X) = shape * scale<sup>2</sup>
     * @param scale 
     * @param numIterations 200 to 1000 generally insures convergence, but 1-5 is often enough to step in the right direction
     * @return the new alpha values
     */ 
    public static double[] findPriorFromCounts(int[][] counts, double[] alphaStart, double shape, double scale, int numIterations) {
        double parametersSum = 0;
        //	Initialize the parameters
        double[] alpha = new double[alphaStart.length];
        for (int k=0; k < alphaStart.length; k++) {
            parametersSum += alphaStart[k];
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
        return alpha;
    }

    /*
    The code below is based on Python code for MLE of Dirichlet copyright Max Sklar.
    // Source: https://github.com/maxsklar/research/tree/master/2014_05_Dirichlet/python
    // Paper: http://arxiv.org/pdf/1405.0099.pdf
    # Copyright 2013 Max Sklar
    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
    */
    /**
     * Find the "sufficient statistic" for a group of multinomials.
     * Essentially, it's the average of the log probabilities.
     * @param multinomials instances of enumerable distributions
     * @param k number of values in distribution
     */ 
    public static double[] getSufficientStatistic(EnumDistrib[] multinomials, int k) {
        int N = multinomials.length;
        int K = k;
        double[] retVal = new double[K];
        for (EnumDistrib m : multinomials) {
            for (int cnt = 0; cnt < K; cnt ++) {
                double prob = m.get(cnt);
                double logProb = Math.log(prob);
                retVal[cnt] += logProb;
            }
        }
        for (int cnt = 0; cnt < K; cnt ++)
            retVal[cnt] /= N;
        return retVal;
    }

    /**
     * Find the log probability of the data for a given Dirichlet.
     * This is equal to the log probability of the data. up to a linear transform.
     * @param alphaTrials initial alpha
     * @param ss sufficient statistics for the data
     * @return the log probability
     */
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

    /**
     * Gives the derivative with respect to the log of prior.  
     * This will be used to adjust the loss.
     * @param alphaTrials initial alpha
     * @param ss sufficient statistics for the data
     * @return the gradient
     */
    private static double[] getGradientForMultinomials(double[] alphaTrials, double[] ss) {
        int K = alphaTrials.length;
        double alpha_sum = 0;
        for (int i = 0; i < alphaTrials.length; i ++) {
            alpha_sum += alphaTrials[i];
        }
        double C = GammaDistrib.digamma(alpha_sum);
        double[] retVal = new double[K];
        for (int i = 0; i < alphaTrials.length; i ++) {
            retVal[i] += (C + ss[i] - GammaDistrib.digamma(alphaTrials[i]));
        }
        return retVal;
    }

    /**
     * The Hessian is the sum of two matrices: 
     * a diagonal matrix and a constant-value matrix.
     * Two functions define them. This is the constant.
     * @param alpha alpha parameters
     * @return Hessian constant
     */
    private static double priorHessianConst(double[] alpha) { 
        double alpha_sum = 0;
        for (int i = 0; i < alpha.length; i ++) {
            alpha_sum += alpha[i];
        }
        return -GammaDistrib.trigamma(alpha_sum);
    }

    /**
     * The Hessian is the sum of two matrices: 
     * a diagonal matrix and a constant-value matrix.
     * Two functions define them. This is the diagonal.
     * @param alpha alpha parameters
     * @return Hessian diagonal
     */
    private static double[] priorHessianDiag(double[] alpha) { 
        double[] retVal = new double[alpha.length];
        for (int i = 0; i < alpha.length; i ++) 
            retVal[i] = GammaDistrib.trigamma(alpha[i]);
        return retVal;
    }
    
    /**
     * Compute the next value to try here.
     * http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf (eq 18)
     * @param alpha alpha parameters
     * @param hConst constant of Hessian
     * @param hDiag diagonal of Hessian
     * @param gradient the gradient 
     * @return step to change alpha
     */
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

    /**
     * Use the Hessian on the log-alpha value to compute the next step.
     * @param alpha alpha parameters
     * @param hConst constant of Hessian
     * @param hDiag diagonal of Hessian
     * @param gradient the gradient 
     * @return step to change alpha
     */
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

    /**
     * Get the loss for current parameters (alpha-values) with the data represented by "sufficient statistics".
     * @param alpha alpha (parameters defining prior)
     * @param ss sufficient statistics for the data
     * @return loss
     */
    private static double getTotalLoss(double[] alpha, double[] ss) {
        return -1.0*logProbForMultinomials(alpha, ss);
    }

    /**
     * Use Hessian to predict step from current parameters (alpha-values).
     * @param alpha initial alpha
     * @param gradient gradient
     * @return step (change to alpha)
     */
    private static double[] predictStepUsingHessian(double[] alpha, double[] gradient) {
        double totalHConst = priorHessianConst(alpha);
        double[] totalHDiag = priorHessianDiag(alpha);		
        return getPredictedStep(alpha, totalHConst, totalHDiag, gradient);
    }
 
    /**
     * Use Hessian to predict step from current parameters (alpha-values).
     * @param alpha initial alpha
     * @param gradient gradient
     * @return step (change to alpha)
     */
    private static double[] predictStepLogSpace(double[] alpha, double[] gradient) {
        double totalHConst = priorHessianConst(alpha);
        double[] totalHDiag = priorHessianDiag(alpha);
        return getPredictedStepAlt(alpha, totalHConst, totalHDiag, gradient);
    }

    /**
     * Returns whether it's a good alpha, and the loss.
     * @param alpha alpha (parameters defining prior)
     * @param ss sufficient statistics for the data
     * @return loss (+infinite if invalid (negative or zero alpha)
     */
    private static double testTrialPriors(double[] alpha, double[] ss) {
        for (double a : alpha) {
            if (a <= 0)
                return Double.POSITIVE_INFINITY;
        }
        return getTotalLoss(alpha, ss);
    }

    private static double sqVectorSize(double[] v) {
        double s = 0;
        for (int doc = 0; doc < v.length; doc ++) 
            s += Math.pow(v[doc], 2);
        return s;
    }

    /**
     * Find the optimal alpha parameters for a Dirichlet based on given instances of enumerable distribution.
     * @param dists array of enumerable distributions compatible with this Dirichlet
     * @return the alpha values
     */
    public double[] findPrior(EnumDistrib[] dists) {
        if (dists.length > 0) {
            double[] ss = DirichletDistrib.getSufficientStatistic(dists, alpha.length);
            return DirichletDistrib.findPrior(alpha, ss);
        }
        throw new RuntimeException("Invalid data for estimation of Dirichlet");
    }
    
    /**
     * Find the optimal alpha parameters for a Dirichlet based on given instances of enumerable distribution.
     * @param dists array of enumerable distributions compatible with this Dirichlet
     * @param alphaStart initial alpha values
     * @return the optimal alpha values
     */
    public static double[] findPrior(EnumDistrib[] dists, double[] alphaStart) {
        double[] ss = DirichletDistrib.getSufficientStatistic(dists, alphaStart.length);
        return DirichletDistrib.findPrior(alphaStart, ss);
    }
    
    /**
     * Find the optimal alpha parameters for a Dirichlet based on sufficient statistics extracted for
     * data set.
     * @param alphaStart initial alpha values
     * @param ss sufficient statistics for the data
     * @return the optimal alpha values
     */
    public static double[] findPrior(double[] alphaStart, double[] ss) {
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
                // System.out.println("Exited after " + count + " rounds with loss = " + currentLoss + " gradient = " + gradient);
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
            trialStep = predictStepLogSpace(currentPriors, gradient);
            for (int a = 0; a < currentPriors.length; a ++)
                trialPriors[a] = currentPriors[a] * Math.exp(trialStep[a]);
            loss = testTrialPriors(trialPriors, ss);
            if (loss < currentLoss) {
                currentLoss = loss;
                // replace with new-found priors
                System.arraycopy(trialPriors, 0, currentPriors, 0, currentPriors.length);
                continue;
            }
            
            // Step in the direction of the gradient until there is a loss improvement
            loss = Double.POSITIVE_INFINITY;
            double learnRate = 1.0;
            while (loss > currentLoss) {
                for (int a = 0; a < currentPriors.length; a ++)
                    trialPriors[a] = currentPriors[a] + gradient[a] * learnRate;
                loss = testTrialPriors(trialPriors, ss);
                learnRate *= 0.9;
            }
            if (learnRate < learnRateTolerance) {
                // Converged with small learn rate
                // System.out.println("Exited after " + count + " rounds with loss = " + loss + " learning rate = " + learnRate);
                return trialPriors;
            }
            currentLoss = loss;
            for (int a = 0; a < currentPriors.length; a ++)
                currentPriors[a] = trialPriors[a];
        }
        // Reached max iterations
        // System.out.println("Exited after " + count + " rounds with loss = " + currentLoss);
        return currentPriors;
    }

    
    public static void main(String[] args) {
        java.util.Random rand = new java.util.Random(1);
        Enumerable dom = new Enumerable(2);
        EnumDistrib[] samples = {
            new EnumDistrib(dom, 0.3, 0.7),
            new EnumDistrib(dom, 0.2, 0.8),
            new EnumDistrib(dom, 0.1, 0.9),
            new EnumDistrib(dom, 0.2, 0.8),
            new EnumDistrib(dom, 0.1, 0.9),
            new EnumDistrib(dom, 0.4, 0.6),
        };
        double[] ss = DirichletDistrib.getSufficientStatistic(samples, 2);
        double[] alpha = DirichletDistrib.findPrior(new double[] {.5, 0.5}, ss);
        DirichletDistrib d = new DirichletDistrib(dom, alpha);
        System.out.println(d);

//        int N = 10;
//        int K = 3;
//        Enumerable dom = new Enumerable(K);
//        DirichletDistrib d0 = new DirichletDistrib(dom, 12, 2, 15);
//        System.out.println(d0);
//        int[][] counts = new int[N][K];
//        EnumDistrib[] samples = new EnumDistrib[N];
//        for (int i = 0; i < N; i ++) {
//            EnumDistrib d = (EnumDistrib) d0.sample();
//            samples[i] = d;
//            double p = d0.get(d);
//            System.out.println(i + ": " + d + " p = " + String.format("%4.2f;", p));
//        }
//        double[] ss = DirichletDistrib.getSufficientStatistic(samples, K);
//        double[] alpha = DirichletDistrib.findPrior(new double[] {0.33,0.33,0.33}, ss);
//        d0.setPrior(alpha);
//        System.out.println("Dirichlet Alpha: " + d0);
//        for (int i = 0; i < 30; i ++) {
//            EnumDistrib d = (EnumDistrib) d0.sample();
//            double p = d0.get(d);
//            System.out.println(i + ": " + d + " p = " + String.format("%4.2f;", p));
//        }
        System.out.println("Final loss = " +DirichletDistrib.getTotalLoss(alpha, ss));
        System.out.println("Best loss = " + DirichletDistrib.getTotalLoss(new double[] {1,2}, ss));
        
//        DirichletDistrib d1 = new DirichletDistrib(Enumerable.nacid, 1.0);
//        EnumDistrib[] dd = new EnumDistrib[dna.length];
//        for (int j = 0; j < dna.length; j ++) {
//            dd[j] = new EnumDistrib(Enumerable.nacid);
//            for (Object sym : Enumerable.nacid.getValues()) {
//                int cnt = 0;
//                for (int i = 0; i < dna[j].length(); i ++)
//                    cnt += (dna[j].charAt(i) == (Character)sym ? 1 : 0);
//                dd[j].set(sym, (double) cnt / dna[j].length());
//            }
//            System.out.println(dd[j]);
//        }
//        d1.setPrior(dd);
//        System.out.println(d1);
//        for (int j = 0; j < 10; j ++)
//            System.out.println(d1.sample());
    }

    static String[] dna = {
            "TTCGGCACGAGTCTCGGGCGGGAGAAGAAGAAGAATTAGTAAAGTGTGATCATAATGTCTGCTAGCGGCG",
            "GCACCGGAGATGAAGATAAGAAGCCTAATGATCAGATGGTTCATATCAATCTCAAGGTTAAGGGTCAGGA",
            "TGGGAATGAAGTTTTTTTCAGGATCAAACGTAGCACACAGATGCGCAAGCTCATGAATGCTTATTGTGAC",
            "CGGCAGTCAGTGGACATGAACTCAATTGCATTCTTATTTGATGGGCGCAGGCTTAGGGCAGAGCAAACTC",
            "CTGATGAGCTGGAGATGGAGGAGGGTGATGAAATCGATGCAATGCTACATCAAACTGGAGGCAGTTGCTG",
            "CACTTGTTTCTCTAATTTTTAACTTGGTTTATGTTAGTAGATTGTTTAGGGTAATACTTTCAACTCCCTC",
            "ATCTGCTCTAAGATGGGTAAATTTATGAATGTTTAGTTTTCAGTATTAGATGATGACACTACTAAATGGT",
            "TCAATTTTCATGGCATTTGTAAAAGTTTACTCTTAATATGGTTAAAAA"};
    
}
