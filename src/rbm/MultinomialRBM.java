/*
    bnkit -- software for building and using probabilistic models
    including Bayesian networks.
    Copyright (C) 2014-2016  M. Boden et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rbm;

import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A restricted Boltzmann Machine with multinomial visible and Boolean hidden nodes.
 * @author mikael
 */
public class MultinomialRBM extends AbstractRBM {

    private  double[][] w;     // weights
    private  double[]   a;     // bias on visible
    private  double[]   b;     // bias on hidden

    private Enumerable domain; // domain of visible variables (set of possible values)

    /**
     * Set the SD for random weights (around mean 0)
     */
    public static double WEIGHT_STDEV = 0.01;

    public int getK() {
        return domain.size();
    }

    public int getIndex(Object value) {
        return domain.getIndex(value);
    }

    /**
     * Create an RBM with a specified numbers of visible and hidden variables.
     * Multinomial visible and Boolean hidden variables will be generated.
     *
     * @param nVis number of visible variables
     * @param nHid number of hidden variables
     * @param domain the values that all visible variables can take
     */
    public MultinomialRBM(int nVis, int nHid, Enumerable domain) {
        this.domain = domain;
        int K = getK(); // degree of each visible variable
        rand = new Random(RANDOM_SEED);
        this.v = new EnumVariable[nVis];
        this.Pv = new EnumDistrib[nVis];
        this.h = new EnumVariable[nHid];
        this.Ph = new EnumDistrib[nHid];
        this.a = new double[v.length * K];
        this.b = new double[h.length];
        this.w = new double[h.length][v.length * K];
        this.lnk = new boolean[h.length][v.length];
        for (int j = 0; j < h.length; j ++) {
            this.h[j] = new EnumVariable(Enumerable.bool, "hid_" + j);
            this.h[j].setPredef("Boolean");
            this.Ph[j] = new EnumDistrib(Enumerable.bool);
            this.Ph[j].setSeed(rand.nextLong());
            for (int i = 0; i < v.length; i ++) {
                for (int k = 0; k < K; k ++) {
                    this.w[j][i * K + k] = rand.nextGaussian() * WEIGHT_STDEV;
                    this.lnk[j][i * K + k] = true; // linked, by default
                }
                if (j == 0) {
                    for (int k = 0; k < K; k ++)
                        this.a[i * K + k] = rand.nextGaussian() * WEIGHT_STDEV;
                    this.v[i] = new EnumVariable(domain, "vis_" + i);
                    this.Pv[i] = new EnumDistrib(domain);
                    this.Pv[i].setSeed(rand.nextLong());
                }
                if (i == 0)
                    this.b[j] = rand.nextGaussian() * WEIGHT_STDEV;
            }
        }
    }


    public double softmax_denom(double[] x) {
        double denom = 0;
        for (int i = 0; i < x.length; i ++)
            denom += Math.exp(x[i]);
        return denom;
    }

    public double[] softmax_distrib(double[] x) {
        double[] d = new double[x.length];
        double denom = softmax_denom(x);
        for (int i = 0; i < d.length; i ++)
            d[i] = Math.exp(x[i]) / denom;
        return d;
    }

    public void setVisible(Object[] input) {
        for (int i = 0; i < v.length; i ++) {
            if (input[i] != null)
                assignVisible(i, input[i]);
        }
    }

    @Override
    public Object[] encode(Object[] input, Object[] hinst) {
        for (int j = 0; j < h.length; j ++) {
            double net = 0.0;
            for (int i = 0; i < v.length; i ++) {
                if (input[i] != null && lnk[j][i]) {
                    int index = getIndex(input[i]);
                    net += w[j][i * getK() + index];
                }
            }
            double p = logistic(net + b[j]);
            Ph[j].set(new double[] {p, 1.0 - p}) ;
            hinst[j] = Ph[j].sample();
        }
        return hinst;
    }

    @Override
    public Object[] decode(Object[] hidden, Object[] vinst) {
        int K = getK(); // degree of each visible variable
        for (int i = 0; i < v.length; i ++) {
            double[] net = new double[getK()];
            for (int k = 0; k < K; k ++) {
                int ik = i * K + k;
                net[ik] = a[ik];
                for (int j = 0; j < h.length; j ++) {
                    if (hidden[j] != null && lnk[j][ik]) {
                        if ((Boolean)hidden[j])
                            net[k] +=  w[j][ik];
                    }
                }
            }
            Pv[i].set(softmax_distrib(net)) ;
            vinst[i] = Pv[i].sample();
        }
        return vinst;
    }

    @Override
    public Double[][] getCDGradient(Object[][] minibatch, int niter) {
        int K = getK(); // degree of each visible variable
        Double[][] sum = new Double[getNHidden() + 1][getNVisible() * K + K]; // plus one to include biases
        int[][] cnt = new int[getNHidden() + 1][getNVisible() + 1];
        err = 0;
        // TODO: multi-thread the following iteration
        for (int p = 0; p < minibatch.length; p ++) {
            Object[] input0 = minibatch[p];
            Object[] hidd0 = encode(input0);
            // determine the positive signal for the gradient
            double[][] pos = new double[getNHidden() + 1][getNVisible() * K + K];
            for (int i = 0; i < getNVisible(); i ++) {
                int index = getIndex(input0[i]);
                if (minibatch[p][i] != null) {
                    for (int j = 0; j < getNHidden(); j ++) {
                        if (i == 0) { // biases; do this only for one round
                            cnt[j][getNVisible()] ++;
                            for (int k = 0; k < K; k ++) {
                                sum[j][getNVisible() * K + k] = 0.0;
                                pos[j][getNVisible() * K + k] = Ph[j].get(0);
                            }
                        }
                        cnt[j][i] ++;
                        for (int k = 0; k < K; k ++) {
                            int ik = i * K + k;
                            sum[j][ik] = 0.0;
                            pos[j][ik] = ((k == index) ? 1 /*true*/ * Ph[j].get(0) : 0);
                        }
                    }
                    cnt[getNHidden()][i] ++;
                    for (int k = 0; k < K; k ++) {
                        int ik = i * K + k;
                        sum[getNHidden()][ik] = 0.0;
                        pos[getNHidden()][ik] = (k == index) ? 1 : 0;
                    }
                }
            }
            // do the back-and-forth, Gibbs sampling, in preparation for the negative signal
            for (int n = 0; n < niter; n ++) {
                Object[] input1 = decode_restricted(hidd0, input0);
                for (int i = 0; i < input0.length; i ++) {
                    int index = getIndex(input0[i]);
                    for (int k = 0; k < K; k ++) {
                        int ik = i * K + k;
                        double x0 = (k == index) ? 1 : 0;
                        double p1 = this.Pv[i].get(k);
                        err += Math.sqrt((x0 - p1) * (x0 - p1));
                    }
                }
                Object[] hidd1 = encode(input1);
            }
            // determine the gradient by computing the negative signal
            for (int i = 0; i < getNVisible(); i ++) {
                if (minibatch[p][i] != null) {
                    int index = getIndex(input0[i]);
                    for (int k = 0; k < K; k ++) {
                        int ik = i * K + k;
                        for (int j = 0; j < getNHidden(); j ++) {
                            if (i == 0) { // do this one round only
                                double neg = Ph[j].get(0);
                                sum[j][getNVisible() * K + k] += pos[j][getNVisible() * K + k] - neg;
                            }
                            double neg = Pv[i].get(k) * Ph[j].get(0);
                            sum[j][ik] += pos[j][ik] - neg;
                        }
                        double neg = Pv[i].get(k);
                        sum[getNHidden()][ik] += pos[getNHidden()][ik] - neg;
                    }
                }
            }
        }
        // find average, given counts
        for (int i = 0; i < getNVisible() + 1; i ++) {
            for (int j = 0; j < getNHidden() + 1; j ++) {
                for (int k = 0; k < K; k++) {
                    if (cnt[j][i] > 0)
                        sum[j][i * K + k] /= cnt[j][i];
                    else
                        sum[j][i * K + k] = null;
                }
            }
        }
        return sum;
    }

    @Override
    public void setCDGradient(Double[][] delta) {
        int K = getK(); // degree of each visible variable
        for (int j = 0; j < getNHidden() + 1; j ++) {
            for (int i = 0; i < getNVisible() + 1; i ++) {
                if (delta[j][i] != null) {
                    if (j == getNHidden())
                        for (int k = 0; k < getK(); k++)
                            a[i * K + k] += delta[j][i * getK() + k];
                    else if (i == getNVisible())
                        b[j] += delta[j][i * K];
                    else
                        for (int k = 0; k < K; k++)
                            w[j][i * K + k] += delta[j][i * K + k];
                }
            }
        }
    }


}
