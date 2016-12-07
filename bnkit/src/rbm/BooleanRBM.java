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

import bn.Predef;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Enumerable;

import java.io.*;
import java.util.Random;

/**
 * A restricted Boltzmann Machine with Boolean visible and hidden nodes.
 * @author mikael
 */
public class BooleanRBM extends AbstractRBM {

    private final double[][] w;     // weights
    private final double[]   a;     // bias on visible
    private final double[]   b;     // bias on hidden
//    private final double[][] w;     // weights
//    private final double[]   a;     // bias on visible
//    private final double[]   b;     // bias on hidden

    /**
     * Set the SD for random weights (around mean 0)
     */
    public static double WEIGHT_STDEV = 0.01;
    
    /**
     * Create an RBM with a specified numbers of visible and hidden Boolean variables.
     * Visible and hidden variables will be generated.
     * @param nVis number of Boolean, visible variables
     * @param nHid number of Boolean, hidden variables
     */
    public BooleanRBM(int nVis, int nHid) {
        rand = new Random(RANDOM_SEED);
        this.v = new EnumVariable[nVis];
        this.Pv = new EnumDistrib[nVis];
        this.h = new EnumVariable[nHid];
        this.Ph = new EnumDistrib[nHid];
        this.a = new double[v.length];
        this.b = new double[h.length];
        this.w = new double[h.length][v.length];
        this.lnk = new boolean[h.length][v.length];
        for (int j = 0; j < h.length; j ++) {
            this.h[j] = new EnumVariable(Enumerable.bool, "hid_" + j);
            this.h[j].setPredef("Boolean");
            this.Ph[j] = new EnumDistrib(Enumerable.bool);
            this.Ph[j].setSeed(rand.nextLong());
            for (int i = 0; i < v.length; i ++) {
                this.w[j][i] = rand.nextGaussian() * WEIGHT_STDEV;
                this.lnk[j][i] = true; // linked, by default
                if (j == 0) {
                    this.a[i] = rand.nextGaussian() * WEIGHT_STDEV;
                    this.v[i] = new EnumVariable(Enumerable.bool, "vis_" + i);
                    this.v[i].setPredef("Boolean");
                    this.Pv[i] = new EnumDistrib(Enumerable.bool);
                    this.Pv[i].setSeed(rand.nextLong());
                }
                if (i == 0) 
                    this.b[j] = rand.nextGaussian() * WEIGHT_STDEV;
            }
        }
    }

    public void save(String filename) throws IOException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");
        writer.printf("w\t");
        for (int i = 0; i < getNVisible(); i ++)
            writer.printf("%s\t", v[i].getName());
        writer.println();
        for (int j = 0; j < getNHidden(); j ++) {
            writer.printf("%s\t", h[j].getName());
            for (int i = 0; i < getNVisible(); i ++) {
                writer.printf("%f\t", w[j][i]);
            }
            writer.println();
        }
        writer.printf("a\t");
        for (int i = 0; i < getNVisible(); i ++)
            writer.printf("%f\t", a[i]);
        writer.println();
        writer.printf("b\n");
        for (int j = 0; j < getNHidden(); j ++)
            writer.printf("%s\t%f\n", h[j].getName(), b[j]);
        writer.close();
    }

    public double logistic(double x) {
        double y = 1.0 / (1.0 + Math.exp(-x));
        return y;
    }

    public void setVisible(Object[] input) {
        for (int i = 0; i < v.length; i ++) {
            if (input[i] != null)
                assignVisible(i, input[i]);
        }
    }

    @Override
    public Object[] encode(Object[] input) {
        Object[] hinst = new Object[this.h.length];
        for (int j = 0; j < h.length; j ++) {
            double net = 0.0;
            for (int i = 0; i < v.length; i ++) {
                if (input[i] != null && lnk[j][i]) {
                    if ((Boolean)input[i])
                        net += w[j][i];
                }
            }
            double p = logistic(net + b[j]);
            Ph[j].set(new double[] {p, 1.0 - p}) ;
            hinst[j] = Ph[j].sample();
        }
        return hinst;
    }

    @Override
    public Object[] decode(Object[] hidden) {
        Object[] vinst = new Object[this.v.length];
        for (int i = 0; i < v.length; i ++) {
            double net = 0.0;
            for (int j = 0; j < h.length; j ++) {
                if (hidden[j] != null && lnk[j][i]) {
                    if ((Boolean)hidden[j])
                        net +=  w[j][i];
                }
            }
            double p = logistic(net + a[i]);
            Pv[i].set(new double[] {p, 1.0 - p}) ;
            vinst[i] = Pv[i].sample();
        }
        return vinst;
    }

    @Override
    public Double[][] getCDGradient(Object[][] minibatch, int niter) {
        Double[][] sum = new Double[getNHidden() + 1][getNVisible() + 1]; // plus one to include biases
        int[][] cnt = new int[getNHidden() + 1][getNVisible() + 1];
        err = 0;
        // TODO: multi-thread the following iteration
        for (int p = 0; p < minibatch.length; p ++) {
            Object[] input0 = minibatch[p];
            Object[] hidd0 = encode(input0);
            // positive signal
            double[][] pos = new double[getNHidden() + 1][getNVisible() + 1];
            for (int i = 0; i < getNVisible(); i ++) {
                boolean doneHidBias = false; // do it once only
                if (minibatch[p][i] != null) {
                    for (int j = 0; j < getNHidden(); j ++) {
                        if (!doneHidBias) {
                            cnt[j][getNVisible()] ++;
                            sum[j][getNVisible()] = 0.0;
                            pos[j][getNVisible()] = Ph[j].get(0);
                        }
                        cnt[j][i] ++;
                        sum[j][i] = 0.0;
                        pos[j][i] = ((Boolean)input0[i]) ? Ph[j].get(0) : 0.0;
                    }
                    if (!doneHidBias)
                        doneHidBias = true;
                    cnt[getNHidden()][i] ++;
                    sum[getNHidden()][i] = 0.0;
                    pos[getNHidden()][i] = ((Boolean)input0[i]) ? 1 : 0;
                }
            }
            for (int n = 0; n < niter; n ++) {
                Object[] input1 = decode_restricted(hidd0, input0);
                for (int i = 0; i < input0.length; i ++) {
                    double x0 = (Boolean)input0[i] ? 1 : 0;
                    double p1 = this.Pv[i].get(0);
                    err += Math.sqrt((x0 - p1)*(x0 - p1));
                }
                Object[] hidd1 = encode(input1);
            }
            for (int i = 0; i < getNVisible(); i ++) {
                boolean doneHidBias = false; // do it once only
                if (minibatch[p][i] != null) {
                    for (int j = 0; j < getNHidden(); j ++) {
                        if (!doneHidBias) {
                            double neg = Ph[j].get(0);
                            sum[j][getNVisible()] += pos[j][getNVisible()] - neg;
                        }
                        double neg = Pv[i].get(0) * Ph[j].get(0);
                        sum[j][i] += pos[j][i] - neg;
                    }
                    if (!doneHidBias)
                        doneHidBias = true;
                    double neg = Pv[i].get(0);
                    sum[getNHidden()][i] += pos[getNHidden()][i] - neg;
                }
            }
        }
        // find average, given counts
        for (int i = 0; i < getNVisible() + 1; i ++) {
            for (int j = 0; j < getNHidden() + 1; j ++) {
                if (cnt[j][i] > 0)
                    sum[j][i] /= cnt[j][i];
                else
                    sum[j][i] = null;
            }
        }
        //System.out.println("\t" + err);
        return sum;
    }

    @Override
    public void setCDGradient(Double[][] delta) {
        for (int j = 0; j < getNHidden() + 1; j ++) {
            for (int i = 0; i < getNVisible() + 1; i ++) {
                if (delta[j][i] != null) {
                    if (j == getNHidden())
                        a[i] += delta[j][i];
                    else if (i == getNVisible())
                        b[j] += delta[j][i];
                    else
                        w[j][i] += delta[j][i];
                }
            }
        }
    }

    public Object[] encode_decode_clamped(Object[] clamp) {
        return encode_decode_clamped(clamp, 0);
    }

    public Object[] encode_decode_clamped(Object[] clamp, int niter) {
        Object[] hinst = encode(clamp);
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i++)
            decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(decoded);
            vinst = decode(hinst);
            for (int i = 0; i < vinst.length; i++)
                decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        }
        return decoded;
    }

    public Object[] encode_decode_restricted(Object[] input) {
        return encode_decode_restricted(input, 0);
    }

    private Object[] decode_restricted(Object[] hinst, Object[] input) {
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (input[i] == null ? null : vinst[i]);
        return decoded;
    }

    public Object[] encode_decode_restricted(Object[] input, int niter) {
        Object[] hinst = encode(input);
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (input[i] == null ? null : vinst[i]);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(decoded);
            vinst = decode(hinst);
            for (int i = 0; i < vinst.length; i++)
                decoded[i] = (input[i] == null ? null : vinst[i]);
        }
        return decoded;
    }

    public Object[] encode_decode_full(Object[] input) {
        return encode_decode_full(input, 0);
    }

    public Object[] encode_decode_full(Object[] input, int niter) {
        Object[] hinst = encode(input);
        Object[] vinst = decode(hinst);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(vinst);
            vinst = decode(hinst);
        }
        return vinst;
    }


}
