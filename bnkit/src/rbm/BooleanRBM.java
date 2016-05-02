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
import java.util.Random;

/**
 * A restricted Boltzmann Machine with Boolean visible and hidden nodes.
 * @author mikael
 */
public class BooleanRBM extends AbstractRBM {
    
    private final double[][] w;     // weights
    private final double[]   a;     // bias on visible
    private final double[]   b;     // bias on hidden
        
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

    public double logistic(double x) {
        double y = 1.0 / (1.0 + Math.exp(-x));
        return y;
    }
    
    @Override
    public Object[] encode(Object[] input) {
        Object[] hinst = new Object[this.h.length];
        for (int j = 0; j < h.length; j ++) {
            double net = 0.0;
            for (int i = 0; i < v.length; i ++) {
                if (input[i] != null && lnk[j][i]) {
                    if ((Boolean)input[i])
                        net +=  w[j][i];
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
    
    public Object[] encodedecode(Object[] clamp) {
        Object[] hinst = encode(clamp);
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        return decoded;
    }
}
