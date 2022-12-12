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
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A restricted Boltzmann Machine with Boolean visible and hidden nodes.
 * @author mikael
 */
public class BooleanRBM extends AbstractRBM {

    private  double[][] w;     // weights
    private  double[]   a;     // bias on visible
    private  double[]   b;     // bias on hidden

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

    private void setWeights(List<Double[]> weights) {
        this.a = new double[v.length];
        this.b = new double[h.length];
        this.w = new double[h.length][v.length];
        this.lnk = new boolean[h.length][v.length];
        for (int j = 0; j < weights.size(); j ++) {
            for (int i = 0; i < weights.get(j).length; i ++) {
                if (weights.get(j)[i] == null)
                    this.lnk[j][i] = false;
                else {
                    this.w[j][i] = weights.get(j)[i];
                    this.lnk[j][i] = true;
                }
            }
        }
    }

    public BooleanRBM(File file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line = br.readLine();
        Character collecting = '?';
        int collectCount = 0;
        int nVis = 0;
        int nHid = 0;
        List<Double[]> weights = new ArrayList<>();
        List<String> hidvars = new ArrayList<>();
        while (line != null) {
            String[] words = line.split("\t");
            if (words[0].startsWith("*w")) {
                collecting = 'w';
                collectCount = 0;
            } else if (words[0].startsWith("*a")) {
                if (collecting == 'w') {
                    nHid = weights.size();
                    this.h = new EnumVariable[nHid];
                    this.Ph = new EnumDistrib[nHid];
                    for (int j = 0; j < nHid; j ++) {
                        this.h[j] = new EnumVariable(Enumerable.bool, hidvars.get(j));
                        this.h[j].setPredef("Boolean");
                        this.Ph[j] = new EnumDistrib(Enumerable.bool);
                    }
                    setWeights(weights);
                }
                collecting = 'a';
                collectCount = 0;
            } else if (words[0].startsWith("*b")) {
                if (collecting == 'w') {
                    nHid = weights.size();
                    this.h = new EnumVariable[nHid];
                    this.Ph = new EnumDistrib[nHid];
                    for (int j = 0; j < nHid; j ++) {
                        this.h[j] = new EnumVariable(Enumerable.bool, hidvars.get(j));
                        this.h[j].setPredef("Boolean");
                        this.Ph[j] = new EnumDistrib(Enumerable.bool);
                    }
                    setWeights(weights);
                }
                collecting = 'b';
                collectCount = 0;
            }

            if (collecting.equals('w')) {
                if (collectCount == 0) { // header: create variables
                    nVis = words.length - 1;
                    this.v = new EnumVariable[nVis];
                    this.Pv = new EnumDistrib[nVis];
                    for (int i = 0; i < nVis; i ++) {
                        this.v[i] = new EnumVariable(Enumerable.bool, words[1 + i]);
                        this.v[i].setPredef("Boolean");
                        this.Pv[i] = new EnumDistrib(Enumerable.bool);
                    }
                    collectCount = 1;
                } else { // instantiate weights
                    Double[] hidw = new Double[nVis]; // the weights for one hidden node
                    hidvars.add(words[0]);
                    for (int i = 0; i < nVis; i ++) {
                        try {
                            hidw[i] = Double.parseDouble(words[1 + i]);
                        } catch (NumberFormatException e) {
                            hidw[i] = null;
                        }
                    }
                    weights.add(hidw);
                    collectCount += 1;
                }
            } else if (collecting.equals('a')) {
                if (collectCount == 0 && this.a != null) { // header: instantiate biases
                    for (int i = 0; i < nVis; i++) {
                        this.a[i] = Double.parseDouble(words[1 + i]);
                    }
                    collectCount = 1;
                }
            } else if (collecting.equals('b')) {
                if (collectCount == 0 && this.b != null) { // header: do nothing
                    collectCount = 1;
                } else if (words[0].equals(this.h[collectCount - 1].getName())) {
                    this.b[collectCount - 1] = Double.parseDouble(words[1]);
                    collectCount += 1;
                }
            }
            line = br.readLine();
        }
        rand = new Random(RANDOM_SEED);
    }

    public BooleanRBM(String filename) throws IOException {
        this(new File(filename));
    }

    public void save(String filename) throws IOException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");
        writer.printf("*w\t");
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
        writer.printf("*a\t");
        for (int i = 0; i < getNVisible(); i ++)
            writer.printf("%f\t", a[i]);
        writer.println();
        writer.printf("*b\n");
        for (int j = 0; j < getNHidden(); j ++)
            writer.printf("%s\t%f\n", h[j].getName(), b[j]);
        writer.close();
    }

    public void setVisible(Object[] input) {
        for (int i = 0; i < v.length; i ++) {
            if (input[i] != null)
                assignVisible(i, input[i]);
        }
    }

    /**
     * Calculate a hidden state vector from a given visible/input state vector;
     * note that this function updates the (internal) probabilities of the hidden state.
     * @param input input state (array)
     * @param hinst hidden state (array); over-written as a result of calling this function
     * @return hidden state (array)
     */
    @Override
    public Object[] encode(Object[] input, Object[] hinst) {
        for (int j = 0; j < h.length; j ++) {
            double net = b[j];
            for (int i = 0; i < v.length; i ++) {
                if (input[i] != null && lnk[j][i]) {
                    if ((Boolean)input[i])
                        net += w[j][i];
                }
            }
            double p = logistic(net);
            Ph[j].set(new double[] {p, 1.0 - p}) ;
            hinst[j] = Ph[j].sample();
        }
        return hinst;
    }

    /**
     * Calculate the visible state from the hidden;
     * note that this function updates the (internal) probabilities of the visible state.
     * @param hidden hidden state (array)
     * @param vinst visible state array; over-written when called to enable the calling function to re-use memory
     * @return visible state (array)
     */
    @Override
    public Object[] decode(Object[] hidden, Object[] vinst) {
        for (int i = 0; i < v.length; i ++) {
            double net = a[i];
            for (int j = 0; j < h.length; j ++) {
                if (hidden[j] != null && lnk[j][i]) {
                    if ((Boolean)hidden[j])
                        net +=  w[j][i];
                }
            }
            double p = logistic(net);
            Pv[i].set(new double[] {p, 1.0 - p}) ;
            vinst[i] = Pv[i].sample();
        }
        return vinst;
    }

    @Override
    public Double[][] getCDGradient(Object[][] minibatch, int niter) {
        Double[][] sum = new Double[getNHidden() + 1][getNVisible() + 1]; // plus one to include bias weights
        int[][] cnt = new int[getNHidden() + 1][getNVisible() + 1];
        err = 0;
        // TODO: multi-thread the following iteration
        for (int p = 0; p < minibatch.length; p ++) {
            Object[] input0 = minibatch[p];
            Object[] hidd0 = encode(input0);
            // positive signal
            double[][] pos = new double[getNHidden() + 1][getNVisible() + 1];
            for (int i = 0; i < getNVisible(); i ++) {
                if (minibatch[p][i] != null) {
                    for (int j = 0; j < getNHidden(); j ++) {
                        if (i == 0) { // do this only for one round (bias)
                            cnt[j][getNVisible()] ++;
                            sum[j][getNVisible()] = 0.0;
                            pos[j][getNVisible()] = Ph[j].get(0);
                        }
                        if (lnk[j][i]) {
                            cnt[j][i] ++;
                            pos[j][i] = ((Boolean)input0[i]) ? 1 /*true*/ * Ph[j].get(0) : 0.0;
                        }
                        sum[j][i] = 0.0;
                    }
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
                if (minibatch[p][i] != null) {
                    for (int j = 0; j < getNHidden(); j ++) {
                        if (i == 0) { // do this one round only
                            double neg = Ph[j].get(0);
                            sum[j][getNVisible()] += pos[j][getNVisible()] - neg;
                        }
                        if (lnk[j][i]) {
                            double neg = Pv[i].get(0) * Ph[j].get(0);
                            sum[j][i] += pos[j][i] - neg;
                        }
                    }
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


}
