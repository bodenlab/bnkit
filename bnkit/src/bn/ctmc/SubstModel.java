/*
 bnkit -- software for building and using Bayesian networks
 Copyright (C) 2014  M. Boden et al.

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

package bn.ctmc;

import asr.ASRRuntimeException;
import bn.prob.EnumDistrib;
import dat.EnumTable;
import dat.EnumVariable;
import dat.Enumerable;
import bn.ctmc.matrix.*;
import bn.math.Matrix.Exp;

import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;

/**
 * Conditional probability table for CTMC based on discrete alphabets
 * @author mikael
 *
 * TODO: optimise the access to conditional probabilities @ time,
 * thread-aware cache and possibly compromise on precision; also be mindful of alphabet size.
 * ALternatively, clone the model for each node.
 *
 */
public abstract class SubstModel {

    final double[][] R; // This is the IRM, sometimes referred to as Q
    final double[] F;   // This is the frequencies of the character states
    final Exp Rexp;     // exp(IRM)
    final Enumerable alpha;
    private EnumTable<EnumDistrib> table = null;

    /**
     * Create time reversible evolutionary model.
     * @param F stationary base frequencies
     * @param S Symmetric, un-scaled version of Q matrix Q_ij = s_ij*pi_j as defined by PAML
     * @param alphabet the values that substitutable variables can take, listed strictly in the order of the array and matrix
     */
    public SubstModel(double[] F, double[][] S, Enumerable alphabet) {
        this(F, S, alphabet, true);
    }

    /**
     * Create evolutionary model.
     * @param F stationary base frequencies
     * @param IRM Instantaneous rate matrix
     * @param alphabet the values that substitutable variables can take, listed strictly in the order of the array and matrix
     * @param symmetric set to true if the IRM is unscaled and symmetric (Sij == Sji) as per PAML format, or
     *                  false if the IRM is the "Q" matrix as often provided in papers
     */
    public SubstModel(double[] F, double[][] IRM, Enumerable alphabet, boolean symmetric) {
        if (IRM.length != F.length)
            throw new IllegalArgumentException("Invalid size of either IRM or F");
        if (alphabet.size() != F.length)
            throw new IllegalArgumentException("Invalid size of alphabet");
        this.F = F;
        R = new double[IRM.length][IRM.length];
        for (int i = 0; i < IRM.length; i ++)  {
            if (IRM[i].length != F.length)
                throw new IllegalArgumentException("IRM must be a square matrix");
            if (symmetric) {
                for (int j = i + 1; j < IRM[i].length; j ++) {
                    double s = IRM[i][j];
                    R[i][j] = s*F[j];
                    R[j][i] = s*F[i];
                }
            } else { // the supplied IRM is Q, hence no conversion necessary
                for (int j = 0; j < IRM[i].length; j ++) {
                    if (i == j)
                        continue;
                    R[i][j] = IRM[i][j];
                }
            }
        }
        this.alpha = alphabet;
        SubstModel.makeValid(R);
        SubstModel.normalize(F, R);
        Rexp = new Exp(R);
    }

    /**
     * Get the name of the evolutionary model
     * @return the name as a text string
     */
    public abstract String getName();

    /**
     * The the frequencies of the character states
     * @return a priori probability of each character
     */
    public double[] getF() {
        return F;
    }

    /**
     * Get the IRM a.k.a. "Q"
     * @return the IRM
     */
    public double[][] getR() {
        return R;
    }

    /**
     * Get the domain that defines the character states.
     * @return the domain
     */
    public Enumerable getDomain() {
        return alpha;
    }
    /**
     * Make it a valid rate matrix (make sum of rows = 0) "in place"
     * @param R the potentially invalid R, to be modified in place
     */
    private static void makeValid(double[][] R) {
        int dim = R.length;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) {
                if (i != j)
                    sum += R[i][j];
            }
            R[i][i] = -sum;
        }
    }

    /**
     * Normalize rate matrix "in place" to one expected substitution per unit time.
     * @param F stationary base frequencies
     * @param R the potentially un-normalised R-matrix, to be modified in place
     */
    private static void normalize(double[] F, double[][] R) {
        int dim = R.length;
        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
            sum += -R[i][i]*F[i];
        }
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                R[i][j] = R[i][j]/sum;
        }
    }

    private double[][] probs = null;
    private double time = 0.0;
    // private boolean updateRequired = true;

    // To speed up calculation, store recent probability matrices
    private Map<Double, double[][]> probscache = new HashMap<>();

    // size of cache
    public int CACHE_SIZE = 10000;

    /**
     * Set another cached prob matrix
     */
    private void setCache(double time, double[][] p) {
        if (probscache.size() < CACHE_SIZE)
            probscache.put(time, p);
    }

    private double[][] getCache(double time) {
        double[][] p = probscache.get(time);
        if (p == null) {
            synchronized (this) {
                p = getProbs(time);
                setCache(time, p);
            }
        }
        return p;
    }
    /**
     * Get conditional probability P(X=x|Y=y,time)
     * @param X
     * @param Y
     * @param time
     * @return
     */
    public double getProb(Object X, Object Y, double time) {
        if (this.time != time || probs == null) // only re-compute matrix if time has changed
            probs = getCache(time);
        int index_X = alpha.getIndex(X);
        int index_Y = alpha.getIndex(Y);
        return probs[index_Y][index_X];
    }

    /**
     * Helper method. Returns the corresponding entry from a user supplied
     * probability matrix using the model alphabet.
     */
    public double getProb(Object X, Object Y, double[][] probMatrix) {
        int index_X = alpha.getIndex(X);
        int index_Y = alpha.getIndex(Y);
        return probMatrix[index_Y][index_X];
    }

    /**
     * Get probability P(X=x)
     * @param X
     * @return
     */
    public double getProb(Object X) {
        int index_X = alpha.getIndex(X);
        return F[index_X];
    }

    public EnumDistrib getDistrib(Object Y, double time) {
        if (this.time != time || probs == null || table == null) { // only re-compute matrix if time has changed
            probs = getCache(time);
            table = new EnumTable<>(new EnumVariable(alpha));
            for (int i = 0; i < probs.length; i ++) {
                EnumDistrib d = new EnumDistrib(alpha, probs[i]);
                table.setValue(i, d);
            }
        }
        return table.getValue(new Object[] {Y});
    }

    /**
     * Compute the transition probabilities for an expected distance
     * using the pre-specified rate matrix
     *
     * @param time expected distance
     * @return the conditional probabilities of a symbol at time t+time GIVEN a symbol at time t  [row: X(t)][col: X(t+time)]
     */
    public final double[][] getProbs(double time) {
        this.time = time;
        int i, j, k;
        double temp;
        double[] eval = Rexp.getEigval();
        double[][] ievec = Rexp.getInvEigvec();
        double[][] evec = Rexp.getEigvec();
        double[][] iexp = new double[R.length][R.length];
        double[][] prob = new double[R.length][R.length];
        for (k = 0; k < R.length; k++) {
            temp = Math.exp(time * eval[k]);
            for (j = 0; j < R.length; j++) {
                iexp[k][j] = ievec[k][j] * temp;
            }
        }
        for (i = 0; i < R.length; i++) {
            for (j = 0; j < R.length; j++) {
                temp = 0.0;
                for (k = 0; k < R.length; k++) {
                    temp += evec[i][k] * iexp[k][j];
                }
                prob[i][j] = Math.abs(temp);
            }
        }
        return prob;
    }

    public static SubstModel createModel(String name) {
        if (name.equalsIgnoreCase("Gap")) {
            return new Gap();
        } else if (name.equalsIgnoreCase("Yang")) {
            return new Yang();
        } else if (name.equalsIgnoreCase("JC")) {
            return new JC(1);
        } else if (name.equalsIgnoreCase("GLOOME1")) {
            return new GLOOME1();
        } else if (name.equalsIgnoreCase("WAG")) {
            return new WAG();
        } else if (name.equalsIgnoreCase("LG")) {
            return new LG();
        }   else if (name.equalsIgnoreCase("JTT")) {
            return new JTT();
        } else if (name.equalsIgnoreCase("Dayhoff")) {
            return new Dayhoff();
        } else if (name.equalsIgnoreCase("SIMPLE_3")) {
            return new SIMPLE_3();
        } else
            return null;
    }

    public static class ModelCache {
        final Map<Object, SubstModel> cache = new HashMap<>();
        final Map<Object, Integer> count = new HashMap<>();
        final int maxsize;
        Object replaceme = null; // pointer to tag that should be replaced next
        public ModelCache(int maxsize) {
            this.maxsize = maxsize;
        }
        public int size() {
            return cache.size();
        }
        public void add(SubstModel model, Object tag) {
            if (cache.size() >= maxsize) {
                cache.remove(replaceme);
                count.remove(replaceme);
            }
            cache.put(tag, model);
            count.put(tag, 0);
            replaceme = tag;
        }
        public SubstModel get(Object tag) {
            SubstModel m = cache.get(tag);
            if (m == null)
                return null;
            int cnt = count.get(tag);
            if (tag == replaceme) { // may need to update which tag is at the bottom of the list
                for (Map.Entry<Object, Integer> entry : count.entrySet()) {
                    if (entry.getValue() < cnt + 1) {
                        replaceme = entry.getKey();
                        break;
                    }
                }
            }
            count.put(tag, cnt + 1);
            return m;
        }
    }

    public static void main(String[] argv) {
        SubstModel sm_gap = new Gap();
        SubstModel sm_gloome1 = new GLOOME1();
        SubstModel sm_yang = new Yang();
        SubstModel sm_wag = new WAG();
        SubstModel sm_lg = new LG();
        SubstModel sm_jtt = new JTT();
        SubstModel sm_dh = new Dayhoff();

        System.out.println("R (Gap)");
        bn.math.Matrix.print(sm_gap.getR());
        System.out.println("R (GLOOME1)");
        bn.math.Matrix.print(sm_gloome1.getR());
        System.out.println("R (Yang)");
        bn.math.Matrix.print(sm_yang.getR());

        System.out.println("R (Dayhoff)");
        bn.math.Matrix.print(sm_dh.getR());
        System.out.println("R (WAG)");
        bn.math.Matrix.print(sm_wag.getR());
        System.out.println("R (LG)");
        bn.math.Matrix.print(sm_lg.getR());

        double time = .1;
        System.out.println("\n\nTransition probabilities of R (Gap) @ time = " + time);
        double[][] prob = sm_gap.getProbs(time);
        bn.math.Matrix.print(prob);

        System.out.println("\n\nTransition probabilities of R (GLOOME1) @ time = " + time);
        prob = sm_gloome1.getProbs(time);
        bn.math.Matrix.print(prob);

        System.out.println("\n\nTransition probabilities of R (Yang) @ time = " + time);
        prob = sm_yang.getProbs(time);
        bn.math.Matrix.print(prob);

        System.out.println("\n\nTransition probabilities of R (WAG) @ time = " + time);
        prob = sm_wag.getProbs(time);
        bn.math.Matrix.print(prob);

        System.out.println("\nTransition probabilities of R (LG) @ time = " + time);
        prob = sm_lg.getProbs(time);
        bn.math.Matrix.print(prob);
        bn.math.Matrix.printLaTeX(prob, sm_lg.getDomain().getValues(), sm_lg.getDomain().getValues());

        System.out.println("\nTransition probabilities of R (JTT) @ time = " + time);
        prob = sm_jtt.getProbs(time);
        bn.math.Matrix.print(prob);

        System.out.println("\nTransition probabilities of R (Dayhoff) @ time = " + time);
        prob = sm_dh.getProbs(time);
        bn.math.Matrix.print(prob);
        double p = sm_dh.getProb('K', 'R', time);
        System.out.println("P(K|R) = " + p);
    }
}
