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

import bn.prob.EnumDistrib;
import dat.EnumTable;
import dat.EnumVariable;
import dat.Enumerable;
import bn.ctmc.matrix.*;
import bn.math.Matrix.Exp;

/**
 *
 * @author mikael
 */
public abstract class SubstModel {
    
    final double[][] R;
    final double[] F;
    final Exp Rexp;
    final Enumerable alpha;
    private EnumTable<EnumDistrib> table = null;

    /**
     * 
     * @param F stationary base frequencies
     * @param Q Q matrix (general reversible model)
     * @param alphabet the values that substitutable variables can take, listed strictly in the order of the array and matrix
     */
    public SubstModel(double[] F, double[][] Q, Enumerable alphabet) {
        if (Q.length != F.length)
            throw new IllegalArgumentException("Invalid size of either Q or F");
        if (alphabet.size() != F.length)
            throw new IllegalArgumentException("Invalid size of alphabet");
        this.F = F;
        R = new double[Q.length][Q.length];
        for (int i = 0; i < Q.length; i ++)  {
            if (Q[i].length != F.length)
                throw new IllegalArgumentException("Q must be a square matrix");
            for (int j = i + 1; j < Q[i].length; j ++) {
                double q = Q[i][j];
                R[i][j] = q*F[j];
                R[j][i] = q*F[i];
            }
        }
        this.alpha = alphabet;
        SubstModel.makeValid(R);
        SubstModel.normalize(F, R);
        Rexp = new Exp(R);
    }

    public abstract String getName();
    
    public double[] getF() {
        return F;
    }
    
    public double[][] getR() {
        return R;
    }
        
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
    
    public double getProb(Object X, Object Y, double time) {
        if (this.time != time || probs == null) // only re-compute matrix if time has changed
            probs = getProbs(time);
        int index_X = alpha.getIndex(X);
        int index_Y = alpha.getIndex(Y);
        return probs[index_Y][index_X];
    }

    /**
     * Helper method. Returns the corresponding entry from an user supplied
     * probability matrix using the model alphabet.
     */
    public double getProb(Object X, Object Y, double[][] probMatrix) {
        int index_X = alpha.getIndex(X);
        int index_Y = alpha.getIndex(Y);
        return probMatrix[index_Y][index_X];
    }
    
    public double getProb(Object X) {
        int index_X = alpha.getIndex(X);
        return F[index_X];
    }
    
    public EnumDistrib getDistrib(Object Y, double time) {
        if (this.time != time || probs == null || table == null) { // only re-compute matrix if time has changed
            probs = getProbs(time);
            table = new EnumTable<>(new EnumVariable(alpha), new EnumVariable(alpha));
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
    private final double[][] getProbs(double time) {
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
    
    public static void main(String[] argv) {
        SubstModel sm_wag = new WAG();
        SubstModel sm_lg = new LG();
        SubstModel sm_jtt = new JTT();
        SubstModel sm_dh = new Dayhoff();

        System.out.println("R (Dayhoff)");
        bn.math.Matrix.print(sm_dh.getR());
        System.out.println("R (WAG)");
        bn.math.Matrix.print(sm_wag.getR());
        System.out.println("\nR (LG)");
        bn.math.Matrix.print(sm_lg.getR());

        double time = 1.0;
        System.out.println("\n\nTransition probabilities of R (WAG) @ time = " + time);
        double[][] prob = sm_wag.getProbs(time);
        bn.math.Matrix.print(prob);
        System.out.println("\nTransition probabilities of R (LG) @ time = " + time);
        prob = sm_lg.getProbs(time);
        bn.math.Matrix.print(prob);
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
