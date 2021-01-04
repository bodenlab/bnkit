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

package bn.ctmc.matrix;

import bn.ctmc.SubstModel;
import dat.Enumerable;

import java.util.Arrays;

/**
 *
 * @author mikael
 */
public class JC extends SubstModel {
    // The original theoretical model of DNA evolution (single nucleotide substitution) but here provided in a more flexible, any-alphabet form.

    static double[] F(int n) {
        double[] freqs = new double[n];
        Arrays.fill(freqs, 1.0/n);
        return freqs;
    }

    static double[][] Q(double gamma, int n) {
        double[][] q = new double[n][n];
        for (int row = 0; row < n; row ++) {
            for (int col = 0; col < n; col ++) {
                if (row != col)
                    q[row][col] = gamma / n;
            }
        }
        return q;
    }

    public static Character[] S = {'A','C','G','T'};

    public JC(double mu) {
        this(mu, S);
    }
    public JC(double mu, Object[] alphabet) {
        super(F(alphabet.length), Q(mu, alphabet.length), new Enumerable((alphabet)), false);
    }
    public String getName() {
        return "JC";
    }

}
