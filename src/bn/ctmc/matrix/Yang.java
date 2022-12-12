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

/**
 *
 * @author mikael
 */
public class Yang extends SubstModel {
    // Empirical model of DNA evolution (single nucleotide substitution)
    // Yang Z (1994) Estimating the pattern of nucleotide substitution. J Mol Evol 39:105–111
    // (Table 1: Using the best tree and the REV model. Note column order.)
    public static Character[] S = {'A','C','G','T'};

    public static double[] F = {
//            0.25, 0.25, 0.25, 0.25};
            // A      C   	 G      T
              0.308, 0.185, 0.199, 0.308};
    public static double[][] Q = {
      //A      C      G      T
    {0.000, 0.139, 0.566, 0.115}, //
    {0.232, 0.000, 0.223, 0.882},
    {0.877, 0.207, 0.000, 0.215},
    {0.115, 0.530, 0.139, 0.000}};

    public Yang() {
        super(F, Q, new Enumerable((S)), false);
    }
    public String getName() {
        return "Yang";
    }

}
