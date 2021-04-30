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
 * @author gabe
 */
public class SIMPLE_3 extends SubstModel {

    public static Character[] S = {'V', 'A', 'N'};

    public static double[] F = {
            // V    A       N
            0.3333, 0.3333, 0.3333};
    public static double[][] Q = {
            //V    A       N
            {0.0, 0.2, 0.4, }, //
            {0.2, 0.0, 0.6 },
            {0.2, 0.4, 0.0}};

    public SIMPLE_3() {
        super(F, Q, new Enumerable((S)), false);
    }
    public String getName() {
        return "SIMPLE_3";
    }

}
