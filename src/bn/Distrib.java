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
package bn;

import bn.prob.GaussianDistrib;
import stats.IndelModel;

/**
 * Probability distribution defined. 
 * A variable in a BNode will be assigned an instance of a distribution in every row.
 * @author mikael
 *
 */
public interface Distrib {

    /**
     * Retrieve the probability of the value
     *
     * @param value
     * @return the probability
     */
    public double get(Object value);

    /**
     * Sample from the distribution
     *
     * @return a sample value, drawn from the distribution in proportion to its
     * probability
     */
    public Object sample();

    public static double[] parseParams(String str) throws RuntimeException {
        try {
            String[] parts = str.split(",");
            double[] result = new double[parts.length];
            for (int i = 0; i < parts.length; i++) {
                result[i] = Double.parseDouble(parts[i].trim());
            }
            return result;
        } catch (NumberFormatException e) {
            throw new RuntimeException("Failed to parse doubles: " + str, e);
        }
    }

    public static Distrib create(String distrib_name, String params) {
        double[] params_arr = parseParams(params);
        switch (distrib_name) {
            case "GaussianDistrib":
            case "gaussiandistrib":
            case "Gaussian":
            case "gaussian":
            case "GDT":
            case "GDF":
            case "Normal":
            case "normal":
                return new GaussianDistrib(params_arr[0], params_arr[1]);
        }
        return null;
    }
}
