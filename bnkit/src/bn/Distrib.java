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


}
