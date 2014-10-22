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

package dat;

/**
 * Continuous domain definition. 
 * For checking validity of values for variables that belong to this domain.
 * @author mikael
 */
public class Continuous implements Domain {

    public boolean isValid(Object value) {
        try {
            Double x = (Double) value;
            return true;
        } catch (ClassCastException e) {
            return false;
        }
    }
}
