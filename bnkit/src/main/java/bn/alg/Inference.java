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
package bn.alg;

import bn.BNet;
import dat.Variable;

/**
 * Interface for inference engines.
 * @author mikael
 */
public interface Inference {

    /**
     * Use instantiation as per specified Bayesian network
     *
     * @param bn
     */
    public void instantiate(BNet bn);

    /**
     * Construct a query handle
     *
     * @param vars the variables that are queried
     * @return the handle
     */
    public Query makeQuery(Variable[] vars);

    /**
     * Infer the probabilities based on the query handle
     *
     * @param q query handle
     * @return the joint probability table with all variables in query handle
     */
    public QueryResult infer(Query q);

}
