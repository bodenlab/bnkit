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

import bn.Distrib;
import dat.EnumTable;
import bn.JPT;
import dat.Variable;
import java.util.Map;

/**
 * Interface for Bayesian network query results.
 * @author mikael
 */
public interface QueryResult {
    
    JPT getJPT();
    Map<Variable, EnumTable<Distrib>> getNonEnum();
    Map<Variable, Distrib> getNonEnumDistrib();
}
