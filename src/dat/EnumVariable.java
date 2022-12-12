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

import java.util.ArrayList;
import java.util.List;

/**
 * Enumerable variable class.
 *
 * @author mikael
 */
public class EnumVariable extends Variable<Enumerable> {

    public EnumVariable(Enumerable domain) {
        super(domain);
    }

    public EnumVariable(Enumerable domain, String name) {
        super(domain, name);
    }

    public int size() {
        return super.getDomain().order;
    }

    public int getIndex(Object value) {
        return super.getDomain().getIndex(value);
    }

    public static List<EnumVariable> toList(EnumVariable[] variables) {
        if (variables == null) {
            return null;
        }
        List<EnumVariable> list = new ArrayList<EnumVariable>(variables.length);
        for (EnumVariable var : variables) {
            list.add(var);
        }
        return list;
    }

}
