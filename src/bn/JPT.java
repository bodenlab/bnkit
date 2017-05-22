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

import bn.factor.FactorTable;
import dat.EnumVariable;
import dat.EnumTable;
import java.util.List;
import java.util.Map;

/**
 * Data structure for Joint Probability Table (JPT).
 *
 * @author mikael
 */
public class JPT {

    // the JPT keeps track of probabilities (summing to 1.0) for entries specified a list of value-assigned variables.
    public final EnumTable<Double> table;

    /**
     * Construct a JPT from raw counts.
     * Useful after training.
     * This is simply a normalisation exercise.
     * @param counts the table with counts
     */
    public JPT(EnumTable<Double> counts) {
        table = new EnumTable<>(counts.getParents());
        double sum = 0.0;
        for (Map.Entry<Integer, Double> entry : counts.getMapEntries()) {
            sum += entry.getValue();
        }
        for (Map.Entry<Integer, Double> entry : counts.getMapEntries()) {
            table.setValue(entry.getKey(), entry.getValue() / sum);
        }
    }

    /**
     * Construct a JPT from a FactorTable.
     * Useful after training.
     * This is simply a normalisation exercise.
     * @param ft the table with counts
     */
    public JPT(FactorTable ft) {
        table = new EnumTable<>(ft.getEnumVariables());
        double sum = ft.getSum();
        for (Map.Entry<Integer, Double> entry : ft.getMapEntries()) {
            table.setValue(entry.getKey(), entry.getValue() / sum);
        }
    }

    /**
     * Retrieve a single probability from a JPT.
     * If multiple entries match, probabilities will be summed.
     * @param key values assigned to variables, where individual variables can be left unspecified by using the value null
     * @return a single probability
     */
    public double get(Object... key) {
        // check first if the key has wildcards, and if so return the sum of matching entries
        for (int i = 0; i < key.length; i++) {
            if (key[i] == null) {
                return getSum(key);
            }
        }
        Double y = table.getValue(key);
        if (y == null) {
            return 0.0;
        } else {
            return y;
        }
    }

    /**
     * Retrieve a list of probabilities corresponding to a partial key.
     * @param key values assigned to variables, where individual variables can be left unspecified by using the value null
     * @return list of probabilities
     */
    public List<Double> getAll(Object[] key) {
        List<Double> all = table.getValues(key);
        return all;
    }

    /**
     * Retrieve the sum of probabilities for entries matching the partial key.
     * Note that the normal get method is quicker for some cases (and just as flexible)
     * @param key values assigned to variables, where individual variables can be left unspecified by using the value null
     * @return sum of probabilities
     * @see bn.JPT#get(Object[])
     */
    public double getSum(Object[] key) {
        List<Double> all = getAll(key);
        double sum = 0.0;
        for (Double y : all) {
            if (y != null) {
                sum += y;
            }
        }
        return sum;
    }

    /**
     * Transform a table with counts to a JPT.
     * @param ctable table of counts
     * @return a new JPT
     */
    public static JPT toJPT(CountTable ctable) {
        return new JPT(ctable.table);
    }

    /**
     * Retrieve the variables that define the JPT.
     * @return variables
     */
    public List<EnumVariable> getParents() {
        return table.getParents();
    }

    /**
     * Pretty-print the JPT
     */
    public void display() {
        table.display();
    }

    /**
     * Dense print-representation of the JPT.
     * @return string representation useful for debugging
     */
    @Override
    public String toString() {
        String[] parents = table.getLabels();
        StringBuilder sbuf = new StringBuilder();
        for (int i = 0; i < parents.length; i++) {
            sbuf.append(parents[i]).append(i < parents.length - 1 ? "," : "");
        }
        return "JPT(" + sbuf.toString() + ")";
    }
}
