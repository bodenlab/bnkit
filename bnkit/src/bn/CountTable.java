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

import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import dat.EnumTable;
import dat.EnumVariable;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * A table with counts for entries of Enumerable variables.
 * @see JPT for a joint probability table
 * @author m.boden
 */
public class CountTable implements Serializable {

    private static final long serialVersionUID = 1L;
    private boolean totalNeedsUpdate = true;
    private double totalCount = 0.0;

    public final EnumTable<Double> table; // table of counts
    private AbstractFactor ftable = null; // factor table

    public CountTable(EnumVariable[] variables) {
        ftable = new DenseFactor(variables);
        List<EnumVariable> list = new ArrayList<>(variables.length);
        list.addAll(Arrays.asList(variables));
        table = new EnumTable<>(list);
    }

    public CountTable(Collection<EnumVariable> variables) {
        table = new EnumTable<>(variables);
    }

    public List<EnumVariable> getParents() {
        return table.getParents();
    }
    
    public double get(Object[] key) {
        int index = table.getIndex(key);
        Double cnt = table.getValue(index);
        if (cnt == null) {
            return 0;
        } else {
            return cnt;
        }
    }

    public double get(int key_index) {
        Double cnt = table.getValue(key_index);
        if (cnt == null) {
            return 0;
        } else {
            return cnt;
        }
    }

    public boolean hasParents() {
        return (table != null);
    }

    public AbstractFactor getFactor() {
        double total = getTotal();
        for (int idx : table) {
            Object[] key = table.getKey(idx);
            ftable.setValue(key, table.getValue(idx) / total);
        }
        return ftable;
    }

    public void put(Object[] key, double count) {
        table.setValue(key, count);
        this.totalNeedsUpdate = true;
    }

    public void put(int key_index, double count) {
        table.setValue(key_index, count);
        this.totalNeedsUpdate = true;
    }

    public int[] getIndices(Object[] key) {
        return table.getIndices(key);
    }
    
    public int getIndex(Object[] key) {
        return table.getIndex(key);
    }

    public double sum(int[] indices) {
        double sum = 0;
        for (int idx : indices)
            sum += get(idx);
        return sum;
    }

    public double sum(Object[] key) {
        return sum(getIndices(key));
    }

    synchronized public void count(Object[] key, double count) {
        int index = table.getIndex(key);
        count(index, count);
    }

    public void count(Object[] key) {
        count(key, 1.0);
        this.totalNeedsUpdate = true;
    }
    
    synchronized public void count(int key_index, double count) {
        Double cnt = table.getValue(key_index);
        if (cnt == null) {
            table.setValue(key_index, count);
        } else {
            table.setValue(key_index, cnt + count);
        }
        this.totalNeedsUpdate = true;
    }

    public void count(int key_index) {
        count(key_index, 1.0);
        this.totalNeedsUpdate = true;
    }
    
    /**
     * Calculate the sum of all entries.
     * @return summed count
     */
    public double getTotal() {
        if (totalNeedsUpdate) {
            for (Double val : table.getValues()) {
                totalCount += val;
            }
            this.totalNeedsUpdate = false;
        }
        return totalCount;
    }

    public void display() {
        table.display();
    }

    @Override
    public String toString() {
        String[] parents = table.getLabels();
        StringBuilder sbuf = new StringBuilder();
        for (int i = 0; i < parents.length; i++) {
            sbuf.append(parents[i]).append(i < parents.length - 1 ? "," : "");
        }
        return "CountTable(" + sbuf.toString() + ")";
    }

}
