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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * A factor table is a data structure for efficiently processing products and
 * sums over nodes. It can represent a CPT or GDT (when the continuous variable
 * is instantiated). The class implements basic operations over factor tables.
 *
 * @author Mikael Boden
 */
public class FactorTable extends EnumTable<Double> {

    public boolean evidenced = false; // flag to indicate if the Factor was reduced due to evidence (used by inference)
    public boolean function = false;  // flag to indicate if the factor was the result of a product or marginalization

    /**
     * Data structure for holding entries associated with keys. Keys are defined
     * by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). An entry is always a double.
     *
     * @param useParents the variables that are included in the table
     */
    public FactorTable(Collection<EnumVariable> useParents) {
        super(useParents);
    }

    /**
     * Data structure for holding entries associated with keys. Keys are defined
     * by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). An entry is always a double.
     *
     * @param useParents the variables that are included in the table
     */
    public FactorTable(EnumVariable[] useParents) {
        super(useParents);
    }

    @Override
    public int setValue(Object[] key, Double value) {
        return super.setValue(key, value);
    }

    @Override
    public Double getValue(Object[] key) {
        Double y = super.getValue(key);
        if (y == null) {
            return 0.0;
        } else {
            return y;
        }
    }

    @Override
    public Double getValue(int index) {
        Double y = super.getValue(index);
        if (y == null) {
            return 0.0;
        } else {
            return y;
        }
    }

    @Override
    public List<Double> getValues(Object[] key) {
        return super.getValues(key);
    }

    /**
     * Add a value to a specified entry identified by index
     *
     * @param index the entry index
     * @param value the value to be added to the current instance (assumed to be
     * 0 if not instantiated)
     * @return the index
     */
    public int addValue(int index, double value) {
        Double prev = map.get(index);
        if (prev == null) {
            map.put(index, value);
        } else {
            map.put(index, prev + value);
        }
        return index;
    }

    /**
     * Add a value to a specified entry identified by instantiated key
     *
     * @param key the entry key
     * @param value the value to be added to the current instance (assumed to be
     * 0 if not instantiated)
     * @return the index for the modified entry
     */
    public int addValue(Object[] key, double value) {
        int index = this.getIndex(key);
        return addValue(index, value);
    }

    /**
     * Sums-out (marginalizes) one or more parents from the current table. Time
     * complexity is O(n+m*3n) where n is the number of parents and m is the
     * specified number of entries in the current table.
     *
     * @param parentsToSumOut variables that will be summed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public FactorTable marginalize(Collection<EnumVariable> parentsToSumOut) {
        Collection<EnumVariable> newparents = new ArrayList<EnumVariable>(nParents - parentsToSumOut.size());
        for (int i = 0; i < nParents; i++) {
            if (!parentsToSumOut.contains(this.parents.get(i))) {
                newparents.add(this.parents.get(i));
            }
        }
        FactorTable ft = new FactorTable(newparents);
        for (Map.Entry<Integer, Double> entry : this.map.entrySet()) {
            int newindex = maskIndex(entry.getKey(), parentsToSumOut);
            ft.addValue(newindex, entry.getValue());
        }
        ft.function = true;
        return ft;
    }

    /**
     * Sums-out (marginalizes) one or more parents from the current table. Time
     * complexity is O(n+m*3n) where n is the number of parents and m is the
     * specified number of entries in the current table.
     *
     * @param parentsToSumOut variables that will be summed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public FactorTable marginalize(EnumVariable[] parentsToSumOut) {
        return this.marginalize(EnumVariable.toList(parentsToSumOut));
    }

    /**
     * Determine the number of entries that (potentially) need to be computed
     * for a factorisation involving involving a given number of variables
     *
     * @param ft the variables of table
     * @return the (maximum) number of products to compute
     * @see bn.FactorTable#union(java.util.Collection, java.util.Collection) to create the union of variables from two collections of
     * variables
     */
    public static int getComplexity(Collection<EnumVariable> ft) {
        int order = 1;
        for (EnumVariable v1 : ft) {
            order *= v1.getDomain().size();
        }
        return order;
    }

    public static Collection<EnumVariable> union(Collection<EnumVariable> ft1, Collection<EnumVariable> ft2) {
        Collection<EnumVariable> isect = new HashSet<EnumVariable>();
        isect.addAll(ft1);
        isect.addAll(ft2);
        return isect;
    }

    /**
     * Takes the product between the two factors provided as arguments.
     * Variables to be "joined" from the factors must have the same name and
     * implement the same domain for the product to be determined.
     *
     * @param ft1 factor one
     * @param ft2 factor two
     * @return the factor table that is the product of the specified factors
     */
    public static FactorTable product(FactorTable ft1, FactorTable ft2) {
        Map<Integer, Integer> overlap = new HashMap<Integer, Integer>();
        List<EnumVariable> unique_ft2 = new ArrayList<EnumVariable>(); // in ft2 but not in ft1
        List<Integer> unique_ft2_idx = new ArrayList<Integer>(); // in ft2 but not in ft1
        for (int i = 0; i < ft2.nParents; i++) {
            boolean match = false;
            for (int j = 0; j < ft1.nParents; j++) {
                if (ft1.parents.get(j) == ft2.parents.get(i)) {
                    overlap.put(j, i); // link index in ft1 to index in ft2, resulting index in ft3 is the same as in ft1
                    match = true;
                    break;
                }
            }
            if (!match) {
                // ft2.parent[i] did not match any var in ft1
                unique_ft2.add(ft2.parents.get(i));
                unique_ft2_idx.add(i);
            }
        }
        Collection<EnumVariable> newparents = new ArrayList<EnumVariable>(ft1.nParents);
        for (int i = 0; i < ft1.nParents; i++) {
            newparents.add(ft1.parents.get(i));
        }
        newparents.addAll(unique_ft2);
        FactorTable ft3 = new FactorTable(newparents);
        for (Map.Entry<Integer, Double> entry1 : ft1.map.entrySet()) {
            if (entry1.getValue() == 0.0) // if the value associated with the entry is 0, the product will be zero, 
            {
                continue; // ... so we abort
            }
            for (Map.Entry<Integer, Double> entry2 : ft2.map.entrySet()) {
                if (entry2.getValue() == 0.0) // if the value associated with the entry is 0, the product will be zero, 
                {
                    continue; // ... so we abort
                }
                Object[] ft1_key = ft1.getKey(entry1.getKey().intValue());
                Object[] ft2_key = ft2.getKey(entry2.getKey().intValue());
                boolean match = true;
                for (Map.Entry<Integer, Integer> link : overlap.entrySet()) {
                    int index_ft1 = link.getKey();
                    int index_ft2 = link.getValue();
                    if (ft1_key[index_ft1] != ft2_key[index_ft2]) {
                        match = false;
                        break;
                    }
                }
                if (match) { // all linked variables have the same value
                    Object[] ft3_key = new Object[ft3.nParents];
                    for (int i = 0; i < ft1.nParents; i++) {
                        ft3_key[i] = ft1_key[i];
                    }
                    int i = ft1.nParents;
                    for (int ft2_idx : unique_ft2_idx) {
                        ft3_key[i] = ft2_key[ft2_idx];
                        i++;
                    }
                    ft3.addValue(ft3_key, entry1.getValue() * entry2.getValue());
                }
            }
        }
        ft3.function = true;
        return ft3;
    }

    public double getLogLikelihood() {
        double sum = 0.0;
        for (Map.Entry<Integer, Double> entry : this.getMapEntries()) {
            sum += entry.getValue().doubleValue();
        }
        return Math.log(sum);
    }

    public static void main(String[] args) {
        /*
         EnumVariable a = EnumVariable.Boolean();
         EnumVariable d = EnumVariable.Boolean();
         EnumVariable b = EnumVariable.Boolean();
         EnumVariable c = EnumVariable.Number(30);
		
         FactorTable ft1 = new FactorTable(new EnumVariable[] {a, b, d});
         ft1.setValue(new Object[] {true,  false, false}, 0.1);
         ft1.setValue(new Object[] {false,  false, false}, 0.3);
         ft1.setValue(new Object[] {true,  true, false},  0.6);
		
         FactorTable ft2 = new FactorTable(new EnumVariable[] {d, b, c});
         ft2.setValue(new Object[] {true, true, 0}, 0.3);
         ft2.setValue(new Object[] {true, false, 1}, 0.2);
         ft2.setValue(new Object[] {false, true,  1}, 0.4);
         ft2.setValue(new Object[] {false, false, 0}, 0.1);
		
         ft1.display();
         ft2.display();
		
         FactorTable ft3 = FactorTable.product(ft1, ft2);
		
         ft3.display();
         */

        EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Number(4);
        EnumVariable v3 = Predef.Nominal(new String[]{"Yes", "No", "Maybe"});

        System.out.println(EnumVariable.pool.get(v1));
        System.out.println(EnumVariable.pool.get(v2));
        System.out.println(EnumVariable.pool.get(v3));

        FactorTable ft3 = new FactorTable(new EnumVariable[]{v1, v2, v3});
        ft3.setValue(new Object[]{true, 1, "No"}, 0.05);
        ft3.setValue(new Object[]{true, 0, "No"}, 0.02);
        ft3.setValue(new Object[]{false, 3, "No"}, 0.17);
        ft3.setValue(new Object[]{false, 0, "Yes"}, 0.07);
        ft3.setValue(new Object[]{false, 1, "Yes"}, 0.04);
        ft3.setValue(new Object[]{true, 2, "Yes"}, 0.13);
        ft3.setValue(new Object[]{true, 0, "Yes"}, 0.18);
        ft3.setValue(new Object[]{true, 3, "Yes"}, 0.08);
        ft3.display();
        FactorTable ft2 = ft3.marginalize(new EnumVariable[]{v2});
        ft2.display();
        FactorTable ft1 = ft3.marginalize(new EnumVariable[]{v1, v3});
        ft1.display();
        System.out.println("Size of ft3: " + ft3.getSize());
    }

    public String toString() {
        StringBuffer sbuf = new StringBuffer("F(");
        for (Variable v : this.parents) {
            sbuf.append(v.toString() + ";");
        }
        return sbuf.toString() + ")";
    }

}

class FactorTableRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public FactorTableRuntimeException(String string) {
        message = string;
    }
}
