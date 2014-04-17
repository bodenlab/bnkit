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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A factor table is a data structure for efficiently processing products and
 * sums over nodes--sometimes referred to as "local" computations. 
 * It can represent a CPT or GDT (when the continuous variable
 * is instantiated or when it occurs as a child of an enumerable variable).
 * 
 * @author Mikael Boden
 */
public class FactorTable extends EnumTable<Double> {

    public boolean evidenced = false; // flag to indicate if the Factor was reduced due to evidence (used by inference)
    public boolean function = false;  // flag to indicate if the factor was the result of a product or marginalization
    /**
     * To manage variables that cannot be directly factorized, their distributions corresponding to 
     * the keyed/instantiated variables are stored as separate tables, for "lazy evaluation".
     * These are empty to begin with.
     */
    private Map<Variable, EnumTable<Distrib>> nonEnumTables = null;
    
    /**
     * Data structure for holding entries associated with keys. Keys are defined
     * by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). An entry is always a double.
     *
     * @param enumVars the enumerable variables that are included in the table
     */
    public FactorTable(Collection<EnumVariable> enumVars) {
        super(enumVars);
    }

    /**
     * Data structure for holding entries associated with keys. Keys are defined
     * by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). An entry is always a double.
     *
     * @param enumVars the enumerable variables that are included in the table
     */
    public FactorTable(EnumVariable[] enumVars) {
        super(enumVars);
    }

    /**
     * FactorTable is a data structure for holding entries associated with keys. 
     * Keys are defined by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). 
     * 
     * Use this constructor for FactorTable with non-Enumerable variables to be evaluated lazily.
     *
     * @param enumVars the variables that are included in the table
     * @param nonEnumVars the variables that are included but not evaluated explicitly as part of the standard operations
     */
    public FactorTable(Collection<EnumVariable> enumVars, Collection<Variable> nonEnumVars) {
        super(enumVars);
        nonEnumTables = new HashMap<>();
        for (Variable v : nonEnumVars)
            nonEnumTables.put(v, new EnumTable<Distrib>(enumVars));
    }

    /**
     * FactorTable is a data structure for holding entries associated with keys. 
     * Keys are defined by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). 
     * 
     * Use this constructor for FactorTable with non-Enumerable variables to be evaluated lazily.
     *
     * @param enumVars the variables that are included in the table
     * @param nonEnumVars the variables that are included but not evaluated explicitly as part of the standard operations
     */
    public FactorTable(EnumVariable[] enumVars, Variable[] nonEnumVars) {
        super(enumVars);
        nonEnumTables = new HashMap<>();
        for (Variable v : nonEnumVars)
            nonEnumTables.put(v, new EnumTable<Distrib>(enumVars));
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key the key identifying the setting of the enumerable variables
     * @param value the value to assign the entry
     * @return the index of the entry
     */
    @Override
    public int setValue(Object[] key, Double value) {
        return super.setValue(key, value);
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key the key identifying the setting of the enumerable variables
     * @param value the value to assign the entry
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setValue(Object[] key, Double value, Variable nonenum, Distrib d) {
        int index = super.setValue(key, value);
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        return nonEnumTable.setValue(index, d);
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(Object[] key, Variable nonenum, Distrib d) {
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        int key_index = nonEnumTable.getIndex(key);
        return nonEnumTable.setValue(key_index, d);
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(int key_index, Variable nonenum, Distrib d) {
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        return nonEnumTable.setValue(key_index, d);
    }

    /**
     * Add a non-enumerable distribution as required by marginalisation of enumerable variables.
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @param weight the weight assigned to the distribution
     * @return the index for the key
     */
    public int addDistrib(Object[] key, Variable nonenum, Distrib d, double weight) {
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        int key_index = nonEnumTable.getIndex(key);
        Distrib prev_d = nonEnumTable.getValue(key_index);
        if (prev_d == null) {
            MixtureDistrib md = new MixtureDistrib(d, weight);
            return nonEnumTable.setValue(key_index, md);
        } else {
            ((MixtureDistrib)prev_d).addDistrib(d, weight);
            return key_index;
        }
    }

    /**
     * Add a non-enumerable distribution as required by marginalisation of enumerable variables.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @param weight the weight assigned to the distribution
     * @return the index for the key
     */
    public int addDistrib(int key_index, Variable nonenum, Distrib d, double weight) {
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        Distrib prev_d = nonEnumTable.getValue(key_index);
        if (prev_d == null) {
            MixtureDistrib md = new MixtureDistrib(d, weight);
            return nonEnumTable.setValue(key_index, md);
        } else {
            ((MixtureDistrib)prev_d).addDistrib(d, weight);
            return key_index;
        }
    }

    public Distrib getDistrib(Object[] key, Variable nonenum) {
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        return nonEnumTable.getValue(key);
    }
    
    public Distrib getDistrib(int key_index, Variable nonenum) {
        EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
        if (nonEnumTable == null)
            throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
        return nonEnumTable.getValue(key_index);
    }
    
    public Set<Variable> getNonEnumVariables() {
        if (nonEnumTables == null)
            return Collections.EMPTY_SET;
        else 
            return nonEnumTables.keySet();
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
     * This method has been partially adapted to accommodate continuous (non-enumerable) variables
     * that are associated with summed-out variables:
     * densities are combined into mixtures, which need to be resolved later.
     *
     * @param parentsToSumOut variables that will be summed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public FactorTable marginalize(Collection<EnumVariable> parentsToSumOut) {
        Collection<EnumVariable> newparents = new ArrayList<>(nParents - parentsToSumOut.size());
        for (int i = 0; i < nParents; i++) {
            if (!parentsToSumOut.contains(this.parents.get(i))) {
                newparents.add(this.parents.get(i));
            }
        }
        FactorTable ft; 
        if (nonEnumTables == null)
            ft = new FactorTable(newparents);
        else
            ft = new FactorTable(newparents, this.getNonEnumVariables());
        for (Map.Entry<Integer, Double> entry : this.map.entrySet()) {
            int oldindex = entry.getKey().intValue();
            double weight = entry.getValue().doubleValue();
            int newindex = maskIndex(oldindex, parentsToSumOut);
            ft.addValue(newindex, weight);
            for (Variable nonenum : this.getNonEnumVariables())
                ft.addDistrib(newindex, nonenum, getDistrib(oldindex, nonenum), weight);
        }
        ft.function = true;
        return ft;
    }

    /**
     * Sums-out (marginalizes) one or more parents from the current table. Time
     * complexity is O(n+m*3n) where n is the number of parents and m is the
     * specified number of entries in the current table.
     * 
     * This method has been partially adapted for continuous (non-enumerable) variables.
     * Specifically, it is able to sum out enumerable variables even if continuous variables 
     * are uninstantiated, by convolution of densities into mixtures.
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
        Collection<EnumVariable> isect = new HashSet<>();
        isect.addAll(ft1);
        isect.addAll(ft2);
        return isect;
    }

    /**
     * Takes the product between the two factors provided as arguments.
     * Variables to be "joined" from the factors must have the same name and
     * implement the same domain for the product to be determined.
     *
     * This method has been partially adapted for continuous variables, but it assumes
     * that continuous variables cannot be reached by multiple paths. For that
     * to work, they have to be joined the same way as discrete variables. Currently,
     * continuous variables are just added without checking if the other factors has them.
     *
     * @param ft1 factor one
     * @param ft2 factor two
     * @return the factor table that is the product of the specified factors
     */
    public static FactorTable product(FactorTable ft1, FactorTable ft2) {
        Map<Integer, Integer> overlap = new HashMap<>();
        List<EnumVariable> unique_ft2 = new ArrayList<>(); // in ft2 but not in ft1
        List<Integer> unique_ft2_idx = new ArrayList<>(); // in ft2 but not in ft1
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
        Collection<EnumVariable> newparents = new ArrayList<>(ft1.nParents);
        for (int i = 0; i < ft1.nParents; i++) {
            newparents.add(ft1.parents.get(i));
        }
        newparents.addAll(unique_ft2);
        Collection<Variable> newEnumParents = new ArrayList<>(ft1.getNonEnumVariables());
        newEnumParents.addAll(ft2.getNonEnumVariables());
        FactorTable ft3 = null;
        if (newEnumParents.size() > 0)
            ft3 = new FactorTable(newparents, newEnumParents);
        else
            ft3 = new FactorTable(newparents);
        for (Map.Entry<Integer, Double> entry1 : ft1.map.entrySet()) {
            if (entry1.getValue() == 0.0) // if the value associated with the entry is 0, the product will be zero, 
            {
                continue; // ... so we abort, to save time as there is no point in doing any calculations
            }
            for (Map.Entry<Integer, Double> entry2 : ft2.map.entrySet()) {
                if (entry2.getValue() == 0.0) // if the value associated with the entry is 0, the product will be zero, 
                {
                    continue; // ... so we abort as there is no point in doing any more calculcations
                }
                int ft1_index = entry1.getKey().intValue();
                int ft2_index = entry2.getKey().intValue();
                Object[] ft1_key = ft1.getKey(ft1_index);
                Object[] ft2_key = ft2.getKey(ft2_index);
                boolean match = true;
                for (Map.Entry<Integer, Integer> link : overlap.entrySet()) {
                    if (ft1_key[link.getKey()] != ft2_key[link.getValue()]) {
                        match = false;
                        break;
                    }
                }
                if (match) { // all linked variables have the same value
                    Object[] ft3_key = new Object[ft3.nParents];
                    System.arraycopy(ft1_key, 0, ft3_key, 0, ft1.nParents);
                    int i = ft1.nParents;
                    for (int ft2_idx : unique_ft2_idx) {
                        ft3_key[i] = ft2_key[ft2_idx];
                        i++;
                    }
                    int ft3_index = ft3.addValue(ft3_key, entry1.getValue() * entry2.getValue());
                    for (Variable ft1_nonenum : ft1.getNonEnumVariables())
                        ft3.setDistrib(ft3_index, ft1_nonenum, ft1.getDistrib(ft1_index, ft1_nonenum));
                    for (Variable ft2_nonenum : ft2.getNonEnumVariables())
                        ft3.setDistrib(ft3_index, ft2_nonenum, ft2.getDistrib(ft2_index, ft2_nonenum));
                }
            }
        }
        ft3.function = true;
        return ft3;
    }

    public double getLogLikelihood() {
        double sum = 0.0;
        if (map.isEmpty()) 
            return 0.0;
        for (Map.Entry<Integer, Double> entry : this.getMapEntries()) {
            sum += entry.getValue().doubleValue();
        }
        if (sum == 0.0)
            throw new RuntimeException("All outcomes are impossible: Log-likelihood is " + Math.log(sum));
        return Math.log(sum);
    }

    public static void main(String[] args) {
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

    @Override
    public String toString() {
        StringBuilder sbuf = new StringBuilder("F(");
        for (Variable v : this.parents) {
            sbuf.append(v.toString()).append(";");
        }
        sbuf.append(";");
        for (Variable v : this.getNonEnumVariables()) {
            sbuf.append(v.toString()).append(";");
        }
        return sbuf.toString() + ")";
    }
    
    private String constantLength(String s, int len) {
        if (s.length() > len)
            return s.substring(0, len);
        else return String.format("%-10s", s);
        
    }
    
    @Override
    public void display() {
        System.out.print("Idx ");
        for (int j = 0; j < this.nParents; j++)
            System.out.print(String.format("[%10s]", constantLength(this.parents.get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" F");
        for (int i = 0; i < this.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = this.getKey(i);
            for (Object key1 : key)
                System.out.print(String.format(" %-10s ", constantLength(key1.toString(), 10)));
            for (Variable nonenum : nonenums) {
                Distrib d = this.getDistrib(i, nonenum);
                System.out.print(String.format(" %s ", d.toString()));
                //System.out.print(String.format(" %-10s ", constantLength(d.toString(), 10)));
            }
            Object val = this.getValue(i);
            if (val != null) 
                System.out.println(String.format(" %7.5f", this.getValue(i)));
            else
                System.out.println(" null ");
        }
    }

    public void displaySampled() {
        System.out.print("Idx ");
        for (int j = 0; j < this.nParents; j++)
            System.out.print(String.format("[%10s]", constantLength(this.parents.get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" F");
        for (int i = 0; i < this.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = this.getKey(i);
            for (Object key1 : key)
                System.out.print(String.format(" %-10s ", constantLength(key1.toString(), 10)));
            for (Variable nonenum : nonenums) {
                Distrib d = this.getDistrib(i, nonenum);
                double sum = 0;
                for (int j = 0; j < 1000; j ++) {
                    Double sample = (Double)d.sample();
                    sum += sample.doubleValue();
                }
                System.out.print(String.format(" %5.3f ", sum/1000.0));
            }
            Object val = this.getValue(i);
            if (val != null) 
                System.out.println(String.format(" %7.5f", this.getValue(i)));
            else
                System.out.println(" null ");
        }
    }
}

class FactorTableRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public FactorTableRuntimeException(String string) {
        message = string;
    }
}
