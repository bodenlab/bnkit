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
package bn.factor;

import bn.Distrib;
import bn.Predef;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import bn.prob.MixtureDistrib;
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
public class FactorTable {

    private EnumTable<Double> enumTable = null;
    private Double atomic = null; // if there are no variables, here's where the only factor is stored
    
    private boolean normalized = false; // flag to indicate if this factor table is normalized or not
    public boolean evidenced = false; // flag to indicate if the Factor was reduced due to evidence (used by inference)
    public boolean function = false;  // flag to indicate if the factor was the result of a product or marginalization
    
    /**
     * To manage variables that cannot be directly factorized, their distributions corresponding to 
     * the keyed/instantiated variables are stored as separate tables, currently for "lazy evaluation",
     * later marginalization. These are empty to begin with.
     */
    private Map<Variable, EnumTable<Distrib>> nonEnumTables = null;
    
    private Map<Variable, Distrib> atomicNonEnumDistribs = null; // where non-enumerable variables are accessed if no enumerable variables to condition with
    
    public FactorTable() {
        atomic = new Double(1.0);
    }
    
    /**
     * Data structure for holding entries associated with keys. Keys are defined
     * by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). An entry is always a double.
     *
     * @param enumVars the enumerable variables that are included in the table
     */
    public FactorTable(Collection<EnumVariable> enumVars) {
        if (enumVars.size() > 0)
            enumTable = new EnumTable<>(enumVars);
        else
            atomic = new Double(1.0);
    }

    /**
     * Data structure for holding entries associated with keys. Keys are defined
     * by a size (of elements), and for each element, a valence (a number ]2,n[
     * where n>2). An entry is always a double.
     *
     * @param enumVars the enumerable variables that are included in the table
     */
    public FactorTable(EnumVariable[] enumVars) {
        if (enumVars.length > 0)
            enumTable = new EnumTable<>(enumVars);
        else
            atomic = new Double(1.0);
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
        if (enumVars.size() > 0)
            enumTable = new EnumTable<>(enumVars);
        else
            atomic = new Double(1.0);
        if (nonEnumVars.size() > 0) {
            if (enumTable == null)
                atomicNonEnumDistribs = new HashMap<>();
            else
                nonEnumTables = new HashMap<>();
            for (Variable v : nonEnumVars) {
                if (enumTable == null)
                    atomicNonEnumDistribs.put(v, null);
                else
                    nonEnumTables.put(v, new EnumTable<Distrib>(enumVars));
            }
        }
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
        if (enumVars.length > 0)
            enumTable = new EnumTable<>(enumVars);
        else
            atomic = new Double(1.0);
        if (nonEnumVars.length > 0) {
            nonEnumTables = new HashMap<>();
            for (Variable v : nonEnumVars) {
                if (enumTable == null)
                    atomicNonEnumDistribs.put(v, null);
                else
                    nonEnumTables.put(v, new EnumTable<Distrib>(enumVars));
            }
        }
    }

    public boolean hasNonEnumVariables() {
        if (this.nonEnumTables != null) {
            return this.nonEnumTables.size() > 0;
        } 
        if (atomicNonEnumDistribs != null) {
            return this.atomicNonEnumDistribs.size() > 0;
        }
        return false;
    }
    
    /**
     * Set the value in a FactorTable. Use with care.
     * @param key the key identifying the setting of the enumerable variables
     * @param value the value to assign the entry
     * @return the index of the entry
     */
    public int setValue(Object[] key, Double value) {
        normalized = false;
        if (enumTable != null)
            return enumTable.setValue(key, value);
        else {
            atomic = value;
            return -1; // no index since there are no enum variables
        }
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key the key identifying the setting of the enumerable variables
     * @param value the value to assign the entry
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry (-1 if no enumerable variables)
     */
    public int setValue(Object[] key, Double value, Variable nonenum, Distrib d) {
        normalized = false;
        if (enumTable != null) {
            int index = enumTable.setValue(key, value);
            EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
            if (nonEnumTable == null)
                throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
            return nonEnumTable.setValue(index, d);
        } else {
            atomic = value;
            atomicNonEnumDistribs.put(nonenum, d);
            return -1; // no index since there are no enum variables
        }
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(Object[] key, Variable nonenum, Distrib d) {
        if (enumTable != null) {
            EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
            if (nonEnumTable == null)
                throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
            int key_index = nonEnumTable.getIndex(key);
            return nonEnumTable.setValue(key_index, d);
        } else {
            atomicNonEnumDistribs.put(nonenum, d);
            return -1; // no index since there are no enum variables
        }
    }

    /**
     * Set the value in a FactorTable. Use with care.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(int key_index, Variable nonenum, Distrib d) {
        if (enumTable != null) {
            EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
            if (nonEnumTable == null)
                throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
            return nonEnumTable.setValue(key_index, d);
        } else {
            atomicNonEnumDistribs.put(nonenum, d);
            return -1; // no index since there are no enum variables
        }
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
        if (enumTable != null) {
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
            }
            return key_index;
        } else {
            Distrib prev_d = atomicNonEnumDistribs.get(nonenum);
            if (prev_d == null) {
                MixtureDistrib md = new MixtureDistrib(d, weight);
                atomicNonEnumDistribs.put(nonenum, md);
            } else {
                ((MixtureDistrib)prev_d).addDistrib(d, weight);
            }
            return -1;
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
        if (enumTable != null) {
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
        } else {
            Distrib prev_d = atomicNonEnumDistribs.get(nonenum);
            if (prev_d == null) {
                MixtureDistrib md = new MixtureDistrib(d, weight);
                atomicNonEnumDistribs.put(nonenum, md);
            } else {
                ((MixtureDistrib)prev_d).addDistrib(d, weight);
            }
            return -1;
        }
    }

    public Distrib getDistrib(Object[] key, Variable nonenum) {
        if (enumTable != null) {
            EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
            if (nonEnumTable == null)
                throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
            return nonEnumTable.getValue(key);
        } else {
            return atomicNonEnumDistribs.get(nonenum);
        }
    }
    
    public Distrib getDistrib(int key_index, Variable nonenum) {
        if (enumTable != null) {
            EnumTable<Distrib> nonEnumTable = nonEnumTables.get(nonenum);
            if (nonEnumTable == null)
                throw new FactorTableRuntimeException("Non-enumerable variable " + nonenum.getName() + " is not valid for FactorTable");
            return nonEnumTable.getValue(key_index);
        } else {
            return atomicNonEnumDistribs.get(nonenum);
        }
    }
    
    public List<EnumVariable> getEnumVariables() {
        if (enumTable != null)
            return enumTable.getParents();
        else
            return Collections.EMPTY_LIST;
    }
    
    public Set<Variable> getNonEnumVariables() {
        if (enumTable != null) {
            if (nonEnumTables == null)
                return Collections.EMPTY_SET;
            else 
                return nonEnumTables.keySet();
        } else {
            if (atomicNonEnumDistribs == null)
                return Collections.EMPTY_SET;
            else 
                return atomicNonEnumDistribs.keySet();
        }
    }
    
    public Set<Map.Entry<Integer, Double>> getMapEntries() {
    	if (enumTable == null){
    		return Collections.EMPTY_SET;
    	} else {
    		return enumTable.getMapEntries();
    	}
    }
    
    public Object[] getKey(int index) {
        return enumTable.getKey(index);
    }
    
    public boolean isAtomic() {
        return enumTable == null;
    }
    
    public Collection<Double> getValues() {
        if (enumTable != null) {
            return enumTable.getValues();
        } else {
            ArrayList<Double> v = new ArrayList<>();
            v.add(atomic);
            return v;
        }
    }
    
    public Double getValue(Object[] key) {
        if (enumTable != null) {
            Double y = enumTable.getValue(key);
            if (y == null) {
                return 0.0;
            } else {
                return y;
            }
        } else {
            if (atomic == null)
                return 0.0;
            else
                return atomic;
        }
    }

    public Double getValue(int index) {
        if (enumTable != null) {
            Double y = enumTable.getValue(index);
            if (y == null) {
                return 0.0;
            } else {
                return y;
            }
        } else {
            if (atomic == null)
                return 0.0;
            else
                return atomic;
        }
    }

    public List<Double> getValues(Object[] key) {
        if (enumTable != null) {
            return enumTable.getValues(key);
        } else {
            return Collections.EMPTY_LIST;
        }
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
        normalized = false;
        if (enumTable != null) {
            Double prev = enumTable.getValue(index);
            if (prev == null) {
                enumTable.setValue(index, value);
            } else {
                enumTable.setValue(index, prev + value);
            }
            return index;
        } else {
            atomic += value;
            return -1;
        }
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
        if (enumTable != null) {
            int index = enumTable.getIndex(key);
            return addValue(index, value);
        } else {
            atomic += value;
            return -1;
        }
    }

    public void normalize() {
        if (normalized)
            return;
        if (enumTable != null) {
            double sum = 0.0;
            for (Map.Entry<Integer, Double> entry : enumTable.getMapEntries()) {
                sum += entry.getValue();
            }
            for (Map.Entry<Integer, Double> entry : enumTable.getMapEntries()) {
                enumTable.setValue(entry.getKey(), entry.getValue() / sum);
            }
        }
        normalized = true;
    }
    
    public double getSum() {
        if (enumTable != null) {
            double sum = 0.0;
            for (Map.Entry<Integer, Double> entry : enumTable.getMapEntries()) 
                sum += entry.getValue();
            return sum;
        } else {
            return atomic;
        }
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
        if (enumTable != null) {
            int nParents = enumTable.nParents;
            Collection<EnumVariable> newparents = new ArrayList<>(nParents - parentsToSumOut.size());
            for (int i = 0; i < nParents; i++) {
                if (!parentsToSumOut.contains(enumTable.getParents().get(i))) {
                    newparents.add(enumTable.getParents().get(i));
                }
            }
            double sum = getSum();
            FactorTable ft; 
            if (nonEnumTables == null)
                ft = new FactorTable(newparents);
            else
                ft = new FactorTable(newparents, this.getNonEnumVariables());
            for (Map.Entry<Integer, Double> entry : enumTable.getMapEntries()) {
                int oldindex = entry.getKey();
                double weight = entry.getValue(); // NOT normalized
                int newindex = enumTable.maskIndex(oldindex, parentsToSumOut);
                ft.addValue(newindex, weight);
                for (Variable nonenum : this.getNonEnumVariables())
                    ft.addDistrib(newindex, nonenum, getDistrib(oldindex, nonenum), weight / sum); // HERE normalized weighting
            }
            ft.function = true;
            return ft;
        } else {
            return null;
        }
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
        if (ft1.enumTable == null && ft2.enumTable == null) {
            List<Variable> nonEnums = new ArrayList<>();
            nonEnums.addAll(ft1.getNonEnumVariables());
            nonEnums.addAll(ft2.getNonEnumVariables());
            FactorTable ft3 = new FactorTable(Collections.EMPTY_LIST, nonEnums);
            ft3.atomic = ft1.atomic * ft2.atomic;
            for (Variable ft1_nonenum : ft1.getNonEnumVariables())
                ft3.setDistrib(-1, ft1_nonenum, ft1.getDistrib(-1, ft1_nonenum));
            for (Variable ft2_nonenum : ft2.getNonEnumVariables())
                ft3.setDistrib(-1, ft2_nonenum, ft2.getDistrib(-1, ft2_nonenum));
            return ft3;
        } else if (ft1.enumTable == null) {
            List<Variable> nonEnums = new ArrayList<>();
            nonEnums.addAll(ft1.getNonEnumVariables());
            nonEnums.addAll(ft2.getNonEnumVariables());
            FactorTable ft3 = new FactorTable(ft2.enumTable.getParents(), nonEnums);
            for (Map.Entry<Integer, Double> entry2 : ft2.enumTable.getMapEntries()) {
                int ft3_index = ft3.enumTable.setValue(entry2.getKey(), entry2.getValue() * ft1.atomic);
                for (Variable ft1_nonenum : ft1.getNonEnumVariables())
                    ft3.setDistrib(ft3_index, ft1_nonenum, ft1.getDistrib(-1, ft1_nonenum));
                for (Variable ft2_nonenum : ft2.getNonEnumVariables())
                    ft3.setDistrib(ft3_index, ft2_nonenum, ft2.getDistrib(entry2.getKey(), ft2_nonenum));
            }
            ft3.function = true;
            return ft3;
        } else if (ft2.enumTable == null) {
            List<Variable> nonEnums = new ArrayList<>();
            nonEnums.addAll(ft1.getNonEnumVariables());
            nonEnums.addAll(ft2.getNonEnumVariables());
            FactorTable ft3 = new FactorTable(ft1.enumTable.getParents(), nonEnums);
            for (Map.Entry<Integer, Double> entry : ft1.enumTable.getMapEntries()) {
                int ft3_index = ft3.enumTable.setValue(entry.getKey(), entry.getValue() * ft2.atomic);
                for (Variable ft1_nonenum : ft1.getNonEnumVariables())
                    ft3.setDistrib(ft3_index, ft1_nonenum, ft1.getDistrib(entry.getKey(), ft1_nonenum));
                for (Variable ft2_nonenum : ft2.getNonEnumVariables())
                    ft3.setDistrib(ft3_index, ft2_nonenum, ft2.getDistrib(-1, ft2_nonenum));
            }
            ft3.function = true;
            return ft3;
        }
        Map<Integer, Integer> overlap = new HashMap<>();
        List<EnumVariable> unique_ft2 = new ArrayList<>(); // in ft2 but not in ft1
        List<Integer> unique_ft2_idx = new ArrayList<>(); // in ft2 but not in ft1
        for (int i = 0; i < ft2.enumTable.nParents; i++) {
            boolean match = false;
            for (int j = 0; j < ft1.enumTable.nParents; j++) {
                if (ft1.enumTable.getParents().get(j) == ft2.enumTable.getParents().get(i)) {
                    overlap.put(j, i); // link index in ft1 to index in ft2, resulting index in ft3 is the same as in ft1
                    match = true;
                    break;
                }
            }
            if (!match) {
                // ft2.parent[i] did not match any var in ft1
                unique_ft2.add(ft2.enumTable.getParents().get(i));
                unique_ft2_idx.add(i);
            }
        }
        Collection<EnumVariable> newparents = new ArrayList<>(ft1.enumTable.nParents);
        for (int i = 0; i < ft1.enumTable.nParents; i++) {
            newparents.add(ft1.enumTable.getParents().get(i));
        }
        newparents.addAll(unique_ft2);
        Collection<Variable> newNonEnumParents = new ArrayList<>(ft1.getNonEnumVariables());
        newNonEnumParents.addAll(ft2.getNonEnumVariables());
        FactorTable ft3 = null;
        if (newNonEnumParents.size() > 0)
            ft3 = new FactorTable(newparents, newNonEnumParents);
        else
            ft3 = new FactorTable(newparents);
        for (Map.Entry<Integer, Double> entry1 : ft1.enumTable.getMapEntries()) {
            if (entry1.getValue() == 0.0) // if the value associated with the entry is 0, the product will be zero, 
            {
                continue; // ... so we abort, to save time as there is no point in doing any calculations
            }
            for (Map.Entry<Integer, Double> entry2 : ft2.enumTable.getMapEntries()) {
                if (entry2.getValue() == 0.0) // if the value associated with the entry is 0, the product will be zero, 
                {
                    continue; // ... so we abort as there is no point in doing any more calculcations
                }
                int ft1_index = entry1.getKey();
                int ft2_index = entry2.getKey();
                Object[] ft1_key = ft1.enumTable.getKey(ft1_index);
                Object[] ft2_key = ft2.enumTable.getKey(ft2_index);
                boolean match = true;
                for (Map.Entry<Integer, Integer> link : overlap.entrySet()) {
                    if (ft1_key[link.getKey()] != ft2_key[link.getValue()]) {
                        match = false;
                        break;
                    }
                }
                if (match) { // all linked variables have the same value
                    Object[] ft3_key = new Object[ft3.enumTable.nParents];
                    System.arraycopy(ft1_key, 0, ft3_key, 0, ft1.enumTable.nParents);
                    int i = ft1.enumTable.nParents;
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

        /**
     * Change the order of defining parents.
     * @param ordered the variables in new order
     * @return 
     */
    private FactorTable rehash(List<Variable> ordered) {
        List<EnumVariable> f_parents = getEnumVariables();
        List<EnumVariable> q_parents = new ArrayList<>();
        List<Variable> c_parents = new ArrayList<>();
        for (Variable var : ordered) {
            try {
                q_parents.add((EnumVariable)var);
            } catch (ClassCastException e) {
                c_parents.add(var);
                ; // non-enumerable variables do not define indices in the resulting table
            }
        }
        int[] map2q = new int[q_parents.size()];
        for (int jf = 0; jf < map2q.length; jf ++) {
            map2q[jf] = -1;
            for (int jq = 0; jq < map2q.length; jq ++) {
                if (f_parents.get(jf).getName().equals(q_parents.get(jq).getName()))  {
                    map2q[jf] = jq;
                    break;
                }
            }
        }
        for (int j = 0; j < map2q.length; j ++) {
            if (map2q[j] == -1)
                throw new FactorTableRuntimeException("Invalid list of variables to extract");
        }
        FactorTable et = new FactorTable(q_parents, c_parents);

        if (this.hasNonEnumVariables()) {
            if (this.isAtomic()) {
                for (Variable nonenum : this.getNonEnumVariables())
                    et.setDistrib(-1, nonenum, this.getDistrib(-1, nonenum));
            } 
        }
        et.atomic = this.atomic;
        for (Map.Entry<Integer, Double> entry : this.getMapEntries()) {
            Object[] fkey = this.getKey(entry.getKey().intValue());
            Object[] qkey = new Object[fkey.length];
            for (int j = 0; j < fkey.length; j ++)
                qkey[map2q[j]] = fkey[j];
            et.setValue(qkey, entry.getValue());
            if (nonEnumTables != null) {
                for (Variable nonenum : this.getNonEnumVariables()) {
                    Distrib d = this.getDistrib(fkey, nonenum);
                    et.setDistrib(qkey, nonenum, d);
                }
            }
        }
        return et;
    }

    public double getLogLikelihood() {
        double sum = 0.0;
        if (enumTable != null) {
            int nParents = enumTable.nParents;
            for (Map.Entry<Integer, Double> entry : enumTable.getMapEntries()) {
                sum += entry.getValue().doubleValue();
            }
        } else
            sum = atomic;
        if (sum == 0.0)
            throw new RuntimeException("All outcomes are impossible: Log-likelihood is " + Math.log(sum));
        return Math.log(sum);
    }

    public static void main(String[] args) {
        EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Number(4);
        EnumVariable v3 = Predef.Nominal(new String[]{"Yes", "No", "Maybe"});

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
        System.out.println("Size of ft3: " + ft3.enumTable.getSize());
    }

    @Override
    public String toString() {
        StringBuilder sbuf = new StringBuilder("F(");
        for (Variable v : enumTable.getParents()) {
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
    
    public void display() {
        System.out.print("Idx ");
        for (int j = 0; j < enumTable.nParents; j++)
            System.out.print(String.format("[%10s]", constantLength(enumTable.getParents().get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" F");
        for (int i = 0; i < enumTable.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = enumTable.getKey(i);
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
        for (int j = 0; j < enumTable.nParents; j++)
            System.out.print(String.format("[%10s]", constantLength(enumTable.getParents().get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" F");
        for (int i = 0; i < enumTable.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = enumTable.getKey(i);
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
    
    public FactorTable integrateOver(Variable varToIntegrateOver, double start, double end, double delta) {
        throw new FactorTableRuntimeException("Not yet implemented");
    }

}

class FactorTableRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public FactorTableRuntimeException(String string) {
        message = string;
    }
}
