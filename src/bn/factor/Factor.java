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
import bn.JDF;
import bn.Predef;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A factor is a data structure for efficiently processing products and
 * sums over nodes--sometimes referred to as "local" computations. 
 * It can represent enumerable and non-enumerable variables but does so
 * by breaking them up in two parts, where instances of enumerable variables 
 * form keys, to which instances of non-enumerables are linked, 
 * reminiscent of so-called conditional Gaussian distributions 
 * (Lauritzen and Jensen, Statistics and Computing, 2001, 11:191-203).
 * 
 * Keys are defined by a size (of elements), and for each element, a valence (a number ]2,n[
 * where n>2). 
 * 
 * @author Mikael Boden
 */
public class Factor {

    private EnumTable<Double> factorTable = null;  // the factors for each permutation of the enumerable variables
    private EnumTable<JDF> densityTable = null;    // the densities for each permutation of the enumerables
    private List<EnumVariable> evars = null;       // all enumerable variables associated with this factor table
    private List<Variable> nvars = null;           // all non-enumerable variables
    
    private Double atomicFactor = null; // if there are no enumerable variables, here's where the constant factor
    private JDF atomicDensity = null;   // if there are no enumerable but some non-enumerable variables, here's where the only joint density is stored
    private boolean normalized = false; // flag to indicate if this factor table is normalized or not
    public boolean evidenced = false;   // flag to indicate if the Factor was reduced due to evidence (used by inference)
    public boolean function = false;    // flag to indicate if the factor was the result of a product or marginalization
    
    /**
     * CG_SAFETY forces factors to store entries that are associated with 0 (impossible) to save densities that are associated with them.
     * Set CG_SAFETY to false if time and space efficiency is critical.
     */
    public static boolean CG_SAFETY = true;    // save conditional densities, accept space and time overhead. Must be static/simulation-wide.
    
    /**
     * Construct a factor with no variables.
     */
    public Factor() {
        atomicFactor = new Double(1.0);
    }
    
    /**
     * Construct a factor from a collection of variables.
     * @param vars the variables that are included in the table (enumerable or not)
     */
    public Factor(List<Variable> vars) {
        if (vars == null)
            atomicFactor = new Double(1.0);
        else if (vars.isEmpty())
            atomicFactor = new Double(1.0);
        else { // proper factor with one or more variables
            // First, bin variables in lists for enumerables and non-enumerables, resp.
            evars = new ArrayList<>();
            nvars = new ArrayList<>();
            for (Variable var : vars) {
                try {
                    EnumVariable evar = (EnumVariable)var;
                    evars.add(evar);
                } catch (ClassCastException e) { // not enumerable
                    nvars.add(var);
                }
            }
            if (!evars.isEmpty()) { // there are enumerable variables
                factorTable = new EnumTable<>(evars);
                if (!nvars.isEmpty()) // there are also non-enumerables
                    densityTable = new EnumTable<>(evars);
            } else { // there are NO enumerables, but some non-enumerables
                atomicDensity = new JDF(nvars);
                atomicFactor = new Double(1.0);
            }
        }
    }

    /**
     * Construct a factor from an array of variables.
     * @param vars the variables that are included in the table (enumerable or not)
     */
    public Factor(Variable... vars) {
        if (vars == null)
            atomicFactor = new Double(1.0);
        else if (vars.length == 0)
            atomicFactor = new Double(1.0);
        else { // proper factor with one or more variables
            // First, bin variables in lists for enumerables and non-enumerables, resp.
            evars = new ArrayList<>();
            nvars = new ArrayList<>();
            for (Variable var : vars) {
                try {
                    EnumVariable evar = (EnumVariable)var;
                    evars.add(evar);
                } catch (ClassCastException e) { // not enumerable
                    nvars.add(var);
                }
            }
            if (!evars.isEmpty()) { // there are enumerable variables
                factorTable = new EnumTable<>(evars);
                if (!nvars.isEmpty()) // there are also non-enumerables
                    densityTable = new EnumTable<>(evars);
            } else { // there are NO enumerables, but some non-enumerables
                atomicDensity = new JDF(nvars);
            }
        }
    }

    /**
     * Find if factor indexes enumerable variables.
     * @return true if factor is associated with enumerable variables
     */
    public boolean hasEnumVariables() {
        if (evars == null)
            return false;
        else 
            return (!evars.isEmpty());
    }
    
    /**
     * Find if factor stores distributions of non-enumerable variables.
     * @return true if factor is associated with non-enumerable variables
     */
    public boolean hasNonEnumVariables() {
        if (nvars == null)
            return false;
        else 
            return (!nvars.isEmpty());
    }
    
    /**
     * Set the value in a factor. 
     * Use with care, as this does not automatically checks the integrity of table.
     * @param key the index key identifying the setting of the enumerable variables
     * @param value the value to assign the entry
     * @return the index of the entry
     */
    public int setFactor(int key, Double value) {
        normalized = false;
        if (factorTable != null)
            return factorTable.setValue(key, value);
        else {
            atomicFactor = value;
            return -1; // no index since there are no enum variables
        }
    }

    /**
     * Set the value in a factor. 
     * Use with care, as this does not automatically checks the integrity of factor table.
     * @param key the instance key identifying the setting of the enumerable variables
     * @param value the value to assign the entry
     * @return the index of the entry
     */
    public int setFactor(Object[] key, Double value) {
        normalized = false;
        if (factorTable != null)
            return factorTable.setValue(key, value);
        else {
            atomicFactor = value;
            return -1; // no index since there are no enum variables
        }
    }

    /**
     * Set the value in a factor without any enumerable variables.
     * @param value the value to assign the entry
     * @return the index of the entry, always -1
     * @exception FactorRuntimeException if the factor has enumerable variables
     */
    public int setFactor(Double value) {
        if (hasEnumVariables())
            throw new FactorRuntimeException("Cannot set factor: " + this);
        normalized = false;
        atomicFactor = value;
        return -1; // no index since there are no enum variables
    }

    /**
     * Set the Joint Density Function for an entry.
     * @param key_index the condition as encoded by a key index
     * @param jdf the JDF
     * @return the index at which the JDF was placed, -1 if atomic
     */
    public int setJDF(int key_index, JDF jdf) {
        if (densityTable != null)
            return densityTable.setValue(key_index, jdf);
        else {
            atomicDensity = jdf;
            return -1; // no index since there are no enum variables
        }
    }
    
    /**
     * Set the Joint Density Function for an entry.
     * @param key the condition as encoded by a key
     * @param jdf the JDF
     * @return the index at which the JDF was placed, -1 if atomic
     */
    public int setJDF(Object[] key, JDF jdf) {
        if (densityTable != null)
            return densityTable.setValue(key, jdf);
        else {
            atomicDensity = jdf;
            return -1; // no index since there are no enum variables
        }
    }
    
    /**
     * Set the Joint Density Function a factor without enumerable variables.
     * @param jdf the JDF
     * @return the index at which the JDF was placed, -1 if atomic
     */
    public int setJDF(JDF jdf) {
        if (hasEnumVariables())
            throw new FactorRuntimeException("Cannot set JDF in factor: " + this);
        atomicDensity = jdf;
        return -1; // no index since there are no enum variables
    }

    /**
     * Get the JDF for an entry.
     * @param key_index the entry as identified by the condition encoded by a key index
     * @return the JDF, null if not available/inapplicable
     */
    public JDF getJDF(int key_index) {
        if (densityTable != null)
            return densityTable.getValue(key_index);
        return null;
    }
    
    /**
     * Get the JDF for an entry.
     * @param key the entry as identified by the condition encoded by a key
     * @return the JDF
     */
    public JDF getJDF(Object[] key) {
        if (densityTable != null)
            return densityTable.getValue(key);
        return null;
    }
    
    /**
     * Get the JDF for a factor without enumerable variables.
     * @return the JDF
     */
    public JDF getJDF() {
        if (atomicDensity != null)
            return atomicDensity;
        return null;
    }
    
    /**
     * Set the conditional distribution of a non-enumerable variable. 
     * Use with care as the method may not do all integrity checks.
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(Object[] key, Variable nonenum, Distrib d) {
        if (densityTable != null) {
            int index = densityTable.getIndex(key);
            return setDistrib(index, nonenum, d);
        } else if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            atomicDensity.setDistrib(d, nonenum);
            return -1; // no index since there are no enum variables
        }
        throw new FactorRuntimeException("Cannot set distribution of factor: " + this);
    }

    /**
     * Set the conditional distribution of a non-enumerable variable. 
     * Use with care as the method may not do all integrity checks.
     * @param key_index the key index identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry
     */
    public int setDistrib(int key_index, Variable nonenum, Distrib d) {
        if (densityTable != null) {
            JDF current = densityTable.getValue(key_index);
            if (current == null) {
                JDF replacement = new JDF(nvars);
                replacement.setDistrib(d, nonenum);
                densityTable.setValue(key_index, replacement);
            } else {
                current.setDistrib(d, nonenum);
            }
            return key_index;
        } else if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            atomicDensity.setDistrib(d, nonenum);
            return -1; // no index since there are no enum variables
        }
        throw new FactorRuntimeException("Cannot set distribution of factor: " + this);
    }

    /**
     * Set the distribution of a non-enumerable variable for a factor without enumerable variables. 
     * Use with care as the method may not do all integrity checks.
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @return the index of the entry, always -1
     */
    public int setDistrib(Variable nonenum, Distrib d) {
        if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            atomicDensity.setDistrib(d, nonenum);
            return -1; // no index since there are no enum variables
        }
        throw new FactorRuntimeException("Cannot set distribution of factor: " + this);
    }

    /**
     * Add a non-enumerable distribution as required by marginalization (removal) of enumerable variables.
     * @param key_index the key index identifying the setting of the enumerable variables ("condition")
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @param weight the weight assigned to the distribution
     * @return the index for the key
     */
    public int addDistrib(int key_index, Variable nonenum, Distrib d, double weight) {
        if (densityTable != null) {
            JDF current = densityTable.getValue(key_index);
            if (current == null) { // no previous density for this "condition", so we create a new one
                JDF replacement = new JDF(nvars);
                replacement.mixDistrib(nonenum, d, weight); 
                densityTable.setValue(key_index, replacement);
            } else {
                current.mixDistrib(nonenum, d, weight);
            }
            return key_index;
        } else if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            atomicDensity.mixDistrib(nonenum, d, weight);
            return -1; // no index since there are no enum variables
        }
        throw new FactorRuntimeException("Cannot add distribution of factor: " + this);
    }

    /**
     * Add a non-enumerable distribution as required by marginalization (removal) of enumerable variables.
     * @param key the key identifying the setting of the enumerable variables
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @param weight the weight assigned to the distribution
     * @return the index for the key
     */
    public int addDistrib(Object[] key, Variable nonenum, Distrib d, double weight) {
        if (densityTable != null) {
            int index = densityTable.getIndex(key);
            return addDistrib(index, nonenum, d, weight);
        } else if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            atomicDensity.mixDistrib(nonenum, d, weight);
            return -1; // no index since there are no enum variables
        }
        throw new FactorRuntimeException("Cannot add distribution of factor: " + this);
    }

    /**
     * Add a non-enumerable distribution as required by marginalization (removal) of enumerable variables.
     * @param nonenum the non-enumerable variable to be set
     * @param d the distribution that applies to the non-enumerable variable
     * @param weight the weight assigned to the distribution
     * @return the index for the key
     */
    public int addDistrib(Variable nonenum, Distrib d, double weight) {
        if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            atomicDensity.mixDistrib(nonenum, d, weight);
            return -1; // no index since there are no enum variables
        }
        throw new FactorRuntimeException("Cannot add distribution of factor: " + this);
    }

    /**
     * Get a conditional distribution for a non-enumerable variable.
     * @param key the key representing the instantiation of enumerable variables
     * @param nonenum the non-enumerable variable
     * @return the distribution, null if unavailable
     */
    public Distrib getDistrib(Object[] key, Variable nonenum) {
        if (densityTable != null) {
            int index = densityTable.getIndex(key);
            return getDistrib(index, nonenum);
        } else  if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            return atomicDensity.getDistrib(nonenum);
        }
        return null;
    }
    
    /**
     * Get a conditional distribution for a non-enumerable variable.
     * @param key_index the index key representing the instantiation of enumerable variables
     * @param nonenum the non-enumerable variable
     * @return the distribution, null if unavailable
     */
    public Distrib getDistrib(int key_index, Variable nonenum) {
        if (densityTable != null) {
            JDF current = densityTable.getValue(key_index);
            if (current != null) {
                return current.getDistrib(nonenum);
            } 
        } else if (atomicDensity != null) { // no enumerable, but some non-enumerable variables
            return atomicDensity.getDistrib(nonenum);
        }
        return null;
    }
    
    /**
     * Get a distribution for a non-enumerable variable in a factor without enumerable variables.
     * @param nonenum the non-enumerable variable
     * @return the distribution, null if unavailable
     */
    public Distrib getDistrib(Variable nonenum) {
        if (atomicDensity != null) // no enumerable, but some non-enumerable variables
            return atomicDensity.getDistrib(nonenum);
        return null;
    }
    
    /**
     * Get the number of enumerable variables in factor.
     * @return the number
     */
    public int getNEnum() {
        if (factorTable == null)
            return 0;
        return factorTable.nParents;
    }
    
    /**
     * Get the number of non-enumerable variables in factor.
     * @return the number
     */
    public int getNNonEnum() {
        if (nvars == null)
            return 0;
        return nvars.size();
    }
    
    /**
     * Retrieve the list of enumerable variables
     * @return the list of enumerable variables associated with this factor, empty list if none
     */
    public List<EnumVariable> getEnumVariables() {
        if (evars == null)
            return Collections.EMPTY_LIST;
        return evars;
    }
    
    /**
     * Retrieve the list of non-enumerable variables
     * @return the list of non-enumerable variables associated with this factor, empty list if none
     */
    public List<Variable> getNonEnumVariables() {
        if (nvars == null)
            return Collections.EMPTY_LIST;
        return nvars;
    }
    
    /**
     * Retrieve the entries for manual iteration of factor table entries. 
     * This list contains all instantiations of the enumerable variables.
     * @return the entries of the factor table
     */
    public Set<Map.Entry<Integer, Double>> getMapEntries() {
    	if (factorTable == null){
            return Collections.EMPTY_SET;
    	} else {
            return factorTable.getMapEntries();
    	}
    }
    
    /**
     * Decode an index to an array of instantiations of enumerable variables.
     * The order of values is the same as that when enumerable variables are retrieved.
     * @param index
     * @return 
     */
    public Object[] getKey(int index) {
        return factorTable.getKey(index);
    }
    
    /**
     * Decode a series of assignments of enumerable variables into a key.
     * @param assignment assignment of values to enumerable variables
     * @return the key which may contain null values
     */
    public Object[] getKey(Variable.Assignment[] assignment) {
        return factorTable.getKey(assignment);
    }
    
    /**
     * Get the indices for entries that match the provided, potentially incomplete key 
     * (instantiation of variables).
     * @param key the key
     * @return the indices
     */
    public int[] getIndices(Object[] key) {
        return factorTable.getIndices(key);
    }
    
    /**
     * Check if the factor has no association with enumerable variables.
     * @return true if there are no enumerable variables indexed by this factor, false otherwise.
     */
    public boolean isAtomic() {
        return factorTable == null;
    }
    
    /**
     * Retrieve all factor values. 
     * Note that order of values is undefined.
     * @return all factor values
     **/
    public Collection<Double> getFactors() {
        if (factorTable != null) {
            return factorTable.getValues();
        } else {
            return Collections.singletonList(atomicFactor);
        }
    }
    
    /**
     * Retrieve all factor values that match a key. 
     * Note that order of returned values is undefined.
     * @param key the instance key, with null representing unspecified
     * @return matching factor values
     **/
    public List<Double> getFactors(Object[] key) {
        if (factorTable != null) {
            return factorTable.getValues(key);
        } else {
            return Collections.EMPTY_LIST;
        }
    }

    /**
     * Retrieve a factor value that match a fully specified key. 
     * @param key the instance key
     * @return factor value for key, 0 if not found
     **/
    public Double getFactor(Object[] key) {
        if (factorTable != null) {
            Double y = factorTable.getValue(key);
            if (y == null) {
                return 0.0;
            } else {
                return y;
            }
        } else if (atomicFactor != null) {
            return atomicFactor;
        }
        throw new FactorRuntimeException("Operation on invalid factor: " + this);
    }

    /**
     * Retrieve a factor value that match a fully specified key. 
     * @param index the index that specifies the instance key
     * @return factor value for key, 0 if not found
     **/
    public Double getFactor(int index) {
        if (factorTable != null) {
            Double y = factorTable.getValue(index);
            if (y == null) {
                return 0.0;
            } else {
                return y;
            }
        } else if (atomicFactor != null) {
            return atomicFactor;
        }
        throw new FactorRuntimeException("Operation on invalid factor: " + this);
    }

    /**
     * Retrieve a factor value from a factor without enumerable variables. 
     * @return factor value for atomic factor, 0 if not found
     **/
    public Double getFactor() {
        if (atomicFactor != null)
            return atomicFactor;
        throw new FactorRuntimeException("Operation on invalid factor: " + this);
    }

    /**
     * Add a value to a specified entry identified by index
     * @param index the entry index or -1 for setting an atomic factor
     * @param value the value to be mixed to the current instance (assumed to be
 0 if not instantiated)
     * @return the index
     */
    synchronized public int addFactor(int index, double value) {
        normalized = false;
        if (factorTable != null) {
            Double prev = factorTable.getValue(index);
            if (prev == null) {
                factorTable.setValue(index, value);
            } else {
                factorTable.setValue(index, prev + value);
            }
            return index;
        } else if (atomicFactor != null) {
            atomicFactor += value;
            return -1;
        }
        throw new FactorRuntimeException("Operation on invalid factor: ");
    }

    /**
     * Add a value to a specified entry identified by instantiated key
     * @param key the entry key
     * @param value the value to be mixed to the current instance (assumed to be
 0 if not instantiated)
     * @return the index for the modified entry
     */
    public int addFactor(Object[] key, double value) {
        if (factorTable != null) {
            int index = factorTable.getIndex(key);
            return addFactor(index, value);
        } else if (atomicFactor != null) {
            atomicFactor += value;
            return -1;
        }
        throw new FactorRuntimeException("Operation on invalid factor: " + this);
    }

    /**
     * Add a value to an atomic factor, i.e. one with no enumerable variables
     * @param value the value to be mixed
     * @return the index for the modified entry, always -1
     */
    public int addFactor(double value) {
        if (atomicFactor != null) {
            atomicFactor += value;
            return -1;
        }
        throw new FactorRuntimeException("Operation on invalid factor: " + this);
    }

    /**
     * Normalise the factors so they mix up to 1.0.
     * Use with caution as sometimes the raw factors are more useful.
     */
    public void normalize() {
        if (normalized)
            return;
        if (factorTable != null) {
            double sum = 0.0;
            for (Map.Entry<Integer, Double> entry : factorTable.getMapEntries()) {
                sum += entry.getValue();
            }
            for (Map.Entry<Integer, Double> entry : factorTable.getMapEntries()) {
                factorTable.setValue(entry.getKey(), entry.getValue() / sum);
            }
        }
        normalized = true;
    }
    
    /**
     * Calculate the sum of factor values.
     * @return the sum
     */
    public double getSum() {
        if (factorTable != null) {
            double sum = 0.0;
            for (Map.Entry<Integer, Double> entry : factorTable.getMapEntries()) 
                sum += entry.getValue().doubleValue();
            return sum;
        } else {
            return atomicFactor.doubleValue();
        }
    }

    /**
     * Calculate the sum of factor values.
     * @param key the representation of the instantiation of enumerable variables, 
     * could be incomplete (null representing unspecified)
     * @return the sum
     */
    public double getSum(Object[] key) {
        if (factorTable != null) {
            return getSum(factorTable.getIndices(key));
        } else {
            return atomicFactor.doubleValue();
        }
    }
    
    /**
     * Calculate the sum of factor values.
     * @param indices the specific indices of the table
     * @return the sum
     */
    public double getSum(int[] indices) {
        if (factorTable != null) {
            double sum = 0.0;
            for (int index : indices) 
                sum += factorTable.getValue(index);
            return sum;
        } else {
            return atomicFactor.doubleValue();
        }
    }
    
    /**
     * Determine the mixture of distributions for all entries that match key
     * @param key the specific entries as identified by a (partial) key
     * @return the joint mixture distribution
     */
    public JDF getJDFSum(Object[] key) {
        if (densityTable != null) {
            return getJDFSum(densityTable.getIndices(key));
        } else {
            return atomicDensity;
        }
    }
    
    /**
     * Determine the mixture of distributions for all entries that match key
     * @param indices the specific entries as identified by their indices
     * @return the joint mixture distribution
     */
    public JDF getJDFSum(int[] indices) {
        double sum = getSum(indices);
        if (densityTable != null) {
            JDF jdf = new JDF(getNonEnumVariables());
            for (int index : indices) {
                JDF current = densityTable.getValue(index);
                double weight = factorTable.getValue(index) / sum;
                jdf = JDF.mix(jdf, current, weight); 
            }
            return jdf;
        } else {
            return atomicDensity;
        }
    }
    
    /**
     * Sums-out (marginalizes) one or more enumerable variables from the current table. Time
     * complexity is O(n+m*3n) where n is the number of parents and m is the
     * specified number of entries in the current table.
     * 
     * This method has been partially adapted to accommodate continuous (non-enumerable) variables
     * that are associated with summed-out variables:
     * densities are combined into mixtures, which need to be resolved later.
     *
     * @param varsToSumOut variables that will be summed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public Factor marginalize(Collection<EnumVariable> varsToSumOut) {
        if (factorTable != null) {
            int nVars = factorTable.nParents;
            List<Variable> newvars = new ArrayList<>(nVars - varsToSumOut.size());
            for (int i = 0; i < nVars; i++) {
                if (!varsToSumOut.contains(factorTable.getParents().get(i))) {
                    newvars.add(factorTable.getParents().get(i));
                }
            }
            double sum = 0; // sum for normalising weight of distributions for non-enumerables
            Factor ft; 
            if (hasNonEnumVariables()) {
                newvars.addAll(this.getNonEnumVariables());
                sum = getSum(); // only need to compute this if we have non-enumerables, used below...
            }
            ft = new Factor(newvars);
            for (Map.Entry<Integer, Double> entry : this.getMapEntries()) {
                int oldindex = entry.getKey();
                double weight = entry.getValue(); // NOT normalized
                int newindex = factorTable.maskIndex(oldindex, varsToSumOut);
                ft.addFactor(newindex, weight);
                for (Variable nonenum : this.getNonEnumVariables()) {
                    double myWeight = 1.0; // if all factors are zero, sum is 0, so we apply a constant, uniform weight 
                    if (sum != 0) // not all factors are zero
                        myWeight = weight / sum; // HERE normalized weighting
                    if (myWeight != 0.0) // do not add densities with zero weight
                        ft.addDistrib(newindex, nonenum, getDistrib(oldindex, nonenum), myWeight); 
                }
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
     * @param varsToSumOut variables that will be summed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public Factor marginalize(EnumVariable[] varsToSumOut) {
        return this.marginalize(EnumVariable.toList(varsToSumOut));
    }

    /** 
     * Max-out one or more enumerable variables from the current table.
     * This is a way to remove variables from the table in a way that tracks the maximally
     * probably explanation of the evidence.
     * 
     * TODO: This method has not yet been tested extensively, esp. with non-enums. 
     * TODO: Variables that are max:ed out can be tracked so as to provide all assignments that form
     * part of the assignments (and not only those part of the query). This tracking probably needs
     * to be for every entry in the factor, and must be transferred and aggregate across products and 
     * sum-outs too. Since this may consume memory, it may need to be an option, or turned on for MPE 
     * only scenarios.
     * 
     * @param varsToMaxOut variables that will be max:ed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public Factor maximize(Collection<EnumVariable> varsToMaxOut) {
        if (factorTable != null) {
            int nVars = factorTable.nParents;
            List<Variable> newvars = new ArrayList<>(nVars - varsToMaxOut.size());
            for (int i = 0; i < nVars; i++) {
                if (!varsToMaxOut.contains(factorTable.getParents().get(i))) {
                    newvars.add(factorTable.getParents().get(i));
                }
            }
            double sum = 0; // sum for normalising weight of distributions for non-enumerables
            Factor ft; 
            if (hasNonEnumVariables()) {
                newvars.addAll(this.getNonEnumVariables());
                sum = getSum(); // only need to compute this if we have non-enumerables, used below...
            }
            ft = new Factor(newvars);
            Map<Integer, Integer> indexMap = new HashMap<>();
            for (Map.Entry<Integer, Double> entry : this.getMapEntries()) {
                int oldindex = entry.getKey();
                double weight = entry.getValue(); // NOT normalized
                int newindex = factorTable.maskIndex(oldindex, varsToMaxOut);
                Integer maxindex = indexMap.get(newindex);
                if (maxindex == null) // nothing to compare with
                    indexMap.put(newindex, oldindex); 
                else {
                    double oldweight = this.getFactor(maxindex);
                    if (weight > oldweight)
                        indexMap.put(newindex, oldindex);
                }
            }
            for (Map.Entry<Integer, Integer> entry : indexMap.entrySet()) {
                int newindex = entry.getKey();
                int oldindex = entry.getValue();
                double weight = this.getFactor(oldindex);
                ft.setFactor(newindex, weight);
                for (Variable nonenum : this.getNonEnumVariables()) {
                    double myWeight = 1.0; // if all factors are zero, sum is 0, so we apply a constant, uniform weight 
                    if (sum != 0) // not all factors are zero
                        myWeight = weight / sum; // HERE semi-normalized weighting (disregarding max:ed out entries)
                    if (myWeight != 0.0) // do not add densities with zero weight
                        ft.addDistrib(newindex, nonenum, getDistrib(oldindex, nonenum), myWeight); 
                }
            }
            ft.function = true;
            return ft;
        } else {
            return null;
        }
    }
    
    /** 
     * Max-out one or more enumerable variables from the current table.
     * This is a way to remove variables from the table in a way that tracks the maximally
     * probably explanation of the evidence.
     * 
     * @param varsToMaxOut variables that will be max:ed out from the current
     * table
     * @return the new FactorTable containing the remaining variables
     */
    public Factor maximize(EnumVariable[] varsToMaxOut) {
        return this.maximize(EnumVariable.toList(varsToMaxOut));
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

    /**
     * Takes the product between the two factors provided as arguments.
     * Variables to be "joined" from the factors must have the same name and
     * implement the same domain for the product to be determined.
     * This method has been partially adapted for continuous variables, but it assumes
     * that continuous variables cannot be reached by multiple paths. For that
     * to work, they have to be joined the same way as discrete variables. Currently,
     * continuous variables are just mixed without checking if the other factors has them.
     *
     * @param ft1 factor one
     * @param ft2 factor two
     * @return the factor table that is the product of the specified factors
     */
    public static Factor product(Factor ft1, Factor ft2) {
        List<Variable> ft3vars = new ArrayList<>();
        if (ft1.factorTable == null && ft2.factorTable == null) {
            ft3vars.addAll(ft1.getNonEnumVariables());
            ft3vars.addAll(ft2.getNonEnumVariables());
            Factor ft3 = new Factor(ft3vars);
            ft3.atomicFactor = ft1.atomicFactor * ft2.atomicFactor;
            if (ft1.hasNonEnumVariables() || ft2.hasNonEnumVariables())
                ft3.setJDF(JDF.combine(ft1.getJDF(), ft2.getJDF()));
            return ft3;
        } else if (ft1.factorTable == null) {
            ft3vars.addAll(ft2.getEnumVariables());
            ft3vars.addAll(ft1.getNonEnumVariables());
            ft3vars.addAll(ft2.getNonEnumVariables());
            Factor ft3 = new Factor(ft3vars);
            for (Map.Entry<Integer, Double> entry2 : ft2.getMapEntries()) {
                int ft2_index = entry2.getKey();
                int ft3_index = ft3.factorTable.setValue(entry2.getKey(), entry2.getValue() * ft1.atomicFactor);
                if (ft2.hasNonEnumVariables())
                    ft3.setJDF(ft3_index, ft2.getJDF(ft2_index));
            }
            ft3.function = true;
            return ft3;
        } else if (ft2.factorTable == null) {
            ft3vars.addAll(ft1.getEnumVariables());
            ft3vars.addAll(ft1.getNonEnumVariables());
            ft3vars.addAll(ft2.getNonEnumVariables());
            Factor ft3 = new Factor(ft3vars);
            for (Map.Entry<Integer, Double> entry1 : ft1.factorTable.getMapEntries()) {
                int ft1_index = entry1.getKey();
                int ft3_index = ft3.factorTable.setValue(entry1.getKey(), entry1.getValue() * ft2.atomicFactor);
                if (ft1.hasNonEnumVariables())
                    ft3.setJDF(ft3_index, ft1.getJDF(ft1_index));
            }
            ft3.function = true;
            return ft3;
        }
        // Since both ft1 and ft2 have enumerable variables, we need to work out their overlap
        // Note that we do not check non-enumerables since currently they can only occur as children.
        // This restriction means that they can only be members of one factor table, hence two FTs
        // will never share non-enumerable variables.
        Map<Integer, Integer> overlap = new HashMap<>();
        List<EnumVariable> unique_ft2 = new ArrayList<>(); // in ft2 but not in ft1
        List<Integer> unique_ft2_idx = new ArrayList<>();  // in ft2 but not in ft1
        for (int i = 0; i < ft2.factorTable.nParents; i++) {
            boolean match = false;
            for (int j = 0; j < ft1.factorTable.nParents; j++) {
                if (ft1.factorTable.getParents().get(j) == ft2.factorTable.getParents().get(i)) {
                    overlap.put(j, i); // link index in ft1 to index in ft2, resulting index in ft3 is the same as in ft1
                    match = true;
                    break;
                }
            }
            if (!match) {
                // ft2.parent[i] did not match any var in ft1
                unique_ft2.add(ft2.factorTable.getParents().get(i));
                unique_ft2_idx.add(i);
            }
        }
        // now we know the overlap so can define the result FT
        ft3vars.addAll(ft1.getEnumVariables());
        ft3vars.addAll(unique_ft2);
        List<Variable> ft1_nonenum = ft1.getNonEnumVariables();
        ft3vars.addAll(ft1_nonenum);
        List<Variable> ft2_nonenum = ft2.getNonEnumVariables();
        ft3vars.addAll(ft2_nonenum); // as explained above, FT1 and FT2 do not overlap in terms of non-enums
        Factor ft3 = new Factor(ft3vars);
        boolean CG_SAFETY_REQ = /* (ft1_nonenum.size() + ft2_nonenum.size() > 0) && */ Factor.CG_SAFETY;
//        if (CG_SAFETY_REQ) {
//            // TODO?
//        }
        for (Map.Entry<Integer, Double> entry1 : ft1.factorTable.getMapEntries()) {
            double f1 = entry1.getValue();
            if (f1 == 0.0 && !CG_SAFETY_REQ) // if the value associated with the entry is 0, the product will always be zero, 
                continue; // ... so without the CG_SAFETY active, we abort, to save time as there is no point in doing any calculations
            for (Map.Entry<Integer, Double> entry2 : ft2.factorTable.getMapEntries()) {
                double f2 = entry2.getValue();
                if ((f2 == 0.0 || f1 == 0.0) && !CG_SAFETY_REQ) // if the value associated with the entry is 0, the product will be zero, 
                    continue; // unless CG safety is active, there is no point in doing any more calculcations
                int ft1_index = entry1.getKey();
                int ft2_index = entry2.getKey();
                Object[] ft1_key = ft1.getKey(ft1_index);
                Object[] ft2_key = ft2.getKey(ft2_index);
                boolean match = true;
                // check that all "linked" variables between FT1 and FT2 match
                for (Map.Entry<Integer, Integer> link : overlap.entrySet()) {
                    if (ft1_key[link.getKey()] != ft2_key[link.getValue()]) {
                        match = false;
                        break;
                    }
                }
                if (match) { // all linked variables have indeed the same value
                    Object[] ft3_key = new Object[ft3.factorTable.nParents];
                    System.arraycopy(ft1_key, 0, ft3_key, 0, ft1.factorTable.nParents);
                    int i = ft1.factorTable.nParents;
                    for (int ft2_idx : unique_ft2_idx) {
                        ft3_key[i] = ft2_key[ft2_idx];
                        i++;
                    }
//                    if ((f2 == 0.0 || f1 == 0.0) && CG_SAFETY_REQ) { // if the value associated with the entry is 0, the product will be zero, but CG safety is active
//                        // TODO?
//                    } else {
                    int ft3_index = ft3.addFactor(ft3_key, entry1.getValue() * entry2.getValue());
                    if (ft1.hasNonEnumVariables() || ft2.hasNonEnumVariables())
                        ft3.setJDF(ft3_index, JDF.combine(ft1.getJDF(ft1_index), ft2.getJDF(ft2_index)));
//                    }
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
     * 
     * TODO Needs rewriting to use JDF etc
     */
    public Factor rehash(List<Variable> ordered) {
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
        // construct a map from current to new parent index in factor
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
                throw new FactorRuntimeException("Invalid list of variables to extract");
        }
        Factor et = new Factor(ordered);
        et.atomicDensity = this.atomicDensity;
        et.atomicFactor = this.atomicFactor;
        for (Map.Entry<Integer, Double> entry : this.getMapEntries()) {
            int fkey_index = entry.getKey().intValue();
            Object[] fkey = this.getKey(fkey_index);
            Object[] qkey = new Object[fkey.length];
            for (int j = 0; j < fkey.length; j ++)
                qkey[map2q[j]] = fkey[j];
            int qkey_index = et.setFactor(qkey, entry.getValue());
            if (hasNonEnumVariables()) 
                et.setJDF(qkey_index, this.getJDF(fkey_index));
        }
        return et;
    }

    public static void main(String[] args) {
        EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Number(4);
        EnumVariable v3 = Predef.Nominal(new String[]{"Yes", "No", "Maybe"});

        Factor ft3 = new Factor(new EnumVariable[]{v1, v2, v3});
        ft3.setFactor(new Object[]{true, 1, "No"}, 0.05);
        ft3.setFactor(new Object[]{true, 0, "No"}, 0.02);
        ft3.setFactor(new Object[]{false, 3, "No"}, 0.17);
        ft3.setFactor(new Object[]{false, 0, "Yes"}, 0.07);
        ft3.setFactor(new Object[]{false, 1, "Yes"}, 0.04);
        ft3.setFactor(new Object[]{true, 2, "Yes"}, 0.13);
        ft3.setFactor(new Object[]{true, 0, "Yes"}, 0.18);
        ft3.setFactor(new Object[]{true, 3, "Yes"}, 0.08);
        ft3.display();
        Factor ft2 = ft3.marginalize(new EnumVariable[]{v2});
        ft2.display();
        Factor ft1 = ft3.marginalize(new EnumVariable[]{v1, v3});
        ft1.display();
        System.out.println("Size of ft3: " + ft3.factorTable.getSize());
    }

    @Override
    public String toString() {
        StringBuilder sbuf = new StringBuilder("F(");
        for (Variable v : factorTable.getParents()) {
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
        for (int j = 0; j < factorTable.nParents; j++)
            System.out.print(String.format("[%10s]", constantLength(factorTable.getParents().get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" F");
        for (int i = 0; i < factorTable.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = factorTable.getKey(i);
            for (Object key1 : key)
                System.out.print(String.format(" %-10s ", constantLength(key1.toString(), 10)));
            for (Variable nonenum : nonenums) {
                Distrib d = this.getDistrib(i, nonenum);
                System.out.print(String.format(" %s ", d.toString()));
                //System.out.print(String.format(" %-10s ", constantLength(d.toString(), 10)));
            }
            Object val = this.getFactor(i);
            if (val != null) 
                System.out.println(String.format(" %7.5f", this.getFactor(i)));
            else
                System.out.println(" null ");
        }
    }

    public void displaySampled() {
        System.out.print("Idx ");
        for (int j = 0; j < factorTable.nParents; j++)
            System.out.print(String.format("[%10s]", constantLength(factorTable.getParents().get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" F");
        for (int i = 0; i < factorTable.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = factorTable.getKey(i);
            for (Object key1 : key)
                System.out.print(String.format(" %-10s ", constantLength(key1.toString(), 10)));
            JDF jdf = this.getJDF(i);
            for (Variable nonenum : nonenums) {
                double sum = 0;
                for (int j = 0; j < 1000; j ++) {
                    Double sample = (Double)jdf.sample(nonenum);
                    sum += sample.doubleValue();
                }
                System.out.print(String.format(" %5.3f ", sum/1000.0));
            }
            Object val = this.getFactor(i);
            if (val != null) 
                System.out.println(String.format(" %7.5f", this.getFactor(i)));
            else
                System.out.println(" null ");
        }
    }
    
    public Factor integrateOver(Variable varToIntegrateOver, double start, double end, double delta) {
        throw new FactorRuntimeException("Not yet implemented");
    }

}

class FactorRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public FactorRuntimeException(String string) {
        message = string;
    }
}
