/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bn.alg;

import bn.Distrib;
import bn.EnumDistrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.Factor;
import bn.JDF;
import bn.JPT;
import bn.MixtureDistrib;
import bn.Variable;
import bn.alg.QueryResult;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Conditional Gaussian Table.
 * 
 * Data structure to hold information about a select set of variables,
 * both enumerable and non-enumerable, to allow single variable distributions 
 * to be pulled out, based on no or some observations.
 * 
 * The data structure is used to store the result of a query over a complete
 * Bayesian network. It is largely based on the structure of a factor, but is
 * properly normalised and immutable.
 * 
 * Developers of Bayesian network inference engines should use this structure 
 * and potentially add a new constructor to suit the inference algorithm.
 * 
 * @author mikael
 */
public class CGTable implements QueryResult {
   
    private final EnumTable<Double> factorTable;  // the factors for each permutation of the enumerable variables
    private final EnumTable<JDF> densityTable;    // the densities for each permutation of the enumerables
    private final List<EnumVariable> evars;       // all enumerable variables associated with this factor table
    private final List<Variable> nvars;           // all non-enumerable variables
    
    private final Double atomicFactor; // if there are no enumerable variables, here's where the constant factor
    private final JDF    atomicDensity;   // if there are no enumerable but some non-enumerable variables, here's where the only joint density is stored
    
    public CGTable(Factor f) {
        evars = f.getEnumVariables();
        nvars = f.getNonEnumVariables();
        boolean hasNonEnum = f.hasNonEnumVariables();
        if (f.hasEnumVariables()) {
            atomicDensity = null;
            if (hasNonEnum) {
                densityTable = new EnumTable<>(evars);
            } else
                densityTable = null;
            factorTable = new EnumTable(evars);
            atomicFactor = null;
            double sum = f.getSum();
            for (Map.Entry<Integer, Double> entry : f.getMapEntries()) {
                int key_index = entry.getKey().intValue();
                double p = entry.getValue().doubleValue() / sum;
                factorTable.setValue(key_index, p);
                if (hasNonEnum) {
                    JDF f_jdf = f.getJDF(key_index);
                    JDF cg_jdf = new JDF(nvars);
                    for (Variable nvar : nvars) {
                        Distrib fd = f_jdf.getDistrib(nvar);
                        try {
                            MixtureDistrib md = (MixtureDistrib) fd;
                            cg_jdf.setDistrib(md.getNormalizedClone(), nvar);
                        } catch (ClassCastException e) {
                            cg_jdf.setDistrib(fd, nvar);
                        }
                    }
                    densityTable.setValue(key_index, cg_jdf);
                }
            }
        } else { // no enumerable variables
            atomicFactor = f.getFactor();
            factorTable = null;
            densityTable = null;
            if (hasNonEnum) {
                JDF f_jdf = f.getJDF();
                JDF cg_jdf = new JDF(nvars);
                for (Variable nvar : nvars) {
                    Distrib fd = f_jdf.getDistrib(nvar);
                    try {
                        MixtureDistrib md = (MixtureDistrib) fd;
                        cg_jdf.setDistrib(md, nvar);
                    } catch (ClassCastException e) {
                        cg_jdf.setDistrib(fd, nvar);
                    }
                }
                atomicDensity = cg_jdf;
            } else
                atomicDensity = null;
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
        return evars.size();
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
     * Get indices for all non-null entries.
     * @return the indices
     */
    public int[] getIndices() {
        return factorTable.getIndices();
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
        return null;
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
        return null;
    }

    /**
     * Retrieve a factor value from a factor without enumerable variables. 
     * @return factor value for atomic factor, 0 if not found
     **/
    public Double getFactor() {
        if (atomicFactor != null)
            return atomicFactor;
        return null;
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
     * Calculate the sum of factor values, accounting for assignments (of non-enumerable variables).
     * @param evid evidence
     * @return the sum
     */
    public double getSum(Variable.Assignment... evid) {
        if (factorTable != null) {
            double sum = 0.0;
            for (Map.Entry<Integer, Double> entry : factorTable.getMapEntries()) {
                int index = entry.getKey().intValue();
                double p1 = entry.getValue().doubleValue(); //  joint prob of enumerables
                double p2 = 1.0;
                if (densityTable != null) {
                    JDF jdf = densityTable.getValue(index);
                    p2 = jdf.density(evid); // joint prob of non-enumerable evidence
                }
                sum += p1 * p2;
            }
            return sum;
        } else {
            if (atomicDensity != null)
                return atomicFactor.doubleValue() * atomicDensity.density(evid); // joint prob of non-enumerable evidence
            return atomicFactor.doubleValue();
        }
    }

    /**
     * Calculate the sum of factor values.
     * @param key the representation of the instantiation of enumerable variables, 
     * could be incomplete (null representing unspecified)
     * @param evid evidence
     * @return the sum
     */
    public double getSum(Object[] key, Variable.Assignment... evid) {
        if (factorTable != null) {
            return getSum(factorTable.getIndices(key), evid);
        } else {
            if (atomicDensity != null)
                return atomicFactor.doubleValue() * atomicDensity.density(evid); // joint prob of non-enumerable evidence
            return atomicFactor.doubleValue();
        }
    }
    
    /**
     * Calculate the sum of factor values.
     * @param indices the specific indices of the table
     * @param evid evidence
     * @return the sum
     */
    public double getSum(int[] indices, Variable.Assignment... evid) {
        if (factorTable != null) {
            double sum = 0.0;
            for (int index : indices) {
                double p1 = factorTable.getValue(index);
                double p2 = 1.0;
                if (densityTable != null) {
                    JDF jdf = densityTable.getValue(index);
                    p2 = jdf.density(evid); // joint prob of non-enumerable evidence
                }
                sum += p1 * p2;
            }
            return sum;
        } else {
            if (atomicDensity != null)
                return atomicFactor.doubleValue() * atomicDensity.density(evid); // joint prob of non-enumerable evidence
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
     * Determine the mixture of distributions for all entries that match key
     * @param key the specific entries as identified by a (partial) key
     * @param evid evidence
     * @return the joint mixture distribution
     */
    public JDF getJDFSum(Object[] key, Variable.Assignment... evid) {
        if (densityTable != null) {
            return getJDFSum(densityTable.getIndices(key), evid);
        } else {
            return atomicDensity;
        }
    }
    
    /**
     * Determine the mixture of distributions for all entries that match key
     * @param indices the specific entries as identified by their indices
     * @param evid evidence
     * @return the joint mixture distribution
     */
    public JDF getJDFSum(int[] indices, Variable.Assignment... evid) {
        double sum = getSum(indices, evid);
        if (densityTable != null) {
            JDF jdf = new JDF(getNonEnumVariables());
            for (int index : indices) {
                JDF current = densityTable.getValue(index);
                double p2 = current.density(evid); // joint prob of non-enumerable evidence
                double weight = (factorTable.getValue(index) * p2) / sum;
                jdf = JDF.mix(jdf, current, weight); 
            }
            return jdf;
        } else {
            return atomicDensity;
        }
    }
    
    @Override
    public String toString() {
        StringBuilder sbuf = new StringBuilder("C(");
        for (Variable v : getEnumVariables()) {
            sbuf.append(v.toString()).append(";");
        }
        sbuf.append(";");
        for (Variable v : getNonEnumVariables()) {
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
        for (int j = 0; j < factorTable.getParents().size(); j++)
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
        for (int j = 0; j < factorTable.getParents().size(); j++)
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

    @Override
    public JPT getJPT() {
        return new JPT(factorTable);
    }

    @Override
    public Map<Variable, EnumTable<Distrib>> getNonEnum() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Map<Variable, Distrib> getNonEnumDistrib() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Method to retrieve a single distribution with all other variables in
     * original query unspecified.
     *
     * @param var the query variable
     * @return the distribution of the query variable GIVEN only the evidence
     * that was provided to the original query, all other variables that were
     * designated query variables in the original query are marginalised.
     */
    public Distrib query(Variable var) {
        return query(var, (Variable.Assignment[]) null);
    }

    /**
     * Method to retrieve a single distribution GIVEN some optional evidence.
     *
     * @param query the query variable
     * @param evid evidence, i.e. assignments of values to variables
     * @return the distribution of the query variable GIVEN the evidence, GIVEN
     * the evidence that was also provided to the original query
     */
    public Distrib query(Variable query, Variable.Assignment... evid) {
        if (evid == null) {
            evid = new Variable.Assignment[]{};
        }
        if (factorTable != null) {
            Object[] key = factorTable.getKey(evid);
            // Two cases: enumerable queries, or non-enumerable
            try { // try enumerable first
                EnumVariable equery = (EnumVariable) query;
                int query_index = getEnumVariables().indexOf(equery);
                if (query_index == -1) 
                    throw new RuntimeException("Error in query variable: " + equery);
                double[] p = new double[equery.size()];
                double sum = 0.0;
                for (int i = 0; i < equery.size(); i++) {
                    Object instance = equery.getDomain().get(i);
                    key[query_index] = instance;
                    p[i] = getSum(key, evid);
                    sum += p[i];
                }
                EnumDistrib d = new EnumDistrib(equery.getDomain(), p);
                return d;
            } catch (ClassCastException e) { // so it must be non-enumerable
                JDF jdf = getJDFSum(key, evid);
                Distrib d = jdf.getDistrib(query);
                return d;
            }
        } else { // atomic, can only be non-enumerables
            JDF jdf = atomicDensity;
            Distrib d = jdf.getDistrib(query);
            return d;
        }
    }
    
}