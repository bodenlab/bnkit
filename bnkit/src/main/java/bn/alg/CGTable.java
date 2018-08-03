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
import bn.prob.EnumDistrib;
import dat.EnumTable;
import dat.EnumVariable;
import bn.factor.Factor;
import bn.JDF;
import bn.JPT;
import bn.prob.MixtureDistrib;
import dat.Variable;
import bn.factor.AbstractFactor;
import bn.factor.Factorize;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
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
    private final EnumTable<Set<Variable.Assignment>> assignTable;    // the traced assignments for each permutation of the enumerables
    private final List<EnumVariable> evars;       // all enumerable variables associated with this factor table
    private final List<Variable> nvars;           // all non-enumerable variables
    
    private final Double atomicFactor;    // if there are no enumerable variables, here's where the constant factor is
    private final JDF    atomicDensity;   // if there are no enumerable but some non-enumerable variables, here's where the only joint density is stored
    private final Set<Variable.Assignment>    atomicAssign;   // if there are no enumerable but some non-enumerable variables, here's where the only assignment set is stored
    
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
                int key_index = entry.getKey();
                double p = entry.getValue() / sum;
                if (p == 0)
                	continue;
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
        assignTable = null; // not used for this constructor, requires AbstractFactor
        atomicAssign = null;
    } 
    
    public CGTable(AbstractFactor f, List<Variable> qvars) {
        f = Factorize.getNormal(f); // normalise to make sure that the factor table is ok to operate on in terms of probability
        evars = new ArrayList<>();
        nvars = new ArrayList<>();
        for (Variable var : qvars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                evars.add(evar);
            } catch (ClassCastException e) {
                nvars.add(var);
            }
        }
        EnumVariable[] q_evars_arr = new EnumVariable[evars.size()];
        evars.toArray(q_evars_arr);
        EnumVariable[] evars_arr = f.getEnumVars();
        Variable[] nvars_arr = f.getNonEnumVars() ;
        int[] fcross2q = new int[evars_arr.length]; // factor index to query
        int[] qcross2f = new int[evars_arr.length]; // query index to factor
        Factorize.getCrossref(evars_arr, fcross2q, q_evars_arr, qcross2f);
        if (f.hasEnumVars()) {
            atomicDensity = null;
            if (f.hasNonEnumVars()) {
                densityTable = new EnumTable<>(evars);
            } else
                densityTable = null;
            if (f.isTraced()) {
                assignTable = new EnumTable<>(evars);
                atomicAssign = null;
            } else {
                assignTable = null;
                atomicAssign = null;
            }
            factorTable = new EnumTable<>(evars);
            atomicFactor = null;
            for (int i = 0; i < f.getSize(); i ++) {
                Object[] fkey = f.getKey(i);
                Object[] qkey = new Object[fkey.length];
                for (int j = 0; j < fkey.length; j ++) 
                    qkey[fcross2q[j]] = fkey[j];
                int key_index = EnumTable.getIndex(qkey, q_evars_arr);
                double p = f.getValue(i);
                if (p == 0 || Double.isNaN(p))
                    continue;
                factorTable.setValue(key_index, p);
                if (f.hasNonEnumVars()) {
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
                if (f.isTraced()) {
                    Set<Variable.Assignment> a = f.getAssign(i);
                    if (a != null)
                        assignTable.setValue(key_index, a);
                }
            }
        } else { // no enumerable variables
            atomicFactor = f.getValue();
            factorTable = null;
            densityTable = null;
            assignTable = null;
            if (f.hasNonEnumVars()) {
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
            if (f.isTraced()) {
                Set<Variable.Assignment> a = f.getAssign();
                if (a != null)
                    atomicAssign = a;
                else
                    atomicAssign = null;
            } else
                atomicAssign = null;
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
                sum += entry.getValue();
            return sum;
        } else {
            return atomicFactor;
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
            return atomicFactor;
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
            return atomicFactor;
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
                int index = entry.getKey();
                double p1 = entry.getValue(); //  joint prob of enumerables
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
                return atomicFactor * atomicDensity.density(evid); // joint prob of non-enumerable evidence
            return atomicFactor;
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
                return atomicFactor * atomicDensity.density(evid); // joint prob of non-enumerable evidence
            return atomicFactor;
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
                return atomicFactor * atomicDensity.density(evid); // joint prob of non-enumerable evidence
            return atomicFactor;
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
    
    /**
     * The usual string representation of a CGTable instance.
     * @return a string representation of CGTable
     */
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
        if (this.isAtomic()) {
            sbuf.append(") = ").append(this.atomicFactor.toString());
            if (this.atomicDensity != null)
                sbuf.append("|").append(this.atomicDensity.toString());
        } else
            sbuf.append(")");
        return sbuf.toString();
    }
    
    private String constantLength(String s, int len) {
        if (s.length() > len)
            return s.substring(0, len);
        else return String.format("%-10s", s);
    }
    
    /**
     * Pretty-print of a CGTable instance. 
     * Prints enumerable variable assignments as table, associated with a probability, then 
     * non-enumerable variable assignments in their original form, conditioned on enumerables.
     */
    public void display() {
        if (this.isAtomic()) { // does not have enumerable variables
            System.out.println(this.atomicFactor.toString());
            List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
            for (Variable nonenum : nonenums) 
                System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
            System.out.println();
            for (Variable nonenum : nonenums) {
                Distrib d = this.getDistrib(nonenum);
                System.out.print(String.format(" %s ", d.toString()));
            }
            return;
        } 
        // has enumerable variables
        System.out.print("Idx ");
        for (int j = 0; j < factorTable.getParents().size(); j++)
            System.out.print(String.format("[%10s]", constantLength(factorTable.getParents().get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" P");
        for (int i = 0; i < factorTable.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = factorTable.getKey(i);
            for (Object key1 : key)
                System.out.print(String.format(" %-10s ", constantLength(key1.toString(), 10)));
            JDF jdf = this.getJDF(i);
            if (jdf != null) {
                for (Variable nonenum : nonenums) {
                    Distrib d = jdf.getDistrib(nonenum);
                    System.out.print(String.format(" %s ", d));
                    //System.out.print(String.format(" %-10s ", constantLength(d.toString(), 10)));
                }
            }
            Object val = this.getFactor(i);
            if (val != null) 
                System.out.println(String.format(" %7.5f", this.getFactor(i)));
            else
                System.out.println(" null ");
        }
    }

    /**
     * Pretty-print of a CGTable instance. 
     * Prints enumerable variable assignments as table, associated with a probability, then 
     * non-enumerable variable assignments as means of samples, conditioned on enumerables.
     */
    public void displaySampled() {
        if (this.isAtomic()) { // does not have enumerable variables
            System.out.println(this.atomicFactor.toString());
            List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
            for (Variable nonenum : nonenums) 
                System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
            System.out.println();
            JDF jdf = this.getJDF();
            if (jdf != null) {
                for (Variable nonenum : nonenums) {
                    double sum = 0;
                    if (jdf.getDistrib(nonenum) != null) {
                        for (int j = 0; j < 1000; j ++) {
                            Double sample = (Double)jdf.sample(nonenum);
                            sum += sample;
                        }
                        System.out.print(String.format(" %10.3f ", sum/1000.0));
                    }
                }
            }
            return;
        } 
        // has enumerable variables
        System.out.print("Idx ");
        for (int j = 0; j < factorTable.getParents().size(); j++)
            System.out.print(String.format("[%10s]", constantLength(factorTable.getParents().get(j).toString(), 10)));
        List<Variable> nonenums = new ArrayList<>(this.getNonEnumVariables());
        for (Variable nonenum : nonenums) 
            System.out.print(String.format("[%10s]", constantLength(nonenum.toString(), 10)));
        System.out.println(" P");
        for (int i = 0; i < factorTable.getSize(); i++) {
            System.out.print(String.format("%3d ", i));
            Object[] key = factorTable.getKey(i);
            for (Object key1 : key)
                System.out.print(String.format(" %-10s ", constantLength(key1.toString(), 10)));
            JDF jdf = this.getJDF(i);
            if (jdf != null) {
                for (Variable nonenum : nonenums) {
                    double sum = 0;
                    for (int j = 0; j < 1000; j ++) {
                        Double sample = (Double)jdf.sample(nonenum);
                        sum += sample;
                    }
                    System.out.print(String.format(" %10.3f ", sum/1000.0));
                }
            }
            Object val = this.getFactor(i);
            if (val != null) 
                System.out.println(String.format(" %7.5f", this.getFactor(i)));
            else
                System.out.println(" null ");
        }
    }

    /**
     * Retrieve JPT from results.
     * 
     * @return JPT of result (result for enumerables only)
     * @deprecated This method is implemented but only because QueryResult requires it. Do not use.
     */
    @Deprecated
    @Override
    public JPT getJPT() {
        return new JPT(factorTable);
    }

    @Deprecated
    @Override
    public Map<Variable, EnumTable<Distrib>> getNonEnum() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Deprecated
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
                for (int i = 0; i < equery.size(); i++) {
                    Object instance = equery.getDomain().get(i);
                    key[query_index] = instance;
                    p[i] = getSum(key, evid);
                }
                EnumDistrib d = new EnumDistrib(equery.getDomain(), p);
                return d;
            } catch (ClassCastException e) { // so it must be non-enumerable
                JDF jdf = getJDFSum(key, evid);
                if (jdf != null) {
                    Distrib d = jdf.getDistrib(query);
                    return d;
                }
                return null;
            }
        } else { // atomic, can only be non-enumerables
            JDF jdf = atomicDensity;
            Distrib d = jdf.getDistrib(query);
            return d;
        }
    }
    
    /**
     * Method to retrieve the most probable explanation to the evidence, in terms of query variables. 
     * @return the most probable assignment of query variables
     */
    public Variable.Assignment[] getMPE() {
        
        if (assignTable == null && atomicAssign == null) { // old version of MPE
            Variable.Assignment[] assign = new Variable.Assignment[evars.size() + nvars.size()];
            int mostProbKey = 0;
            double mostProbValue = -1;
            if (this.isAtomic()) {
                    for (int i = 0; i < nvars.size(); i ++){
                        assign[i] = Variable.assign(nvars.get(i), atomicDensity.getDistrib(nvars.get(i)));
                    }
                    return assign;
            } else {
                    for (Map.Entry<Integer, Double> entry : factorTable.getMapEntries()) {
                        if (entry.getValue() > mostProbValue) {
                            mostProbKey = entry.getKey();
                            mostProbValue = entry.getValue();
                        }
                    }
            }
            //values.length will be no. enumVars
            Object[] values = factorTable.getKey(mostProbKey);
            //add all enum assignments to list
            for (int i = 0; i < values.length; i ++)
                assign[i] = Variable.assign(evars.get(i), values[i]);
            //add all nonEnum assignments to list
            if (nvars.size() > 0) {
                JDF jdf = densityTable.getValue(mostProbKey);
                for (int i = 0; i < nvars.size(); i ++){
                    assign[i + values.length] = Variable.assign(nvars.get(i), jdf.getDistrib(nvars.get(i)));
                }
            }
            return assign;
        } else { // new version based on AbstractFactor
            Set<Variable.Assignment> result = new HashSet<>();

            if (this.isAtomic()) {
                result.addAll(atomicAssign);
                if (this.hasNonEnumVariables() && atomicDensity != null) {
                    for (int i = 0; i < nvars.size(); i ++)
                        result.add(Variable.assign(nvars.get(i), atomicDensity.getDistrib(nvars.get(i))));
                }
            } else {
                int mostProbKey = 0;
                double mostProbValue = -1;
                for (Map.Entry<Integer, Double> entry : factorTable.getMapEntries()) {
                    if (entry.getValue() > mostProbValue) {
                        mostProbKey = entry.getKey();
                        mostProbValue = entry.getValue();
                    }
                }
                result.addAll(assignTable.getValue(mostProbKey));
                if (this.hasNonEnumVariables()) {
                    JDF jdf = densityTable.getValue(mostProbKey);
                    if (jdf != null) {
                        for (int i = 0; i < nvars.size(); i ++)
                            result.add(Variable.assign(nvars.get(i), jdf.getDistrib(nvars.get(i))));
                    }
                }
            }
            Variable.Assignment[] assign = new Variable.Assignment[result.size()];
            result.toArray(assign);
            return assign;
        }
    }
}
