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
package bn.node;

import bn.BNode;
import bn.CountTable;
import bn.Distrib;
import bn.prior.Prior;
import bn.prob.EnumDistrib;
import bn.factor.Factor;
import bn.JPT;
import bn.Predef;
import bn.Sample;
import bn.TiedNode;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import dat.Enumerable;
import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import bn.factor.Factorize;
import static bn.factor.Factorize.exitIfInvalid;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;

/**
 * Class for Conditional Probability Table (CPT). This is a table for a
 * conditioned enumerable variable, that has any number (incl 0) enumerable
 * parents.
 *
 * @author mikael
 */
public class CPT implements BNode, TiedNode<CPT>, Serializable{

    private static final long serialVersionUID = 1L;
    final private EnumVariable var;
    protected EnumTable<EnumDistrib> table; // table of (enumerable) probability distributions
    private EnumDistrib prior; // one (enumerable) probability distribution that is used if this variable is NOT conditioned
    final protected int nParents;
    protected CountTable count = null; // keep counts when learning/observing; first "parent" is the conditioned variable, then same order as in CPT
    private boolean relevant = false; //for inference, track whether the node is relevant to the query
    private BNode tieSource = null;

    /**
     * Create a conditional probability table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public CPT(EnumVariable var, List<EnumVariable> parents) {
        this.var = var;
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<>(parents);
                this.prior = null;
                this.nParents = parents.size();
                List<EnumVariable> cond = new ArrayList<>();
                cond.add(var); // first variable is always the conditioned variable
                cond.addAll(parents);
                this.count = new CountTable(cond);
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
    }

    /**
     * Create a conditional probability table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public CPT(EnumVariable var, EnumVariable... parents) {
        this.var = var;
        if (parents != null) {
            if (parents.length > 0) {
                this.table = new EnumTable<>(EnumVariable.toList(parents));
                this.prior = null;
                this.nParents = parents.length;
                List<EnumVariable> cond = new ArrayList<>();
                cond.add(var); // first variable is always the conditioned variable
                cond.addAll(this.getParents());
                this.count = new CountTable(cond);
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
    }

    /**
     * Create a prior (a CPT without conditioning variables) for a variable.
     *
     * @param var variable
     */
    public CPT(EnumVariable var) {
        this.var = var;
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
    }

    /**
     * Create a CPT from a JPT. The variable in the JPT var is the variable
     * conditioned on in the CPT.
     *
     * @param jpt
     * @param var
     */
    public CPT(JPT jpt, EnumVariable var) {
        this.var = var;
        List<EnumVariable> cptParents = new ArrayList<>(jpt.getParents().size() - 1);
        int index = -1;
        for (int i = 0; i < jpt.getParents().size(); i++) {
            EnumVariable jptParent = jpt.getParents().get(i);
            if (jptParent == var) {
                index = i;
            } else {
                cptParents.add(jptParent);
            }
        }
        if (index == -1) {
            throw new RuntimeException("Invalid variable " + var + " for creating CPT");
        }
        if (cptParents.isEmpty()) { // no parents in CPT
            this.table = null;
            this.prior = new EnumDistrib(var.getDomain());
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                this.prior.set(jptkey[0], entry.getValue());
            }
            this.prior.normalise();
            this.nParents = 0;
        } else { // there are parents
            this.table = new EnumTable<>(cptParents);
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                Object[] cptkey = new Object[jptkey.length - 1];
                int j = 0;
                for (int i = 0; i < jptkey.length; i++) {
                    if (i != index) // selected variable
                    {
                        cptkey[j++] = jptkey[i];
                    }
                }
                int cpt_index = this.table.getIndex(cptkey);
                EnumDistrib d = this.table.getValue(cpt_index);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                }
                d.set(jptkey[index], entry.getValue());
                this.table.setValue(cpt_index, d);
            }
            for (Map.Entry<Integer, EnumDistrib> cpt_entry : this.table.getMapEntries()) {
                cpt_entry.getValue().normalise();
            }
            this.nParents = cptParents.size();
        }
    }

    /**
     * Assign entries in CPT according to a given JPT. Not bullet-proof, watch
     * out for CPT - JPT incompatibilities, e.g. variable order which is not
     * checked currently
     *
     * @param jpt
     */
    public void put(JPT jpt) {
        int index = -1;
        for (int i = 0; i < jpt.getParents().size(); i++) {
            EnumVariable jptParent = jpt.getParents().get(i);
            if (jptParent == var) {
                index = i;
            } else {
                if (!this.getParents().contains(jptParent)) {
                    throw new RuntimeException("No variable " + jptParent + " in CPT");
                }
            }
        }
        if (index == -1) {
            throw new RuntimeException("No variable " + var + " in JPT");
        }
        if (this.table == null && jpt.getParents().size() == 1) { // no parents in CPT
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                this.prior.set(jptkey[0], entry.getValue());
            }
            this.prior.normalise();
        } else if (jpt.getParents().size() == this.getParents().size() + 1) { // there are parents
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                Object[] cptkey = new Object[jptkey.length - 1];
                int j = 0;
                for (int i = 0; i < jptkey.length; i++) {
                    if (i != index) // selected variable
                    {
                        cptkey[j++] = jptkey[i];
                    }
                }
                int cpt_index = this.table.getIndex(cptkey);
                EnumDistrib d = this.table.getValue(cpt_index);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                }
                d.set(jptkey[index], entry.getValue());
                this.table.setValue(cpt_index, d);
            }
            for (Map.Entry<Integer, EnumDistrib> cpt_entry : this.table.getMapEntries()) {
                cpt_entry.getValue().normalise();
            }
        } else {
            throw new RuntimeException("Cannot set CPT from given JPT: " + jpt);
        }
    }

    /**
     * Retrieve the distribution for this node that applies GIVEN the parents' instantiations.
     * Requires all parent nodes to be instantiated.
     * @param key the parent values
     * @return the distribution of the variable for this node
     */
    @Override
    public Distrib getDistrib(Object[] key) {
        if (this.table == null || key == null)
            return this.getDistrib();
        try {
            return this.table.getValue(key);
        } catch (RuntimeException e) {
            throw new RuntimeException("Evaluation of CPT " + this.toString() + " failed since condition was not fully specified: " + e.getMessage());
        }
    }
    
    @Override
    public Distrib makeDistrib(Collection<Sample> samples) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }
    
    /**
     * Make a Factor out of this CPT. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * @param relevant only include relevant nodes, with instantiations if available
     * @return factor of CPT considering if parents are relevant (rel)
     */
    @Override
    public Factor makeFactor(Map<Variable, Object> relevant) {
        List<EnumVariable> vars_old = this.getParents();
        EnumVariable myvar = this.getVariable();
        // get value of this node if any assigned
        Object varinstance = relevant.get(myvar); 
        Enumerable dom = myvar.getDomain();
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<Variable> vars_new = new ArrayList<>(vars_old.size() + 1);
            List<EnumVariable> irrel_pars = new ArrayList<>(); //irrelevant parents
            for (int i = 0; i < vars_old.size(); i++) {
                EnumVariable parent = vars_old.get(i);
                boolean parent_is_relevant = relevant.containsKey(parent);
                // Record irrelevant parents to sum out
                //FIXME when should a parent be removed? Allow it to influence factor table then remove it?
                // If parent is evidenced it will not be included in factor table; it is removed later through marginalization
                if (!parent_is_relevant) 
                    irrel_pars.add(parent);
                else
                    searchkey[i] = relevant.get(parent);
                if (searchkey[i] == null) // new factor needs to include this variable (to be summed out before returned if irrelevant)
                    vars_new.add(parent);
            }
            if (varinstance == null) {
                vars_new.add(myvar);
            }
            Factor ft = new Factor(vars_new);
            if (varinstance != null) {
                ft.evidenced = true;
            } else {
                for (Object searchkey1 : searchkey) {
                    if (searchkey != null) {
                        ft.evidenced = true;
                        break;
                    }
                }
            }
            int[] indices = table.getIndices(searchkey);
            Object[] newkey = new Object[vars_new.size()];
            for (int index : indices) {
                EnumDistrib d = table.getValue(index);
                if (d != null) {
                    Object[] key = table.getKey(index);
                    int newcnt = 0;
                    for (int i = 0; i < key.length; i++) {
                        if (searchkey[i] == null) {
                            newkey[newcnt++] = key[i];
                        }
                    }
                    if (varinstance != null) { // the variable for this CPT is instantiated
                        if (newkey.length == 0) // atomic FactorTable
                            ft.setFactor(null, d.get(varinstance));
                        else
                            ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this CPT is NOT instantiated so we add one entry for each possible instantiation
                        for (int j = 0; j < dom.size(); j++) {
                            newkey[newkey.length - 1] = dom.get(j);
                            Double p = d.get(j);
                            if (p != null) {
                                ft.addFactor(newkey, p);
                            }
                        }
                    }
                } else { // this entry is null
                    //
                }
            }
            if (!irrel_pars.isEmpty()) {
            	ft = ft.marginalize(irrel_pars);
            }
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                Factor ft = new Factor();
                ft.setFactor(this.prior.get(varinstance));
                return ft;
            }
            List<Variable> vars_new = new ArrayList<>(1);
            vars_new.add(myvar);
            Factor ft = new Factor(vars_new);
            Object[] newkey = new Object[1];
            EnumDistrib d = this.prior;
            for (int j = 0; j < dom.size(); j++) {
                newkey[0] = dom.get(j);
                ft.addFactor(newkey, d.get(j));
            }
            return ft;
        }
    }
    /**
     * Make a Factor out of this CPT. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor.
     *
     * Note: In some instances, the CPT is poorly populated with distributions.
     * This introduces the possibility that NO entries are selected for the FT.
     * There appears no principled way to resolve this scenario.
     * Currently, this is managed by creating a reduced FT with only the current variable,
     * ASSUMING that a priori the assignment of values to the variable is described by
     * a uniform distribution.
     *
     * TODO: consider the case when the above is true, and there are other, not-selected
     * entries; use other distributions to create a more realistic distribution for the
     * null case.
     *
     * TODO: Fix GDT and other node types in an analogous manner.
     *
     * @param relevant Only include relevant nodes, with instantiations if available
     * @return factor of CPT considering if parents are relevant (rel)
     */
    @Override
    public AbstractFactor makeDenseFactor(Map<Variable, Object> relevant) {
        List<EnumVariable> parents = this.getParents();
        EnumVariable myvar = this.getVariable();
        // get value of this node if any assigned
        Object varinstance = relevant.get(myvar); 
        Enumerable dom = myvar.getDomain();
        if (parents != null) { // there are parent variables
            Object[] searchcpt = new Object[parents.size()];
            List<Variable> fvars = new ArrayList<>(parents.size() + 1); // factor variables
            List<EnumVariable> sumout = new ArrayList<>();  // irrelevant variables to be summed out later
            for (int i = 0; i < parents.size(); i++) {
                EnumVariable parent = parents.get(i);
                // If parent is evidenced it will not be included in factor table; 
                // Record irrelevant parents to sum out: removed later through marginalization
                if (!relevant.containsKey(parent)) 
                    sumout.add(parent);
                else
                    searchcpt[i] = relevant.get(parent);
                if (searchcpt[i] == null) // new factor will include this variable
                    fvars.add(parent);
            }
            if (varinstance == null) {
                fvars.add(myvar);
            }
            Variable[] vars_arr = new Variable[fvars.size()];
            fvars.toArray(vars_arr);
            AbstractFactor ft = new DenseFactor(vars_arr);
            EnumVariable[] evars = ft.getEnumVars(); // the order may have changed
            int[] xcross = new int[parents.size()];
            int[] ycross = new int[evars.length];
            table.crossReference(xcross, evars, ycross);
            int missing = -1; // the position of *myvar* in the new factor (if applicable)
            for (int i = 0; i < ycross.length; i ++)
                if (ycross[i] == -1) {
                    missing = i;
                    break;
                }
            // set factor to be "evidenced" if there was evidence used
            if (varinstance != null) {
                ft.evidenced = true;
            } else {
                for (Object instcpt : searchcpt) {
                    if (instcpt != null) {
                        ft.evidenced = true;
                        break;
                    }
                }
            }
            int[] indices = table.getIndices(searchcpt);
            Object[] fkey = new Object[evars.length];
            if (indices.length == 0) { // The are no entries in CPT that match the search key (which itself indicates the state of parents)
                // create a new, reduced FT with only myvar, since there is no relationship with parents evident from the entries
                ft = new DenseFactor(new Variable[] {myvar});
                fkey = new Object[1];
                if (varinstance != null) { // the variable for this CPT is instantiated
                    fkey[0] = varinstance;
                    ft.setValue(fkey, 1.0 / dom.size());
                } else {
                    for (int j = 0; j < dom.size(); j++) {
                        fkey[0] = dom.get(j);
                        Double p = 1.0 / dom.size();
                        ft.setValue(fkey, p); //add one entry for each possible instantiation to create uniform outcome
                    }
                }
                return ft;
            }
            for (int index : indices) {
                EnumDistrib d = table.getValue(index);
                if (d != null) { // there is a distribution associated with this entry in the CPT
                    Object[] cptkey = table.getKey(index); // work out the condition for this entry
                    for (int i = 0; i < cptkey.length; i++) {
                        if (xcross[i] != -1) 
                            fkey[xcross[i]] = cptkey[i];
                    }
                    if (varinstance != null) { // the variable for this CPT is instantiated
                        if (fkey.length == 0) // atomic factor
                            ft.setValue(d.get(varinstance));
                        else
                            ft.setValue(fkey, d.get(varinstance));
                    } else { // the variable for this CPT is NOT instantiated so we add one entry for each possible instantiation
                        for (int j = 0; j < dom.size(); j++) {
                            fkey[missing] = dom.get(j);
                            Double p = d.get(j);
                            ft.setValue(fkey, p);
                        }
                    }
                } 
            }
            if (!sumout.isEmpty()) {
                Variable[] sumout_arr = new Variable[sumout.size()];
                sumout.toArray(sumout_arr);
            	ft = Factorize.getMargin(ft, sumout_arr);
            }
            return ft;
        } else { // no parents, just a prior
            EnumDistrib d = this.prior;
            if (varinstance != null) { // instantiated prior
                AbstractFactor ft = new DenseFactor();
                if (d == null)
                    ft.setValue(1.0/dom.size());
                else
                    ft.setValue(d.get(varinstance));
                return ft;
            }
            AbstractFactor ft = new DenseFactor(myvar);
            Object[] newkey = new Object[1];
            for (int j = 0; j < dom.size(); j++) {
                newkey[0] = dom.get(j);
                if (d == null)
                    ft.setValue(newkey, 1.0/dom.size());
                else
                    ft.setValue(newkey, d.get(j));
            }
            return ft;
        }
    }
    
    /**
     * Get the name of the CPT
     * @return 
     */
    @Override
    public String getName() {
        return getVariable().toString();
    }

    /**
     * Get the variable of the CPT.
     *
     * @return the variable of the CPT
     */
    @Override
    public EnumVariable getVariable() {
        return var;
    }

    /**
     * Retrieve the names of all parent variables (that is all variables that
     * are conditioning the CPT variable)
     *
     * @return the variables of the parent variables
     */
    @Override
    public List<EnumVariable> getParents() {
        if (table == null) {
            return null;
        }
        List<EnumVariable> parents = table.getParents();
        return parents;
    }

    /**
     * Check if this CPT is a "root" CPT, i.e. has no parents.
     *
     * @return true if the CPT has no parents, false if it has
     */
    @Override
    public boolean isRoot() {
        return table == null;
    }

    /**
     * Get the conditional probability of the variable (represented by this CPT)
     *
     * @param key parent cptkey (condition); if null the CPT is assumed to be a
 prior.
     * @return the probability of the variable
     */
    @Override
    public Double get(Object[] key, Object value) {
        if (key == null) {
            return prior.get(value);
        }
        if (key.length == 0) {
            return prior.get(value);
        }
        return table.getValue(key).get(value);
    }

    /**
     * Get the conditional probability of the variable (represented by this CPT)
     *
     * @param key parent cptkey (condition); if null the CPT is assumed to be a
 prior.
     * @return the probability of the variable
     */
    @Override
    public Double get(Object value, Object... key) {
        if (key == null) {
            return prior.get(value);
        }
        if (key.length == 0) {
            return prior.get(value);
        }
        return table.getValue(key).get(value);
    }

    @Override
    public Double get(Object value) {
        if (isRoot()) {
            return prior.get(value);
        } else {
            throw new RuntimeException("Not a prior");
        }
    }

    /**
     * @deprecated Do NOT use, will be removed in the future
     */
    @Override
    public EnumTable getTable() {
        return table;
    }

    /**
     *
     * @return
     */
    @Override
    public EnumDistrib getDistrib() {
        return prior;
    }

    /**
     * Set entry (or entries) of the CPT to the specified probability value
     * (variable is true).
     *
     * @param key the boolean cptkey (probabilistic condition)
     * @param prob the probability value (must be >=0 and <=1)
     */
    public void put(Object[] key, Distrib prob) {
        if (key == null) {
            put(prob);
        } else if (key.length == 0) {
            put(prob);
        } else {
            table.setValue(key, (EnumDistrib)prob);
        }
    }
    
    /**
     * Set entry (or entries) of the CPT to the specified probability value index
     * (variable is true).
     */
    public void put(int index, Distrib prob) {
    	table.setValue(index, (EnumDistrib)prob);
    }

    /**
     * Set entry (or entries) of the CPT to the specified probability value
     * (variable is true).
     *
     * @param prob the probability value (must be >=0 and <=1)
     * @param key the cptkey (the condition)
     */
    public void put(Distrib prob, Object... key) {
        if (key == null) {
            put(prob);
        } else if (key.length == 0) {
            put(prob);
        } else {
            table.setValue(key, (EnumDistrib) prob);
        }
    }

    /**
     * Set the prior probability of this CPT that has no parents.
     *
     * @param prob
     */
    public void put(Distrib prob) {
        if (!isPrior()) {
            throw new RuntimeException("Unable to set prior. CPT " + var + " is conditioned.");
        }
        if (!((EnumDistrib)prob).isNormalised()) {
            throw new RuntimeException("Probability value is invalid: " + prob);
        }
        prior = (EnumDistrib)prob;
    }

    /**
     * Checks if this CPT has no parents.
     * @return 
     */
    public boolean isPrior() {
        return table == null;
    }

    /**
     * Provide a non-unique string representation of this CPT.
     */
    @Override
    public String toString() {
        if (isPrior()) {
            return "CPT(" + getVariable().getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            StringBuilder sbuf = new StringBuilder();
            for (int i = 0; i < table.nParents; i++) {
                sbuf.append(table.getParents().get(i).toString()).append(i < table.nParents - 1 ? "," : "");
            }
            return "CPT(" + getVariable().getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        }
    }

    /**
     * Just a pretty-print of the title (can be modified for sub-classes so the
     * tables look nice)
     * @return 
     */
    protected String formatTitle() {
        return String.format(" %10s", var.getName());
    }

    /**
     * Just a pretty-print of the value (can be modified for sub-classes so the
     * tables look nice)
     * @param x
     * @return 
     */
    protected String formatValue(EnumDistrib x) {
        StringBuilder sbuf = new StringBuilder("<");
        double[] distrib = x.get();
        for (int i = 0; i < distrib.length; i++) {
            sbuf.append(String.format("%4.2f ", distrib[i]));
        }
        sbuf.replace(sbuf.length() - 1, sbuf.length() - 1, ">");
        return sbuf.toString();
    }

    /**
     * Pretty-print of whole table
     */
    @Override
    public void print() {
        System.out.println(formatTitle());
        if (!isPrior()) { // variables in condition
            table.display();
        } else { // prior
            if (prior != null) {
                System.out.println(formatValue(prior));
            }
        }
    }

    private Object instance = null;

    /**
     * Set the variable of this CPT to a constant value. This means that parent
     * variables will NOT influence inference.
     *
     * @param value the value that is assigned to this instantiated CPT
     */
    @Override
    public void setInstance(Object value) {
        instance = value;
    }

    /**
     * Set the variable of this CPT to unspecified, or NOT instantiated.
     */
    @Override
    public void resetInstance() {
        instance = null;
    }

    /**
     * Retrieve the instantiated value of this CPT.
     *
     * @return the value of this CPT if instantiated, null if the CPT is not
     * instantiated.
     */
    @Override
    public Object getInstance() {
        return instance;
    }

    /**
     * Count this observation. Note that for it (E-step in EM) to affect the
     * CPT, {@link CPT#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the CPT variable
     * @param prob the expectation of seeing this observation (1 if we actually
     * see it, otherwise the probability)
     * @see CPT#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
    	if (prob == 0.0 || Double.isNaN(prob)) {
            return;
    	}
        if (key == null) {
            key = new Object[0];
        }
        Object[] mykey = new Object[key.length + 1];
        mykey[0] = value;
        System.arraycopy(key, 0, mykey, 1, key.length);
        count.count(mykey, prob);
    }
    
    /**
     * Count this observation. Note that for it (E-step in EM) to affect the
     * CPT, {@link CPT#maximizeInstance()} must be called.
     * Probability of observation is 1.0 (the value is definitely being observed, 
     * as opposed to expected with a probability).
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the CPT variable
     * @see CPT#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value) {
        if (key == null) {
            key = new Object[0];
        }
        Object[] mykey = new Object[key.length + 1];
        mykey[0] = value;
        System.arraycopy(key, 0, mykey, 1, key.length);
        count.count(mykey, 1.0);
    }

    /**
     * Take stock of all observations counted via
     * {@link CPT#countInstance(Object[], Object, Double)}, ie implement the
     * M-step locally.
     */
    @Override
    public void maximizeInstance() {
        if (count.table.isEmpty()) {
            return;
        }
        if (table != null) { // there are parents in the CPT
            // Set all 'old' distributions in the CPT to valid = false, i.e.
            // we are marking entries so we can remove 'ghosts' after counting
            for (EnumDistrib d : this.table.getValues()) {
                d.setValid(false);
            }

            // add the counts to the CPT
            for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
                double nobserv = entry.getValue();
                Object[] cntkey = count.table.getKey(entry.getKey().intValue());
                Object[] cptkey = new Object[cntkey.length - 1];
                for (int i = 0; i < cptkey.length; i++) {
                    cptkey[i] = cntkey[i + 1];
                }
                EnumDistrib d = this.table.getValue(cptkey);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain()); // this will automatically set entry to "valid" 
                    d.set(cntkey[0], nobserv);
                    this.put(cptkey, d);
                } else {
                    d.set(cntkey[0], nobserv);
                }
            } // normalisation happens internally when values are required	
    	
            //Remove 'old' (or 'ghost' entries from CPT (for which no counts
            for (Iterator<Entry<Integer, EnumDistrib>> it = table.getMapEntries().iterator(); it.hasNext(); ) {
                Entry<Integer, EnumDistrib> entry = it.next();
            	EnumDistrib obs = entry.getValue();
            	if (!obs.isValid()) 
                    it.remove();
            }
            
        } else { // there are no parents
            Object[] cntkey = new Object[1];
            double[] cnts = new double[var.size()];
            for (int i = 0; i < var.size(); i++) {
                cntkey[0] = var.getDomain().get(i);
                cnts[i] = count.get(cntkey);
            }
            prior = new EnumDistrib(this.var.getDomain(), cnts);	// EnumDistrib normalises the counts internally
        }
        count.table.setEmpty(); // reset counts
    }

    protected CountTable getCount() {
        return count;
    }

    /**
     * Put random entries in the CPT if not already set.
     * @param seed
     */
    @Override
    public void randomize(long seed) {
        Random rand = new Random(seed);
        if (table == null) {
            if (prior == null)
                prior = EnumDistrib.random(var.getDomain());
        } else {
            int nrows = table.getSize();
            for (int i = 0; i < nrows; i++) {
                if (!table.hasValue(i))
                    table.setValue(i, EnumDistrib.random(var.getDomain(), rand.nextLong()));
            }
        }
    }

    /**
     * Put random entries in the CPT.
     */
//	public void randomize(Object[] observations, int seed) {
//		// TODO: use "observations" to come up with good initial probabilities
//		// we do not have access to parents so this is only a rough estimate
//		Random rand=new Random(seed);
//		int nrows=getMaxIndex();
//		if (nrows<2)
//			put(rand.nextDouble());
//		else {
//			for (int i=0; i<nrows; i++)
//				put(entry(i), rand.nextDouble());
//		}
//	}
    /**
     * Set this CPT to be trained when the Bayesian network it is part of is
     * trained. A CPT is trainable (true) by default.
     *
     * @param status true if trainable, false otherwise
     */
    @Override
    public void setTrainable(boolean status) {
        trainable = status;
    }

    protected boolean trainable = true;

    /**
     * Check if this CPT should be trained or not
     */
    @Override
    public boolean isTrainable() {
        return trainable;
    }

    @Override
    public String getStateAsText() {
        StringBuilder sbuf = new StringBuilder("\n");
        if (isPrior()) {
            EnumDistrib d = prior;
            if (d != null) {
                double[] distrib = d.get();
                for (int j = 0; j < distrib.length; j++) {
                    sbuf.append("").append(distrib[j]);
                    if (j < distrib.length - 1) {
                        sbuf.append(", ");
                    }
                }
                sbuf.append(";\n");
            }
        } else {
            for (int i = 0; i < table.getSize(); i++) {
                EnumDistrib d = table.getValue(i);
                if (d != null) {
                    double[] distrib = d.get();
                    sbuf.append(i).append(": ");	// use index as cptkey because values above can be of different non-printable types
                    for (int j = 0; j < distrib.length; j++) {
                        sbuf.append("").append(distrib[j]);
                        if (j < distrib.length - 1) {
                            sbuf.append(", ");
                        }
                    }
                    sbuf.append("; (");
                    // If we want to *see* the cptkey, may not work well for some non-printable types
                    Object[] key = table.getKey(i);
                    for (int j = 0; j < key.length; j++) {
                        if (j < key.length - 1) {
                            sbuf.append(key[j]).append(", ");
                        } else {
                            sbuf.append(key[j]).append(")\n");
                        }
                    }
                }
            }
        }
        return sbuf.toString();
    }

    @Override
    public boolean setState(String dump) {
        if (isPrior()) {
            String[] line = dump.split(";");
            if (line.length >= 1) {
                String[] y = line[0].split(",");
                if (y.length == var.size()) {
                    double[] distrib = new double[y.length];
                    try {
                        for (int i = 0; i < distrib.length; i++) {
                            distrib[i] = Double.parseDouble(y[i]);
                        }
                    } catch (NumberFormatException e) {
                        return false;
                    }
                    this.put(new EnumDistrib(var.getDomain(), distrib));
                    return true;
                }
            }
        } else {
            for (String line : dump.split("\n")) {
                // 0: 0.4, 0.6; (true, true)
                String[] specline = line.split(";");
                if (specline.length >= 1) {
                    String[] parts = specline[0].split(":");
                    if (parts.length >= 2) {
                        try {
                            int index = Integer.parseInt(parts[0]);
                            String[] y = parts[1].split(",");
                            if (y.length == var.size()) {
                                double[] distrib = new double[y.length];
                                try {
                                    for (int i = 0; i < distrib.length; i++) {
                                        distrib[i] = Double.parseDouble(y[i]);
                                    }
                                } catch (NumberFormatException e) {
                                    return false;
                                }
                                this.put(table.getKey(index), new EnumDistrib(var.getDomain(), distrib));
                            }
                        } catch (NumberFormatException e) {
                            System.err.println("Number format wrong and ignored: " + line);
                        }
                    }
                }
            }
        }
        return false;
    }

    @Override
    public String getType() {
        return "CPT";
    }
    
    @Override
    public boolean isRelevant() {
        return relevant;
    }

    @Override
    public void setRelevant(boolean relevant) {
        this.relevant = relevant;
    }

    public BNode getTieSource(){
        return this.tieSource;
    }

    /**
     * Tie all parameters essential to inference and training for this CPT to those of another CPT.
     * Variables should be separate but they are required to (1) be of the same type/domain, and (2) be listed in the same order.
     * @param source the CPT from which parameters will be copied and held fixed.
     */
    @Override
    public void tieTo(CPT source) {
        CPT src = (CPT)source;
        if (!this.var.getDomain().equals(source.getVariable().getDomain()))
            throw new RuntimeException("Invalid sharing: " + var.getName() + " does not share domain with " + source.getVariable().getName());
        if (this.nParents != src.nParents)
            throw new RuntimeException("Invalid sharing: " + var.getName() + " has different number of parents from " + source.getVariable().getName());
        for (int i = 0; i < this.nParents; i ++) {
            Variable p1 = this.getParents().get(i);
            Variable p2 = src.getParents().get(i);
            if (!p1.getDomain().equals(p2.getDomain()))
                throw new RuntimeException("Invalid sharing: " + p1.getName() + " does not share domain with " + p2.getName());
        }
        this.tieSource = source;
        // need to tie:
        // - count (used during learning)
        // - prior (if applicable)
        // - table (if applicable)
        this.prior = source.prior;
        if (this.nParents > 0)
            this.table = source.table.retrofit(this.getParents());
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Boolean();
        EnumVariable v3 = Predef.Boolean();

        CPT cpt1 = new CPT(v1, v2, v3);
        cpt1.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{1, 0}));
        cpt1.print();
        CPT cpt3 = new CPT(v1, v2, v3);
        cpt3.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{0.4, 0.6}));
        cpt3.print();
        cpt3.tieTo(cpt1);
        cpt3.print();
        cpt1.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{0.3, 0.7}));
        cpt1.print();
        cpt3.print();
        CPT cpt4 = new CPT(v1, v2, v3);
        CPT cpt5 = new CPT(v1, v2, v3);
        cpt5.print();
        cpt5.tieTo(cpt4);
        cpt4.countInstance(new Object[]{true, false}, true);
        cpt4.countInstance(new Object[]{false, false}, true);
        cpt4.countInstance(new Object[]{false, false}, false);
        cpt4.countInstance(new Object[]{true, false}, false);
        cpt4.countInstance(new Object[]{true, false}, true);
        cpt4.maximizeInstance();
        cpt4.print();
        cpt5.print();
        System.out.println("-----");
        CountTable cnt = new CountTable(new EnumVariable[]{v1, v2, v3});
        cnt.count(new Boolean[]{false, true, false}, 2);
        cnt.count(new Boolean[]{true, true, false}, 1);
        cnt.count(new Boolean[]{true, false, false}, 3);
        cnt.count(new Boolean[]{false, false, false}, 1);
        cnt.count(new Boolean[]{true, true, true}, 1);
        JPT jpt = new JPT(cnt.table);
        jpt.display();
        CPT cpt2 = new CPT(jpt, v2);
        cpt2.print();
        System.out.println("-----");
        
    }

	@Override
	public List<Sample> getConditionDataset(int conditionIndex) {
		Object[] parentValues = null;
		if(conditionIndex >= 0) { // non-root node
			parentValues = table.getKey(conditionIndex);
		} else { // root node
			parentValues = new Object[0];
		}
		
		List<Sample> data = new LinkedList<Sample>();
		for(int i = 0; i < var.getDomain().size(); i++) {
			Object[] key = new Object[parentValues.length + 1];
			key[0] = var.getDomain().get(i);
			for(int j = 0; j < parentValues.length; j++) {
				key[1 + j] = parentValues[j];
			}
			double sampleProb = count.get(key);
			Sample sample = new Sample(key[0], sampleProb);
			data.add(sample);
		}
		return data;
	}

	@Override
	public Distrib getlikelihoodDistrib() {
		return new EnumDistrib(this.var.getDomain());
	}


}
