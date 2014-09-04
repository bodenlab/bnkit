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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;

/**
 * Class for Dirichlet Density Table (DirDT). This is a table for a variable that can take enumerable 
 * distributions as values, with each row specifying a single Dirichlet density, specifying what 
 * distributions it produces. The node has one or more enumerable parents.
 *
 * @author m.boden
 */
public class DirDT implements BNode, Serializable {

    private static final long serialVersionUID = 1L;
    
    final private Variable<EnumDistrib> var;
    private DirichletDistrib prior = null;
    private EnumTable<DirichletDistrib> table = null;

    // Parameters used for training. Some of which are allocated prior to training, and then re-used to save time.
    private SampleTable<EnumDistrib> count = null; // the table that will contain all samples of the type "EnumDistrib" during learning
    
    private boolean relevant = false;
    private EnumDistrib instance = null; // the value this node takes, null if unspecified
    private String tag = null;

    /**
     * Create a Dirichlet density table for a variable. The variable is
     * conditioned on a set of enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public DirDT(Variable<EnumDistrib> var, List<EnumVariable> parents) {
        this.var = var;
        if (parents != null) {
            if (parents.size() > 0) 
                this.table = new EnumTable<>(parents);
        }
    }

    /**
     * Create a Dirichlet density table for a variable. The variable is
     * conditioned on a set of enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public DirDT(Variable<EnumDistrib> var, EnumVariable... parents) {
        this(var, EnumVariable.toList(parents));
    }

    /**
     * Create a Dirichlet prior for a variable. Note: the prior distribution is
     * not set (@see bn.DirDT#put(DirichletDistrib)).
     * @param var variable
     */
    public DirDT(Variable<EnumDistrib> var) {
        this.var = var;
    }

    /**
     * Assign a tag name for this node.
     * @param name
     */
    public void setTag(String name){
        this.tag = name;
    }

    /**
     * Get the tag name for this node
     * @return tag name
     */
    public String getTag(){
        return this.tag;
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
        } catch (EnumTableRuntimeException e) {
            throw new RuntimeException("Evaluation of DirDT " + this.toString() + " failed since condition was not fully specified: " + e.getMessage());
        }
    }
    
    /**
     * Retrieve the distribution for this node with no parents.
     * @return the distribution of the variable for this root node
     */
    @Override
    public DirichletDistrib getDistrib() {
        return prior;
    }

    /**
     * Make a FactorTable out of this DirDT. If a variable is instantiated it will
     * be factored out.
     *
     * @param bn the BNet instance that can be used to check the status of nodes
     * so that factoring can be done (instantiation of variables are done for a
     * BNet node).
     * @return the FactorTable created from the DirDT, provided instantiations of BN
     */
    @Override
    public Factor makeFactor(BNet bn) {
        List<EnumVariable> vars_old = this.getParents();
        Object varinstance = this.getInstance();
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<Variable> vars_new = new ArrayList<>(vars_old.size() + 1);
            for (int i = 0; i < vars_old.size(); i++) {
                EnumVariable parent = vars_old.get(i);
                BNode pnode = bn.getNode(parent);
                if (pnode != null)
                    searchkey[i] = pnode.getInstance();
                if (searchkey[i] == null)
                    vars_new.add(parent);
            }
            Factor ft;
            if (varinstance == null)
                vars_new.add(this.var);
            ft = new Factor(vars_new);
            if (varinstance != null)
                ft.evidenced = true;
            int[] indices = table.getIndices(searchkey);
            Object[] newkey = new Object[ft.getNEnum()];
            for (int index : indices) {
                DirichletDistrib d = table.getValue(index);
                if (d != null) {
                    Object[] key = table.getKey(index);
                    int newcnt = 0;
                    for (int i = 0; i < key.length; i++) {
                        if (searchkey[i] == null) {
                            newkey[newcnt++] = key[i];
                        }
                    }
                    if (varinstance != null) { // the variable for this DirDT is instantiated
                        ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this DirDT is NOT instantiated...
                        ft.addFactor(newkey, 1.0);
                        ft.setDistrib(newkey, this.var, d);
                    }
                } else { // this entry is null
                    //
                }
            }
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) // instantiated prior is not possible to factorise
                return null;
            throw new RuntimeException("DirDT can not be factorised unless it has enumerable parent variables");
        }
    }

    /**
     * Make a FactorTable out of this DirDT. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * Marginalization technique requires updating.
     *
     * @param bn the BNet instance that can be used to check the status of nodes
     * so that factoring can be done (instantiation of variables are done for a
     * BNet node).
     * @param relevant
     * @return the FactorTable created from the DirDT, provided instantiations of BN
     */
    @Override
    public Factor makeFactor(BNet bn, boolean relevant) {
        List<EnumVariable> vars_old = this.getParents();
        Object varinstance = this.getInstance();
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<Variable> vars_new = new ArrayList<>(vars_old.size() + 1);
            List<EnumVariable> irrel_pars = new ArrayList<>(); //irrelevant parents
            for (int i = 0; i < vars_old.size(); i++) {
                EnumVariable parent = vars_old.get(i);
                // Record irrelevant parents to sum out
                // FIXME when should a parent be removed? Allow it to influence factor table then remove it?
                // If parent is evidenced it will not be included in factor table
                if (!bn.getNode(parent).isRelevant() && bn.getNode(parent).getInstance() == null) {
                	irrel_pars.add(parent);
                }
                BNode pnode = bn.getNode(parent);
                if (pnode != null)
                    searchkey[i] = pnode.getInstance();
                if (searchkey[i] == null)
                    vars_new.add(parent);
            }
            Factor ft;
            if (varinstance == null)
                vars_new.add(this.var);
            ft = new Factor(vars_new);
            if (varinstance != null)
                ft.evidenced = true;
            int[] indices = table.getIndices(searchkey);
            Object[] newkey = new Object[ft.getNEnum()];
            for (int index : indices) {
                DirichletDistrib d = table.getValue(index);
                if (d != null) {
                    Object[] key = table.getKey(index);
                    int newcnt = 0;
                    for (int i = 0; i < key.length; i++) {
                        if (searchkey[i] == null) {
                            newkey[newcnt++] = key[i];
                        }
                    }
                    if (varinstance != null) { // the variable for this DirDT is instantiated
                        ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this DirDT is NOT instantiated...
                        ft.addFactor(newkey, 1.0);
                        ft.setDistrib(newkey, this.var, d);
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
            if (varinstance != null) // instantiated prior is not possible to factorise
                return null;
            throw new RuntimeException("DirDT can not be factorised unless it has enumerable parent variables");
        }
    }
    
    /**
     * Get the conditional probability of the variable (represented by this DirDT)
     * when set to a specified value.
     *
     * @param key condition
     * @param value the value of the variable represented by this DirDT
     * @return the probability density of the variable set to specified value
     */
    public Double get(Object[] key, Object value) {
        if (key == null) {
            if (prior != null) {
                return prior.get(value);
            }
            return null;
        }
        DirichletDistrib d = table.getValue(key);
        if (d != null) {
            Double p = d.get(value);
            return p;
        }
        return null;
    }

    /**
     * Get the conditional probability of the variable (represented by this DirDT)
     * when set to a specified value.
     *
     * @param value the value of the variable represented by this DirDT
     * @param key condition
     * @return the probability density of the variable set to specified value
     */
    public Double get(Object value, Object... key) {
        if (key == null) {
            if (prior != null) {
                return prior.get(value);
            }
            return null;
        }
        DirichletDistrib d = table.getValue(key);
        if (d != null) {
            Double p = d.get(value);
            return p;
        }
        return null;
    }

    /**
     * Get the prior probability of the variable (represented by this DirDT)
     *
     * @param value the value of the variable represented by this DirDT
     * @return the probability of the variable set to specified value
     */
    @Override
    public Double get(Object value) {
        if (prior != null) {
            return prior.get(value);
        }
        return null;
    }

    @Override
    public EnumTable getTable() {
        return table;
    }

    /**
     * Set entry (or entries) of the DirDT to the specified probability
     * distribution
     *
     * @param key the boolean key (probabilistic condition)
     * @param distr the distribution
     */
    public void put(Object[] key, DirichletDistrib distr) {
        table.setValue(key, distr);
    }

    /**
     * Set entry (or entries) of the DirDT to the specified probability
     * distribution
     *
     * @param index the index for the key (probabilistic condition)
     * @param distr the distribution
     */
    public void put(int index, DirichletDistrib distr) {
        table.setValue(index, distr);
    }

    /**
     * Set entry (or entries) of the DirDT to the specified probability
     * distribution
     *
     * @param distr the distribution
     * @param key the boolean key (probabilistic condition)
     */
    public void put(DirichletDistrib distr, Object... key) {
        table.setValue(key, distr);
    }

    public void put(DirichletDistrib distr) {
        prior = distr;
    }

    @Override
    public String toString() {
        List<EnumVariable> parents = table.getParents();
        StringBuilder sbuf = new StringBuilder();
        for (int i = 0; i < parents.size(); i++) {
            sbuf.append(parents.get(i).getName()).append(i < parents.size() - 1 ? "," : "");
        }
        return "DirDT(" + getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
    }

    protected String formatTitle() {
        return String.format(" %10s", var.getName());
    }

    protected String formatValue(DirichletDistrib x) {
        return String.format("<%s>", x.toString());
    }

    @Override
    public void print() {
        if (table.nParents > 0) // variables in condition
        {
            table.display();
        } else { // prior
            System.out.println(formatTitle());
            if (prior != null) {
                System.out.println(formatValue(prior));
            }
        }
    }

    @Override
    public boolean isRoot() {
        return table == null;
    }

    @Override
    public void setInstance(Object value) {
        try {
            instance = (EnumDistrib) value;
        } catch (ClassCastException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void resetInstance() {
        instance = null;
    }

    @Override
    public EnumDistrib getInstance() {
        return instance;
    }
    
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        if (this.isRoot()) {
            throw new RuntimeException("DirDT can not be trained as root");
            // same process as for entries with parents, just a single queue of observations...
        } else {
            if (count == null) // create count table if none exists
                count = new SampleTable<>(this.getParents());
            count.count(key, (EnumDistrib)value, prob);
        }
    }
    
    @Override
    public void countInstance(Object[] key, Object value) {
        countInstance(key, value, 1.0);
    }

    /**
     * Set this DirDT to be trained when the Bayesian network it is part of is
     * trained. A DirDT is trainable (true) by default.
     *
     * @param status true if trainable, false otherwise
     */
    public void setTrainable(boolean status) {
        trainable = status;
    }

    protected boolean trainable = true;

    /**
     * Check if this DirDT should be trained or not
     */
    @Override
    public boolean isTrainable() {
        return trainable;
    }

    /**
     * Put random entries in the DirDT if not already set.
     * @param seed
     */
    @Override
    public void randomize(long seed) {
        Random rand = new Random(seed);
        if (table == null) {
            if (prior == null)
                prior = new DirichletDistrib((Enumerable)prior.getDomain(), (1+rand.nextGaussian()));
        } else {
            int nrows = table.getSize();
            for (int i = 0; i < nrows; i++) {
                if (!table.hasValue(i))
                    throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
        }
    }

    /**
     * Find the best parameter setting for the observed data.
     * Note that this uses observed Double:s which are looked at directly, and 
     * observed Distrib:s which are used to stochastically generate samples. 
     * We cannot use distributions directly since they can be of any kind, as 
     * long as defined over a single continuous variable.
     */
    @Override
    public void maximizeInstance() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String getName() {
        return getVariable().getName();
    }

    @Override
    public Variable<EnumDistrib> getVariable() {
        return var;
    }

    @Override
    public List<EnumVariable> getParents() {
        if (table == null) {
            return null;
        }
        List<EnumVariable> parents = table.getParents();
        return parents;
    }

    @Override
    public String getType() {
        return "DirDT";
    }

    @Override
    public String getStateAsText() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean setState(String dump) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public boolean isRelevant() {
        return relevant;
    }

    public void setRelevant(boolean relevant) {
        this.relevant = relevant;
    }

    @Override
    public Distrib makeDistrib(Collection<Sample> samples) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}