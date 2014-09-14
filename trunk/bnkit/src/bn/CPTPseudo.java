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
import java.util.*;
import java.util.Map.Entry;


/**
 * Class for Pseudo Conditional Probability Table (CPT). This is a table for a
 * conditioned enumerable variable, that has any number (incl 0) enumerable
 * parents.
 *
 * This class is an alternative to the current CPT class, allowing the
 * addition of pseudo counts for the node. Essentially identical to te current CPT implementation,
 * this class also allows users to provide pseudocounts in the form of a PseudoMatrix object.
 * These pseudocounts are added to the CountTable `count' when countinstance is first called.
 * Currently this only applies to the main countinstance method.
 *
 * Given this is a rather `hacked' together implementation, there are several issues that need to be accounted
 * for when using this:
 *      - If this node has many parents (or even a few with many classes), searching and
 *      populating the count table will be incredibly computationally expensive. Additionally, it is unlikely users will
 *      be able to provide a full pseudocounts matrix accounting for every possibility. Hence, the user is required to
 *      set a main parent if the node has more than one parent. Essentially this is equivalent to saying "This parent
 *      is the main parent/variable which is influencing this child node and I only want to apply pseudo counts
 *      specifically in respect to this relationship". This is done with - setMainParentIndex(String name).
 *      This parent should be the one whose domain entries correspond respectively to each row in a provided PseudoMatrix.
 *      (Where columns should correspond to each possible observation across the child's domain).
 *
 *
 * @author julian
 */

public class CPTPseudo implements BNode, TiedNode<CPTPseudo>, Serializable {
    private static final long serialVersionUID = 1L;
    final private EnumVariable var;
    private EnumTable<EnumDistrib> table; // table of (enumerable) probability distributions
    private EnumDistrib prior; // one (enumerable) probability distribution that is used if this variable is NOT conditioned
    final private int nParents;
    private CountTable count = null; // keep counts when learning/observing; first "parent" is the conditioned variable, then same order as in CPT
    // pseudo count matrix; each row correlates to the parent key state and each column the state of this node
    private PseudoMatrix pseudoMatrix; //stores the pseudo counts
    private Integer mainparent_index;
    private boolean relevant = false; //for inference, track whether the node is relevant to the query

    /**
     * Create a conditional probability table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public CPTPseudo(EnumVariable var, List<EnumVariable> parents, PseudoMatrix pseudo) {
        this.var = var;
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<EnumDistrib>(parents);
                this.prior = null;
                this.nParents = parents.size();
                List<EnumVariable> cond = new ArrayList<>();
                cond.add(var); // first variable is always the conditioned variable
                cond.addAll(parents);
                this.count = new CountTable(cond);
                this.pseudoMatrix = pseudo;
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
        this.pseudoMatrix = pseudo;
    }

    public CPTPseudo(EnumVariable var, List<EnumVariable> parents) {
        this.var = var;
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<EnumDistrib>(parents);
                this.prior = null;
                this.nParents = parents.size();
                List<EnumVariable> cond = new ArrayList<>();
                cond.add(var); // first variable is always the conditioned variable
                cond.addAll(this.getParents());
                this.count = new CountTable(cond);
                this.pseudoMatrix = null;
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
        this.pseudoMatrix = null;
    }

    /**
     * Create a conditional probability table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public CPTPseudo(EnumVariable var, PseudoMatrix pseudo, EnumVariable... parents) {
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
                this.pseudoMatrix = pseudo;
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
        this.pseudoMatrix = pseudo;
    }

    public CPTPseudo(EnumVariable var, EnumVariable... parents){
        this(var, null, parents);
    }

    /**
     * Create a prior (a CPT without conditioning variables) for a variable.
     *
     * @param var variable
     */
    public CPTPseudo(EnumVariable var, PseudoMatrix pseudo) {
        this.var = var;
        this.table = null;
        this.nParents = 0;
        List<EnumVariable> cond = new ArrayList<>();
        cond.add(var); // first variable is always the conditioned variable
        this.count = new CountTable(cond);
        this.pseudoMatrix = pseudo;
    }

    /**
     * Create a CPT from a JPT. The variable in the JPT var is the variable
     * conditioned on in the CPT.
     *
     * @param jpt
     * @param var
     */
    public CPTPseudo(JPT jpt, EnumVariable var, PseudoMatrix pseudo) {
        this.var = var;
        this.pseudoMatrix = pseudo;
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
                this.prior.set(jptkey[0], entry.getValue().doubleValue());
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
                d.set(jptkey[index], entry.getValue().doubleValue());
                this.table.setValue(cpt_index, d);
            }
            for (Map.Entry<Integer, EnumDistrib> cpt_entry : this.table.getMapEntries()) {
                cpt_entry.getValue().normalise();
            }
            this.nParents = cptParents.size();
        }
    }

    /**
     * Create a CPTPseudo from a JPT. The variable in the JPT var is the variable
     * conditioned on in the CPTPseudo.
     *
     * @param jpt
     * @param var
     */
    public CPTPseudo(JPT jpt, EnumVariable var) {
        this.var = var;
        List<EnumVariable> CPTPseudoParents = new ArrayList<>(jpt.getParents().size() - 1);
        int index = -1;
        for (int i = 0; i < jpt.getParents().size(); i++) {
            EnumVariable jptParent = jpt.getParents().get(i);
            if (jptParent == var) {
                index = i;
            } else {
                CPTPseudoParents.add(jptParent);
            }
        }
        if (index == -1) {
            throw new RuntimeException("Invalid variable " + var + " for creating CPTPseudo");
        }
        if (CPTPseudoParents.isEmpty()) { // no parents in CPTPseudo
            this.table = null;
            this.prior = new EnumDistrib(var.getDomain());
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                this.prior.set(jptkey[0], entry.getValue());
            }
            this.prior.normalise();
            this.nParents = 0;
        } else { // there are parents
            this.table = new EnumTable<>(CPTPseudoParents);
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                Object[] CPTPseudokey = new Object[jptkey.length - 1];
                int j = 0;
                for (int i = 0; i < jptkey.length; i++) {
                    if (i != index) // selected variable
                    {
                        CPTPseudokey[j++] = jptkey[i];
                    }
                }
                int CPTPseudo_index = this.table.getIndex(CPTPseudokey);
                EnumDistrib d = this.table.getValue(CPTPseudo_index);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                }
                d.set(jptkey[index], entry.getValue());
                this.table.setValue(CPTPseudo_index, d);
            }
            for (Map.Entry<Integer, EnumDistrib> CPTPseudo_entry : this.table.getMapEntries()) {
                CPTPseudo_entry.getValue().normalise();
            }
            this.nParents = CPTPseudoParents.size();
        }
    }

    /**
     * Assign entries in CPTPseudo according to a given JPT. Not bullet-proof, watch
     * out for CPTPseudo - JPT incompatibilities, e.g. variable order which is not
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
                    throw new RuntimeException("No variable " + jptParent + " in CPTPseudo");
                }
            }
        }
        if (index == -1) {
            throw new RuntimeException("No variable " + var + " in JPT");
        }
        if (this.table == null && jpt.getParents().size() == 1) { // no parents in CPTPseudo
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                this.prior.set(jptkey[0], entry.getValue());
            }
            this.prior.normalise();
        } else if (jpt.getParents().size() == this.getParents().size() + 1) { // there are parents
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                Object[] CPTPseudokey = new Object[jptkey.length - 1];
                int j = 0;
                for (int i = 0; i < jptkey.length; i++) {
                    if (i != index) // selected variable
                    {
                        CPTPseudokey[j++] = jptkey[i];
                    }
                }
                int CPTPseudo_index = this.table.getIndex(CPTPseudokey);
                EnumDistrib d = this.table.getValue(CPTPseudo_index);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                }
                d.set(jptkey[index], entry.getValue());
                this.table.setValue(CPTPseudo_index, d);
            }
            for (Map.Entry<Integer, EnumDistrib> CPTPseudo_entry : this.table.getMapEntries()) {
                CPTPseudo_entry.getValue().normalise();
            }
        } else {
            throw new RuntimeException("Cannot set CPTPseudo from given JPT: " + jpt);
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
        } catch (EnumTableRuntimeException e) {
            throw new RuntimeException("Evaluation of CPTPseudo " + this.toString() + " failed since condition was not fully specified: " + e.getMessage());
        }
    }

    @Override
    public Distrib makeDistrib(Collection<Sample> samples) {
        throw new UnsupportedOperationException("Not supported yet.");
    }


    /**
     * Make a Factor out of this CPTPseudo. If a variable is instantiated it will
     * be factored out.
     *
     * @param bn the BNet instance that can be used to check the status of nodes
     * so that factoring can be done (instantiation of variables are done for a
     * BNet node).
     * @return factor of CPTPseudo wrt instantiation of bn
     */
    @Override
    public Factor makeFactor(BNet bn) {
        List<EnumVariable> vars_old = this.getParents();
        EnumVariable myvar = this.getVariable();
        Object varinstance = null;
        BNode cnode = bn.getNode(myvar);
        if (cnode != null) {
            varinstance = cnode.getInstance();
        }
        Enumerable dom = myvar.getDomain();
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<Variable> vars_new = new ArrayList<>(vars_old.size() + 1);
            for (int i = 0; i < vars_old.size(); i++) {
                EnumVariable parent = vars_old.get(i);
                BNode pnode = bn.getNode(parent);
                if (pnode != null) {
                    searchkey[i] = pnode.getInstance();
                }
                if (searchkey[i] == null) {
                    vars_new.add(parent);
                }
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
                    if (varinstance != null) { // the variable for this CPTPseudo is instantiated
                        if (newkey.length == 0) // atomic FactorTable
                            ft.setFactor(null, d.get(varinstance));
                        else
                            ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this CPTPseudo is NOT instantiated so we add one entry for each possible instantiation
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
     * Make a Factor out of this CPTPseudo. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * @param bn the BNet instance that can be used to check the status of nodes
     * so that factoring can be done (instantiation of variables are done for a
     * BNet node).
     * @param rel true if you want to only include relevant nodes
     * @return factor of CPTPseudo wrt bn, considering if parents are relevant (rel)
     */
    @Override
    public Factor makeFactor(BNet bn, boolean rel) {
        List<EnumVariable> vars_old = this.getParents();
        EnumVariable myvar = this.getVariable();
        Object varinstance = null;
        BNode cnode = bn.getNode(myvar);
        if (cnode != null) {
            varinstance = cnode.getInstance();
        }
        Enumerable dom = myvar.getDomain();
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<Variable> vars_new = new ArrayList<>(vars_old.size() + 1);
            List<EnumVariable> irrel_pars = new ArrayList<>(); //irrelevant parents
            for (int i = 0; i < vars_old.size(); i++) {
                EnumVariable parent = vars_old.get(i);
                // Record irrelevant parents to sum out
                //FIXME when should a parent be removed? Allow it to influence factor table then remove it?
                // If parent is evidenced it will not be included in factor table
                if (!bn.getNode(parent).isRelevant() && bn.getNode(parent).getInstance() == null) {
                    irrel_pars.add(parent);
                }
                BNode pnode = bn.getNode(parent);
                if (pnode != null) {
                    searchkey[i] = pnode.getInstance();
                }
                if (searchkey[i] == null) {
                    vars_new.add(parent);
                }
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
                    if (varinstance != null) { // the variable for this CPTPseudo is instantiated
                        if (newkey.length == 0) // atomic FactorTable
                            ft.setFactor(null, d.get(varinstance));
                        else
                            ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this CPTPseudo is NOT instantiated so we add one entry for each possible instantiation
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
     * Get the name of the CPTPseudo
     * @return
     */
    @Override
    public String getName() {
        return getVariable().getName();
    }

    /**
     * Get the variable of the CPTPseudo.
     *
     * @return the variable of the CPTPseudo
     */
    @Override
    public EnumVariable getVariable() {
        return var;
    }

    /**
     * Retrieve the names of all parent variables (that is all variables that
     * are conditioning the CPTPseudo variable)
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
     * Check if this CPTPseudo is a "root" CPTPseudo, i.e. has no parents.
     *
     * @return true if the CPTPseudo has no parents, false if it has
     */
    @Override
    public boolean isRoot() {
        return table == null;
    }

    /**
     * Get the conditional probability of the variable (represented by this CPTPseudo)
     *
     * @param key parent key (condition); if null the CPTPseudo is assumed to be a
     * prior.
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
     * Get the conditional probability of the variable (represented by this CPTPseudo)
     *
     * @param key parent key (condition); if null the CPTPseudo is assumed to be a
     * prior.
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
     * Set entry (or entries) of the CPTPseudo to the specified probability value
     * (variable is true).
     *
     * @param key the boolean key (probabilistic condition)
     * @param prob the probability value (must be >=0 and <=1)
     */
    public void put(Object[] key, EnumDistrib prob) {
        if (key == null) {
            put(prob);
        } else if (key.length == 0) {
            put(prob);
        } else {
            table.setValue(key, prob);
        }
    }

    /**
     * Set entry (or entries) of the CPTPseudo to the specified probability value
     * (variable is true).
     *
     * @param prob the probability value (must be >=0 and <=1)
     * @param key the key (the condition)
     */
    public void put(EnumDistrib prob, Object... key) {
        if (key == null) {
            put(prob);
        } else if (key.length == 0) {
            put(prob);
        } else {
            table.setValue(key, prob);
        }
    }

    /**
     * Set the prior probability of this CPTPseudo that has no parents.
     *
     * @param prob
     */
    public void put(EnumDistrib prob) {
        if (!isPrior()) {
            throw new RuntimeException("Unable to set prior. CPTPseudo " + var + " is conditioned.");
        }
        if (!prob.isNormalised()) {
            throw new RuntimeException("Probability value is invalid: " + prob);
        }
        prior = prob;
    }

    /**
     * Checks if this CPTPseudo has no parents.
     * @return
     */
    public boolean isPrior() {
        return table == null;
    }

    /**
     * Provide a non-unique string representation of this CPTPseudo.
     */
    @Override
    public String toString() {
        if (isPrior()) {
            return "CPTPseudo(" + getVariable().getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            StringBuilder sbuf = new StringBuilder();
            for (int i = 0; i < table.nParents; i++) {
                sbuf.append(table.getParents().get(i).toString()).append(i < table.nParents - 1 ? "," : "");
            }
            return "CPTPseudo(" + getVariable().getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
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
     * Set the variable of this CPTPseudo to a constant value. This means that parent
     * variables will NOT influence inference.
     *
     * @param value the value that is assigned to this instantiated CPTPseudo
     */
    @Override
    public void setInstance(Object value) {
        instance = value;
    }

    /**
     * Set the variable of this CPTPseudo to unspecified, or NOT instantiated.
     */
    @Override
    public void resetInstance() {
        instance = null;
    }

    /**
     * Retrieve the instantiated value of this CPTPseudo.
     *
     * @return the value of this CPTPseudo if instantiated, null if the CPTPseudo is not
     * instantiated.
     */
    @Override
    public Object getInstance() {
        return instance;
    }

    /**
     * Set which parent is the main parent; Stores the index.
     * We rely on the indexing to be consistent
     * @param mainParentName the name of the main parent for this node
     */
    public void setMainParentIndex(String mainParentName) {
        List<EnumVariable> parents = this.getParents();
        if (parents == null){
            throw new RuntimeException("Can't set a main parent, this node has no parents!");
        } else{
            for (int i = 0; i < parents.size(); i++){
                EnumVariable aparent = parents.get(i);
                if (mainParentName.equals(aparent.getName())) {
                    this.mainparent_index = i;
                    break;
                }
            }
        }
    }

    public Integer getMainParentIndex(){
        if (this.mainparent_index != null || isPrior()) {
            //if main parent has been set, do nothing. Indexing will be correct.
        }
        //not set but one parent, so just automatically set
        else if (this.mainparent_index == null && this.nParents == 1) {
            this.mainparent_index = 0;
        } else {
            throw new RuntimeException("No main parent has been specified and number of" +
                    " parents is more than one; use set MainParentIndex");
        }
        return this.mainparent_index;
    }

    /**
     * Count this observation. Note that for it (E-step in EM) to affect the
     * CPT, {@link bn.CPT#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the CPT variable
     * @param prob the expectation of seeing this observation (1 if we actually
     * see it, otherwise the probability)
     * @see bn.CPT#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        if (prob == 0.0) {
            return;
        }
        if (count.table.map.isEmpty()) {
            //        System.out.println();
            //CPTPseudo specific. Here the count table is initialized with pseudo counts
            //Domain lengths determine the matrix[i][j]...up to user to supply correctly formatted pseudo matrix
            Integer p_idx = getMainParentIndex();
            Enumerable cdom = this.getVariable().getDomain(); //child (this node's) domain
            if (key == null) { // if the node is a root
                //then create new key of length 1
                Object[] newKey = new Object[1];
                for (int j = 0; j < cdom.size(); j++) { //go through child domain
                    Object co = new Object[]{cdom.get(j)}; // child observation
                    double obsCount = this.pseudoMatrix.getValue(0, j); //the count for the child observation
                    newKey[0] = co;
                    count.count(newKey, obsCount);
                }
            } else {
                //get the parent domain. p_idx will have to be != 0 if more than one parent
                Enumerable pdom = this.getParents().get(p_idx).getDomain(); //parent domain
                for (int i = 0; i < pdom.size(); i++) {
                    Object po = pdom.get(i); //parent observation
                    for (int j = 0; j < cdom.size(); j++) {
                        Object co = cdom.get(j); // child observation
                        double obsCount = this.pseudoMatrix.getValue(i, j); //the count
                        //                        double obsCount = this.pseudoMatrix.getValue(i, co); //Don't use this...specific use only
                        // add one as count table key includes child observation
                        Object[] newKey = new Object[key.length + 1];
                        newKey[0] = co;
                        for (int x = 1; x < newKey.length; x++) {
                            if (x == p_idx + 1) {
                                newKey[x] = po; //We use the observation for the main parent
                            } else {
                                newKey[x] = null; //All other parent's keys are set to null
                            }
                        }
                        //get all possible indices for the marginalised key
                        int[] countable_idxs = count.table.getTheoreticalIndices(newKey);
                        //get all possible indices for the marginalised key
                        //for each index, add the corresponding count
                        for (int x = 0; x < countable_idxs.length; x++) {
                            count.count(countable_idxs[x], obsCount);
                        }
                    }
                }
            }
        }
        if (key == null) {
            key = new Object[0];
        }
        Object[] mykey = new Object[key.length + 1];
        mykey[0] = value;
        for (int i = 0; i < key.length; i++) {
            mykey[i + 1] = key[i];
        }
        count.count(mykey, prob);
    }

    /**
     * Prob can be set to 1.0 because when counted the value is being observed??
     * Count this observation. Note that for it (E-step in EM) to affect the
     * CPT, {@link bn.CPT#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the CPT variable
     * @see bn.CPT#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value) {
        if (count.table.map.isEmpty()) { // create count table if none exists
            //CPTPseudo specific. Here the count table is initialized with pseudo counts
            //Domain lengths determine the matrix[i][j]...up to user to supply correctly formatted pseudo matrix
            Integer p_idx = getMainParentIndex();
            Enumerable cdom = this.getVariable().getDomain(); //child (this node's) domain
            if (key == null) { // if the node is a root
                //then create new key of length 1
                Object[] newKey = new Object[1];
                for (int j = 0; j < cdom.size(); j++){ //go through child domain
                    Object co = cdom.get(j); // child observation
                    double obsCount = this.pseudoMatrix.getValue(0, j); //the count for the child observation
                    newKey[0] = co;
                    count.count(newKey, obsCount);
                }
            } else {
                //get the parent domain. p_idx will have to be != 0 if more than one parent
                Enumerable pdom = this.getParents().get(p_idx).getDomain(); //parent domain
                for (int i = 0; i < pdom.size(); i++) {
                    Object po = pdom.get(i); //parent observation
                    for (int j = 0; j < cdom.size(); j++) {
                        Object co = cdom.get(j); // child observation
                        double obsCount = this.pseudoMatrix.getValue(i, j); //the count
//                        double obsCount = this.pseudoMatrix.getValue(i, co);
                        // add one as count table key includes child observation
                        Object[] newKey = new Object[key.length + 1];
                        newKey[0] = co;
                        for (int x = 1; x < newKey.length; x++) {
                            if (x == p_idx + 1) {
                                newKey[x] = po; //We use the observation for the main parent
                            } else {
                                newKey[x] = null; //All other parent's keys are set to null
                            }
                        }
                        //get all possible indices for the marginalised key
                        int[] countable_idxs = count.table.getTheoreticalIndices(newKey);
                        //get all possible indices for the marginalised key
                        //for each index, add the corresponding count
                        for (int x = 0; x < countable_idxs.length; x++) {
                            count.count(countable_idxs[x], obsCount);
                        }
                    }
                }
            }
        }
        if (key == null) {
            key = new Object[0];
        }
        Object[] mykey = new Object[key.length + 1];
        mykey[0] = value;
        for (int i = 0; i < key.length; i++) {
            mykey[i + 1] = key[i];
        }
        //FIXME - is prob = 1.0 for observed instance accurate?
        count.count(mykey, 1.0);
    }

    /**
     * Take stock of all observations counted via
     * {@link bn.CPTPseudo#countInstance(Object[], Object, Double)}, ie implement the
     * M-step locally.
     */
    @Override
    public void maximizeInstance() {
        if (count.table.map.isEmpty()) {
            return;
        }
        if (table != null) { // there are parents in the CPTPseudo
            //Set all 'old' distributions in the CPTPseudo to valid = false
            for (EnumDistrib d : this.table.getValues()) {
                d.setValid(false);
            }
            // add the counts to the CPTPseudo
            for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
                double nobserv = entry.getValue();
                Object[] cntkey = count.table.getKey(entry.getKey().intValue());
                Object[] CPTPseudokey = new Object[cntkey.length - 1];
                for (int i = 0; i < CPTPseudokey.length; i++) {
                    CPTPseudokey[i] = cntkey[i + 1];
                }
                EnumDistrib d = this.table.getValue(CPTPseudokey);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                    d.set(cntkey[0], nobserv);
                    this.put(CPTPseudokey, d);
                } else {
                    d.set(cntkey[0], nobserv);
                }
            } // normalisation happens internally when values are required

            //Remove 'old' entries from CPTPseudo
            for (Entry<Integer, EnumDistrib> entry : table.getMapEntries()) {
                EnumDistrib obs = entry.getValue();
                Object[] CPTPseudokey = table.getKey(entry.getKey().intValue());
                if (!obs.isValid()) {
                    table.map.remove(CPTPseudokey);
                }
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
     * Put random entries in the CPTPseudo if not already set.
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
                    table.setValue(i, EnumDistrib.random(var.getDomain()));
            }
        }
    }

    /**
     * Put random entries in the CPTPseudo.
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
     * Set this CPTPseudo to be trained when the Bayesian network it is part of is
     * trained. A CPTPseudo is trainable (true) by default.
     *
     * @param status true if trainable, false otherwise
     */
    @Override
    public void setTrainable(boolean status) {
        trainable = status;
    }

    protected boolean trainable = true;

    /**
     * Check if this CPTPseudo should be trained or not
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
                EnumDistrib d = table.map.get(new Integer(i));
                if (d != null) {
                    double[] distrib = d.get();
                    sbuf.append(i).append(": ");	// use index as key because values above can be of different non-printable types
                    for (int j = 0; j < distrib.length; j++) {
                        sbuf.append("").append(distrib[j]);
                        if (j < distrib.length - 1) {
                            sbuf.append(", ");
                        }
                    }
                    sbuf.append("; (");
                    // If we want to *see* the key, may not work well for some non-printable types
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
        return "CPTPseudo";
    }

    @Override
    public boolean isRelevant() {
        return relevant;
    }

    @Override
    public void setRelevant(boolean relevant) {
        this.relevant = relevant;
    }

    /**
     * Tie all parameters essential to inference and training for this CPTPseudo to those of another CPTPseudo.
     * Variables should be separate but they are required to (1) be of the same type/domain, and (2) be listed in the same order.
     * @param source the CPTPseudo from which parameters will be copied and held fixed.
     */
    @Override
    public void tieTo(CPTPseudo source) {
        CPTPseudo src = (CPTPseudo)source;
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
        // need to tie:
        // - count (used during learning)
        // - prior (if applicable)
        // - table (if applicable)
        this.prior = source.prior;
        if (this.nParents > 0)
            this.table = source.table.retrofit(this.getParents());
    }

}

