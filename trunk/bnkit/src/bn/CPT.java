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
import java.util.Map;
import java.util.Random;

/**
 * Class for Conditional Probability Table (CPT). This is a table for a
 * conditioned enumerable variable, that has any number (incl 0) enumerable
 * parents.
 *
 * @author mikael
 */
public class CPT implements BNode, Serializable {

    private static final long serialVersionUID = 1L;
    final private EnumVariable var;
    final private EnumTable<EnumDistrib> table; // table of (enumerable) probability distributions
    private EnumDistrib prior; // one (enumerable) probability distribution that is used if this variable is NOT conditioned
    final private int nParents;
    private CountTable count = null; // keep counts when learning/observing; first "parent" is the conditioned variable, then same order as in CPT

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
                this.table = new EnumTable<EnumDistrib>(parents);
                this.prior = null;
                this.nParents = parents.size();
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
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
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
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
                this.prior.set(jptkey[0], entry.getValue().doubleValue());
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
                d.set(jptkey[index], entry.getValue().doubleValue());
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
    public Distrib getDistrib(Object[] key) {
        if (this.table == null || key == null)
            return this.getDistrib();
        try {
            return this.table.getValue(key);
        } catch (EnumTableRuntimeException e) {
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
     *
     * @param bn the BNet instance that can be used to check the status of nodes
     * so that factoring can be done (instantiation of variables are done for a
     * BNet node).
     */
    @Override
    public Factor makeFactor(BNet bn) {
        List<EnumVariable> vars_old = this.getParents();
        EnumVariable var = this.getVariable();
        Object varinstance = null;
        BNode cnode = bn.getNode(var);
        if (cnode != null) {
            varinstance = cnode.getInstance();
        }
        Enumerable dom = var.getDomain();
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
                vars_new.add(var);
            }
            Factor ft = new Factor(vars_new);
            if (varinstance != null) {
                ft.evidenced = true;
            } else {
                for (int i = 0; i < searchkey.length; i++) {
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
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                Factor ft = new Factor();
                ft.setFactor(this.prior.get(varinstance));
                return ft;
            }
            List<Variable> vars_new = new ArrayList<>(1);
            vars_new.add(var);
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
     * Get the name of the CPT
     */
    public String getName() {
        return getVariable().getName();
    }

    /**
     * Get the variable of the CPT.
     *
     * @return the variable of the CPT
     */
    public EnumVariable getVariable() {
        return var;
    }

    /**
     * Retrieve the names of all parent variables (that is all variables that
     * are conditioning the CPT variable)
     *
     * @return the variables of the parent variables
     */
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
    public boolean isRoot() {
        if (table == null) {
            return true;
        }
        return false;
    }

    /**
     * Get the conditional probability of the variable (represented by this CPT)
     *
     * @param key parent key (condition); if null the CPT is assumed to be a
     * prior.
     * @return the probability of the variable
     */
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
     * @param key parent key (condition); if null the CPT is assumed to be a
     * prior.
     * @return the probability of the variable
     */
    public Double get(Object value, Object... key) {
        if (key == null) {
            return prior.get(value);
        }
        if (key.length == 0) {
            return prior.get(value);
        }
        return table.getValue(key).get(value);
    }

    public Double get(Object value) {
        if (isRoot()) {
            return prior.get(value);
        } else {
            throw new RuntimeException("Not a prior");
        }
    }

    public EnumTable getTable() {
        return table;
    }

    public EnumDistrib getDistrib() {
        return prior;
    }

    /**
     * Set entry (or entries) of the CPT to the specified probability value
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
     * Set entry (or entries) of the CPT to the specified probability value
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
     * Set the prior probability of this CPT that has no parents.
     *
     * @param prob
     */
    public void put(EnumDistrib prob) {
        if (!isPrior()) {
            throw new RuntimeException("Unable to set prior. CPT " + var + " is conditioned.");
        }
        if (!prob.isValid()) {
            throw new RuntimeException("Probability value is invalid: " + prob);
        }
        prior = prob;
    }

    /**
     * Checks if this CPT has no parents.
     */
    public boolean isPrior() {
        return table == null;
    }

    /**
     * Provide a non-unique string representation of this CPT.
     */
    public String toString() {
        if (isPrior()) {
            return "CPT(" + getVariable().getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            StringBuffer sbuf = new StringBuffer();
            for (int i = 0; i < table.nParents; i++) {
                sbuf.append(table.getParents().get(i).toString() + (i < table.nParents - 1 ? "," : ""));
            }
            return "CPT(" + getVariable().getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        }
    }

    /**
     * Just a pretty-print of the title (can be modified for sub-classes so the
     * tables look nice)
     */
    protected String formatTitle() {
        return String.format(" %10s", var.getName());
    }

    /**
     * Just a pretty-print of the value (can be modified for sub-classes so the
     * tables look nice)
     */
    protected String formatValue(EnumDistrib x) {
        StringBuffer sbuf = new StringBuffer("<");
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
        if (count == null) { // create count table if none exists
            List<EnumVariable> cond = new ArrayList<EnumVariable>();
            cond.add(var); // first variable is always the conditioned variable
            if (table != null) { // then add parents, if any
                cond.addAll(table.getParents());
            }
            count = new CountTable(cond);
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
        if (count == null) { // create count table if none exists
            List<EnumVariable> cond = new ArrayList<EnumVariable>();
            cond.add(var); // first variable is always the conditioned variable
            if (table != null) { // then add parents, if any
                cond.addAll(table.getParents());
            }
            count = new CountTable(cond);
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
     * {@link bn.CPT#countInstance(Object[], Object, Double)}, ie implement the
     * M-step locally.
     */
    @Override
    public void maximizeInstance() {
        if (count == null) {
            return;
        }
        if (table != null) { // there are parents in the CPT
            // add the counts to the CPT
            for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
                double nobserv = entry.getValue().doubleValue();
                Object[] cntkey = count.table.getKey(entry.getKey().intValue());
                Object[] cptkey = new Object[cntkey.length - 1];
                for (int i = 0; i < cptkey.length; i++) {
                    cptkey[i] = cntkey[i + 1];
                }
                EnumDistrib d = this.table.getValue(cptkey);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                    d.set(cntkey[0], nobserv);
                    this.put(cptkey, d);
                } else {
                    d.set(cntkey[0], nobserv);
                }
            } // normalisation happens internally when values are required			
        } else { // there are no parents
            Object[] cntkey = new Object[1];
            double[] cnts = new double[var.size()];
            for (int i = 0; i < var.size(); i++) {
                cntkey[0] = var.getDomain().get(i);
                cnts[i] = count.get(cntkey);
            }
            prior = new EnumDistrib(this.var.getDomain(), cnts);	// EnumDistrib normalises the counts internally
        }
        count = null; // reset counts
    }

    protected CountTable getCount() {
        return count;
    }

    /**
     * Put random entries in the CPT if not already set.
     */
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
    public void setTrainable(boolean status) {
        trainable = status;
    }

    protected boolean trainable = true;

    /**
     * Check if this CPT should be trained or not
     */
    public boolean isTrainable() {
        return trainable;
    }

    @Override
    public String getStateAsText() {
        StringBuffer sbuf = new StringBuffer("\n");
        if (isPrior()) {
            EnumDistrib d = prior;
            if (d != null) {
                double[] distrib = d.get();
                for (int j = 0; j < distrib.length; j++) {
                    sbuf.append("" + distrib[j]);
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
                    sbuf.append(i + ": ");	// use index as key because values above can be of different non-printable types
                    for (int j = 0; j < distrib.length; j++) {
                        sbuf.append("" + distrib[j]);
                        if (j < distrib.length - 1) {
                            sbuf.append(", ");
                        }
                    }
                    sbuf.append("; (");
                    // If we want to *see* the key, may not work well for some non-printable types
                    Object[] key = table.getKey(i);
                    for (int j = 0; j < key.length; j++) {
                        if (j < key.length - 1) {
                            sbuf.append(key[j] + ", ");
                        } else {
                            sbuf.append(key[j] + ")\n");
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
                        e.printStackTrace();
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
                                    e.printStackTrace();
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

    /**
     * @param args
     */
    public static void main(String[] args) {
        EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Boolean();
        EnumVariable v3 = Predef.Boolean();

        CPT cpt1 = new CPT(v1, new EnumVariable[]{v2, v3});
        cpt1.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{1, 0}));
        cpt1.print();
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
    }

}
