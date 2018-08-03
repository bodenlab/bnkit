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
import bn.Distrib;
import bn.prob.EnumDistrib;
import bn.factor.Factor;
import bn.Sample;
import bn.SampleTable;
import bn.TiedNode;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import dat.Enumerable;
import bn.prob.DirichletDistrib;
import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import bn.factor.Factorize;
import dat.Domain;
import dat.IntegerSeq;
import java.io.Serializable;
import java.util.*;

/**
 * Class for Dirichlet Density Table (DirDT). This is a table for a variable that can take enumerable 
 * distributions as values, with each row specifying a single Dirichlet density, specifying what 
 * distributions it produces. The node has one or more enumerable parents.
 *
 * @author m.boden
 */
public class DirDT implements BNode, TiedNode, Serializable {

    private static final long serialVersionUID = 1L;
    
    final private Variable<EnumDistrib> var;
    private DirichletDistrib prior = null;
    private EnumTable<DirichletDistrib> table = null;

    // Parameters used for training. Some of which are allocated prior to training, and then re-used to save time.
    private SampleTable<IntegerSeq> count = null; // the table that will contain all samples of the type "IntegerSeq" during learning
    
    private boolean relevant = false;
    private Domain instance = null; // the value this node takes, null if unspecified, needs to implement a domain, currently EnumDistrib and IntegerSeq are supported

    private DirDT tieSource;
    
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
            if (parents.size() > 0) {
                this.table = new EnumTable<>(parents);
                this.prior = null;
                this.count = new SampleTable(parents);
            }
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
     * If a parent is not relevant, it will not be included in the factor
     *
     * Marginalization technique requires updating.
     *
     * @param relevant only include relevant nodes, with instantiations if available
     * @return the FactorTable created from the DirDT, provided instantiations of BN
     */
    @Override
    public Factor makeFactor(Map<Variable, Object> relevant) {
        List<EnumVariable> vars_old = this.getParents();
        Variable myvar = this.getVariable();
        // get value of this node if any assigned
        Object varinstance = relevant.get(myvar); 
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
                        try {
                            IntegerSeq iseq = (IntegerSeq) varinstance;
                            int[] counts = IntegerSeq.intArray(iseq.get());
                            ft.addFactor(newkey, Math.exp(d.logLikelihood(counts)));
                        } catch (ClassCastException e) {
                            // Assume the instance is an EnumDistrib
                            ft.addFactor(newkey, d.get(varinstance));
                        }
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
     * Make factor out of this DirDT.
     * This method will deconstruct the condition, accounting for which variables that are relevant 
     * (as nominated by the caller of this function, e.g. inference algorithm) and those that are
     * instantiated. Instantiated and irrelevant variables are removed before the factor is returned.
     * @param relevant a map containing all relevant variables as keys, and values which are non-null if the variable
     * is instantiated
     * @return the factor
     */
    @Override
    public AbstractFactor makeDenseFactor(Map<Variable, Object> relevant) {
        List<EnumVariable> parents = this.getParents();
        Variable myvar = this.getVariable();
        // get value of this node if any assigned
        Object varinstance = relevant.get(myvar); 
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
                fvars.add(myvar); // add to factor, but note this is a non-enumerable variable so handled differently inside
            }
            Variable[] vars_arr = new Variable[fvars.size()];
            fvars.toArray(vars_arr);
            AbstractFactor ft = new DenseFactor(vars_arr);
            EnumVariable[] evars = ft.getEnumVars(); // the order may have changed; excludes non-enums
            int[] xcross = new int[parents.size()];
            int[] ycross = new int[evars.length];
            table.crossReference(xcross, evars, ycross);
            // set factor to be "evidenced" is there was an evidence used
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
            for (int index : indices) {
                DirichletDistrib d = table.getValue(index);
                if (d != null) { // there is a distribution associated with this entry in the CPT
                    Object[] cptkey = table.getKey(index); // work out the condition for this entry
                    for (int i = 0; i < cptkey.length; i++) {
                        if (xcross[i] != -1)
                            fkey[xcross[i]] = cptkey[i];
                    }
                    if (varinstance != null) { // the variable for this DirDT is instantiated
                        try {
                            IntegerSeq iseq = (IntegerSeq) varinstance;
                            int[] counts = IntegerSeq.intArray(iseq.get());
                            if (fkey.length == 0) // and the parents are too
                                ft.setLogValue(d.logLikelihood(counts)); // these probs can be very small, stay in log-space...
                            else
                                ft.setLogValue(fkey, d.logLikelihood(counts));
                        } catch (ClassCastException e) {
                            // Assume the instance is an EnumDistrib
                            if (fkey.length == 0) // and the parents are too
                                ft.setValue(d.get(varinstance));
                            else
                                ft.setValue(fkey, d.get(varinstance));
                        }
                    } else { // the variable for this DirDT is NOT instantiated so we put it in the JDF
                        if (fkey.length == 0) { // but the parents are instantiated
                            ft.setLogValue(0.0); 
                            ft.setDistrib(myvar, d);
                        } else {
                            ft.setLogValue(fkey, 0.0);
                            ft.setDistrib(fkey, myvar, d);
                        }
                    }
                } 
            }
            //ft = Factorize.getNormal(ft);
            if (!sumout.isEmpty()) {
                Variable[] sumout_arr = new Variable[sumout.size()];
                sumout.toArray(sumout_arr);
            	ft = Factorize.getMargin(ft, sumout_arr);
            }
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) // instantiated prior is not possible to factorise
                return null;
            throw new RuntimeException("DirDTs can not be factorised unless it has enumerable parent variables");
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

    /**
     * @deprecated Do NOT use, will be removed in the future
     */
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
    public void put(Object[] key, Distrib distr) {
        table.setValue(key, (DirichletDistrib) distr);
    }

    /**
     * Set entry (or entries) of the DirDT to the specified probability
     * distribution
     *
     * @param index the index for the key (probabilistic condition)
     * @param distr the distribution
     */
    public void put(int index, Distrib distr) {
        table.setValue(index, (DirichletDistrib) distr);
    }

    /**
     * Set entry (or entries) of the DirDT to the specified probability
     * distribution
     *
     * @param distr the distribution
     * @param key the boolean key (probabilistic condition)
     */
    public void put(Distrib distr, Object... key) {
        table.setValue(key, (DirichletDistrib)distr);
    }

    public void put(Distrib distr) {
        prior = (DirichletDistrib)distr;
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
        if (table.nParents > 0) { // variables in condition
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
        } catch (ClassCastException e1) {
            try {
                instance = (IntegerSeq) value;
            } catch (ClassCastException e2) {
                System.err.println("Invalid setInstance: " + this.getName() + " = " + value);
            }
        }
    }

    @Override
    public void resetInstance() {
        instance = null;
    }

    @Override
    public Object getInstance() {
        return instance;
    }
    
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        if (prob == 0)
            return;
        if (this.isRoot()) {
            throw new RuntimeException("DirDT can not be trained as root");
            // same process as for entries with parents, just a single queue of observations...
        } else {
            try {
                count.count(key, (IntegerSeq)value, prob);
            } catch (ClassCastException e) {
                System.err.println("Invalid instance, must implement IntegerSeq: " + this.getName() + " = " + value);
            }
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
                prior = new DirichletDistrib((Enumerable)prior.getDomain(), Math.max(1+rand.nextInt(100), 100*Math.abs(rand.nextGaussian())));
        } else {
            int nrows = table.getSize();
            for (int i = 0; i < nrows; i++) {
                if (!table.hasValue(i))
                    table.setValue(i, new DirichletDistrib((Enumerable)this.var.getDomain().getDomain(), Math.max(1+rand.nextInt(100), 100*Math.abs(rand.nextGaussian()))));
            }
        }
    }

    private class DnIpair { 
        Double toss;
        Integer i;
        DnIpair(Double toss, Integer i) { this.toss = toss; this.i = i; }
    }
    
    /**
     * Find the best parameter setting for the observed data.
     * This uses a method to find the ML estimate for the Dirichlet from
     * a set of IntegerSeq, counted in countInstance.
     * The current version doe not weight samples but collects them similar to
     * how k-means operate.
     * 
     * Note that this method is NOT used currently.
     */
    //@Override
    public void maximizeInstance_Gibbs() {
        if (count.isEmpty()) {
            return;
        }
        Random rand = new Random();
        Enumerable e = this.var.getDomain().getDomain();
        EnumTable<List<Sample<IntegerSeq>>> samples = count.getTable();
        if (samples != null) {
            Map<IntegerSeq, DnIpair> sample2index = new HashMap<>(); 
            for (Map.Entry<Integer, List<Sample<IntegerSeq>>> entry : samples.getMapEntries()) {
                int index = entry.getKey();
                for (Sample<IntegerSeq> sample : entry.getValue()) {
                    double p = sample.prob;
                    IntegerSeq ed = sample.instance;
                    DnIpair prev = sample2index.get(sample.instance);
                    if (prev != null) { // found the sample already associated with an index, re-consider?
                        if (prev.toss <= p) // pick this one? Yes... so set others to impossible
                            sample2index.put(ed, new DnIpair(1.0, index));
                        else 
                            sample2index.put(ed, new DnIpair(prev.toss - p, index));
                    } else { // first appearance
                        double toss = rand.nextDouble(); // the value that will decide which index to use for this sample
                        if (toss <= p) // use this one, set prob of forthcoming to impossible
                            sample2index.put(ed, new DnIpair(1.0, index));
                        else // do not use this index, set the prob of next: toss - toss
                            sample2index.put(ed, new DnIpair(toss - p, index));
                    }
                }
            }
            List[] select = new ArrayList[this.table.getSize()];
            for (int index = 0; index < this.table.getSize(); index ++)
                select[index] = new ArrayList();
            int cnt = sample2index.size();
            for (Map.Entry<IntegerSeq, DnIpair> entry : sample2index.entrySet()) {
                IntegerSeq ed = entry.getKey();
                int i = entry.getValue().i;
                select[i].add(ed);
            }
            for (int index = 0; index < this.table.getSize(); index ++) {
                List selected = select[index];
                if (selected.isEmpty()) {
                    for (int attempt = 0; attempt < cnt / this.table.getSize() * 2; attempt ++) {
                        IntegerSeq stash = null;
                        int pick = rand.nextInt(cnt);
                        int accum = 0;
                        for (int j = 0; j < select.length; j ++) {
                            if (pick < accum + select[j].size()) {
                                stash = (IntegerSeq)select[j].remove(pick - accum);
                                selected.add(stash);
                                break;
                            } 
                            accum += select[j].size();
                        }
                    }
                }
            }
            for (int index = 0; index < this.table.getSize(); index ++) {
                List selected = select[index];
                if (selected.size() > 0) {
                    IntegerSeq[] dists = new IntegerSeq[selected.size()];
                    selected.toArray(dists);
                    int[][] hists = new int[selected.size()][];
                    for (int j = 0; j < hists.length; j ++) {
                        hists[j] = IntegerSeq.intArray(((IntegerSeq)selected.get(j)).get());
                    }
                    double[] location = DirichletDistrib.getAlpha(hists);
                    DirichletDistrib dd = new DirichletDistrib(e, location);
                    table.setValue(index, dd);
                } else {
                    System.err.println("Cannot happen");
                }
            }
        }
        count.setEmpty();
    }

    /**
     * Find the best parameter setting for the observed data.
     * This uses a method to find the ML estimate for the Dirichlet from
     * a set of IntegerSeq, counted in countInstance.
     */
    @Override
    public void maximizeInstance() {
        if (count.isEmpty()) {
            return;
        }
        Random rand = new Random();
        Enumerable e = this.var.getDomain().getDomain();
        for (int index = 0; index < this.table.getSize(); index ++) {
            List<Sample<IntegerSeq>> samples = count.getAll(index);
            if (samples != null) {
                int[][] hists = new int[samples.size()][];
                double[] probs = new double[samples.size()];
                // do the calculations
                for (int j = 0; j < samples.size(); j ++) {
                    Sample<IntegerSeq> sample = samples.get(j);
                    IntegerSeq d = (IntegerSeq)sample.instance;
                    probs[j] = sample.prob;
                    hists[j] = IntegerSeq.intArray(sample.instance.get());
                }
                try {
                    double[] alpha = DirichletDistrib.getAlpha(hists, probs);
                    DirichletDistrib dd = new DirichletDistrib(e, alpha);
                    this.put(index, dd);
                } catch (StackOverflowError ex) {
                    System.err.println("Stack overflow in node " + this);
                }
            } else { // no counts
                this.table.removeValue(index);
            }
        }
        count.setEmpty();
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
        StringBuilder sbuf = new StringBuilder("\n");
        if (isRoot()) {
            DirichletDistrib d = prior;
            if (d != null) {
                double[] alpha = d.getAlpha();
                for (int a = 0; a < alpha.length; a ++)
                    sbuf.append(alpha[a]).append((a == alpha.length - 1 ? ";" : ", "));
                sbuf.append("\n");
            }
        } else {
            for (int i = 0; i < table.getSize(); i++) {
                DirichletDistrib d = table.getValue(i);
                if (d != null) {
                    sbuf.append(i).append(": ");	// use index as key because values above can be of different non-printable types
                    double[] alpha = d.getAlpha();
                    for (int a = 0; a < alpha.length; a ++)
                        sbuf.append(alpha[a]).append((a == alpha.length - 1 ? ";" : ", "));
                    // If we want to *see* the key, may not work well for some non-printable types
                    sbuf.append(" (");
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
        Enumerable e = this.var.getDomain().getDomain();
        if (isRoot()) {
            String[] line = dump.split(";");
            if (line.length >= 1) {
                String[] y = line[0].split(",");
                if (y.length == e.size()) {
                    double[] alpha = new double[y.length];
                    try {
                        for (int i = 0; i < alpha.length; i++) {
                            alpha[i] = Double.parseDouble(y[i]);
                        }
                    } catch (NumberFormatException ex) {
                        ex.printStackTrace();
                        return false;
                    }
                    this.put(new DirichletDistrib(e, alpha));
                    return true;
                }
            }
        } else {
            for (String line : dump.split("\n")) {
                line = line.trim();
                // 0: 0.4, 0.6; (true, true)
                String[] specline = line.split(";");
                if (specline.length >= 1) {
                    String[] parts = specline[0].split(":");
                    if (parts.length >= 2) {
                        try {
                            int index = Integer.parseInt(parts[0]);
                            String[] y = parts[1].split(",");
                            if (y.length == e.size()) {
                                double[] alpha = new double[y.length];
                                try {
                                    for (int i = 0; i < alpha.length; i++) {
                                        alpha[i] = Double.parseDouble(y[i]);
                                    }
                                } catch (NumberFormatException ex) {
                                    ex.printStackTrace();
                                    return false;
                                }
                                this.put(table.getKey(index), new DirichletDistrib(e, alpha));
                            }
                        } catch (NumberFormatException ex) {
                            System.err.println("Number format wrong and ignored: " + line);
                        }
                    }
                }
            }
        }
        return false;
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

    public void tieTo(BNode source) {
        DirDT src = (DirDT)source;
        if (!this.var.getDomain().equals(source.getVariable().getDomain()))
            throw new RuntimeException("Invalid sharing: " + var.getName() + " does not share domain with " + source.getVariable().getName());
        if (this.table.nParents != src.table.nParents)
            throw new RuntimeException("Invalid sharing: " + var.getName() + " has different number of parents from " + source.getVariable().getName());
        for (int i = 0; i < this.table.nParents; i ++) {
            Variable p1 = this.getParents().get(i);
            Variable p2 = src.getParents().get(i);
            if (!p1.getDomain().equals(p2.getDomain()))
                throw new RuntimeException("Invalid sharing: " + p1.getName() + " does not share domain with " + p2.getName());
        }
        this.tieSource = src;
        // need to tie:
        // - count (used during learning)
        // - prior (if applicable)
        // - table (if applicable)
        this.prior = src.prior;
        if (this.table.nParents > 0)
            this.table = src.table.retrofit(this.getParents());
    }
    public BNode getTieSource(){
        return this.tieSource;
    }

	@Override
	public List<Sample> getConditionDataset(int conditionIndex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Distrib getlikelihoodDistrib() {
		// TODO Auto-generated method stub
		return null;
	}


}