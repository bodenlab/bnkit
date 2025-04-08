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

import api.JSONUtils;
import bn.*;
import bn.alg.EM;
import bn.factor.Factor;
import dat.*;
import bn.prob.GaussianDistrib;
import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import bn.factor.Factorize;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.phylo.IdxTree;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;

/**
 * Class for Gaussian Density Table (GDT). This is a table for a continuous
 * variable, with each row specifying a single Gaussian density. The node has
 * one or more enumerable parents.
 *
 * @author m.boden
 */
public class GDT implements BNode, TiedNode<GDT>, Serializable {

    private static final long serialVersionUID = 1L;

    public static final int VARIANCE_UNTIED = 0;
    public static final int VARIANCE_TIED_MAX = 1;
    public static final int VARIANCE_TIED_POOLED = 2;

//    private int tieVariances = VARIANCE_TIED_POOLED;
    private int tieVariances = VARIANCE_UNTIED;

    final private Variable<Continuous> var;
    private GaussianDistrib prior = null;
    private EnumTable<GaussianDistrib> table = null;    // the table with entries for each "component" distribution

    // Parameters used for training. Some of which are allocated prior to training, and then re-used to save time.
    
    private SampleTable<Double> countDouble = null;     // the table that will contain all samples of the type "Double" during learning
    private SampleTable<Distrib> countDistrib = null;   // the table that will contain all samples of the type "Distrib" during learning

    private double[] observed = new double[0];  // observations are stored here AFTER sample collecting, some caching is possible (values and numbers there-of)
    private double[] weight = new double[0];      // probability of component/combination of parent values
    private int[] row = new int[0];             // identifies to which "row" the observation belongs
    private double[] responsibilities;          // assigned to each sample, so that it sums to 1.0 across components


    final private double[] means;       // save the means 
    final private double[] vars; 	    // save the variances
    final private double[] n;           // save the numbers of samples
    
    private boolean relevant = false;
    private GDT tiedMaster = null; // set to another GDT in case the current GDT is "tied" to it (meaning that the master supplies the parameters for inference and learning)

    /**
     * Create a Gaussian density table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public GDT(Variable<Continuous> var, List<EnumVariable> parents) {
        this.var = var;
        int maxrows = 0;
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<>(parents);
                maxrows = this.table.getSize();
            }
        }
        means = new double[maxrows];
        vars = new double[maxrows];
        n = new double[maxrows];
    }

    /**
     * Create a Gaussian density table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public GDT(Variable<Continuous> var, EnumVariable... parents) {
        this(var, EnumVariable.toList(parents));
    }

    /**
     * Create a Gaussian prior for a variable. Note: the prior distribution is
     * not set (@see bn.GDT#put(GaussianDistrib)).
     * @param var variable
     */
    public GDT(Variable<Continuous> var) {
        this.var = var;
        int maxrows = 0;
        means = new double[maxrows];
        vars = new double[maxrows];
        n = new double[maxrows];
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
            throw new RuntimeException("Evaluation of GDT " + this.toString() + " failed since condition was not fully specified: " + e.getMessage());
        }
    }
    
    /**
     * Retrieve the distribution for this node with no parents.
     * @return the distribution of the variable for this root node
     */
    @Override
    public GaussianDistrib getDistrib() {
        return prior;
    }

    /**
     * Make a native Distrib instance out of a collection of samples.
     * @param samples
     * @return an instance of Distrib that can be used to populate this node.
     */
    public GaussianDistrib makeDistrib(Collection<Sample> samples) {
        if (samples == null)
            return new GaussianDistrib(0,0.01);
        double sum = 0;		// we keep track of the sum of observed values
        double tot = 0;		// we keep track of the total of counts
        double[] infobs = new double[samples.size()];
        double[] infprb = new double[samples.size()];
        // go through observed samples
        int j = 0;
        for (Sample sample : samples) {// look at each distribution
            try {
                Double y = (Double) sample.instance;
                infobs[j] = y;        // actual value (or score)
                infprb[j] = sample.prob;            // p(class=key) i.e. the height of the density for this parent config  
                sum += infobs[j] * infprb[j];            // update the numerator of the mean calc
                tot += infprb[j];                   // update the denominator of the mean calc
            } catch (ClassCastException e1) {
                try {
                    Distrib d = (Distrib) sample.instance;
                    infobs[j] = (double) d.sample();        // actual value (or score)
                    infprb[j] = sample.prob;                   // p(class=key) i.e. the height of the density for this parent config  
                    sum += infobs[j] * infprb[j];				// update the numerator of the mean calc
                    tot += infprb[j];					// update the denominator of the mean calc
                } catch (ClassCastException e2) {
                    throw new RuntimeException("Evaluation of GDT distribution failed since sample was of unknown type: " + sample.instance);
                }
            }
            j++; 
        }                
        // calculate mean
        double mean = 0;
        if (tot > 0) 
            mean = sum / tot;
        // now for calculating the variance
        double diff = 0;
        for (int jj = 0; jj < j; jj++)
            diff += (mean - infobs[jj]) * (mean - infobs[jj]) * infprb[jj];
        double variance = Math.max(diff / tot, 0.01);
        return new GaussianDistrib(mean, variance);
    }

    /**
     * Make a FactorTable out of this GDT. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * Marginalization technique requires updating
     *
     * @param relevant only include relevant nodes, with instantiations if available
     * @return the FactorTable created from the GDT, provided instantiations
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
                GaussianDistrib d = table.getValue(index);
                if (d != null) {
                    Object[] key = table.getKey(index);
                    int newcnt = 0;
                    for (int i = 0; i < key.length; i++) {
                        if (searchkey[i] == null) {
                            newkey[newcnt++] = key[i];
                        }
                    }
                    if (varinstance != null) { // the variable for this GDT is instantiated
                        ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this GDT is NOT instantiated...
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
            throw new RuntimeException("GDTs can not be factorised unless it has enumerable parent variables");
        }
    }

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
                fvars.add(myvar);
            }
            Variable[] vars_arr = new Variable[fvars.size()];
            fvars.toArray(vars_arr);
            AbstractFactor ft = new DenseFactor(vars_arr);
            AbstractFactor.FactorFiller ff = ft.getFiller();
            EnumVariable[] evars = ft.getEnumVars(); // the order may have changed
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
                GaussianDistrib d = table.getValue(index);
                if (d != null) { // there is a distribution associated with this entry in the CPT
                    Object[] cptkey = table.getKey(index); // work out the condition for this entry
                    for (int i = 0; i < cptkey.length; i++) {
                        if (xcross[i] != -1)
                            fkey[xcross[i]] = cptkey[i];
                    }
                    if (varinstance != null) { // the variable for this CPT is instantiated
                        if (fkey.length == 0) // and the parents are too
                            ft.setValue(d.get(varinstance));
                        else
                            ff.setValue(fkey, d.get(varinstance));
                    } else { // the variable for this CPT is NOT instantiated so we put it in the JDF
                        if (fkey.length == 0) { // but the parents are instantiated
                            ft.setValue(1.0);
                            ft.setDistrib(myvar, d);
                        } else {
                            ff.setValue(fkey, 1.0);
                            ft.setDistrib(fkey, myvar, d);
                        }
                    }
                } 
            }
            ff.setNormalised(); // normalise Gaussian factors, ensuring that they sum to 1
            ft.setValuesByFiller(ff);
            if (!sumout.isEmpty()) {
                Variable[] sumout_arr = new Variable[sumout.size()];
                sumout.toArray(sumout_arr);
            	ft = Factorize.getMargin(ft, sumout_arr);
            }
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) // instantiated prior is not possible to factorise
                return null;
            throw new RuntimeException("GDTs can not be factorised unless it has enumerable parent variables");
        }
    }
    
    /**
     * Get the conditional probability of the variable (represented by this GDT)
     * when set to a specified value.
     *
     * @param key condition
     * @param value the value of the variable represented by this GDT
     * @return the probability density of the variable set to specified value
     */
    public Double get(Object[] key, Object value) {
        if (key == null) {
            if (prior != null) {
                return prior.get(value);
            }
            return null;
        }
        GaussianDistrib d = table.getValue(key);
        if (d != null) {
            Double p = d.get(value);
            return p;
        }
        return null;
    }

    /**
     * Get the conditional probability of the variable (represented by this GDT)
     * when set to a specified value.
     *
     * @param value the value of the variable represented by this GDT
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
        GaussianDistrib d = table.getValue(key);
        if (d != null) {
            Double p = d.get(value);
            return p;
        }
        return null;
    }

    /**
     * Get the prior probability of the variable (represented by this GDT)
     *
     * @param value the value of the variable represented by this CDT
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
     * Get the conditional probability (density if continuous) of the variable
     * (represented by this CDT).
     *
     * @return the probability
     */
//	public GaussianDistrib get(Object[] key) {
//		return table.getValue(key);
//	}
    
    /**
     * @deprecated Do NOT use, will be removed in the future
     */
    @Override
    public EnumTable getTable() {
        return table;
    }

    /**
     * Set entry (or entries) of the CDT to the specified probability
     * distribution
     *
     * @param key the boolean key (probabilistic condition)
     * @param distr the distribution
     */
    public void put(Object[] key, Distrib distr) {
        table.setValue(key, (GaussianDistrib)distr);
    }

    /**
     * Set entry (or entries) of the CDT to the specified probability
     * distribution
     *
     * @param index the index for the key (probabilistic condition)
     * @param distr the distribution
     */
    public void put(int index, Distrib distr) {
        table.setValue(index, (GaussianDistrib)distr);
    }

    /**
     * Set entry (or entries) of the CDT to the specified probability
     * distribution
     *
     * @param distr the distribution
     * @param key the boolean key (probabilistic condition)
     */
    public void put(Distrib distr, Object... key) {
        table.setValue(key, (GaussianDistrib) distr);
    }

    public void put(Distrib distr) {
        prior = (GaussianDistrib) distr;
    }

    @Override
    public String toString() {
        List<EnumVariable> parents = table.getParents();
        StringBuilder sbuf = new StringBuilder();
        for (int i = 0; i < parents.size(); i++) {
            sbuf.append(parents.get(i).getName()).append(i < parents.size() - 1 ? "," : "");
        }
        return "GDT(" + getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
    }

    protected String formatTitle() {
        return String.format(" %10s", var.getName());
    }

    protected String formatValue(GaussianDistrib x) {
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

    private Double instance = null;

    @Override
    public void setInstance(Object value) {
        try {
            instance = (Double) value;
        } catch (ClassCastException e) {
            try {
                instance = ((Integer) value)/1.0;
            } catch (ClassCastException e2) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public void resetInstance() {
        instance = null;
    }

    @Override
    public Double getInstance() {
        return instance;
    }

    /**
     * Store the continuous observation, which is associated with a specific condition, and a probability.
     * @param key the boolean key (how conditioning variables are set)
     * @param value the value of the conditioned, continuous variable
     * @param prob the expectation
     */
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        if (this.isRoot()) {
            throw new RuntimeException("GDT can not be trained as root");
            // same process as for entries with parents, just a single queue of observations...
        } else {
            GDT master = getMaster();
            if (master == null) { // this GDT is a master
                if (countDouble == null) // create countDouble table if none exists
                    countDouble = new SampleTable<>(this.getParents());
                if (countDistrib == null)
                    countDistrib = new SampleTable<>(this.getParents());
                try {
                    countDouble.count(key, (Double) value, prob);
                } catch (ClassCastException e) {
                    countDistrib.count(key, (Distrib) value, prob);
                }
            } else { // this GDT is slave to a master
                if (master.countDouble == null) // create countDouble table if none exists
                    master.countDouble = new SampleTable<>(this.getParents());
                if (master.countDistrib == null)
                    master.countDistrib = new SampleTable<>(this.getParents());
                try {
                    master.countDouble.count(key, (Double) value, prob);
                } catch (ClassCastException e) {
                    master.countDistrib.count(key, (Distrib) value, prob);
                }
            }
        }
    }
    
    @Override
    public void countInstance(Object[] key, Object value) {
        countInstance(key, value, 1.0);
    }

    /**
     * Set this CDT to be trained when the Bayesian network it is part of is
     * trained. A CDT is trainable (true) by default.
     *
     * @param status true if trainable, false otherwise
     */
    public void setTrainable(boolean status) {
        trainable = status;
    }

    protected boolean trainable = true;

    /**
     * Check if this GDT should be trained or not
     */
    @Override
    public boolean isTrainable() {
        return trainable;
    }

    /**
     * Put random entries in the GDT if not already set.
     * @param seed
     */
    @Override
    public void randomize(long seed) {
        Random rand = new Random(seed);
        if (table == null) {
            if (prior == null)
                prior = new GaussianDistrib(rand.nextGaussian(), rand.nextDouble());
        } else {
            int nrows = table.getSize();
            for (int i = 0; i < nrows; i++) {
                if (!table.hasValue(i))
                    table.setValue(i, new GaussianDistrib(rand.nextGaussian(), rand.nextDouble()));
            }
            // table.display();
        }
    }

    public void randomize(Object[] observations, int seed) {
        // use observations to initialise this GDT
        // we do not have access to parents so this is only a rough estimate
        double min=0, max=0, sum=0, var=0, cnt=0;
        for (int i=0; i<observations.length; i++) {
            if (observations[i] == null)
                continue;

            Double y;
            if (observations[i] instanceof Integer) {
                // catches ClassCastExceptions for integer values
                y = ((Integer) observations[i]).doubleValue();
            } else {
                y = (Double) observations[i];
            }

            sum+=y;
            if (i==0)
                max=min=y;
            else {
                if (y>max)
                    max=y;
                else if (y<min)
                    min=y;
            }
            cnt += 1;
        }
        double mean=sum/cnt;
        for (int i=0; i<observations.length; i++) {
            if (observations[i] == null)
                continue;

            Double y;
            if (observations[i] instanceof Integer) {
                // catches ClassCastExceptions for integer values
                y = ((Integer) observations[i]).doubleValue();
            } else {
                y = (Double) observations[i];
            }
            var+=(mean-y)*(mean-y);
        }
        var/=cnt;
        int nrows= table.getSize();
        if (nrows<2)
            put(new GaussianDistrib(mean, var));
        else {
            double range=max-min;
            double stepsize=range/nrows;
            for (int i=0; i<nrows; i++) {
                put(i, new GaussianDistrib((max-stepsize*i)-stepsize/2, var));
            }
        }
     }

    /**
     * Couple/uncouple variances.
     *     public final int VARIANCE_UNTIED = 0;
     *     public final int VARIANCE_TIED_MAX = 1;
     *     public final int VARIANCE_TIED_POOLED = 2;
     * (see Hastie & Tibshirani).
     *
     * @param status variance update policy
     */
    public void setTieVariances(int status) {
        tieVariances = status;
    }

    public int getNumberObservedDistrib() {
        int number = 0;
        if (countDistrib == null)
            return 0;
        for (List<Sample<Distrib>> samplesDistrib : countDistrib.getTable().getValues()) {
            number += samplesDistrib.size();
        }
        return number;
    }
    
    public int getNumberObservedSample() {
        int number = 0;
        if (countDouble == null)
            return 0;
        for (List<Sample<Double>> samplesDouble : countDouble.getTable().getValues()) {
            number += samplesDouble.size();
        }
        return number;
    }

    /**
     * Weights when this GDT is used as a mixture of Gaussians
     */
    public double[] GMM_WEIGHTS = null;

    /**
     * Find the parameter setting to best reproduce the observed data.
     * Note that this uses observed Double:s which are looked at directly, and 
     * observed Distrib:s which are used to stochastically generate samples. 
     * We cannot use distributions directly since they can be of any kind, as 
     * long as defined over a single continuous variable.
     */
    @Override
    public void maximizeInstance() {
        int nComponents = table.getSize();

        GDT master = getMaster();
        if (master == null) {
            double maxVar = 0;                      // the largest variance of any class

            int nTotal = getNumberObservedSample();
            if (nTotal == 0) // no data, no learning
                return;
            if (observed.length != nTotal) { // check if we need to re-allocate memory; if not we will re-use
                observed = new double[nTotal];
                weight = new double[nTotal];
                row = new int[nTotal]; // the row in the table to which the observation belongs
                responsibilities = new double[nTotal];
            }
            // Go through each possible row, each with a unique combination of parent values (no need to know parent values actually)
            int j = 0;                                          // sample count
            double sum = 0;
            Map<Object, Double> respsum = new HashMap<>();
            for (int index = 0; index < nComponents; index ++) {    // go through all "components" (combinations of parent values)
                Distrib d = this.table.getValue(index);
                List<Sample<Double>> samplesDouble = null;
                n[index] = 0;
                if (countDouble != null)
                    samplesDouble = countDouble.get(index);
                int jStart = j;     // start index for samples in the row
                if (samplesDouble == null)
                    continue;                 // no samples of any kind
                // go through actual values... for this index/component
                if (samplesDouble != null) {
                    n[index] = samplesDouble.size();
                    for (Sample<Double> sample : samplesDouble) {   // look at each entry
                        observed[j] = sample.instance;              // actual value (or score)
                        // weight of the parent
                        weight[j] = sample.prob;                      // p(class=key) i.e. the (normalised) height of the density for this parent config
                        // sum += observed[j] * prob[j];              // update the numerator of the mean calc
                        responsibilities[j] = weight[j] * d.get(observed[j]); // tentative assignment of responsibilities (sample-specific)
                        if (!respsum.containsKey(observed[j])) {
                            respsum.put(observed[j], responsibilities[j]);
                        } else {
                            respsum.put(observed[j], respsum.get(observed[j]) + responsibilities[j]);
                        }
                        // FIXME: responsibilities are not normalised across components
                        sum += responsibilities[j];               // keep sum so we can normalise responsibilities
                        // tot += prob[j];                            // update the denominator of the mean calc
                        row[j] = index;
                        j++;
                    }
                }
            }
            // normalise responsibilities
            for (int i = 0; i < j; i ++) {
                responsibilities[i] /= respsum.get(observed[i]);
            }
            // M-step
            double[] Nk = new double[nComponents];
            GMM_WEIGHTS = new double[nComponents];
            j = 0;
            for (int index = 0; index < nComponents; index ++) {    // go through all "components" (combinations of parent values)
                GaussianDistrib gd = this.table.getValue(index);
                List<Sample<Double>> samplesDouble = null;
                if (countDouble != null)
                    samplesDouble = countDouble.get(index);
                int jStart = j;     // start index for samples in the row
                if (samplesDouble == null)
                    continue;                 // no samples of any kind
                double meansum = 0;
                double varsum = 0;
                for (Sample<Double> sample : samplesDouble) {   // look at each entry
                    Nk[index] += responsibilities[j];
                    meansum += responsibilities[j] * observed[j];
                    varsum += responsibilities[j] * Math.pow(observed[j] - gd.getMean(), 2);
                    j ++;
                }
                // calculate mean
                means[index] = meansum / Nk[index];
                // now for calculating the variance
                vars[index] = varsum / Nk[index];
                if (vars[index] < 0.001) {
                    vars[index] = 0.001;
                }
                if (vars[index] > maxVar) {
                    maxVar = vars[index];
                }
                GMM_WEIGHTS[index] = Nk[index] / samplesDouble.size();
            }
            countDistrib = null;    // reset counts
            countDouble = null;     // reset counts

            if (tieVariances == VARIANCE_UNTIED) { // if we use the individual variances
                for (int i = 0; i < nComponents; i ++) {
                    if (n[i] > 0)
                        this.put(i, new GaussianDistrib(means[i], vars[i]));
                }
            } else { // re-compute variances if they need to be tied
                // there are different ways of dealing with this
                if (tieVariances == VARIANCE_TIED_MAX) {
                    // (1) simply use the max of the existing variances
                    for (int i = 0; i < nComponents; i ++) {
                        if (n[i] > 0)
                            this.put(i, new GaussianDistrib(means[i], maxVar));
                    }
                } else if (tieVariances == VARIANCE_TIED_POOLED) {
                    // (2) use the pooled existing variances (http://en.wikipedia.org/wiki/Pooled_variance)
                    double num = 0.0;
                    double denom = 0.0;
                    for (int i = 0; i < nComponents; i ++) {
                        if (n[i] >= 1) {
                            num += (n[i] - 1) * vars[i];
                            denom += (n[i] - 1);
                        }
                    }
                    for (int i = 0; i < nComponents; i ++) {
                        if (n[i] > 0)
                            this.put(i, new GaussianDistrib(means[i], num / denom));
                    }
                } else {
                    // (3) compute the variance of all the values (Hastie and Tibshirani did this--but I have not had great success with this /MB)
                    throw new RuntimeException("This variant of tied variance is not yet implemented");
                    // ...
                }
            }
        } else { // there is a master
            master.maximizeInstance();
            // the alternative is to wait for the master, but that may not happen if not in the update set
        }
    }

    @Override
    public String getName() {
        return getVariable().getName();
    }

    @Override
    public Variable<Continuous> getVariable() {
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
        return "GDT";
    }

    /**
     * Convert GDT to JSON
     * @return
     */
    @Override
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        JSONArray jidxs = new JSONArray();
        JSONArray jcond = new JSONArray();
        JSONArray jdist = new JSONArray();
        for (int i = 0; i < table.getSize(); i++) {
            GaussianDistrib d = table.getValue(i);
            if (d != null) {
                jidxs.put(i);
                Object[] key = table.getKey(i);
                jcond.put(new JSONArray(key));
                jdist.put(d.toJSONArray());
            }
        }
        json.put("Variable", var.toJSON());
        json.put("Index", jidxs);
        json.put("Condition", jcond);
        json.put("Pr", jdist);
        json.put("Nodetype", getType());
        json.put("TieVariance", this.tieVariances);
        return json;
    }

    /**
     * Recover the (conditioned) GDT from JSON, based on a specified variable, and the parent variables
     * @param json JSON representation
     * @param nodevar the non-enumerable variable defining this node
     * @param parvars  the enumerable variables that are the parents of this node
     * @return new GDT
     * @throws JSONException if the JSON specification is invalid
     */
    public static GDT fromJSON(JSONObject json, Variable nodevar, List<EnumVariable> parvars) throws JSONException {
        GDT gdt = new GDT(nodevar, parvars);
        JSONArray jidxs = json.getJSONArray("Index");
        JSONArray jeds = json.getJSONArray("Pr");
        for (int i = 0; i < jidxs.length(); i++) {
            int idx = jidxs.getInt(i);
            JSONArray jed = jeds.getJSONArray(i);
            gdt.put(idx, GaussianDistrib.fromJSONArray(jed));
        }
        Integer tvar = json.optInt("TieVariance");
        if (tvar != null)
            gdt.setTieVariances(tvar);
        return gdt;
    }

    /**
     * Recover the (conditioned) GDT from JSON, based on a specified variable, and the parent variables
     * @param json JSON representation
     * @param nodevar the non-enumerable variable defining this node
     * @param parvars  the enumerable variables that are the parents of this node
     * @return new GDT
     * @throws JSONException if the JSON specification is invalid
     */
    public static GDT fromJSON(JSONObject json, Variable nodevar, EnumVariable ... parvars) throws JSONException {
        return fromJSON(json, nodevar, EnumVariable.toList(parvars));
    }

    /**
     * Recover the (conditioned) GDT from JSON, based on the parent variables
     * @param json JSON representation
     * @param parvars  the enumerable variables that are the parents of this node
     * @return new GDT
     * @throws JSONException if the JSON specification is invalid
     */
    public static GDT fromJSON(JSONObject json, EnumVariable[] parvars) throws JSONException {
        JSONObject jvar = json.getJSONObject("Variable");
        Variable nvar = Variable.fromJSON(jvar);
        GDT gdt = new GDT(nvar, parvars);
        JSONArray jidxs = json.getJSONArray("Index");
        JSONArray jeds = json.getJSONArray("Pr");
        Domain dom = nvar.getDomain();
        for (int i = 0; i < jidxs.length(); i++) {
            int idx = jidxs.getInt(i);
            JSONArray jed = jeds.getJSONArray(i);
            gdt.put(idx, GaussianDistrib.fromJSONArray(jed));
        }
        Integer tvar = json.optInt("TieVariance");
        if (tvar != null)
            gdt.setTieVariances(tvar);
        return gdt;
    }

    @Override
    public String getStateAsText() {
        StringBuilder sbuf = new StringBuilder("\n");
        if (isRoot()) {
            GaussianDistrib d = prior;
            if (d != null) {
                sbuf.append("").append(d.getMean()).append(", ").append(d.getVariance()).append(";\n");
            }
        } else {
            for (int i = 0; i < table.getSize(); i++) {
                GaussianDistrib d = table.getValue(i);
                if (d != null) {
                    sbuf.append(i).append(": ");	// use index as key because values above can be of different non-printable types
                    sbuf.append("").append(d.getMean()).append(", ").append(d.getVariance()).append("; (");
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
        if (isRoot()) {
            String[] line = dump.split(";");
            if (line.length >= 1) {
                String[] y = line[0].split(",");
                if (y.length == 2) {
                    double[] distrib = new double[y.length];
                    try {
                        for (int i = 0; i < distrib.length; i++) {
                            distrib[i] = Double.parseDouble(y[i]);
                        }
                    } catch (NumberFormatException e) {
                        e.printStackTrace();
                        return false;
                    }
                    this.put(new GaussianDistrib(distrib[0], distrib[1]));
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
                            if (y.length == 2) {
                                double[] distrib = new double[y.length];
                                try {
                                    for (int i = 0; i < distrib.length; i++) {
                                        distrib[i] = Double.parseDouble(y[i]);
                                    }
                                } catch (NumberFormatException e) {
                                    e.printStackTrace();
                                    return false;
                                }
                                this.put(table.getKey(index), new GaussianDistrib(distrib[0], distrib[1]));
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

    public boolean isRelevant() {
        return relevant;
    }

    public void setRelevant(boolean relevant) {
        this.relevant = relevant;
    }
    
    /**
     * @param args
     */
    public static void main0(String[] args) {
        Variable<Continuous> v1 = Predef.Real();
        EnumVariable v2 = Predef.Boolean();
        EnumVariable v3 = Predef.Boolean();

        GDT gdt1 = new GDT(v1, new EnumVariable[]{v2, v3});
        gdt1.put(new Object[]{true, false}, new GaussianDistrib(0.0, 1.0));
        gdt1.print();
        System.out.println(gdt1.get(new Object[]{true, false}, 0.5));
    }
    
    /**
     * get the data from both distribution and the double sample
     */
	@Override
	public List<Sample> getConditionDataset(int conditionIndex) {
		int nSample = 20; 
		List<Sample> result = new LinkedList<Sample>();
		// if there is countDistrib, then get..
        if (countDistrib != null) {
        	List<Sample<Distrib>> samplesDistrib = countDistrib.get(conditionIndex);
        	for(Sample<Distrib> sample: samplesDistrib) {
        		for(int j = 0; j < nSample; j++) {
        			Double r = (Double) sample.instance.sample();
        			Sample<Double> newSample = new Sample<Double>(r, 1.0);
        			result.add(newSample);
        		}
        	}
        }
        	
        // if there is countDouble, then get..
        if (countDouble != null) {
        	List<Sample<Double>> samplesDouble = countDouble.get(conditionIndex);
        	for(Sample<Double> s: samplesDouble) {
        		Sample newSample = new Sample(s.instance, s.prob);
        		result.add(newSample);
        	}
        }
       
        return result;
        
	}
	
	/**
	 * pay attention here, we should note that if we assume that
	 * our variable to be fixed, our prior is Guassian distribution 
	 * otherwise, our prior is Gaumma distribution.
	 */
	@Override
	public Distrib getlikelihoodDistrib() {
		return new GaussianDistrib(0.0, 1.0);
	}

    public GDT getMaster(){
        return this.tiedMaster;
    }


    /**
     * Tie parameters from training for this CPT to those of another CPT.
     * Variables should be separate, but they are required to (1) be of the same type/domain, and (2) be listed in the same order.
     * @param master the CPT from which parameters will be copied and held fixed.
     * @return true if successful, false otherwise
     */
    @Override
    public boolean tieTo(GDT master) {
        if (master.getMaster() != null)
            return false;
        if (!this.var.getDomain().equals(master.getVariable().getDomain()))
            throw new RuntimeException("Invalid sharing: " + var.getName() + " does not share domain with " + master.getVariable().getName());
        if (this.getParents() != null && master.getParents() != null) {
            if (this.getParents().size() != master.getParents().size())
                throw new RuntimeException("Invalid sharing: " + var.getName() + " has different number of parents from " + master.getVariable().getName());
            for (int i = 0; i < this.getParents().size(); i++) {
                Variable p1 = this.getParents().get(i);
                Variable p2 = master.getParents().get(i);
                if (!p1.getDomain().equals(p2.getDomain()))
                    throw new RuntimeException("Invalid sharing: " + p1.getName() + " does not share domain with " + p2.getName());
            }
        }
        this.tiedMaster = master;
        // need to tie:
        // - prior (if applicable)
        // - table (if applicable)
        this.prior = master.prior; // copy the reference to the value
        if (this.getParents() != null) {
            // we assume that parents are in the same order, belong to the same domains etc (as checked above)
            this.table.setMapRef(master.table.getMapRef()); // copy the reference to the actual "master" table (not a copy of)
//            this.table = master.table; // copy the reference to the actual "master" table (not a copy of)
            // removed "retro-fitted" table as this makes a shallow copy of the master table
            // this.table = src.table.retrofit(this.getParents());
        }
        // also need to tie:
        // - count (used during learning, i.e. when "counting", so this happens in countInstance, not here
        return true;
    }

    /**
     * Train a Gaussian mixture model using a dataset
     *
     * example
     * {"Items":["WT14","WT04","WT48","WT49","WT38","WT16","WT06","WT07","WT51","WT40","WT31","WT43"],
     * "Features":["Tm"],
     * "Data":[[[null],[null],[null],[0.49],[0.39],[null],[0.5],[0.49],[0.46],[0.39],[0.41],[0.43]]]}
     *
     * @param ds
     * @return
     */
    public static GDT trainGMM(JSONUtils.DataSet ds, int nComponents, int nRounds, int seed) {
        EnumVariable X = Predef.Number(nComponents, "Component");
        String[] features = ds.getFeatures();
        Variable<Continuous> Y = Predef.Real(features[0]);
        GDT gdt = new GDT(Y, X);
        CPT cpt = new CPT(X);
        gdt.randomize(seed);
        cpt.randomize(seed + 1);
        Object[][] data = ds.getItemisedData()[0];
        Double max = null;
        Double min = null;
        for (Object[] val : data) {
            if (max == null)
                max = (Double) val[0];
            if (min == null)
                min = (Double) val[0];
            if (val[0] != null && max !=null)
                max = ((Double) val[0] > max) ? (Double) val[0] : max;
            if (val[0] != null && min !=null)
                min = ((Double) val[0] < min) ? (Double) val[0] : min;
        }
        // assume uniform distribution
        if (max == min || max == null || min == null)
            throw new RuntimeException("Dataset is invalid");
        double step = (max - min) / (nComponents);
        double base = min + step / 2;
        double nData = data.length;
        double[] weights = new double[nComponents];
        int j = 0;
        for (Object component : cpt.getDistrib().getDomain().getValues()) {
            weights[j ++] = (double)nData / (double)nComponents;
            GaussianDistrib gd = (GaussianDistrib) gdt.getDistrib(new Object[]{component});
            gd.setMean(base);
            gd.setVariance(step);
            base += step;
        }
        gdt.setTieVariances(gdt.VARIANCE_UNTIED);
        //System.out.println("Before...");
        //System.out.println("Mixing variable: " + cpt.getDistrib());
        //System.out.println("GDT: ");
        //gdt.print();
        cpt.getDistrib().set(weights);
        for (int i = 0; i < 100; i ++) {
            for (int d = 0; d < data.length; d++) {
                if (data[d][0] == null)
                    continue;
                for (Object component : cpt.getDistrib().getDomain().getValues()) {
                    double weight = cpt.getDistrib().get(component);
                    gdt.countInstance(new Object[] {component}, data[d][0], weight);
                }
            }
            gdt.maximizeInstance();
            cpt.getDistrib().set(gdt.GMM_WEIGHTS);
        }
        System.out.println("After...");
        System.out.println("Mixing variable: " + cpt.getDistrib());
        System.out.println("GDT: ");
        gdt.print();
        return gdt;
    }


    public static void main(String[] args) {
        String folder = "/Users/mikael/simhome/ASR/ReconMode/";
        try {
            Set<String> extants = new HashSet<>();
            Set<String> ancestors = new HashSet<>();
            TSVFile names = new TSVFile(folder + "ERED_names.tsv", true);
            Map<String, String> namemap = new HashMap<>();
            for (Object[] row : names.getRows())
                namemap.put((String) row[1], (String) row[4]);
            TSVFile wt = new TSVFile(folder + "s2-2.tsv", true);
            Map<String, Object[]> propmap = new HashMap<>();
            for (Object[] row : wt.getRows()) {
                propmap.put((String) row[0], row);
                extants.add((String) row[0]);
            }
            TSVFile anc = new TSVFile(folder + "s2-1.tsv", true);
            for (Object[] row : anc.getRows()) {
                propmap.put((String) row[0], row);
                ancestors.add((String) row[0]);
            }
            String[] extants_arr = new String[extants.size()];
            String[] ancestors_arr = new String[ancestors.size()];
            extants.toArray(extants_arr);
            ancestors.toArray(ancestors_arr);
            IdxTree tree = Newick.load(folder + "lewis_tree.nwk");

            int NPROTEINS = 40;
            int NCOMPONENTS = 3;
            int[] percent_ancestors = new int[]{0, 20, 40, 60, 80, 100};
            for (int percent : percent_ancestors) {
                int NANCESTORS = NPROTEINS * percent / 100;
                int NEXTANTS = NPROTEINS - NANCESTORS;
                for (int SEED = 0; SEED < 5; SEED += 1) {
                    Random rand = new Random(SEED);
                    Set<String> select = new HashSet<>();
                    while (select.size() < NEXTANTS)
                        select.add(extants_arr[rand.nextInt(extants_arr.length)]);
                    while (select.size() < NPROTEINS)
                        select.add(ancestors_arr[rand.nextInt(ancestors_arr.length)]);
                    String[] select_arr = new String[select.size()];
                    select.toArray(select_arr);
                    Object[][][] data = new Object[1][select_arr.length][1];
                    for (int i = 0; i < select_arr.length; i++) {
                        data[0][i][0] = propmap.get(select_arr[i])[1];
                        if (data[0][i][0] != null)
                            data[0][i][0] = (((double)((Integer)data[0][i][0]) + i/100.0));

//                        data[0][i][0] = ((double) ((Integer) data[0][i][0]) / (10.0 + Math.abs(rand.nextGaussian())));
                    }
                    JSONUtils.DataSet ds = new JSONUtils.DataSet(select_arr, new String[]{"Tm"}, data);
                    GDT gdt = GDT.trainGMM(ds, NCOMPONENTS, 10, SEED);
                    System.out.println(ds.toJSON());
                    System.out.println(gdt.toJSON());
                    BNet bn = new BNet();
                    //bn.add(new BNode[] {gdt, cpt});
                }
            }


        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
