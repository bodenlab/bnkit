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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Class for Gaussian Density Table (GDT). This is a table for a continuous
 * variable, with each row specifying a single Gaussian density. The node has
 * one or more enumerable parents.
 *
 * @author m.boden
 */
public class GDT implements BNode, Serializable {

    private static final long serialVersionUID = 1L;
    private boolean tieVariances = false;
    final private Variable<Continuous> var;
    private GaussianDistrib prior = null;
    private EnumTable<GaussianDistrib> table = null;
    private SampleTable<Double> count = null; // the table that will contain all samples during learning

    /**
     * Create a Gaussian density table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public GDT(Variable<Continuous> var, List<EnumVariable> parents) {
        this.var = var;
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<GaussianDistrib>(parents);
            }
        }
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
    }

    /**
     * Make a FactorTable out of this GDT. If a variable is instantiated it will
     * be factored out.
     *
     * @param bn the BNet instance that can be used to check the status of nodes
     * so that factoring can be done (instantiation of variables are done for a
     * BNet node).
     */
    @Override
    public FactorTable makeFactor(BNet bn) {
        // TODO: fix factorisation of un-instantiated continuous variables, perhaps by discretizing...
        List<EnumVariable> vars_old = this.getParents();
        Object varinstance = null;
        BNode cnode = bn.getNode(var);
        if (cnode != null) {
            varinstance = cnode.getInstance();
        }
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<EnumVariable> vars_new = new ArrayList<EnumVariable>(vars_old.size() + 1);
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
                throw new RuntimeException("GDTs can not be factorised if continuous variable is un-instantiated");
            }
            FactorTable ft = new FactorTable(vars_new);
            if (varinstance != null) {
                ft.evidenced = true;
            }
            int[] indices = table.getIndices(searchkey);
            Object[] newkey = new Object[vars_new.size()];
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
                        ft.addValue(newkey, d.get(varinstance));
                    } else { // the variable for this GDT is NOT instantiated...
                        throw new RuntimeException("GDTs can not be factorised if continuous variable is un-instantiated");
                    }
                } else { // this entry is null
                    //
                }
            }
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) // instantiated prior is not possible to factorise
            {
                return null;
            }
            throw new RuntimeException("GDTs can not be factorised if continuous variable is un-instantiated");
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
    public Double get(Object value) {
        if (prior != null) {
            return prior.get(value);
        }
        return null;
    }

    /**
     * Get the prior probability of the variable (represented by this CDT)
     *
     * @return the density
     */
    public GaussianDistrib get() {
        return prior;
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
    public EnumTable getTable() {
        return table;
    }

    public GaussianDistrib getDistrib() {
        return prior;
    }

    /**
     * Set entry (or entries) of the CDT to the specified probability
     * distribution
     *
     * @param key the boolean key (probabilistic condition)
     * @param distr the distribution
     */
    public void put(Object[] key, GaussianDistrib distr) {
        table.setValue(key, distr);
    }

    /**
     * Set entry (or entries) of the CDT to the specified probability
     * distribution
     *
     * @param distr the distribution
     * @param key the boolean key (probabilistic condition)
     */
    public void put(GaussianDistrib distr, Object... key) {
        table.setValue(key, distr);
    }

    public void put(GaussianDistrib distr) {
        prior = distr;
    }

    public String toString() {
        List<EnumVariable> parents = table.getParents();
        StringBuffer sbuf = new StringBuffer();
        for (int i = 0; i < parents.size(); i++) {
            sbuf.append(parents.get(i).getName() + (i < parents.size() - 1 ? "," : ""));
        }
        return "GDT(" + getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
    }

    protected String formatTitle() {
        return String.format(" %10s", var.getName());
    }

    protected String formatValue(GaussianDistrib x) {
        return String.format("<%s>", x.toString());
    }

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

    public void setInstance(Object value) {
        try {
            instance = (Double) value;
        } catch (ClassCastException e) {
            e.printStackTrace();
        }
    }

    public void resetInstance() {
        instance = null;
    }

    public Double getInstance() {
        return instance;
    }

    public void countInstance(Object[] key, Object value, Double prob) {
        if (this.isRoot()) {
            throw new RuntimeException("GDT can not be trained as root, should be implemented soon...");
            // same process as for extries with parents, just a single queue of observations...
        } else {
            if (count == null) // create count table if none exists
            {
                count = new SampleTable<Double>(this.getParents());
            }
            count.count(key, prob);
        }
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
    public boolean isTrainable() {
        return trainable;
    }

    /**
     * Put random entries in the GDT if not already set.
     */
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
        }
    }

    /*	public void randomize(Object[] observations, int seed) {
     // use observations to initialise this GDT
     // we do not have access to parents so this is only a rough estimate
     double min=0, max=0, sum=0, var=0;
     for (int i=0; i<observations.length; i++) {
     Double y=(Double)observations[i];
     sum+=y;
     if (i==0)
     max=min=y;
     else {
     if (y>max)
     max=y;
     else if (y<min)
     min=y;
     }
     }
     double mean=sum/observations.length;
     for (int i=0; i<observations.length; i++) {
     Double y=(Double)observations[i];
     var+=(mean-y)*(mean-y);
     }
     var/=observations.length;
     int nrows=getMaxIndex();
     if (nrows<2)
     put(new Gaussian(mean, var));
     else {
     double range=max-min;
     double stepsize=range/nrows;
     for (int i=0; i<nrows; i++) {
     put(entry(i), new Gaussian((max-stepsize*i)-stepsize/2, var));
     }
     }
     }
     */
    /**
     * Couple each entry in this GDT so that the variances are equal (Hastie &
     * Tibshirani).
     *
     * @param status tie variances if true, individual variances if false
     * (default)
     */
    public void setTieVariances(boolean status) {
        tieVariances = status;
    }

    /**
     * Use the maximum of any variance as the global variance (set to false if
     * the variance of all values should be used instead)
     */
    public static boolean USE_MAX_VARIANCE = true;

    public void maximizeInstance() {
        int maxrows = table.getSize();
        double maxVar = 0;  				// the largest variance of any class  
        double[] means = new double[maxrows]; // save the means 
        double[] vars = new double[maxrows]; 	// save the variances
        double middleMean = 0; 				// the mean of all values
        double middleVar = 0;					// the variance of all values
        double middleTot = 0;					// the sum of count for all parent configs
        // go through each possible setting of the parent nodes
        for (Map.Entry<Integer, List<SampleTable<Double>.Sample>> entry : count.table.getMapEntries()) {
            int index = entry.getKey().intValue();
            Object[] key = count.table.getKey(index); // the values of the parent nodes, same as for CPT
            List<SampleTable<Double>.Sample> samples = entry.getValue();
            double sum = 0;		// we keep track of the sum of observed values for a specific parent config weighted by count, e.g. 4x0.5 + 3x1.3 + ... 
            double tot = 0;		// we keep track of the total of counts for a specific parent config so we can compute the mean of values, e.g. there are 23 counted for parent is "true"
            double[] y = new double[samples.size()];	// the observations, i.e. observed values or scores
            double[] p = new double[y.length]; 		// how often do we see this score GIVEN the parents (an absolute count)
            int j = 0;
            for (SampleTable<Double>.Sample sample : samples) {	// look at each value entry
                y[j] = sample.instance.doubleValue();				// actual value (or score)
                p[j] = sample.prob;			// p(class=key) i.e. the height of the density for this parent config  
                sum += y[j] * p[j];					// update the numerator of the mean calc
                tot += p[j];						// update the denominator of the mean calc
                j++;
            }
            means[index] = sum / tot;
            double diff = 0;
            for (int jj = 0; jj < y.length; jj++) {
                diff += (means[index] - y[jj]) * (means[index] - y[jj]) * p[jj];
            }
            vars[index] = diff / tot;
            if (vars[index] < 0.01) {
                vars[index] = 0.01;
            }
            if (vars[index] > maxVar) {
                maxVar = vars[index];
            }
            // note the same key/index for both the CPT and the Sample table
            this.put(key, new GaussianDistrib(means[index], vars[index]));
            middleTot += tot;
            middleMean += sum;
        }
        middleMean /= middleTot;

        if (tieVariances) { // re-compute variances if they need to be tied
            // there are different ways of dealing with this
            if (USE_MAX_VARIANCE) {
                // (1) simply use the max of the existing variances
                for (Map.Entry<Integer, GaussianDistrib> entry : table.getMapEntries()) {
                    int index = entry.getKey().intValue();
                    Object[] key = count.table.getKey(index); // the values of the parent nodes
                    this.put(key, new GaussianDistrib(means[index], maxVar));
                }
            } else {
                // (2) compute the variance of all the values (Hastie and Tibshirani did this--but I have not had great success with this /MB)
                throw new RuntimeException("This variant of tied variance is not yet implemented");
                // ...
            }
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

    @Override
    public String getStateAsText() {
        StringBuffer sbuf = new StringBuffer("\n");
        if (isRoot()) {
            GaussianDistrib d = prior;
            if (d != null) {
                sbuf.append("" + d.getMean() + ", " + d.getVariance() + ";\n");
            }
        } else {
            for (int i = 0; i < table.getSize(); i++) {
                GaussianDistrib d = table.map.get(new Integer(i));
                if (d != null) {
                    sbuf.append(i + ": ");	// use index as key because values above can be of different non-printable types
                    sbuf.append("" + d.getMean() + ", " + d.getVariance() + "; (");
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

    /**
     * @param args
     */
    public static void main(String[] args) {
        Variable<Continuous> v1 = Predef.Real();
        EnumVariable v2 = Predef.Boolean();
        EnumVariable v3 = Predef.Boolean();

        GDT gdt1 = new GDT(v1, new EnumVariable[]{v2, v3});
        gdt1.put(new Object[]{true, false}, new GaussianDistrib(0.0, 1.0));
        gdt1.print();
        System.out.println(gdt1.get(new Object[]{true, false}, 0.5));
    }

}
