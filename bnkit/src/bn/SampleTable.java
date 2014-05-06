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
import java.util.Map.Entry;

/**
 * A table of entries of Enumerable variables each associated with a list of raw
 * samples. Primary purpose is to support learning algorithms. Currently, each
 * sample is a value, count pair. An alternative design is to store samples more
 * efficiently by a value -> count map.
 *
 * @author m.boden
 */
public class SampleTable<T> implements Serializable {

    private static final long serialVersionUID = 1L;

    public class Sample {

        public final T instance;
        public final double prob;

        Sample(T instance, double prob) {
            this.instance = instance;
            this.prob = prob;
        }
    }

    public final EnumTable<List<Sample>> table; // table of counts

    public SampleTable(EnumVariable[] variables) {
        List<EnumVariable> list = new ArrayList<EnumVariable>(variables.length);
        for (EnumVariable var : variables) {
            list.add(var);
        }
        table = new EnumTable<List<Sample>>(list);
    }

    public SampleTable(Collection<EnumVariable> variables) {
        table = new EnumTable<List<Sample>>(variables);
    }

    public List<Sample> get(Object[] key) {
        int index = table.getIndex(key);
        List<Sample> samples = table.getValue(index);
        return samples;
    }

    /**
     * Set the complete sequence of observations, erasing any previous.
     *
     * @param key the condition under which the observations are valid
     * @param observations list of samples
     */
    public void put(Object[] key, List<Sample> observations) {
        table.setValue(key, observations);
    }

    /**
     * Make one valid observation of value, under a specified condition (key),
     * associated with a probability (or count).
     *
     * @param key the condition
     * @param value the observed value
     * @param count the weight of the observation (usually a probability)
     */
    synchronized public void count(Object[] key, T value, double count) {
        int index = table.getIndex(key);
        List<Sample> samples = table.getValue(index);
        if (samples == null) {
            samples = new ArrayList<Sample>();
            samples.add(new Sample(value, count));
            table.setValue(key, samples);
        } else {
            samples.add(new Sample(value, count));
        }
    }

    public void count(Object[] key, T value) {
        count(key, value, 1.0);
    }

    /**
     * Create distribution of samples collected in table
     * @param qTab query table containing sample table and query node
     */
    public void maximizeInstance(QueryTable qTab) {
        int maxrows = qTab.getQuery().getTable().getSize();
        double maxVar = 0;  				// the largest variance of any class  
        double[] means = new double[maxrows]; // save the means 
        double[] vars = new double[maxrows]; 	// save the variances
        double middleMean = 0; 				// the mean of all values
        double middleVar = 0;					// the variance of all values
        double middleTot = 0;					// the sum of count for all parent configs

        // go through each possible setting of the parent nodes
        for (Entry<Integer, List<Sample>> entry : table.getMapEntries()) {
            int index = entry.getKey().intValue();
            Object[] key = table.getKey(index); // the values of the parent nodes, same as for CPT
            List<Sample> samples = entry.getValue();
            double sum = 0;		// we keep track of the sum of observed values for a specific parent config weighted by count, e.g. 4x0.5 + 3x1.3 + ... 
            double tot = 0;		// we keep track of the total of counts for a specific parent config so we can compute the mean of values, e.g. there are 23 counted for parent is "true"
            double[] y = new double[samples.size()];	// the observations, i.e. observed values or scores
            double[] p = new double[y.length]; 		// how often do we see this score GIVEN the parents (an absolute count)
            int j = 0;
            for (Sample sample : samples) {	// look at each value entry
                y[j] = (Double)sample.instance;				// actual value (or score)
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
//            this.put(key, new GaussianDistrib(means[index], vars[index]));
            middleTot += tot;
            middleMean += sum;
            GaussianDistrib d = new GaussianDistrib(means[index], vars[index]);
//            qTab.getNonEnumTable().setValue(key,  d);
            qTab.getQuery().getTable().setValue(key, d);
        }
    }
    
    public void display() {
        table.display();
    }

    public String toString() {
        String[] parents = table.getLabels();
        StringBuffer sbuf = new StringBuffer();
        for (int i = 0; i < parents.length; i++) {
            sbuf.append(parents[i] + (i < parents.length - 1 ? "," : ""));
        }
        return "SampleTable(" + sbuf.toString() + ")";
    }

}
