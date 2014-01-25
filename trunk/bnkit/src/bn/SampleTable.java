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
