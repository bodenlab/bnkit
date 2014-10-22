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

import dat.EnumVariable;
import dat.EnumTable;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * A table of entries of Enumerable variables each associated with a list of raw
 * samples. Primary purpose is to support learning algorithms. Currently, each
 * sample is a value, count pair. An alternative design is to store samples more
 * efficiently by a value -> count map.
 *
 * TODO: Implement more efficient data structure for cases when number of samples is known a priori.
 * 
 * @author m.boden
 * @param <T> The class used for samples, e.g. Double or Distrib (see GDT for an example)
 */


public class SampleTable<T> implements Serializable {

    
    private static final long serialVersionUID = 1L;
    protected final EnumTable<List<Sample<T>>> table; // table of counts
    protected final List<Sample<T>> list; // counts if no enumerable variables as key
    
    public SampleTable(EnumVariable[] variables) {
        if (variables == null) {
            table = null;
            list = new ArrayList<>();
        } else if (variables.length == 0) {
            table = null;
            list = new ArrayList<>(); 
        } else {
            List<EnumVariable> varlist = new ArrayList<>(variables.length);
            varlist.addAll(Arrays.asList(variables));
            table = new EnumTable<>(varlist);
            this.list = null;
        }
    }

    public SampleTable(Collection<EnumVariable> variables) {
        if (variables == null) {
            table = null;
            list = new ArrayList<>(); 
        } else if (variables.isEmpty()) {
            table = null;
            list = new ArrayList<>(); 
        } else {
            table = new EnumTable<>(variables);
            list = null;
        }
    }

    public SampleTable() {
        table = null;
        list = new ArrayList<>(); 
    }
    
    public SampleTable(Collection<EnumVariable> variables, int nSamples) {
        if (variables == null) {
            table = null;
            list = new ArrayList<>(); 
        } else if (variables.isEmpty()) {
            table = null;
            list = new ArrayList<>(); 
        } else {
            table = new EnumTable<>(variables);
            list = null;
        }
    }
    
    public List<Sample<T>> get(Object[] key) {
        if (table == null) 
            throw new RuntimeException("Invalid call to SampleTable: has no key, but one given");
        int index = table.getIndex(key);
        return get(index);
    }

    public List<Sample<T>> get(int index) {
        if (table == null) 
            throw new RuntimeException("Invalid call to SampleTable: has no key, but one given");
        List<Sample<T>> samples = table.getValue(index);
        return samples;
    }

    public List<Sample<T>> get() {
        if (table == null) 
            return list;
        throw new RuntimeException("Invalid call to SampleTable: has key, but none given");
    }
    
    public List<Sample<T>> getAll(Object[] key) {
        if (table == null) 
            throw new RuntimeException("Invalid call to SampleTable: has no key, but one given");
        int index = table.getIndex(key);
        return getAll(index);
    }
    
    public List<Sample<T>> getAll(int index) {
        if (table == null) 
            throw new RuntimeException("Invalid call to SampleTable: has no key, but one given");
        return table.getValue(index);
    }
    
    public List<Sample<T>> getAll() {
        if (table == null)
            return list;
        throw new RuntimeException("Invalid call to SampleTable: has key, but none given");
    }

    public EnumTable<List<Sample<T>>> getTable() {
        return table;
    }
    
    /**
     * Set the complete sequence of observations, erasing any previous.
     *
     * @param key the condition under which the observations are valid
     * @param observations list of samples
     */
    public void put(Object[] key, List<Sample<T>> observations) {
        table.setValue(key, observations);
    }

    /**
     * Set the complete sequence of observations, erasing any previous.
     *
     * @param observations list of samples
     */
    public void put(List<Sample<T>> observations) {
        if (table == null) {
            list.clear();
            list.addAll(observations);
        }
        throw new RuntimeException("Invalid call to SampleTable: has key, but none given");
    }

    /**
     * Make one valid observation of value, under a specified condition (key),
     * associated with a probability (or count).
     *
     * @param key_index the condition encoded by a key index
     * @param value the observed value
     * @param count the weight of the observation (usually a probability)
     */
    synchronized public void count(int key_index, T value, double count) {
        List<Sample<T>> samples = table.getValue(key_index);
        if (samples == null) {
            samples = new ArrayList<>();
            samples.add(new Sample(value, count));
            table.setValue(key_index, samples);
        } else {
            samples.add(new Sample(value, count));
        }
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
        count(index, value, count);
    }

    public void count(Object[] key, T value) {
        count(key, value, 1.0);
    }
    
    public void count(int key_index, T value) {
        count(key_index, value, 1.0);
    }
    
    /**
     * Make one valid observation of value, under no condition,
     * associated with a probability (or count).
     *
     * @param value the observed value
     * @param count the weight of the observation (usually a probability)
     */
    synchronized public void count(T value, double count) {
        List<Sample<T>> samples = list;
        if (list != null) 
            list.add(new Sample(value, count));
        else
            throw new RuntimeException("Invalid call to SampleTable: has key, but none given");
    }

    public void count(T value) {
        count(value, 1.0);
    }
    
    public boolean isEmpty() {
        if (list != null)
            return list.isEmpty();
        else
            return table.isEmpty();
    }
    
    public void setEmpty() {
        if (list != null)
            list.clear();
        else
            table.setEmpty();
    }
    
    public void display() {
        if (table != null)
            table.display();
        else
            System.out.println(this);
    }

    @Override
    public String toString() {
        if (table != null) {
            String[] parents = table.getLabels();
            StringBuilder sbuf = new StringBuilder();
            for (int i = 0; i < parents.length; i++) {
                sbuf.append(parents[i]).append(i < parents.length - 1 ? "," : "");
            }
            return "SampleTable(" + sbuf.toString() + ")";
        } else {
            return "SampleTable()";
        }
    }

}

