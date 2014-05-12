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
 * A table with counts for entries of Enumerable variables.
 * @see JPT for a joint probability table
 * @author m.boden
 */
public class CountTable implements Serializable {

    private static final long serialVersionUID = 1L;
    private boolean totalNeedsUpdate = true;
    private double totalCount = 0.0;

    public final EnumTable<Double> table; // table of counts

    public CountTable(EnumVariable[] variables) {
        List<EnumVariable> list = new ArrayList<EnumVariable>(variables.length);
        for (EnumVariable var : variables) {
            list.add(var);
        }
        table = new EnumTable<Double>(list);
    }

    public CountTable(Collection<EnumVariable> variables) {
        table = new EnumTable<Double>(variables);
    }

    public double get(Object[] key) {
        int index = table.getIndex(key);
        Double cnt = table.getValue(index);
        if (cnt == null) {
            return 0;
        } else {
            return cnt.doubleValue();
        }
    }

    public double get(int key_index) {
        Double cnt = table.getValue(key_index);
        if (cnt == null) {
            return 0;
        } else {
            return cnt.doubleValue();
        }
    }

    public void put(Object[] key, double count) {
        table.setValue(key, count);
        this.totalNeedsUpdate = true;
    }

    public void put(int key_index, double count) {
        table.setValue(key_index, count);
        this.totalNeedsUpdate = true;
    }

    synchronized public void count(Object[] key, double count) {
        int index = table.getIndex(key);
        Double cnt = table.getValue(index);
        if (cnt == null) {
            table.setValue(key, count);
        } else {
            table.setValue(key, cnt.doubleValue() + count);
        }
        this.totalNeedsUpdate = true;
    }

    public void count(Object[] key) {
        count(key, 1.0);
        this.totalNeedsUpdate = true;
    }
    
    /**
     * Take stock of all observations counted via
     * {@link bn.CountTable#count(Object[])}, ie implement the
     * M-step locally.
     * 
     * @param qTab query table containing sample table and query node
     */
    public void maximizeInstance(QueryTable qTab) {
    	BNode query = qTab.getQuery();
    	EnumVariable var = (EnumVariable)query.getVariable();
    	//FIXME - check parents in enumTable!!
    	if (query.getTable() != null){
    		// add the counts to the CPT
    		for (Map.Entry<Integer, Double> entry : table.getMapEntries()) {
    			double nobserv = entry.getValue().doubleValue();
    			Object[] cntkey = table.getKey(entry.getKey().intValue());
    			Object[] cptkey = new Object[cntkey.length - 1];
    			for (int i = 0; i < cptkey.length; i++) {
    				cptkey[i] = cntkey[i + 1];
    			}
    			EnumDistrib d = (EnumDistrib)query.getTable().getValue(cptkey);
    			if (d == null) {
    				d = new EnumDistrib(var.getDomain());
    				d.set(cntkey[0], nobserv);
    				//directly alters the CPT of the node - changed in future versions
    				query.getTable().setValue(cptkey,  d);
    			} else {
    				d.set(cntkey[0], nobserv);
    				//directly alters the CPT of the node - changed in future versions
    				query.getTable().setValue(cptkey, d);
    			}
    		} // normalisation happens internally when values are required			
    	} else { // there are no parents
    		Object[] cntkey = new Object[1];
    		double[] cnts = new double[var.size()];
    		for (int i = 0; i < var.size(); i++) {
    			cntkey[0] = var.getDomain().get(i);
    			cnts[i] = this.get(cntkey);
    		}
    		EnumDistrib prior = (EnumDistrib) query.getDistrib();
    		prior = new EnumDistrib(var.getDomain(), cnts);	// EnumDistrib normalises the counts internally
    	}
    	//RESET COUNTS?
    }

    /**
     * Calculate the sum of all entries.
     * @return summed count
     */
    public double getTotal() {
        if (totalNeedsUpdate) {
            for (Double val : table.getValues()) {
                totalCount += val.doubleValue();
            }
            this.totalNeedsUpdate = false;
        }
        return totalCount;
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
        return "CountTable(" + sbuf.toString() + ")";
    }

}
