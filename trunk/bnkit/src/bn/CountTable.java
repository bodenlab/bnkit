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

    public void put(Object[] key, double count) {
        table.setValue(key, count);
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
