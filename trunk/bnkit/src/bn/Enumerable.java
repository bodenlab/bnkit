/**
 *
 */
package bn;

/**
 * This class represents a (countable) domain of values that discrete variables
 * can take.
 *
 * @author mikael
 */
public class Enumerable implements Domain {

    final int order;
    final Object[] values;

    public Enumerable(int order) {
        if (order < 2) {
            throw new RuntimeException("An Enumerable must have at least two values");
        }
        this.order = order;
        this.values = new Integer[order];
        for (int i = 0; i < order; i++) {
            this.values[i] = i;
        }
    }

    public Enumerable(Object[] values) {
        this.values = values;
        this.order = values.length;
        if (order < 2) {
            throw new RuntimeException("An Enumerable must have at least two values");
        }
    }

    public int size() {
        return order;
    }

    public int getIndex(Object value) {
        if (value instanceof java.lang.String) {
            for (int i = 0; i < values.length; i++) {
                if (((String) value).equals(values[i])) {
                    return i;
                }
            }
            throw new RuntimeException("Value " + value.toString() + " unknown to enumerable domain " + this.toString());
        } else {
            for (int i = 0; i < values.length; i++) {
                if (value == values[i]) {
                    return i;
                }
            }
            throw new RuntimeException("Value " + value.toString() + " unknown to enumerable domain " + this.toString());
        }
    }

    public Object get(int index) {
        return values[index];
    }

    public boolean isValid(Object value) {
        try {
            getIndex(value);
            return true;
        } catch (RuntimeException e) {
            return false;
        }
    }

    public static Enumerable bool = new Enumerable(new Boolean[]{true, false});
    public static Enumerable nacid = new Enumerable(new Character[]{'A', 'C', 'G', 'T'});
    public static Enumerable aacid = new Enumerable(new Character[]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'});

}
