/**
 *
 */
package bn;

/**
 * Probability distribution defined. 
 * A variable in a BNode will be assigned an instance of a distribution in every row.
 * @author mikael
 *
 */
public interface Distrib {

    /**
     * Retrieve the probability of the value
     *
     * @param value
     * @return the probability
     */
    public double get(Object value);

    /**
     * Sample from the distribution
     *
     * @return a sample value, drawn from the distribution in proportion to its
     * probability
     */
    public Object sample();

}
