/**
 *
 */
package bn;

/**
 * Continuous domain definition. 
 * For checking validity of values for variables that belong to this domain.
 * @author mikael
 */
public class Continuous implements Domain {

    public boolean isValid(Object value) {
        try {
            Double x = (Double) value;
            return true;
        } catch (ClassCastException e) {
            return false;
        }
    }
}
