package bn;

/**
 * Interface for classes that implement type checking of variables
 * @author mikael
 */
public interface Domain {

    /**
     * Check if the specified value is valid for the domain
     *
     * @param value
     * @return true if valid, false otherwise
     */
    public boolean isValid(Object value);

}
