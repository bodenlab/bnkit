/**
 *
 */
package bn;

/**
 * Interface for inference engines.
 * @author mikael
 */
public interface Inference {

    /**
     * Use instantiation as per specified Bayesian network
     *
     * @param bn
     */
    public void instantiate(BNet bn);

    /**
     * Construct a query handle
     *
     * @param vars the variables that are queried
     * @return the handle
     */
    public Query makeQuery(Variable[] vars);

    /**
     * Infer the probabilities based on the query handle
     *
     * @param q query handle
     * @return the joint probability table with all variables in query handle
     */
    public JPT infer(Query q);

    public double getLogLikelihood();
}
