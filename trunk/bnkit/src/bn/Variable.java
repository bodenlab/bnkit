/**
 *
 */
package bn;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Class that defines the term (random) variable as it is used for Bayesian networks.
 * @author mikael
 * @param <E>
 *
 */
public class Variable<E extends Domain> {

    // A domain is a type checking entity (e.g. to check if a value is valid)
    private final E domain;
    // Variables are associated with a unique counter (to impose a canonical ordering amongst variables)
    final static Map<Variable<?>, Integer> pool = new HashMap<Variable<?>, Integer>();
    // counter of variables used in the current namespace
    private final int canonicalIndex;
    // default name of variable
    private String name = "";
    // predef name/type, needs to be specified for the variable to be saved with a BN
    private String type = null;
    // parameter-string to be saved to file with variable, if any
    private String params = null;

    /**
     * Constructor intended for subclasses.
     * @param domain the domain of the variable
     * @param name the name of the variable
     */
    protected Variable(E domain, String name) {
        this.domain = domain;
        Variable.pool.put(this, Variable.pool.size());
        this.canonicalIndex = Variable.pool.get(this);
        this.name = name;
    }

    /**
     * Constructor intended for subclasses.
     * @param domain the domain of the variable
     */
    protected Variable(E domain) {
        this.domain = domain;
        Variable.pool.put(this, Variable.pool.size());
        this.canonicalIndex = Variable.pool.get(this);
    }

    /**
     * Set name of variable
     * @param name 
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Get name of variable
     * @return the name (not necessarily unique)
     */
    public String getName() {
        return this.name;
    }

    /**
     * Gives a printable, unique, system name to the variable.
     * @return system name
     */
    public String toString() {
        return this.name + "." + this.canonicalIndex;
    }

    /**
     * Compares variables by their canonical index (not name)
     * @param var
     * @return true of the same (in the current namespace)
     */
    public boolean equals(Object var) {
        return (this.canonicalIndex == Variable.pool.get(var));
    }

    /**
     * Re-package an array of Variable to a list of them
     *
     * @param useParents an array of variables
     * @return a list of variables
     */
    public static Collection<Variable<?>> toList(Variable<?>[] useParents) {
        List<Variable<?>> list = new ArrayList<Variable<?>>(useParents.length);
        for (Variable<?> var : useParents) {
            list.add(var);
        }
        return list;
    }

    /**
     * Intersection between two variable collections.
     * @param c1 collection 1
     * @param c2 collection 2
     * @return the intersection between c1 and c2
     */
    public static Collection<Variable<Domain>> isect(Collection<Variable<Domain>> c1, Collection<Variable<Domain>> c2) {
        throw new RuntimeException("Not implemented");
    }

    /**
     * Get the domain of the variable
     * @return the domain
     */
    public E getDomain() {
        return domain;
    }

    /**
     * Get all variables in the current variable namespace 
     * @return all variables
     */
    public Set<Variable<?>> getAll() {
        return pool.keySet();
    }

    /* -------------------- Functions to define variables outside Java -------------------- */
    /* Includes quick-use variable pre-definitions, and variables defined in XML files etc. */
    /* ------------------------------------------------------------------------------------ */
    
    /**
     * For the variable to be saved it needs to have a predef name. 
     * Must be unique (programmer's responsibility).
     *
     * @param type
     * @see bn.Predef
     */
    public void setPredef(String type) {
        this.type = type;
    }

    /**
     * Get the descriptor for the type to which this variable belongs 
     * @return the type descriptor, e.g. "Boolean", "String", "Real"
     */
    public String getPredef() {
        return this.type;
    }

    /**
     * Set any parameters necessary to fully define the variable, including the domain.
     * Example: for a variable called Gender the params would be "Male;Female"
     * as in the XML BN file format:
     * <var name="Gender" params="Male;Female;" type="String"/>
     * @param params 
     */
    public void setParams(String params) {
        this.params = params;
    }

    /**
     * Get parameters for this variable, to be used when saving to text file etc
     * @return string representation of parameter values
     */
    public String getParams() {
        return this.params;
    }

}
