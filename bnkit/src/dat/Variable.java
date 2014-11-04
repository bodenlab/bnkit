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
package dat;

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
public class Variable<E extends Domain> implements Comparable {

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
    public Variable(E domain, String name) {
        this.domain = domain;
        Variable.pool.put(this, Variable.pool.size());
        this.canonicalIndex = Variable.pool.get(this);
        this.name = name;
    }

    /**
     * Constructor intended for subclasses.
     * @param domain the domain of the variable
     */
    public Variable(E domain) {
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
     * Get the system-wide index of the variable, given at the time of construction.
     * @return the canonical index, which is unique for all variables 
     */
    public int getCanonicalIndex() {
        return canonicalIndex;
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
    public static Set<Variable<?>> getAll() {
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

    @Override
    public int compareTo(Object othvar) {
        try {
            Variable var = (Variable)othvar;
            return this.canonicalIndex - var.canonicalIndex;
        } catch (ClassCastException e) {
            throw new RuntimeException("Not a variable: " + othvar.toString());
        }
    }

    /**
     * A class to associate a variable to a value.
     * Useful for passing variable/value pairs to methods where multiples of these are required.
     */
    public static class Assignment {
        public final Variable var;
        public final Object val;
        public Assignment(Variable var, Object val) {
            this.var = var; 
            this.val = val;
        }
        public static Assignment[] array(Variable[] vars, Object[] vals) {
            Assignment[] ret = new Assignment[vars.length];
            for (int i = 0; i < vars.length && i < vals.length; i ++)
                ret[i] = new Assignment(vars[i], vals[i]);
            return ret;
        }
        public static <T extends Variable> Assignment[] array(List<T> vars, Object[] vals) {
            Assignment[] ret = new Assignment[vars.size()];
            for (int i = 0; i < ret.length && i < vals.length; i ++)
                ret[i] = new Assignment(vars.get(i), vals[i]);
            return ret;
        }
        public static Map<Variable, Object> toMap(Assignment[] assign_array) {
            Map<Variable, Object> map = new HashMap<>();
            for (Assignment assign : assign_array)
                map.put(assign.var, assign.val);
            return map;
        } 
        public static Map<Variable, Object> toMap(Collection<Assignment> assign_list) {
            Map<Variable, Object> map = new HashMap<>();
            for (Assignment assign : assign_list)
                map.put(assign.var, assign.val);
            return map;
        } 
        public String toString() {
            return var.toString() + "=" + val.toString();
        }
    }

    /**
     * A factory method for creating a single variable/value pair.
     * @param var
     * @param val
     * @return the variable paired with a value
     */
    public static Assignment assign(Variable var, Object val) {
        return new Variable.Assignment(var, val);
    }
    

}

