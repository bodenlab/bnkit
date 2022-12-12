/*
 * bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2014 M. Boden et al.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bn;

import dat.Variable;
import bn.prob.MixtureDistrib;
import bn.prob.GaussianDistrib;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * Joint Density Function.
 * Class that represents a collection of continuous variables, and their distribution(s).
 *
 * @author mikael
 */
public class JDF {
  
    private final Map<Variable, Distrib> varmap;
    
    /**
     * Create a uni-variate JDF.
     * @param var the variable over which the density is defined
     */
    public JDF(Variable var) {
        this.varmap = new HashMap<>();
        this.varmap.put(var, null);
    }
    
    /**
     * Create a multi-variate JDF.
     * @param vars the variables over which the JDF is defined
     */
    public JDF(Variable... vars) {
        this.varmap = new HashMap<>();
        for (Variable var : vars)
            this.varmap.put(var, null);
    }

    /**
     * Create a multi-variate JDF.
     * @param vars the variables over which the JDF is defined
     */
    public JDF(Collection<Variable> vars) {
        this.varmap = new HashMap<>();
        for (Variable var : vars)
            this.varmap.put(var, null);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Entry<Variable, Distrib> entry : varmap.entrySet())
            sb.append(entry.getKey().getName() + ":" + entry.getValue().toString() + "|");
        return sb.toString();
    }
    
    public Variable[] getVarArray() {
        Variable[] arr = new Variable[varmap.size()];
        varmap.keySet().toArray(arr);
        return arr;
    }
    
    /**
     * Get the number of variables (and distributions) of this JDF.
     * @return the number of variables
     */
    public int getSize() {
        return varmap.size();
    }
    
    /**
     * Assign a distribution to a variable.
     * Note that variables can be used for vectors for multi-variate distribs.
     * @param d distribution
     * @param var variable
     */
    public void setDistrib(Distrib d, Variable var) {
        this.varmap.put(var, d);
    }
  
    /**
     * Retrieve the distribution for a single variable.
     * Note that a distribution must be defined for the variable, which may be a vector.
     * @param var variable
     * @return the distribution
     */
    public Distrib getDistrib(Variable var) {
        Distrib d = varmap.get(var);
        return d;
    }
    
    /**
     * Mix the distribution currently associated with specified variable,
     * with that provided as a parameter, weighted as specified.
     * The result will always be a new mixture distribution, even if the original distribution 
     * is not. If the latter, we assume the weight of the original distribution is 1.
     * If there is no distribution already assigned to the variable, a new singleton mix will 
     * be created (with weight assigned as provided).
     * @param var the variable that will have a new, or amended mixed distribution
     * @param d the distribution that will form part of the mixture
     * @param weight the (absolute) weight that the new distribution will have, if NaN it is set to 1.0
     * @return the new "mixed" distribution
     */
    public Distrib mixDistrib(Variable var, Distrib d, double weight) {
        if (Double.isNaN(weight))
            weight = 1.0;
        Distrib prev = varmap.get(var);
        if (prev != null) {
            try {
                MixtureDistrib prev_mix = (MixtureDistrib) prev;
                // slot in new distribution, by examining previous distribs first, only mix new if must
                prev_mix.addDistrib(d, weight);
                return prev_mix;
            } catch (ClassCastException e) {
                // The distribution used previously is NOT a mixture
                MixtureDistrib mix = new MixtureDistrib(prev, 1.0); // we have to assume that the weight was 1.0
                mix.addDistrib(d, weight);
                return mix;
            }
        } else { // no previous distribution associated with variable, so we create a new singleton "mixture"
            MixtureDistrib mix = new MixtureDistrib(d, weight);
            varmap.put(var, mix);
            return mix;
        }
    }
    
    /**
     * Sample the assignment of a variable assuming nothing about others.
     * @param var variable
     * @return a sample value
     */
    public Object sample(Variable var) {
        return sample(var, new Variable.Assignment[] {});
    }

    /**
     * Sample the assignment of a variable given the assignments of others.
     * Note that currently the variable is assumed to be independent of the others.
     * @param var variable
     * @param assignments variable value assignments
     * @return a sample assignment of variable, chosen per density
     */
    public Object sample(Variable var, Variable.Assignment... assignments) {
        /* Here, assignments should be considered.
           Currently, all variables are assumed to be independent of one another. */
        Distrib d = varmap.get(var);
        if (d != null)
            return d.sample();
        return null;
    }
    
    /**
     * Calculate the probability of the assignments according to densities of each matching variable.
     * Assumes independence, i.e. works accurately if non-enumerables occur only as children in BN.
     * @param assignments variable assignments
     * @return 
     */
    public double density(Variable.Assignment... assignments) {
        double p = 1.0;
        for (Variable.Assignment assign : assignments) {
            Distrib d = this.getDistrib(assign.var);
            if (d != null) // found the distribution for the variable
                p *= d.get(assign.val);
        }
        return p;
    }
    
    /**
     * Create a new JDF from variables with distributions in the two provided.
     * Non-overlapping variables extend the resulting JDF.
     * Overlapping variables should not occur unless the variables are parents,
     * which is currently not supported, so this possibility is ignored.
     * @param jdf1 joint density function one (may be null)
     * @param jdf2 joint density function two (may be null; note if both are null, null is returned)
     * @return the convolution of the two JDFs
     */
    public static JDF combine(JDF jdf1, JDF jdf2) {
        if (jdf1 == null && jdf2 == null)
            return null;
        if (jdf1 == null) {
            Variable[] vars = jdf2.getVarArray();
            JDF jdf = new JDF(vars);
            for (Variable var : vars)
                jdf.setDistrib(jdf2.getDistrib(var), var);
            return jdf;
        } else if (jdf2 == null) {
            Variable[] vars = jdf1.getVarArray();
            JDF jdf = new JDF(vars);
            for (Variable var : vars)
                jdf.setDistrib(jdf1.getDistrib(var), var);
            return jdf;
        }
        Variable[] all = new Variable[jdf1.getSize() + jdf2.getSize()];
        Variable[] vars1 = jdf1.getVarArray();
        Variable[] vars2 = jdf2.getVarArray();
        System.arraycopy(vars1, 0, all, 0, vars1.length);
        System.arraycopy(vars2, 0, all, vars1.length, vars2.length);
        JDF jdf = new JDF(all);
        for (Variable var : vars1)
            jdf.setDistrib(jdf1.getDistrib(var), var);
        for (Variable var : vars2)
            jdf.setDistrib(jdf2.getDistrib(var), var);
        return jdf;
    }
    
    /**
     * Add (mix distributions of) a (weighted) JDF to the current.
     * This function may be used to combine JDFs defined over compatible sets of variables,
     * as marginalisation of other variables occur.
     * Note that the current JDF changes as a result.
     * @param partner the other JDF
     * @param weight the weight that applies to mix distributions in partner with those in current JDF
     * @return 
     */
    public JDF mix(JDF partner, double weight) {
        // TODO: this needs to create a *new* JDF
        for (Entry<Variable, Distrib> entry : partner.varmap.entrySet()) {
            Variable var = entry.getKey();
            Distrib d = entry.getValue();
            if (d == null)
                throw new JDFRuntimeException("Unassigned variable: " + var.toString());
            mixDistrib(var, d, weight);
        }
        return this;
    }
    
    /**
     * Mix two JDFs, by mixing their distributions. Mix is unweighted, giving both JDFs equal influence.
     * Note that the resulting JDF is new JDF, and the two component JDFs are not changed.
     * @param orig start JDF
     * @param other other JDF to be mixed with the start JDF
     * @return a new JDF containing the mixture distributions
     */
    public static JDF mix(JDF orig, JDF other) {
        return JDF.mix(orig, 1.0, other, 1.0);
    }
    
    /**
     * Mix two JDFs, by mixing their distributions, the second weighted as specified.
     * Note that the resulting JDF is new JDF, and the two component JDFs are not changed.
     * @param orig start JDF
     * @param other other JDF to be mixed with the start JDF
     * @param other_weight the weight with which the other JDF is mixed in
     * @return a new JDF containing the mixture distributions
     */
    public static JDF mix(JDF orig, JDF other, double other_weight) {
        return JDF.mix(orig, 1.0, other, other_weight);
    }
    
    /**
     * Mix two JDFs, by mixing their distributions, weighted as specified.
     * Note that the resulting JDF is new JDF, and the two component JDFs are not changed.
     * @param orig start JDF
     * @param weight_orig weighting for start JDF
     * @param other other JDF to be mixed with the start JDF
     * @param weight_other the weight with which the other JDF is mixed in
     * @return a new JDF containing the mixture distributions
     */
    public static JDF mix(JDF orig, double weight_orig, JDF other, double weight_other) {
        Set<Variable> vars = orig.varmap.keySet();
        JDF mixed = new JDF(vars);
        for (Entry<Variable, Distrib> entry : orig.varmap.entrySet()) {
            Variable orig_var = entry.getKey();
            Distrib orig_d = entry.getValue();
            Distrib other_d = other.getDistrib(orig_var);
            if (orig_d == null && other_d == null) { // if none is initialised before
                //throw new JDFRuntimeException("Cannot mix variable without distributions: " + orig_var.toString());
                mixed.setDistrib(null, orig_var);
                continue;
            }
            MixtureDistrib d;
            if (orig_d != null) {
                d = new MixtureDistrib(orig_d, weight_orig); // BUG fixed 26/9/14 (weight_orig used to be 1.0)
                if (other_d != null)
                    d.addDistrib(other_d, weight_other);
            } else
                d = new MixtureDistrib(other_d, weight_other);
            mixed.setDistrib(d, orig_var);
        }
        return mixed;
    }
    
    public static class JDFRuntimeException extends RuntimeException {
        public JDFRuntimeException(String msg) {
            super(msg);
        }
    }
    
    public static void main(String[] args) {
        
        Variable a = Predef.Real("A");
        Variable b = Predef.Real("B");
        Variable c = Predef.Real("C");
        
        Distrib d1 = new GaussianDistrib(0,1);
        Distrib d2 = new GaussianDistrib(1,1);
        Distrib d3 = new GaussianDistrib(2,1);
        
        JDF jdf1 = new JDF(a,b);
        JDF jdf2 = new JDF(c);
        JDF jdf3 = new JDF(a,b);
        
        jdf1.setDistrib(d1, a);
        jdf1.setDistrib(d2, b);
        jdf2.setDistrib(d3, c);
        jdf3.setDistrib(d1, a);
        jdf3.setDistrib(d2, b);
        JDF jdf4 = JDF.combine(jdf1, jdf2);
        JDF jdf5 = jdf1.mix(jdf3, 3.3);
        System.out.println(jdf1);
        System.out.println(jdf2);
        System.out.println(jdf3);
        System.out.println(jdf4);
        System.out.println(jdf5);
    }
}
