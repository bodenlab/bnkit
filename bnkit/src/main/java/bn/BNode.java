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

import bn.factor.Factor;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import bn.factor.AbstractFactor;
import bn.prior.Prior;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Interface that must be implemented by anything that appears as a 
 * node in a Bayesian network (BNet).
 * A node is a network-aware representation of a variable.
 * This is where the parameters/settings of the variable are stored, 
 * and managed.
 * The node is what is instantiated when a variable takes on a value.
 * The variable decides on what type of values that can be assigned.
 * The node will point to what (enumerable) variables that condition 
 * the assignments of the node variable.
 * @author Mikael Boden
 */
public interface BNode {

    public String getName();

    /**
     * Compute the probability of this variable taking the specified value
     *
     * @param key the condition
     * @param value the value
     * @return the (conditional) probability
     */
    public Double get(Object[] key, Object value);

    /**
     * Compute the probability of this variable taking the specified value
     *
     * @param value the value
     * @param key the condition
     * @return the (conditional) probability
     */
    public Double get(Object value, Object... key);

    /**
     * Get the prior probability of this variable taking the specified value.
     *
     * @param value
     * @return the (prior) probability
     */
    public Double get(Object value);

    /**
     * Retrieve variable of this node
     *
     * @return the variable of the node
     */
    @SuppressWarnings("rawtypes")
    public Variable getVariable();

    /**
     * Retrieve the variables of the parent nodes
     *
     * @return the variables of the parent nodes
     */
    public List<EnumVariable> getParents();

    /**
     * Retrieve the table with probability distributions
     * @deprecated BNodes should not need to implement this. The only place this is used is in the GUIs.
     * @return the table with entries
     */
    public EnumTable getTable();

    /**
     * Retrieve the distribution for this node that applies GIVEN the parents' instantiations.
     * Requires all parent nodes to be instantiated.
     * @param key the parents' values
     * @return the distribution of the variable for this node
     */
    public Distrib getDistrib(Object[] key);
    
    public Distrib getDistrib();

    /**
     * Pretty-print this node
     */
    public void print();

    /**
     * Predefined type of BNode
     *
     * @return text string so that BNode instance can be created after save
     */
    public String getType();

    /**
     * Dump parameters to text (for saving node to file)
     *
     * @see bn.BNode#setState(String)
     */
    public String getStateAsText();

    /**
     * Set parameters from text dump 
     *
     * @param dump
     * @return true if successful, false otherwise
     * @see bn.BNode#getStateAsText()
     */
    public boolean setState(String dump);

    /**
     * Check if the node is parent-less
     *
     * @return true if a root in the BN
     */
    public boolean isRoot();

    /**
     * Instantiate the node to a specific value
     *
     * @param value the value
     */
    public void setInstance(Object value);

    /**
     * Remove the instantiation of the node
     */
    public void resetInstance();

    /**
     * Retrieve the value of the instantiated node
     *
     * @return the value that the node is instantiated to, null if not
     * instantiated
     */
    public Object getInstance();

    /**
     * Make a native Distrib instance out of a collection of samples.
     * This method is currently only required in very special circumstances:
     * - in approximate inference, AND
     * - where the node is non-enumerable (e.g. GDT and DirDT) as used in bn.SampleTrace
     * @param samples
     * @return an instance of Distrib that can be used to populate this node.
     */
    public Distrib makeDistrib(Collection<Sample> samples);
    
    /**
     * Make BNode into a Factor.
     * The method reduces the factor so that only nominated nodes are included.
     * @param rel relevant variables with evidence if available
     * @return factor from BNode, taking evidence and (ir)relevance of parent variables into account
     */
    public Factor makeFactor(Map<Variable, Object> rel);

    /**
     * Make BNode into a Factor.
     * The method reduces the factor so that only nominated nodes are included.
     * @param rel relevant variables with evidence if available
     * @return factor from BNode, taking evidence and (ir)relevance of parent variables into account
     */
    public AbstractFactor makeDenseFactor(Map<Variable, Object> rel);

    /**
     * Method used to modify the CPT/CDT to be modified (EM uses this).
     *
     * @param key the boolean key (how conditioning variables are set)
     * @param value the value of the conditioned variable
     * @param prob the expectation
     */
    public void countInstance(Object[] key, Object value, Double prob);
    
    /**
     * Method used to modify the CPT/CDT to be modified (ApproxInfer uses this).
     *
     * @param key the boolean key (how conditioning variables are set)
     * @param value the value of the conditioned variable
     */
    public void countInstance(Object[] key, Object value);
    
    /**
     * For countInstance to have an effect on the actual CPT, this method needs
     * to be called. It should compute the CPT that maximizes the likelihood of
     * the counts (e.g. as used by EM).
     */
    public void maximizeInstance();

    /**
     * @return true if this node should be included in training, false otherwise
     */
    public boolean isTrainable();

    public void randomize(long seed);
    
    /**
     * Set whether or not a node is relevant to the current query (Inference)
     * @deprecated may interfere when multiple inferences are run multi-threaded
     */
    public void setRelevant(boolean relevant);
    
    /**
     * Set whether or not a node should be trained
     * @param trainable
     */
    public void setTrainable(boolean trainable);
    
    /**
     * @return true if this node is relevant to the current query (Inference)
     * @deprecated setting may interfere when multiple inferences are run multi-threaded, hence checking this status is moot
     */
    public boolean isRelevant();
    
    /**
     * get the condition data given the condition index in the enum table 
     * @param conditionIndex parent value index
     * @return array of data
     */
    public List<Sample> getConditionDataset(int conditionIndex);
    
    /**
     * get a new, empty distribution used by this Bnode
     * @return distribution
     */
    public Distrib getlikelihoodDistrib();
    
    /**
     * set the distribution conditioned on the parent value
     * @param key parent value array
     * @param distr the distribution you would like to set
     */
    public void put(Object[] key, Distrib distr);
    
    /**
     * set the distribution for root node
     * @param prob distribution
     */
    public void put(Distrib prob);
    
    /**
     * set the distribution conditioned on the parent value
     * @param prob
     * @param key
     */
    public void put(Distrib prob, Object... key);
    
    /**
     * set the distribution conditioned on the parent value index
     * @param index parent value index
     * @param distr
     */
    public void put(int index, Distrib distr);

}
