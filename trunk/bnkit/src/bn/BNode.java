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

import java.util.List;

/**
 * This interface must be implemented by a node in a Bayesian network.
 * @author mikael
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
     *
     * @return the table with entries
     */
    public EnumTable getTable();

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
     * Make BNode into a FactorTable
     */
    public FactorTable makeFactor(BNet bn);

    /**
     * Method used to modify the CPT/CDT to be modified (EM uses this).
     *
     * @param key the boolean key (how conditioning variables are set)
     * @param value the value of the conditioned variable
     * @param prob the expectation
     */
    public void countInstance(Object[] key, Object value, Double prob);

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

}