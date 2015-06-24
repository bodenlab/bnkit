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
package bn.node;

import bn.BNode;
import bn.CountTable;
import bn.Distrib;
import bn.prob.EnumDistrib;
import bn.factor.Factor;
import bn.JPT;
import bn.Predef;
import bn.Sample;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import dat.Enumerable;
import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import bn.factor.Factorize;

import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;

/**
 * Class for Noisy-OR version of a conditional probability table with 'smart' training from data.
 * See Russell and Norvig (2002) p. 500 for an informal explanation of how this table works with Boolean parent variables.
 * This implementation of Noisy-OR allows for any set of discrete parent variables, where the value 
 * in each parent to be used in the Noisy-OR assumption is defined.
 * The difference between this class and 'normal' Noisy-OR lies in training. Here maximizeInstance 
 * uses all samples to set the parameters of the table.
 * It can be used as any other CPT, e.g. as a parent variable and as a latent variable in EM-training.
 *
 * @author ralph
 */
public class SmartNoisyOR implements BNode, Serializable{

    private static final long serialVersionUID = 1L;
    final private EnumVariable var;
    final private EnumTable<EnumDistrib> table; // table of (enumerable) probability distributions
    final private List<Object> plabels; //list of values representing variable names in parent variables
    private EnumDistrib prior; // one (enumerable) probability distribution that is used if this variable is NOT conditioned
    final private int nParents;
    private CountTable count = null; // keep counts when learning/observing; first "parent" is the conditioned variable, then same order as in SmartNoisyOR
    private boolean relevant = false; //for inference, track whether the node is relevant to the query

    /**
     * Create a SmartNoisyOR table for a variable. The variable is
     * conditioned on a set of Enumerable variables. Values from
     * each parent variable should be specified to use in the SmartNoisyOR
     *
     * @param var variable
     * @param parents parent variables
     */
    public SmartNoisyOR(EnumVariable var, List<EnumVariable> parents, List<Object> labels) {
        this.var = var;
        if (parents.size() != labels.size()) {
        	throw new RuntimeException("number of labels should equal number of parents.");
        }
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<EnumDistrib>(parents);
                this.prior = null;
                this.nParents = parents.size();
                this.plabels = labels;
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        this.plabels = null;
    }
    
    /**
     * Create a SmartNoisyOR table for a variable. The variable is
     * conditioned on a set of Enumerable variables. Values from
     * each parent variable should be specified to use in the SmartNoisyOR
     *
     * @param var variable
     * @param parents parent variables
     */
    public SmartNoisyOR(EnumVariable var, List<EnumVariable> parents, Object[] labels) {
        this.var = var;
        if (parents.size() != labels.length) {
        	throw new RuntimeException("number of labels should equal number of parents.");
        }
        if (parents != null) {
            if (parents.size() > 0) {
                this.table = new EnumTable<EnumDistrib>(parents);
                this.prior = null;
                this.nParents = parents.size();
                this.plabels = new ArrayList<Object>();
                for (Object lab : labels){
                	this.plabels.add(lab);
                }
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        this.plabels = null;
    }	
    
    /**
     * Create a SmartNoisyOR table for a variable. The variable is
     * conditioned on a set of Enumerable variables. Values from
     * each parent variable should be specified to use in the SmartNoisyOR
     *
     * @param var variable
     * @param parents parent variables
     */
    public SmartNoisyOR(EnumVariable var, EnumVariable[] parents, Object[] labels) {
        this.var = var;
        if (parents.length != labels.length) {
        	throw new RuntimeException("number of labels should equal number of parents.");
        }
        if (parents != null) {
            if (parents.length > 0) {
                this.table = new EnumTable<EnumDistrib>(parents);
                this.prior = null;
                this.nParents = parents.length;
                this.plabels = new ArrayList<Object>();
                for (Object lab : labels){
                	this.plabels.add(lab);
                }
                return;
            }
        }
        this.table = null;
        this.nParents = 0;
        this.plabels = null;
    }

    /**
     * Create a prior (a SmartNoisyOR without conditioning variables) for a variable.
     *
     * @param var variable
     */
    public SmartNoisyOR(EnumVariable var) {
        this.var = var;
        this.table = null;
        this.nParents = 0;
        this.plabels = null;
    }

    /**
     * Create a SmartNoisyOR from a JPT. The variable in the JPT var is the variable
     * conditioned on in the SmartNoisyOR.
     *
     * @param jpt
     * @param var
     */
    public SmartNoisyOR(JPT jpt, EnumVariable var, List<Object> labels) {
    	this.plabels = labels;
        this.var = var;
        List<EnumVariable> SmartNoisyORParents = new ArrayList<>(jpt.getParents().size() - 1);
        int index = -1;
        for (int i = 0; i < jpt.getParents().size(); i++) {
            EnumVariable jptParent = jpt.getParents().get(i);
            if (jptParent == var) {
                index = i;
            } else {
                SmartNoisyORParents.add(jptParent);
            }
        }
        if (index == -1) {
            throw new RuntimeException("Invalid variable " + var + " for creating SmartNoisyOR");
        }
        if (SmartNoisyORParents.isEmpty()) { // no parents in SmartNoisyOR
            this.table = null;
            this.prior = new EnumDistrib(var.getDomain());
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                this.prior.set(jptkey[0], entry.getValue().doubleValue());
            }
            this.prior.normalise();
            this.nParents = 0;
        } else { // there are parents
            this.table = new EnumTable<>(SmartNoisyORParents);
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                System.out.println(jptkey[0]+" "+jptkey[1]+" "+jptkey[2]);
                Object[] SmartNoisyORkey = new Object[jptkey.length - 1];
                int j = 0;
                int nkey=0;
                for (int i=0; i<jptkey.length; i++) {
                	if ((Boolean) jptkey[i]) {
                		nkey++;
                	}
                }
                if (nkey <= 1) {
                	for (int i = 0; i < jptkey.length; i++) {
                		if (i != index) // selected variable
                		{
                			SmartNoisyORkey[j++] = jptkey[i];
                		}
                	}
                	int SmartNoisyOR_index = this.table.getIndex(SmartNoisyORkey);
                	EnumDistrib d = this.table.getValue(SmartNoisyOR_index);
                	if (d == null) {
                		d = new EnumDistrib(var.getDomain());
                	}
                	d.set(jptkey[index], entry.getValue().doubleValue());
                	this.table.setValue(SmartNoisyOR_index, d);
                }
            }
            for (Map.Entry<Integer, EnumDistrib> SmartNoisyOR_entry : this.table.getMapEntries()) {
                SmartNoisyOR_entry.getValue().normalise();
            }
            this.nParents = SmartNoisyORParents.size();
        }
    }

    /**
     * Assign entries in SmartNoisyOR according to a given JPT. Not bullet-proof, watch
     * out for SmartNoisyOR - JPT incompatibilities, e.g. variable order which is not
     * checked currently
     *
     * @param jpt
     */
    public void put(JPT jpt) {
        int index = -1;
        for (int i = 0; i < jpt.getParents().size(); i++) {
            EnumVariable jptParent = jpt.getParents().get(i);
            if (jptParent == var) {
                index = i;
            } else {
                if (!this.getParents().contains(jptParent)) {
                    throw new RuntimeException("No variable " + jptParent + " in SmartNoisyOR");
                }
            }
        }
        if (index == -1) {
            throw new RuntimeException("No variable " + var + " in JPT");
        }
        if (this.table == null && jpt.getParents().size() == 1) { // no parents in SmartNoisyOR
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                this.prior.set(jptkey[0], entry.getValue().doubleValue());
            }
            this.prior.normalise();
        } else if (jpt.getParents().size() == this.getParents().size() + 1) { // there are parents
            for (Map.Entry<Integer, Double> entry : jpt.table.getMapEntries()) {
                Object[] jptkey = jpt.table.getKey(entry.getKey().intValue());
                Object[] SmartNoisyORkey = new Object[jptkey.length - 1];
                int j = 0;
                for (int i = 0; i < jptkey.length; i++) {
                    if (i != index) // selected variable
                    {
                        SmartNoisyORkey[j++] = jptkey[i];
                    }
                }
                int SmartNoisyOR_index = this.table.getIndex(SmartNoisyORkey);
                EnumDistrib d = this.table.getValue(SmartNoisyOR_index);
                if (d == null) {
                    d = new EnumDistrib(var.getDomain());
                }
                d.set(jptkey[index], entry.getValue().doubleValue());
                this.table.setValue(SmartNoisyOR_index, d);
            }
            for (Map.Entry<Integer, EnumDistrib> SmartNoisyOR_entry : this.table.getMapEntries()) {
                SmartNoisyOR_entry.getValue().normalise();
            }
        } else {
            throw new RuntimeException("Cannot set SmartNoisyOR from given JPT: " + jpt);
        }
    }

    /**
     * Retrieve the distribution for this node that applies GIVEN the parents' instantiations.
     * Requires all parent nodes to be instantiated.
     * @param key the parent values
     * @return the distribution of the variable for this node
     */
    public Distrib getDistrib(Object[] key) {
        if (this.table == null || key == null)
            return this.getDistrib();
        try {
            return this.table.getValue(key);
        } catch (RuntimeException e) {
            throw new RuntimeException("Evaluation of SmartNoisyOR " + this.toString() + " failed since condition was not fully specified: " + e.getMessage());
        }
    }
    
    @Override
    public Distrib makeDistrib(Collection<Sample> samples) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    /**
     * Make a Factor out of this SmartNoisyOR. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * @param relevant only include relevant nodes, with instantiations if available
     * @return factor of SmartNoisyOR considering if parents are relevant
     */
    public Factor makeFactor(Map<Variable, Object> relevant) {
        List<EnumVariable> vars_old = this.getParents();
        EnumVariable myvar = this.getVariable();
        // get value of this node if any assigned
        Object varinstance = relevant.get(myvar); 
        Enumerable dom = myvar.getDomain();
        if (vars_old != null) { // there are parent variables
            Object[] searchkey = new Object[vars_old.size()];
            List<Variable> vars_new = new ArrayList<>(vars_old.size() + 1);
            List<EnumVariable> irrel_pars = new ArrayList<>(); //irrelevant parents
            for (int i = 0; i < vars_old.size(); i++) {
                EnumVariable parent = vars_old.get(i);
                boolean parent_is_relevant = relevant.containsKey(parent);
                // Record irrelevant parents to sum out
                //FIXME when should a parent be removed? Allow it to influence factor table then remove it?
                // If parent is evidenced it will not be included in factor table; it is removed later through marginalization
                if (!parent_is_relevant) 
                    irrel_pars.add(parent);
                else
                    searchkey[i] = relevant.get(parent);
                if (searchkey[i] == null) // new factor needs to include this variable (to be summed out before returned if irrelevant)
                    vars_new.add(parent);
            }
            if (varinstance == null) {
                vars_new.add(var);
            }
            Factor ft = new Factor(vars_new);
            if (varinstance != null) {
                ft.evidenced = true;
            } else {
            	//FIXME what is this loop doing? Should it be searchkey[i]??
                for (int i = 0; i < searchkey.length; i++) {
                    if (searchkey != null) {
                        ft.evidenced = true;
                        break;
                    }
                }
            }
            /*
             * if the search key is not 'OR'; i.e. there is more than
             * one parent labels set, then adding factors is not a simple
             * case of looking up the table.
             * The easiest solution for now may be to populate the table with 
             * the calculated values (from get) as required.
             */
            int nkey = 0;
            for (int i=0; i<searchkey.length; i++) {
            	if (searchkey[i] != null) {
            		if (searchkey[i].equals(this.plabels.get(i))) {
            			nkey++;}
            	}
            }
            if (nkey > 1) {
            	insert(searchkey);
            }
            
            //entries that are needed should exist now
            int[] indices = table.getIndices(searchkey);
            Object[] newkey = new Object[vars_new.size()];
            
            for (int index : indices) {
            	
                EnumDistrib d = table.getValue(index);
                if (d != null) {
                    Object[] key = table.getKey(index);
                    
                    int newcnt = 0;
                    for (int i = 0; i < key.length; i++) {
                        if (searchkey[i] == null) {
                            newkey[newcnt++] = key[i];
                        }
                    }
                    if (varinstance != null) { // the variable for this SmartNoisyOR is instantiated
                        if (newkey.length == 0) // atomic FactorTable
                            ft.setFactor(null, d.get(varinstance));
                        else
                            ft.addFactor(newkey, d.get(varinstance));
                    } else { // the variable for this SmartNoisyOR is NOT instantiated so we add one entry for each possible instantiation
                        for (int j = 0; j < dom.size(); j++) {
                            newkey[newkey.length - 1] = dom.get(j);
                            Double p = d.get(j);
                            if (p != null) {
                                ft.addFactor(newkey, p);
                            }
                        }
                    }
                } else { // this entry is null
                    //
                }
            }
            if (!irrel_pars.isEmpty()) {
            	ft = ft.marginalize(irrel_pars);
            }
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                Factor ft = new Factor();
                ft.setFactor(this.prior.get(varinstance));
                return ft;
            }
            List<Variable> vars_new = new ArrayList<>(1);
            vars_new.add(var);
            Factor ft = new Factor(vars_new);
            Object[] newkey = new Object[1];
            EnumDistrib d = this.prior;
            for (int j = 0; j < dom.size(); j++) {
                newkey[0] = dom.get(j);
                ft.addFactor(newkey, d.get(j));
            }
            return ft;
        }
    }

    /**
     * Make a Factor out of this SmartNoisyOR. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * @param relevant only include relevant nodes, with instantiations if available
     * @return factor of SmartNoisyOR considering if parents are relevant (rel)
     */
    @Override
    public AbstractFactor makeDenseFactor(Map<Variable, Object> relevant) {
        List<EnumVariable> parents = this.getParents();
        EnumVariable myvar = this.getVariable();
        // get value of this node if any assigned
        Object varinstance = relevant.get(myvar); 
        Enumerable dom = myvar.getDomain();
        if (parents != null) { // there are parent variables
            Object[] searchcpt = new Object[parents.size()];
            List<Variable> fvars = new ArrayList<>(parents.size() + 1); // factor variables
            List<EnumVariable> sumout = new ArrayList<>();  // irrelevant variables to be summed out later
            for (int i = 0; i < parents.size(); i++) {
                EnumVariable parent = parents.get(i);
                // If parent is evidenced it will not be included in factor table; 
                // Record irrelevant parents to sum out: removed later through marginalization
                if (!relevant.containsKey(parent)) 
                    sumout.add(parent);
                else
                    searchcpt[i] = relevant.get(parent);
                if (searchcpt[i] == null) // new factor will include this variable
                    fvars.add(parent);
            }
            
            /*
             * if the search key is not 'OR'; i.e. there is more than
             * one parent labels set, then adding factors is not a simple
             * case of looking up the table.
             * The easiest solution for now may be to populate the table with 
             * the calculated values (from get) as required.
             */
            int mkey = 0;
            for (int i=0; i<searchcpt.length; i++) {
            	if (searchcpt[i] != null) {
            		if (searchcpt[i].equals(this.plabels.get(i))) {
            			mkey++;}
            	}
            }
            if (mkey > 1) {
            	insert(searchcpt);
            }
            
            if (varinstance == null) {
                fvars.add(myvar);
            }
            Variable[] vars_arr = new Variable[fvars.size()];
            fvars.toArray(vars_arr);
            AbstractFactor ft = new DenseFactor(vars_arr);

            //if (Factorize.exitIfInvalid2(ft, this.toString())){
            //	System.out.println("invalid");
            //}
            EnumVariable[] evars = ft.getEnumVars(); // the order may have changed
            int[] xcross = new int[parents.size()];
            int[] ycross = new int[evars.length];
            table.crossReference(xcross, evars, ycross);
            int missing = -1; // the position of *myvar* in the new factor (if applicable)
            for (int i = 0; i < ycross.length; i ++)
                if (ycross[i] == -1) {
                    missing = i;
                    break;
                }
            // set factor to be "evidenced" is there was an evidence used
            if (varinstance != null) {
                ft.evidenced = true;
            } else {
                for (Object instcpt : searchcpt) {
                    if (instcpt != null) {
                        ft.evidenced = true;
                        break;
                    }
                }
            }
            //if (Factorize.exitIfInvalid2(ft, this.toString())){
            //	System.out.println("invalid");
            //}
            int[] indices = table.getIndices(searchcpt);
            Object[] fkey = new Object[evars.length];
            for (int index : indices) {
                EnumDistrib d = table.getValue(index);
                //if distribution is null, will need to calculate it
                Object[] orkey = table.getKey(index);
                int nkey = 0;
                for (int i=0; i<orkey.length; i++) {
                	if (orkey[i] != null) {
                		if (orkey[i].equals(this.plabels.get(i))) {
                			nkey++;}
                	}
                }
                if (nkey > 1) {
                	insert(orkey);
                }
                
                if (d != null) { // there is a distribution associated with this entry in the CPT
                    Object[] cptkey = table.getKey(index); // work out the condition for this entry
                    for (int i = 0; i < cptkey.length; i++) {
                        if (xcross[i] != -1) 
                            fkey[xcross[i]] = cptkey[i];
                    }
                    if (varinstance != null) { // the variable for this CPT is instantiated
                        if (fkey.length == 0) // atomic factor
                            ft.setValue(d.get(varinstance));
                        else
                            ft.setValue(fkey, d.get(varinstance));
                        //if (Factorize.exitIfInvalid2(ft, this.toString())){
                        //	System.out.println("invalid");
                        //}
                    } else { // the variable for this CPT is NOT instantiated so we add one entry for each possible instantiation
                        for (int j = 0; j < dom.size(); j++) {
                            fkey[missing] = dom.get(j);
                            Double p = d.get(j);
                            ft.setValue(fkey, p);
                        }
                    }
                } 
            }
            if (!sumout.isEmpty()) {
                Variable[] sumout_arr = new Variable[sumout.size()];
                sumout.toArray(sumout_arr);
            	ft = Factorize.getMargin(ft, sumout_arr);
            }
           // if (Factorize.exitIfInvalid2(ft, this.toString())){
           // 	System.out.println("invalid");
            //}
            return ft;
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                AbstractFactor ft = new DenseFactor();
                ft.setValue(this.prior.get(varinstance));
                Factorize.exitIfInvalid(ft, this.toString());
                return ft;
            }
            AbstractFactor ft = new DenseFactor(myvar);
            Object[] newkey = new Object[1];
            EnumDistrib d = this.prior;
            for (int j = 0; j < dom.size(); j++) {
                newkey[0] = dom.get(j);
                Double p = d.get(j);
                ft.setValue(newkey, p);
            }
            Factorize.exitIfInvalid(ft, this.toString());
            return ft;
        }
    }

    /**
     * Get the name of the SmartNoisyOR
     */
    public String getName() {
        return getVariable().getName();
    }

    /**
     * Get the variable of the SmartNoisyOR.
     *
     * @return the variable of the SmartNoisyOR
     */
    public EnumVariable getVariable() {
        return var;
    }

    /**
     * Retrieve the names of all parent variables (that is all variables that
     * are conditioning the SmartNoisyOR variable)
     *
     * @return the variables of the parent variables
     */
    public List<EnumVariable> getParents() {
        if (table == null) {
            return null;
        }
        List<EnumVariable> parents = table.getParents();
        return parents;
    }
    
    /**
     * Retrieve the labels used for defining the 'on'/true case 
     * for the SmartNoisyOR parents
     * @return the labels of the parent variables
     */
    public List<Object> getParentLabels() {
    	return this.plabels;
    }

    /**
     * Check if this SmartNoisyOR is a "root" SmartNoisyOR, i.e. has no parents.
     *
     * @return true if the SmartNoisyOR has no parents, false if it has
     */
    public boolean isRoot() {
        if (table == null) {
            return true;
        }
        return false;
    }

    /**
     * Get the conditional probability of the variable (represented by this SmartNoisyOR)
     *
     * @param key parent key (condition); if null the SmartNoisyOR is assumed to be a
     * prior.
     * @return the probability of the variable
     */
    public Double get(Object[] key, Object value) {
    	if (key == null) {
            return prior.get(value);
        }
        if (key.length == 0) {
            return prior.get(value);
        }
        //first check if we're looking at the false values
        int nkey = 0;
    	for (int i=0; i<key.length; i++) {
    		if (key[i] != null) {
    			if (key[i].equals(this.plabels.get(i))) {
    				nkey++;}
    		}
    	}
    	if (nkey == 0) {
    		//if none of the valid parent labels are set,
    		//just return the entry for this key
    		//possibly should return the product of all
    		//the rows with no parent switched on?
    		return table.getValue(key).get(value);
    	}
        
        Double pNotTrue = 1.0;
        //get the indices in the table
        int [] indices = this.table.getIndices();
        //create a map of parent variables to keys,
        //we'll need to look this information up later on
        HashMap<Object, List<Object []>> keyMap = new HashMap<Object, List<Object []>>();
        for (Object parent : this.table.getParents()) {
        	keyMap.put(parent, new ArrayList<Object []>());
        }
        //go through each row (key) in the table and assign the key to the relevant parent
        for (Integer i : indices) {
        	Object [] thisKey = this.table.getKey(i);
        	nkey = 0;
        	for (int k=0; k<thisKey.length; k++) {
        		if (thisKey[k].equals(this.plabels.get(k))) nkey++;
        	}
        	for (int j=0; j<thisKey.length; j++) {
        		if (thisKey[j].equals(this.plabels.get(j)) && nkey == 1) {
        			List<Object []> thisList = keyMap.get(this.table.getParents().get(j));
        			thisList.add(thisKey);
        			keyMap.put(this.table.getParents().get(j), thisList);
        		}
        	}
        }
        //finally, go through the key and compute probability
        for (int i=0; i<key.length; i++) {
        	if (this.plabels.get(i).equals(key[i])) {
        		/* key item is a valid parent label:
        		 * go through cases where the label is true
        		 * check every index where the parent label is true
        		 */
        		Object thisParent = this.table.getParents().get(i);
        		List<Object []> parentKeys = keyMap.get(thisParent);
        		for (Object [] pkey : parentKeys) {
        			Double p = this.table.getValue(pkey).get(value);
        			if (p==null)
    					return null;
    				pNotTrue*=(1.0-p);
        		}
        	}
        }
        return 1.0-pNotTrue;
    }

    /**
     * Get the conditional probability of the variable (represented by this SmartNoisyOR)
     *
     * @param key parent key (condition); if null the SmartNoisyOR is assumed to be a
     * prior.
     * @return the probability of the variable
     */
    public Double get(Object value, Object... key) {
    	if (key == null) {
            return prior.get(value);
        }
        if (key.length == 0) {
            return prior.get(value);
        }
        //first check if we're looking at the false values
        int nkey = 0;
    	for (int i=0; i<key.length; i++) {
    		if (key[i].equals(this.plabels.get(i))) {
    			nkey++;}
    	}
    	if (nkey == 0) {
    		//if none of the valid parent labels are set,
    		//just return the entry for this key
    		//possibly should return the product of all
    		//the rows with no parent switched on?
    		return table.getValue(key).get(value);
    	}
        
        Double pNotTrue = 1.0;
        //get the indices in the table
        int [] indices = this.table.getIndices();
        //create a map of parent variables to keys,
        //we'll need to look this information up later on
        HashMap<Object, List<Object []>> keyMap = new HashMap<Object, List<Object []>>();
        for (Object parent : this.table.getParents()) {
        	keyMap.put(parent, new ArrayList<Object []>());
        }
        //go through each row (key) in the table and assign the key to the relevant parent
        for (Integer i : indices) {
        	Object [] thisKey = this.table.getKey(i);
        	nkey = 0;
        	for (int k=0; k<thisKey.length; k++) {
        		if (thisKey[k].equals(this.plabels.get(k))) nkey++;
        	}
        	for (int j=0; j<thisKey.length; j++) {
        		if (thisKey[j].equals(this.plabels.get(j)) && nkey == 1) {
        			List<Object []> thisList = keyMap.get(this.table.getParents().get(j));
        			thisList.add(thisKey);
        			keyMap.put(this.table.getParents().get(j), thisList);
        		}
        	}
        }
        //finally, go through the key and compute probability
        for (int i=0; i<key.length; i++) {
        	if (this.plabels.get(i).equals(key[i])) {
        		/* key item is a valid parent label:
        		 * go through cases where the label is true
        		 * check every index where the parent label is true
        		 */
        		Object thisParent = this.table.getParents().get(i);
        		List<Object []> parentKeys = keyMap.get(thisParent);
        		for (Object [] pkey : parentKeys) {
        			Double p = this.table.getValue(pkey).get(value);
        			if (p==null)
    					return null;
    				pNotTrue*=(1.0-p);
        		}
        	}
        }
        return 1.0-pNotTrue;
    }
    
    /**
     * Method for inserting non-OR values into the table
     * when executing makeFactor.
     * @param key
     * @param prob
     */
    public void insert(Object [] key) {
        int nkey = 0;
        for (int i=0; i<key.length; i++) {
        	if (key[i] != null) {
        		if (key[i].equals(this.plabels.get(i))) {
        			nkey++;}
        	}
        }
        Enumerable en = var.getDomain();
        double [] probabilities = new double[en.size()];
        if (nkey > 1) {
            //need to create a new entry
            for (int i=0; i<en.size(); i++) {
            	Double searchkeyprob = get(key, en.get(i));
            	probabilities[i] = searchkeyprob;
            }
            //if the key contains a null value
            Object [] newkey = new Object [key.length];
            for (int i=0; i<key.length; i++) {
            	if (key[i] != null)
            		newkey[i] = key[i];
            	else {
            		//set to a non positive value of the parent
            		List<EnumVariable> parents = this.table.getParents();
					EnumVariable par = parents.get(i);
            		Object [] values = par.getDomain().getValues();
            		Object lab = null;
					for (Object v : values) {
						if (!v.equals(this.plabels.get(i))) {
							lab = v; break;
						}
					}
					newkey[i] = lab;
            	}
            }
            table.setValue(newkey, new EnumDistrib(var.getDomain(), probabilities));
        }
    }
    
    /**
     * Set entry (or entries) of the NoisyOR to the specified probability value index
     * (variable is true).
     */
    public void put(int index, Distrib prob) {
    	table.setValue(index, (EnumDistrib)prob);
    }
    

    /**
     * Set entry (or entries) of the SmartNoisyOR to the specified probability value
     * (variable is true).
     *
     * @param key the boolean key (probabilistic condition)
     * @param prob the probability value (must be >=0 and <=1)
     */
    public void put(Object[] key, Distrib prob) {
    	
        if (key == null) {
            put(prob);
        } else if (key.length == 0) {
            put(prob);
        
        } else {
        	//check that the key is valid for this Noisy-OR
        	int nkey = 0;
        	for (int i=0; i<key.length; i++) {
        		if (key[i].equals(this.plabels.get(i))) {
        			nkey++;}
        	}
        	if (nkey <= 1) {
        		//if no more than 1 of the parent labels is set, the key is valid
        		table.setValue(key, (EnumDistrib)prob);}
        }
    }

    /**
     * Set entry (or entries) of the SmartNoisyOR to the specified probability value
     * (variable is true).
     *
     * @param prob the probability value (must be >=0 and <=1)
     * @param key the key (the condition)
     */
    public void put(Distrib prob, Object... key) {
        if (key == null) {
            put(prob);
        } else if (key.length == 0) {
            put(prob);
        } else {
        	//check that the key is valid for this Noisy-OR
        	int nkey = 0;
        	for (int i=0; i<key.length; i++) {
        		if (key[i].equals(this.plabels.get(i))) {
        			nkey++;}
        	}
        	if (nkey <= 1) {
        		table.setValue(key, (EnumDistrib)prob);
        	}
        }
    }

    /**
     * Set the prior probability of this SmartNoisyOR that has no parents.
     *
     * @param prob
     */
    public void put(Distrib prob) {
        if (!isPrior()) {
            throw new RuntimeException("Unable to set prior. SmartNoisyOR " + var + " is conditioned.");
        }
        if (!((EnumDistrib)prob).isNormalised()) {
            throw new RuntimeException("Probability value is invalid: " + prob);
        }
        prior = (EnumDistrib)prob;
    }

    /**
     * Checks if this SmartNoisyOR has no parents.
     */
    public boolean isPrior() {
        return table == null;
    }

    /**
     * Provide a non-unique string representation of this SmartNoisyOR.
     */
    public String toString() {
        if (isPrior()) {
            return "SmartNoisyOR(" + getVariable().getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            StringBuffer sbuf = new StringBuffer();
            for (int i = 0; i < table.nParents; i++) {
                sbuf.append(table.getParents().get(i).toString() + (i < table.nParents - 1 ? "," : ""));
            }
            return "SmartNoisyOR(" + getVariable().getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        }
    }

    /**
     * Just a pretty-print of the title (can be modified for sub-classes so the
     * tables look nice)
     */
    protected String formatTitle() {
        return String.format(" %10s", var.getName());
    }

    /**
     * Just a pretty-print of the value (can be modified for sub-classes so the
     * tables look nice)
     */
    protected String formatValue(EnumDistrib x) {
        StringBuffer sbuf = new StringBuffer("<");
        double[] distrib = x.get();
        for (int i = 0; i < distrib.length; i++) {
            sbuf.append(String.format("%4.2f ", distrib[i]));
        }
        sbuf.replace(sbuf.length() - 1, sbuf.length() - 1, ">");
        return sbuf.toString();
    }

    /**
     * Pretty-print of whole table
     */
    @Override
    public void print() {
        System.out.println(formatTitle());
        if (!isPrior()) { // variables in condition
            table.display();
        } else { // prior
            if (prior != null) {
                System.out.println(formatValue(prior));
            }
        }
    }

    private Object instance = null;

    /**
     * Set the variable of this SmartNoisyOR to a constant value. This means that parent
     * variables will NOT influence inference.
     *
     * @param value the value that is assigned to this instantiated SmartNoisyOR
     */
    @Override
    public void setInstance(Object value) {
        instance = value;
    }

    /**
     * Set the variable of this SmartNoisyOR to unspecified, or NOT instantiated.
     */
    @Override
    public void resetInstance() {
        instance = null;
    }

    /**
     * Retrieve the instantiated value of this SmartNoisyOR.
     *
     * @return the value of this SmartNoisyOR if instantiated, null if the SmartNoisyOR is not
     * instantiated.
     */
    @Override
    public Object getInstance() {
        return instance;
    }

    /**
     * Count this observation. Note that for it (E-step in EM) to affect the
     * SmartNoisyOR, {@link bn.SmartNoisyOR#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the SmartNoisyOR variable
     * @param prob the expectation of seeing this observation (1 if we actually
     * see it, otherwise the probability)
     * @see bn.SmartNoisyOR#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
    	if (prob == 0.0) {
    		return;
    	}
        if (count == null) { // create count table if none exists
            List<EnumVariable> cond = new ArrayList<EnumVariable>();
            cond.add(var); // first variable is always the conditioned variable
            if (table != null) { // then add parents, if any
                cond.addAll(table.getParents());
            }
            count = new CountTable(cond);
        }
        if (key == null) {
            key = new Object[0];
        }
        Object[] mykey = new Object[key.length + 1];
        mykey[0] = value;
        for (int i = 0; i < key.length; i++) {
            mykey[i + 1] = key[i];
        }
        count.count(mykey, prob);
    }
    
    /**
     * Prob can be set to 1.0 because when counted the value is being observed??
     * Count this observation. Note that for it (E-step in EM) to affect the
     * SmartNoisyOR, {@link bn.SmartNoisyOR#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the SmartNoisyOR variable
     * @see bn.SmartNoisyOR#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value) {
        if (count == null) { // create count table if none exists
            List<EnumVariable> cond = new ArrayList<EnumVariable>();
            cond.add(var); // first variable is always the conditioned variable
            if (table != null) { // then add parents, if any
                cond.addAll(table.getParents());
            }
            count = new CountTable(cond);
        }
        if (key == null) {
            key = new Object[0];
        }
        Object[] mykey = new Object[key.length + 1];
        mykey[0] = value;
        for (int i = 0; i < key.length; i++) {
            mykey[i + 1] = key[i];
        }
        //FIXME - is prob = 1.0 for observed instance accurate?
        count.count(mykey, 1.0);
    }

    /**
     * Take stock of all observations counted via
     * {@link bn.SmartNoisyOR#countInstance(Object[], Object, Double)}, ie implement the
     * M-step locally.
     * Similar to CPT, but only consider observations that match the SmartNoisyOR assumption
     */
    @Override
    public void maximizeInstance() {
        if (count == null) {
            return;
        }
        if (table != null) { // there are parents in the SmartNoisyOR
        	
        	//Set all 'old' distributions in the SmartNoisyOR to valid = false
        	for (EnumDistrib d : this.table.getValues()) {
            	d.setValid(false);
            }      	
        	
            //before calculating counts, need to generate a map of valid keys
        	Hashtable<Object [], Double> keyMap = new Hashtable<Object [], Double>();
        	//will need a 'normal' (i.e. non-weighted) set of counts
        	Hashtable<Object [], Double> trueCountMap = new Hashtable<Object [], Double>();
        	//also need to remember the number of configurations used for calculating each row
        	Hashtable<Object [], Integer> configMap = new Hashtable<Object [], Integer>();
        	
        	//first go through the count table and identify the rows to calculate.
        	//we'll also record the observations for each entry in the count table.
        	for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
        		Object[] cntkey = count.table.getKey(entry.getKey().intValue());
        		double nobserv = entry.getValue().doubleValue();
                Object[] SmartNoisyORkey = new Object[cntkey.length - 1];
                //count the entries in the key
                
                for (int i = 0; i < SmartNoisyORkey.length; i++) {
                	SmartNoisyORkey[i] = cntkey[i + 1];
                }
                int nkey = 0;
                
                for (int i=0; i<SmartNoisyORkey.length; i++) {
                	if (SmartNoisyORkey[i].equals(this.plabels.get(i))) {
                		nkey++;}
                }
                if (nkey <= 1) {
                	//valid key - add it to the map
                	Boolean key_in_map = false;
        			for (Object [] key : keyMap.keySet()) {
        				if (equals(key, cntkey)) key_in_map = true;
        			}
        			if (!key_in_map) {
        				keyMap.put(cntkey, 0.0);
        				configMap.put(cntkey, 0);
        			}
                	//keyMap.put(cntkey, 0.0);
                	trueCountMap.put(cntkey, nobserv);
                	//configMap.put(cntkey, 0);
                }
                else if (nkey > 1) {
                	trueCountMap.put(cntkey, nobserv);
                	//check if the row contains information on a parent variable
                	//that is not already in the maps
                	for (int i=1; i<cntkey.length; i++) {
                		Object [] newkey = new Object [cntkey.length];
                		newkey[0] = cntkey[0];
                		if (cntkey[i].equals(this.plabels.get(i-1))) {
                			newkey[i] = this.plabels.get(i-1);
                			for (int j=1; j<cntkey.length; j++) {
                				if (i!=j) {
                					List<EnumVariable> parents = this.table.getParents();
                					EnumVariable par = parents.get(i-1);
                					Object [] values = par.getDomain().getValues();
                					Object lab = null;
                					for (Object v : values) {
                						if (!v.equals(this.plabels.get(i-1))) {
                							lab = v; break;
                						}
                					}
                					newkey[j] = lab;
                				}
                			}
                			Boolean key_in_map = false;
                			for (Object [] key : keyMap.keySet()) {
                				if (equals(key, newkey)) key_in_map = true;
                			}
                			if (!key_in_map) {
                				keyMap.put(newkey, 0.0);
                				configMap.put(newkey, 0);
                			}
                		}
                	}
                }
        	}
        	
        	//Go through count table again and calculate the weighted counts for each row in the table
            for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
                double nobserv = entry.getValue().doubleValue();
                Object[] cntkey = count.table.getKey(entry.getKey().intValue());
                Object[] SmartNoisyORkey = new Object[cntkey.length - 1];
                //count the entries in the key
                
                for (int i = 0; i < SmartNoisyORkey.length; i++) {
                	SmartNoisyORkey[i] = cntkey[i + 1];
                }
                int nkey = 0;
                
                for (int i=0; i<SmartNoisyORkey.length; i++) {
                	if (SmartNoisyORkey[i].equals(this.plabels.get(i))) {
                		nkey++;}
                }
                if (nkey <= 1) {
                	//case where at most a single parent is switched on
                	int keypos=0;
                	for (int i=1; i<cntkey.length; i++) {
                		if (cntkey[i].equals(this.plabels.get(i-1))) {
                			keypos = i;
                		}
                	}
                	//single parents, add the full observation values
                	Double current_count = null;
                	Integer current_configs = null;
                	Double full_count = 0.0;
                	Object [] keyMapKey = null;
                	Object [] configMapKey = null;
                	for (Object [] key : trueCountMap.keySet()) {
                		//key = null;
                		for (Object [] k : keyMap.keySet()) {
                			if (equals(k, key)) key = k;
                		}
                		int count = 0;
                		for (int i=1; i<key.length; i++) {
                			if (key[i].equals(this.plabels.get(i-1))) count++;
                		}
                		if (key != null && count <=1) {
	                		if (key[0].equals(cntkey[0]) && key[keypos].equals(cntkey[keypos])) {
	                			current_count = keyMap.get(key);
	                			current_configs = configMap.get(key);
	                			configMapKey = key;
	                			keyMapKey = key;
	                			//now we need to count the total number of observations
		                		//of the parent variable
		                		if (key[keypos].equals(cntkey[keypos])) {
		                			for (Object [] k : trueCountMap.keySet()) {
		                				Object [] orkey = new Object [SmartNoisyORkey.length];
		                				for (int i=0; i<orkey.length; i++) {
		                					orkey[i] = k[i + 1];
		                				}
		                				if (equals(orkey, SmartNoisyORkey)) {
		                					 full_count = full_count + trueCountMap.get(k);
		                				}
		                			}
		                		}
	                		}
                		}
                	}
                	if (current_count != null && equals(cntkey, keyMapKey)) {
	                	Double value = nobserv/full_count;
	                	keyMap.put(keyMapKey, current_count + value);
	                	configMap.put(configMapKey, current_configs + 1);
                	}
                }
                else {
                	//more than 1 parent set to true. Divide observations by number of parents.
                	Double current_count = null;
                	Integer current_configs = null;
                	Double full_count = 0.0;
                	for (Object [] key : trueCountMap.keySet()) {
                		Object [] configkey = new Object [SmartNoisyORkey.length];
                		for (int i=0; i<configkey.length; i++) {
                			configkey[i] = key[i + 1];
                		}
                		if (equals(configkey, SmartNoisyORkey)) {
                			full_count = full_count + trueCountMap.get(key);
                		}
                	}
                	
                	for (int i=1; i<cntkey.length; i++) {
                		Object thisLabel = null;
                		for (int j=1; j<cntkey.length; j++) {
                			if (i==j && cntkey[j].equals(this.plabels.get(j-1))) {
                				thisLabel = this.plabels.get(j-1);}
                		}
                		//need to get each valid key
                		if (thisLabel != null) {
                			for (Object [] key : keyMap.keySet()) {
                				if (key[i].equals(thisLabel) && key[0].equals(cntkey[0])) {
                					current_count = keyMap.get(key);
                					Double value = nobserv/(full_count*nkey);
                					keyMap.put(key,  current_count + value);
                					current_configs = configMap.get(key);
                					configMap.put(key,  current_configs + 1);
                				}
                			}
                		}
                	}
                	
                }
            } 
            //now to add the counts to the table
            for (Object [] cntkey : keyMap.keySet()) {
            	Object[] SmartNoisyORkey = new Object[cntkey.length - 1];
                //count the entries in the key
                for (int i = 0; i < SmartNoisyORkey.length; i++) {
                	SmartNoisyORkey[i] = cntkey[i + 1];
                }
                int index = 0;
                for (int i=1; i<cntkey.length; i++) {
                	if (cntkey[i].equals(this.plabels.get(i-1))) {
                		index = i;
                	}
                }
                if (index > 0) {
	                EnumDistrib d = this.table.getValue(SmartNoisyORkey);
	                if (d == null) {
	                	d = new EnumDistrib(var.getDomain());
	                	Double weighted_count = keyMap.get(cntkey);
	                	Integer nconfigs = configMap.get(cntkey);
	                	d.set(cntkey[0], weighted_count/nconfigs);
	                	this.put(SmartNoisyORkey, d);
	                } else {
	                	Double weighted_count = keyMap.get(cntkey);
	                	Integer nconfigs = configMap.get(cntkey);
	                	d.set(cntkey[0], weighted_count/nconfigs);
	                	this.put(SmartNoisyORkey, d);
	                }
                }
                else { //all parent variables set to negatives
                	EnumDistrib d = this.table.getValue(SmartNoisyORkey);
                	if (d == null) {
                		d = new EnumDistrib(var.getDomain());
                		Double obsCount = trueCountMap.get(cntkey);
                		d.set(cntkey[0], obsCount);
	                	this.put(SmartNoisyORkey, d);
                	}
                	else {
                		Double obsCount = null;
	                	for (Object [] k : trueCountMap.keySet()) {
	                		if (equals(cntkey, k)) {
	                			obsCount = trueCountMap.get(k);}
	                	}
                		d.set(cntkey[0], obsCount);
	                	this.put(SmartNoisyORkey, d);
                	}
                }
            }
            
            /**
             * TO BE REMOVED when pseudo counts are implemented.
             * At this point, some rows may not have been assigned values 
             * if a parent was not observed in the training data. To handle this,
             * rows that currently have null values will be set with uniform distributions. 
             */
            
            Enumerable dom = var.getDomain();
            for (int i=0; i<plabels.size(); i++) {
            	Object [] key = new Object [plabels.size()];
            	key[i] = plabels.get(i);
            	for (int j=0; j<plabels.size(); j++) {
            		if (i!=j) {
    		        	ArrayList<EnumVariable> parents = (ArrayList<EnumVariable>) getParents();
    					EnumVariable par = parents.get(i);
    					Object [] domainValues = par.getDomain().getValues();
    					Object lab = null;
    					for (Object v : domainValues) {
    						if (!v.equals(plabels.get(i))) {
    							lab = v; break;
    						}
    					}
    					key[j] = lab; 
            		}
    	        }
            	//try and pull the key out of the table
            	Object value = table.getValue(key);
            	if (value == null) {
            		Map<Object, Double> distMap = new HashMap<Object, Double>();
            		for (int x=0; x<dom.size(); x++) {
            			distMap.put(dom.get(x), 1.0/dom.size());
            		}
            		EnumDistrib dist = new EnumDistrib(distMap);
            		put(key, dist);
            	}
            }
            
            // normalisation happens internally when values are required	

            //Remove 'ghost' entries from CPT (which have not been made valid above)
            for (Iterator<Entry<Integer, EnumDistrib>> it = table.getMapEntries().iterator(); it.hasNext(); ) {
                Entry<Integer, EnumDistrib> entry = it.next();
            	EnumDistrib obs = entry.getValue();
            	if (!obs.isValid()) 
                    it.remove();
            }
        } else { // there are no parents
            Object[] cntkey = new Object[1];
            double[] cnts = new double[var.size()];
            for (int i = 0; i < var.size(); i++) {
                cntkey[0] = var.getDomain().get(i);
                cnts[i] = count.get(cntkey);
            }
            prior = new EnumDistrib(this.var.getDomain(), cnts);	// EnumDistrib normalises the counts internally
        }
        count = null; // reset counts
    }

    protected CountTable getCount() {
        return count;
    }

    /**
     * Put random entries in the SmartNoisyOR if not already set.
     */
    public void randomize(long seed) {
        Random rand = new Random(seed);
        if (table == null) {
            if (prior == null)
                prior = EnumDistrib.random(var.getDomain());
        } else {
            int nrows = table.getSize();
            for (int i = 0; i < nrows; i++) {
            	//want to check that this row is valid for SmartNoisyOR
            	Object [] key = table.getKey(i);
            	int nkey = 0;
            	for (int j=0; j<key.length; j++) {
            		if (key[j].equals(this.plabels.get(j))) {
            			nkey++;}
            	}
            	if (nkey <= 1) {
                if (!table.hasValue(i)) {
                    table.setValue(i, EnumDistrib.random(var.getDomain()));}
            	}
            }
        }
    }

    /**
     * Set this SmartNoisyOR to be trained when the Bayesian network it is part of is
     * trained. A SmartNoisyOR is trainable (true) by default.
     *
     * @param status true if trainable, false otherwise
     */
    public void setTrainable(boolean status) {
        trainable = status;
    }

    protected boolean trainable = true;

    /**
     * Check if this SmartNoisyOR should be trained or not
     */
    public boolean isTrainable() {
        return trainable;
    }

    @Override
    public String getStateAsText() {
        StringBuffer sbuf = new StringBuffer("\n");
        if (isPrior()) {
            EnumDistrib d = prior;
            if (d != null) {
                double[] distrib = d.get();
                for (int j = 0; j < distrib.length; j++) {
                    sbuf.append("" + distrib[j]);
                    if (j < distrib.length - 1) {
                        sbuf.append(", ");
                    }
                }
                sbuf.append(";\n");
            }
        } else {
            for (int i = 0; i < table.getSize(); i++) {
                EnumDistrib d = table.getValue(i);
                if (d != null) {
                    double[] distrib = d.get();
                    sbuf.append(i + ": ");	// use index as key because values above can be of different non-printable types
                    for (int j = 0; j < distrib.length; j++) {
                        sbuf.append("" + distrib[j]);
                        if (j < distrib.length - 1) {
                            sbuf.append(", ");
                        }
                    }
                    sbuf.append("; (");
                    // If we want to *see* the key, may not work well for some non-printable types
                    Object[] key = table.getKey(i);
                    for (int j = 0; j < key.length; j++) {
                        if (j < key.length - 1) {
                            sbuf.append(key[j] + ", ");
                        } else {
                            sbuf.append(key[j] + ")\n");
                        }
                    }
                }
            }
        }
        return sbuf.toString();
    }

    @Override
    public boolean setState(String dump) {
        if (isPrior()) {
            String[] line = dump.split(";");
            if (line.length >= 1) {
                String[] y = line[0].split(",");
                if (y.length == var.size()) {
                    double[] distrib = new double[y.length];
                    try {
                        for (int i = 0; i < distrib.length; i++) {
                            distrib[i] = Double.parseDouble(y[i]);
                        }
                    } catch (NumberFormatException e) {
                        e.printStackTrace();
                        return false;
                    }
                    this.put(new EnumDistrib(var.getDomain(), distrib));
                    return true;
                }
            }
        } else {
            for (String line : dump.split("\n")) {
                // 0: 0.4, 0.6; (true, true)
                String[] specline = line.split(";");
                if (specline.length >= 1) {
                    String[] parts = specline[0].split(":");
                    if (parts.length >= 2) {
                        try {
                            int index = Integer.parseInt(parts[0]);
                            String[] y = parts[1].split(",");
                            if (y.length == var.size()) {
                                double[] distrib = new double[y.length];
                                try {
                                    for (int i = 0; i < distrib.length; i++) {
                                        distrib[i] = Double.parseDouble(y[i]);
                                    }
                                } catch (NumberFormatException e) {
                                    e.printStackTrace();
                                    return false;
                                }
                                this.put(table.getKey(index), new EnumDistrib(var.getDomain(), distrib));
                            }
                        } catch (NumberFormatException e) {
                            System.err.println("Number format wrong and ignored: " + line);
                        }
                    }
                }
            }
        }
        return false;
    }

    @Override
    public String getType() {
        return "SmartNoisyOR";
    }
    
	public boolean isRelevant() {
		return relevant;
	}

	public void setRelevant(boolean relevant) {
		this.relevant = relevant;
	}
	
	public static boolean equals(final Object[] a, final Object[] b){
		
		if (a == b) return true;
		if (a == null || b == null) return false;
		if (a.length != b.length) return false;
			
		for (int i=0; i<a.length; i++) {
			Object x = a[i];
			Object y = b[i];
			
			if (x != y) return false;
			if (x == null || y == null) return false;
			//if (deep) {
			//	if (x instanceof Object[] && y instanceof Object[]) {
			//		if (! equals((Object[])x, (Object[])y, true)) return false;
			//}
			//else {
			//	return false;
			//	}
			//}
			if (! x.equals(y)) return false;
		}
		
		return true;
}

    /**
     * @param args
     */
    public static void main(String[] args) {
    	Object [] key1 = new Object [] {true, false, false};
    	Object [] key2 = new Object [] {true, false, false};
    	System.out.println(equals(key1, key2));
    	
    	
        EnumVariable v1 = Predef.Boolean();
        EnumVariable v2 = Predef.Boolean();
        EnumVariable v3 = Predef.Boolean();

        SmartNoisyOR SmartNoisyOR1 = new SmartNoisyOR(v1, new EnumVariable[]{v2, v3}, new Object [] {true, true});
        SmartNoisyOR1.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{0.78, 0.22}));
        SmartNoisyOR1.put(new Object[]{false, true}, new EnumDistrib(v2.getDomain(), new double[]{0.6, 0.4}));
        SmartNoisyOR1.put(new Object[]{false, false}, new EnumDistrib(v2.getDomain(), new double[]{0.15, 0.85}));
        SmartNoisyOR1.print();
       
    }

	@Override
	public Double get(Object value) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public EnumTable getTable() {
		// TODO Auto-generated method stub
		return table;
	}

	@Override
	public Distrib getDistrib() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Sample> getConditionDataset(int conditionIndex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Distrib getlikelihoodDistrib() {
		// TODO Auto-generated method stub
		return null;
	}

}
