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

import bn.CountTable;
import bn.JPT;
import dat.EnumVariable;
import dat.Enumerable;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;


/**
 * Class for Pseudo Conditional Probability Table (CPT). This is a table for a
 * conditioned enumerable variable, that has any number (incl 0) enumerable
 * parents.
 *
 * This class is an alternative to the current CPT class, allowing the
 * addition of pseudo counts for the node. Essentially identical to te current CPT implementation,
 * this class also allows users to provide pseudocounts in the form of a PseudoMatrix object.
 * These pseudocounts are added to the CountTable `count' when countinstance is first called.
 * Currently this only applies to the main countinstance method.
 *
 * Given this is a rather `hacked' together implementation, there are several issues that need to be accounted
 * for when using this:
 *      - If this node has many parents (or even a few with many classes), searching and
 *      populating the count table will be incredibly computationally expensive. Additionally, it is unlikely users will
 *      be able to provide a full pseudocounts matrix accounting for every possibility. Hence, the user is required to
 *      set a main parent if the node has more than one parent. Essentially this is equivalent to saying "This parent
 *      is the main parent/variable which is influencing this child node and I only want to apply pseudo counts
 *      specifically in respect to this relationship". This is done with - setMainParentIndex(String name).
 *      This parent should be the one whose domain entries correspond respectively to each row in a provided PseudoMatrix.
 *      (Where columns should correspond to each possible observation across the child's domain).
 *
 *
 * @author julian
 */

public class CPTPseudo extends CPT{
    private static final long serialVersionUID = 1L;
    // pseudo count matrix; each row correlates to the parent key state and each column the state of this node
    private PseudoMatrix pseudoMatrix; //stores the pseudo counts
    private Integer mainparent_index;
    private boolean relevant = false; //for inference, track whether the node is relevant to the query

    /**
     * Create a conditional probability table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public CPTPseudo(EnumVariable var, List<EnumVariable> parents, PseudoMatrix pseudo) {
        super(var,parents);
        this.pseudoMatrix = pseudo;
    }

    public CPTPseudo(EnumVariable var, List<EnumVariable> parents) {
        super(var,parents);
        this.pseudoMatrix = null;
    }

    /**
     * Create a conditional probability table for a variable. The variable is
     * conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parents parent variables
     */
    public CPTPseudo(EnumVariable var, PseudoMatrix pseudo, EnumVariable... parents) {
        super(var,parents);
        this.pseudoMatrix = pseudo;
    }

    public CPTPseudo(EnumVariable var, EnumVariable... parents){
        this(var, null, parents);
    }

    /**
     * Create a prior (a CPT without conditioning variables) for a variable.
     *
     * @param var variable
     */
    public CPTPseudo(EnumVariable var, PseudoMatrix pseudo) {
        super(var);
        this.pseudoMatrix = pseudo;
    }

    /**
     * Create a CPT from a JPT. The variable in the JPT var is the variable
     * conditioned on in the CPT.
     *
     * @param jpt
     * @param var
     */
    public CPTPseudo(JPT jpt, EnumVariable var, PseudoMatrix pseudo) {
        super(jpt,var);
        this.pseudoMatrix = pseudo;
    }

    /**
     * Create a CPTPseudo from a JPT. The variable in the JPT var is the variable
     * conditioned on in the CPTPseudo.
     *
     * @param jpt
     * @param var
     */
    public CPTPseudo(JPT jpt, EnumVariable var) {
        this(jpt,var,null);
    }

    /**
     * Set the pseudo matrix
     * @param pseudo
     */
    public void setPseudo(PseudoMatrix pseudo){
        this.pseudoMatrix = pseudo;
    }

    /**
     * Provide a non-unique string representation of this CPTPseudo.
     */
    @Override
    public String toString() {
        if (isPrior()) {
            return "CPTPseudo(" + getVariable().getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            StringBuilder sbuf = new StringBuilder();
            for (int i = 0; i < this.getTable().nParents; i++) {
                sbuf.append(this.getTable().getParents().get(i).toString()).append(i < this.getTable().nParents - 1 ? "," : "");
            }
            return "CPTPseudo(" + getVariable().getName() + "|" + sbuf.toString() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        }
    }


    /**
     * Set which parent is the main parent; Stores the index.
     * We rely on the indexing to be consistent
     * @param mainParentName the name of the main parent for this node
     */
    public void setMainParentIndex(String mainParentName) {
        List<EnumVariable> parents = this.getParents();
        if (parents == null){
            throw new RuntimeException("Can't set a main parent, this node has no parents!");
        } else{
            for (int i = 0; i < parents.size(); i++){
                EnumVariable aparent = parents.get(i);
                if (mainParentName.equals(aparent.getName())) {
                    this.mainparent_index = i;
                    break;
                }
            }
        }
    }

    public Integer getMainParentIndex(){
        if (this.mainparent_index != null || isPrior()) {
            //if main parent has been set, do nothing. Indexing will be correct.
        }
        //not set but one parent, so just automatically set
        else if (this.mainparent_index == null && this.getParents().size() == 1) {
            this.mainparent_index = 0;
        } else {
            throw new RuntimeException("No main parent has been specified and number of" +
                    " parents is more than one; use set MainParentIndex");
        }
        return this.mainparent_index;
    }

    /**
     * Count this observation. Note that for it (E-step in EM) to affect the
     * CPT, {@link CPTPseudo#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the CPT variable
     * @param prob the expectation of seeing this observation (1 if we actually
     * see it, otherwise the probability)
     * @see CPTPseudo#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        if (prob == 0.0) {
            return;
        }
        CountTable count = this.getCount();
        if (this.getCount().table.isEmpty()) {
            //        System.out.println();
            //CPTPseudo specific. Here the count table is initialized with pseudo counts
            //Domain lengths determine the matrix[i][j]...up to user to supply correctly formatted pseudo matrix
            Integer p_idx = getMainParentIndex();
            Enumerable cdom = this.getVariable().getDomain(); //child (this node's) domain
            if (key == null) { // if the node is a root
                //then create new key of length 1
                Object[] newKey = new Object[1];
                for (int j = 0; j < cdom.size(); j++) { //go through child domain
                    Object co = cdom.get(j); // child observation
                    double obsCount = this.pseudoMatrix.getValue(0, j); //the count for the child observation
                    newKey[0] = co;
                    count.count(newKey, obsCount);
                }
            } else {
                //get the parent domain. p_idx will have to be != 0 if more than one parent
                Enumerable pdom = this.getParents().get(p_idx).getDomain(); //parent domain
                for (int i = 0; i < pdom.size(); i++) {
                    Object po = pdom.get(i); //parent observation
                    for (int j = 0; j < cdom.size(); j++) {
                        Object co = cdom.get(j); // child observation
                        double obsCount = this.pseudoMatrix.getValue(i, j); //the count
                        //                        double obsCount = this.pseudoMatrix.getValue(i, co); //Don't use this...specific use only
                        // add one as count table key includes child observation
                        Object[] newKey = new Object[key.length + 1];
                        newKey[0] = co;
                        for (int x = 1; x < newKey.length; x++) {
                            if (x == p_idx + 1) {
                                newKey[x] = po; //We use the observation for the main parent
                            } else {
                                newKey[x] = null; //All other parent's keys are set to null
                            }
                        }
                        //get all possible indices for the marginalised key
                        int[] countable_idxs = count.table.getTheoreticalIndices(newKey);
                        //get all possible indices for the marginalised key
                        //for each index, add the corresponding count
                        for (int x = 0; x < countable_idxs.length; x++) {
                            count.count(countable_idxs[x], obsCount);
                        }
                    }
                }
            }
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
     * CPT, {@link @see CPTPseudo#maximizeInstance()} must be called.
     *
     * @param key the setting of the parent variables in the observation
     * @param value the setting of the CPT variable
     *
     * @see CPTPseudo#maximizeInstance()
     */
    @Override
    public void countInstance(Object[] key, Object value) {
        CountTable count = this.getCount();
        if (count.table.isEmpty()) { // create count table if none exists
            //CPTPseudo specific. Here the count table is initialized with pseudo counts
            //Domain lengths determine the matrix[i][j]...up to user to supply correctly formatted pseudo matrix
            Integer p_idx = getMainParentIndex();
            Enumerable cdom = this.getVariable().getDomain(); //child (this node's) domain
            if (key == null) { // if the node is a root
                //then create new key of length 1
                Object[] newKey = new Object[1];
                for (int j = 0; j < cdom.size(); j++){ //go through child domain
                    Object co = cdom.get(j); // child observation
                    double obsCount = this.pseudoMatrix.getValue(0, j); //the count for the child observation
                    newKey[0] = co;
                    count.count(newKey, obsCount);
                }
            } else {
                //get the parent domain. p_idx will have to be != 0 if more than one parent
                Enumerable pdom = this.getParents().get(p_idx).getDomain(); //parent domain
                for (int i = 0; i < pdom.size(); i++) {
                    Object po = pdom.get(i); //parent observation
                    for (int j = 0; j < cdom.size(); j++) {
                        Object co = cdom.get(j); // child observation
                        double obsCount = this.pseudoMatrix.getValue(i, j); //the count
//                        double obsCount = this.pseudoMatrix.getValue(i, co);
                        // add one as count table key includes child observation
                        Object[] newKey = new Object[key.length + 1];
                        newKey[0] = co;
                        for (int x = 1; x < newKey.length; x++) {
                            if (x == p_idx + 1) {
                                newKey[x] = po; //We use the observation for the main parent
                            } else {
                                newKey[x] = null; //All other parent's keys are set to null
                            }
                        }
                        //get all possible indices for the marginalised key
                        int[] countable_idxs = count.table.getTheoreticalIndices(newKey);
                        //get all possible indices for the marginalised key
                        //for each index, add the corresponding count
                        for (int x = 0; x < countable_idxs.length; x++) {
                            count.count(countable_idxs[x], obsCount);
                        }
                    }
                }
            }
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

    @Override
    public String getType() {
        return "CPTPseudo";
    }
//
//    /**
//     * Tie all parameters essential to inference and training for this CPTPseudo to those of another CPTPseudo.
//     * Variables should be separate but they are required to (1) be of the same type/domain, and (2) be listed in the same order.
//     * @param source the CPTPseudo from which parameters will be copied and held fixed.
//     */
//    @Override
//    public void tieTo(CPTPseudo source) {
//        CPTPseudo src = (CPTPseudo)source;
//        if (!this.var.getDomain().equals(source.getVariable().getDomain()))
//            throw new RuntimeException("Invalid sharing: " + var.getName() + " does not share domain with " + source.getVariable().getName());
//        if (this.nParents != src.nParents)
//            throw new RuntimeException("Invalid sharing: " + var.getName() + " has different number of parents from " + source.getVariable().getName());
//        for (int i = 0; i < this.nParents; i ++) {
//            Variable p1 = this.getParents().get(i);
//            Variable p2 = src.getParents().get(i);
//            if (!p1.getDomain().equals(p2.getDomain()))
//                throw new RuntimeException("Invalid sharing: " + p1.getName() + " does not share domain with " + p2.getName());
//        }
//        // need to tie:
//        // - count (used during learning)
//        // - prior (if applicable)
//        // - table (if applicable)
//        this.prior = source.prior;
//        if (this.nParents > 0)
//            this.table = source.table.retrofit(this.getParents());
//    }
//
}

