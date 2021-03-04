/*
 * bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2014  M. Boden et al.
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
package bn.ctmc;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import dat.EnumVariable;
import dat.Enumerable;
import dat.PhyloTree;
import dat.PhyloTree.Node;
import dat.Variable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Class for a Bayesian network that represent a phylogenetic tree.
 * @author mikael
 * @deprecated
 */
public class PhyloBNet {

    private final BNet bn;
    private final SubstModel model;
    private double rate = 1.0;
    private List<EnumVariable> leaves;
    private SubstNode bnroot;

    private PhyloBNet(SubstModel model) {
        bn = new BNet();
        this.leaves = new ArrayList<>();
        this.model = model;
    }

    public BNet getBN() {
        return bn;
    }

    public BNode getRoot() {
        return bnroot;
    }

    protected void addBNode(BNode node) {
        bn.add(node);
    }

    protected void removeBNode(BNode node){
        bn.remove(node);
    }

    /**
     * Construct a BN for specified phylogenetic tree using supplied model.
     * @param tree phylogenetic tree
     * @param model evolutionary model
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBNet create(PhyloTree tree, SubstModel model) {
        return create(tree, model, 1.0);
    }

    /**
     * Construct a BN for specified phylogenetic tree using supplied model.
     * @param tree phylogenetic tree
     * @param model evolutionary model
     * @param rate the evolutionary rate to be applied
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBNet create(PhyloTree tree, SubstModel model, double rate) {
        PhyloBNet pbn = new PhyloBNet(model);
        pbn.rate = rate;
        Node root = tree.getRoot();
        EnumVariable rvar = new EnumVariable(model.getDomain(), root.getLabel().toString());
        pbn.bnroot = new SubstNode(rvar, model);
        pbn.addBNode(pbn.bnroot);
        pbn.createNodesForSubtree(root, rvar);
        return pbn;
    }



    /**
     * Construct a BN for specified phylogenetic tree using supplied model.
     * @param tree phylogenetic tree
     * @param model evolutionary model
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBNet createGap(PhyloTree tree, SubstModel model) {
        return createGap(tree, model, 1.0);
    }

    /**
     * Construct a BN for specified phylogenetic tree using supplied model and GAP CHARACTER ALPHABET
     * @param tree phylogenetic tree
     * @param model evolutionary model
     * @param rate the evolutionary rate to be applied
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBNet createGap(PhyloTree tree, SubstModel model, double rate) {
        PhyloBNet pbn = new PhyloBNet(model);
        pbn.rate = rate;
        Node root = tree.getRoot();
        EnumVariable rvar = Predef.GapCharacter(root.getLabel().toString());
//        EnumVariable rvar = Predef.AminoAcid(root.toString());
        pbn.bnroot = new SubstNode(rvar, model);
        pbn.addBNode(pbn.bnroot);
        pbn.createNodesForSubtreeGap(root, rvar);
        return pbn;
    }

    /**
     * Get variables that are found at the leaf nodes.
     * @return list of leaf variables
     */
    public List<EnumVariable> getLeaves() {
        return leaves;
    }

    public List<EnumVariable> getInternal() {
        List<Variable> vars = bn.getOrderedVariables();
        List<EnumVariable> internal = new ArrayList<>();
        for (Variable var : vars) {
            try {
                EnumVariable evar = (EnumVariable) var;
                if (!leaves.contains(evar))
                    internal.add(evar);
            } catch (ClassCastException ex) {
            }
        }
        return internal;
    }

    /**
     * Get the rate of change as determined from instantiated nodes.
     * @return rate of change relative evolutionary time
     */
    public double getRate() {
        double rate = 0;
        int nBranches = 0;
        for (BNode node : bn.getNodes()) {
            // only consider leaf nodes
            if (!bn.hasChildren(node)) {
                int nChanges = 0; // no of different instantiations
                double cumTime = 0; // cumulative time according to nodes
                SubstNode snode = null;
                Object prev = node.getInstance();
                try {
                    snode = (SubstNode)node;
                } catch (ClassCastException e) {
                    System.err.println("Not supported node type for PhyloBNet, whole branch ignored");
                }
                do {
                    try {
                        Object inst = snode.getInstance();
                        if (inst != prev && inst != null) {
                            prev = snode.getInstance();
                            nChanges ++;
                        }
                        cumTime += snode.getTime();
                        Set<BNode> parents = bn.getParents(snode);
                        if (parents == null)
                            break;
                        if (parents.isEmpty())
                            break;
                        for (BNode parent : parents) // there can only be one parent
                            snode = (SubstNode)parent;
                    } catch (ClassCastException e) {
                        System.err.println("Not supported node type for PhyloBNet, whole branch ignored");
                    }
                } while (snode != null);
                rate += (nChanges + 0.00001) / cumTime;
                nBranches ++;
            }
        }
        if (nBranches > 0)
            return rate / nBranches;
        else
            return 0;
    }

    private void createNodesForSubtree(Node pnode, EnumVariable evar) {
        Collection<Node> children = pnode.getChildren();
        if (children.isEmpty()) {
            leaves.add(evar);

        } else {
            for (Node child : children) {
                EnumVariable cvar = new EnumVariable(model.getDomain(), child.getLabel().toString());
                SubstNode cnode = new SubstNode(cvar, evar, model, child.getDistance() * this.rate);
                this.addBNode(cnode);
                createNodesForSubtree(child, cvar);
            }
        }
    }

    private void createNodesForSubtreeGap(Node pnode, EnumVariable evar) {
        Collection<Node> children = pnode.getChildren();
        if (children.isEmpty()) {
            leaves.add(evar);
        } else {
            for (Node child : children) {
                EnumVariable cvar = Predef.GapCharacter(child.getLabel().toString());
//                EnumVariable cvar = Predef.AminoAcid(child.toString());
                SubstNode cnode = new SubstNode(cvar, evar, model, child.getDistance() * this.rate);
                this.addBNode(cnode);
                createNodesForSubtreeGap(child, cvar);
            }
        }
    }

    public int purgeGaps() {
        List<BNode> nodes = bn.getOrdered();
        if (nodes == null)
            return 0;
        Set<BNode> purgedAlready = new HashSet<>();
        for (BNode node : nodes) {
            SubstNode snode = (SubstNode)node; //Cast to subst node
//            if (node.getInstance() == null && !bn.hasChildren(node)) { // here's a leaf node which is not instantiated
            if (snode.getGap() && !bn.hasChildren(node)) { // here's a leaf node which is not instantiated
                purgedAlready.add(node);
                purgeMe(node, purgedAlready);
            }
        }
        for (BNode purge : purgedAlready)
            bn.remove(purge);
        bn.compile();
        return purgedAlready.size();
    }

    /**
     * Identify the top-most nodes above the specified, that can be purged safely.
     * This will only be one in the case of phylogenetic trees.
     * @param me
     * @return
     */
    private void purgeMe(BNode me, Set<BNode> purgedAlready) {
        Set<BNode> parents = bn.getParents(me); // should only be one in the case of a phylogenetic tree
        if (parents == null)
            return;
        else {
            for (BNode p : parents) { // only one
                Set<BNode> children = bn.getChildren(p);
                int purged = 0;
                for (BNode c : children) {
                    if (purgedAlready.contains(c))
                        purged ++;
                }
                if (children.size() - purged == 0) { // any other children (than me)?
                    purgedAlready.add(p);              // if not, we can continue up the tree
                    purgeMe(p, purgedAlready);
                }
            }
        }
    }

    /**
     * Remove internal nodes that are singly connected, with exactly one child.
     * @return number of collapsed nodes
     */
    public int collapseSingles() {
        List<BNode> nodes = bn.getOrdered();
        if (nodes == null)
            return 0;
        Set<BNode> toBePurged = new HashSet<>();
        Set<BNode> toBeAdded = new HashSet<>();
        for (BNode node : nodes) {
            SubstNode parent = (SubstNode)node;
            if (parent.getInstance() != null)
                continue;
            Set<BNode> children = bn.getChildren(parent);
            if (children == null)
                continue;
            if (children.size() == 1) { // a "single-child" parent
                Set<BNode> siblings = bn.getSiblings(parent);
                boolean topNode;
                if (siblings == null) // the parent has no siblings, because it has no parents of its own, so do not look above
                    topNode = true;
                else
                    topNode = (siblings.size() > 0); // the parent has siblings, so we do not need to look above
                if (topNode && parent.getGap()) { //only collapse if the parent is recorded as being a gap
                    // we have established that "parent" is at the top level to collapse/remove
                    toBePurged.add(parent);
                    // we need to remember details because we are going to create a replacement node
                    double time = parent.getTime();
                    Object instance = null;
                    // let's look at the (only) child
                    SubstNode child = (SubstNode)(children.toArray()[0]);
                    while (true) {
                        instance = child.getInstance(); // use instance of the most specific child
                        time += child.getTime(); // add to evolutionary time for every node that is purged
                        toBePurged.add(child);
                        Set<BNode> childrensChildren = bn.getChildren(child); // look ahead
                        if (childrensChildren == null) // no more nodes below
                            break;
                        if (childrensChildren.size() != 1) // more nodes below, but they cannot be purged
                            break;
                        child = (SubstNode)(childrensChildren.toArray()[0]); // only one child, let's look at that...
                    }
                    // found the most specific child in a chain of single-child parents
                    // find more details to construct replacement
                    EnumVariable childVar = (EnumVariable)child.getVariable();
                    SubstNode replacement;
                    List<EnumVariable> superParentVars = parent.getParents();
                    if (superParentVars != null) { // the parent had parents of its own
                        EnumVariable superParentVar = superParentVars.get(0); // can only be one
                        replacement = new SubstNode(childVar, superParentVar, child.getModel(), time);
                    } else { // the parent was the root of the tree
                        replacement = new SubstNode(childVar, child.getModel());
                    }
                    replacement.setInstance(instance);
                    toBeAdded.add(replacement);
                }
            }
        }
        for (BNode node : toBePurged)
            bn.remove(node);
        for (BNode node : toBeAdded)
            bn.add(node);
        bn.compile();
        return toBeAdded.size();
    }
}
