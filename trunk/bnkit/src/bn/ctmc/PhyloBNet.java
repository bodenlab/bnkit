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
import dat.PhyloTree;
import dat.PhyloTree.Node;
import dat.Variable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Class for a Bayesian network that represent a phylogenetic tree.
 * @author mikael
 */
public class PhyloBNet {
    
    private final BNet bn;
    private final SubstModel model;
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
    
    /**
     * Construct a BN for specified phylogenetic tree using supplied model.
     * @param tree phylogenetic tree
     * @param model evolutionary model
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBNet create(PhyloTree tree, SubstModel model) {
        PhyloBNet pbn = new PhyloBNet(model);
        Node root = tree.getRoot();
//        EnumVariable rvar = Predef.AminoAcid(root.getContent().toString());
        EnumVariable rvar = Predef.AminoAcid(replacePunct(root.toString()));
        pbn.bnroot = new SubstNode(rvar, model);
        pbn.addBNode(pbn.bnroot);
        pbn.createNodesForSubtree(root, rvar);
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

    private static String replacePunct(String str) {
        return str.replace('.', '_');
    }

    private void createNodesForSubtree(Node pnode, EnumVariable evar) {
        Set<Node> children = pnode.getChildren();
        if (children.isEmpty()) {
            leaves.add(evar);
        } else {
            for (Node child : children) {
//                EnumVariable cvar = Predef.AminoAcid(child.getContent().toString());
                EnumVariable cvar = Predef.AminoAcid(replacePunct(child.toString()));
                SubstNode cnode = new SubstNode(cvar, evar, model, child.getDistance());
                this.addBNode(cnode);
                createNodesForSubtree(child, cvar);
            }
        }
    }
    
    public int purgeGaps() {
        List<BNode> nodes = bn.getOrdered();
        if (nodes == null)
            return 0;
        Set<BNode> purgedAlready = new HashSet<>();
        for (BNode node : nodes) {
            if (node.getInstance() == null && !bn.hasChildren(node)) { // here's a leaf node which is not instantiated
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
                if (siblings == null) // the parent has no siblings either, so will need to defer until we investigate its grandparent
                    topNode = true;
                else
                    topNode = (siblings.size() > 0);
                if (topNode) {
                    // we have established that "parent" is at the top level to collapse/remove
                    toBePurged.add(parent);
                    double time = parent.getTime();
                    Object instance = null;
                    // Let's look at the (only) child
                    SubstNode child = (SubstNode)(children.toArray()[0]);
                    while (true) {
                        instance = child.getInstance();
                        time += child.getTime();
                        toBePurged.add(child);
                        Set<BNode> childrensChildren = bn.getChildren(child);
                        if (childrensChildren == null)
                            break;
                        if (childrensChildren.size() != 1)
                            break;
                        child = (SubstNode)(childrensChildren.toArray()[0]);
                    } 
                    EnumVariable childVar = (EnumVariable)child.getVariable();
                    SubstNode replacement;
                    List<EnumVariable> superParentVars = parent.getParents();
                    if (superParentVars != null) {
                        EnumVariable superParentVar = superParentVars.get(0); // can only be one
                        replacement = new SubstNode(childVar, superParentVar, child.getModel(), time);
                    } else {
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
