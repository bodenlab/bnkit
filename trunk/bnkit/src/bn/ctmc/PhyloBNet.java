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
        EnumVariable rvar = Predef.AminoAcid(root.getContent().toString());
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
    
    private void createNodesForSubtree(Node pnode, EnumVariable evar) {
        Set<Node> children = pnode.getChildren();
        if (children.isEmpty()) {
            leaves.add(evar);
        } else {
            for (Node child : children) {
                EnumVariable cvar = Predef.AminoAcid(child.getContent().toString());
                SubstNode cnode = new SubstNode(cvar, evar, model, child.getDistance());
                this.addBNode(cnode);
                createNodesForSubtree(child, cvar);
            }
        }
    }
    
    public void purgeGaps() {
        List<BNode> nodes = bn.getOrdered();
        if (nodes == null)
            return;
        for (BNode node : nodes) {
            if (node.getInstance() == null && !bn.hasChildren(node)) { // here's a leaf node which is not instantiated
                
                
            }
        }
    }
    
}
