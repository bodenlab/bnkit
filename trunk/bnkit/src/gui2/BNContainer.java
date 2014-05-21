/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui2;

import bn.BNet;
import bn.BNode;
import bn.EnumVariable;
import bn.Predef;
import bn.Variable;
import bn.file.BNBuf;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Container of things related to the BN that is under construction. 
 * The BNet object is constructed only when needed and is not maintained in this class.
 * @author mikael
 */
public class BNContainer {
    private Map<String, Variable> vars = new HashMap<String, Variable>(); // all variables that are considered and can be used to create nodes
    private Map<String, BNode> nodes = new HashMap<String, BNode>(); // all nodes that are considered, will change frequently as the graph is constructed
    
    /**
     * Retrieve the BNet structure for the current set of nodes.
     * Is not maintained internally.
     */
    public BNet getBNet() {
        BNet bnet = new BNet();
        for (Map.Entry<String, BNode> entry : nodes.entrySet())
            bnet.add(entry.getValue());
        return bnet;
    }
    
    /**
     * Add variable to container, replace if the same variable already exists.
     * @param var 
     */
    public void addVariable(Variable var) {
        vars.put(var.getName(), var);
    }
    
    public Variable getVariable(String name) {
        return vars.get(name);
    }
    
    public void addNode(BNode node) {
        System.out.println("in BNContainer.java, adding node:" +
                node.getName() + "::" + node + " to bnc");
        if (!vars.containsKey(node.getVariable().getName()))
            addVariable(node.getVariable());
        nodes.put(node.getName(), node);
    }
    
    public BNode getNode(String name) {
        return nodes.get(name);
    }
    
    public void removeNode(BNode node) {
        // remove actual node (and variable)
        if (node == null) return;
        System.out.println("removing " +node.getName());
        nodes.remove(node.getName());
        Variable parentvar = node.getVariable();
        // also update children to not point to parent
        List<BNode> additions = new ArrayList<BNode>();
        List<BNode> deletions = new ArrayList<BNode>();
        for (BNode child : nodes.values()) {
            List<EnumVariable> parents = child.getParents();
            if (parents != null) {
                if (parents.contains(parentvar)) {
                    deletions.add(child);
                    Set<Variable> newparents = new HashSet<Variable>();
                    for (EnumVariable v:parents)
                        if (!v.getName().equals(parentvar.getName()))
                            newparents.add((Variable)v);
                    BNode newchild = Predef.getBNode(child.getVariable(), new ArrayList<Variable>(newparents), Predef.getType(child));
                    additions.add(newchild);
                }
            }
        }
        for (BNode deleteMe : deletions)
            nodes.remove(deleteMe.getName());
        for (BNode addMe : additions)
            nodes.put(addMe.getName(), addMe);
    }
    
    public void removeParent(BNode child, Variable parent) {
        System.out.println("child is: " + child + " parent is: " + parent);
        List<EnumVariable> parents = child.getParents();
        List<BNode> additions = new ArrayList<BNode>();
        List<BNode> deletions = new ArrayList<BNode>();
        if (parents != null) {
            if (parents.contains(parent)) {
                nodes.remove(child.getName());
                List<Variable> newparents = new ArrayList<Variable>();
                for (EnumVariable v:parents)
                    if (!v.getName().equals(parent.getName()))
                        newparents.add((Variable)v);
                BNode newchild = Predef.getBNode(child.getVariable(), newparents, Predef.getType(child));
                nodes.put(newchild.getName(), newchild);
            }
        }
    }
    
    public void addParent(BNode child, Variable parent) {
        List<EnumVariable> parents = child.getParents();
        Set<Variable> newparents = new HashSet<Variable>();
        if (parents != null)
            newparents.addAll(parents);
        if (newparents.contains(parent)) // already has parent, so nothing needed
            return;
        newparents.add((Variable)parent);
        BNode newchild = Predef.getBNode(child.getVariable(), new ArrayList<Variable>(newparents), Predef.getType(child));
        nodes.remove(child.getName());
        nodes.put(newchild.getName(), newchild);
    }
    
    public void load(String filename, boolean replace) {
        if (replace)
            nodes = new HashMap<String, BNode>();
        BNet bnet = BNBuf.load(filename);
        for (BNode node : bnet.getNodes()) {
            vars.put(node.getVariable().getName(), node.getVariable());
            nodes.put(node.getName(), node);
        }
    }
    
    public void load(String filename) {
        load(filename, true);
    }
    
    public void save(String filename) {
        BNBuf.save(getBNet(), filename);
    }
}
