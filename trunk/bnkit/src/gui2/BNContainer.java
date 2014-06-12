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
 * Container of things related to the BN that is under construction. The BNet
 * object is constructed only when needed and is not maintained in this class.
 *
 * @author mikael 
 *  modified by jun.
 * All instances of BNode have been superceded by NodeModel. NodeModel functions
 * synonymously with BNode but also supports the Observer pattern (as an Observable).
 */
public class BNContainer {

    private Map<String, Variable> vars = new HashMap<String, Variable>(); // all variables that are considered and can be used to create nodes
    private Map<String, NodeModel> nodems = new HashMap<String, NodeModel>(); // analog of nodes using nodemodel objects.

    /**
     * Removes all stored BNodes and NodeModels
     */
    public void clear() {
        nodems = new HashMap<String, NodeModel>();
    }
    
    /**
     * Instantiates and returns BNet.
     * @return BNet
     */
    public BNet getBNet() {
        BNet bnet = new BNet();
        for (Map.Entry m : nodems.entrySet()) {
            BNode bn = ((NodeModel) m.getValue()).getBNode();
            bnet.add(bn);
        }
        return bnet;
    }

    /**
     * Add variable to container, replace if the same variable already exists.
     *
     * @param var
     */
    public void addVariable(Variable var) {
        vars.put(var.getName(), var);
    }

    public Variable getVariable(String name) {
        return vars.get(name);
    }

    /**
     * Add NodeModel to BNContainer's list.
     * @param node 
     */
    public void addNode(NodeModel node) {
        if (!vars.containsKey(node.getVariable().getName())) {
            addVariable(node.getVariable());
        }
        nodems.put(node.getName(), node);
    }

    /**
     * Get NodeModel by name.
     * @param name
     * @return NodeModel 
     */
    public NodeModel getNodeModel(String name) {
        System.out.println("nodems" + nodems.keySet());
        return nodems.get(name);
    }

    /**
     * Return BNContainer's NodeModel array.
     * @return 
     */
    public Map<String, NodeModel> getNodeModelArr() {
        return nodems;
    }

    
    public void removeNode(NodeModel node) {
        // remove actual node (and variable)
        if (node == null) {
            return;
        }
        nodems.remove(node.getName());
        Variable parentvar = node.getVariable();
        // also update children to not point to parent
        List<NodeModel> additions = new ArrayList<NodeModel>();
        List<NodeModel> deletions = new ArrayList<NodeModel>();
        for (NodeModel child : nodems.values()) {
            List<EnumVariable> parents = child.getParents();
            if (parents != null) {
                if (parents.contains(parentvar)) {
                    deletions.add(child);
                    Set<Variable> newparents = new HashSet<Variable>();
                    for (EnumVariable v : parents) {
                        if (!v.getName().equals(parentvar.getName())) {
                            newparents.add((Variable) v);
                        }
                    }
                    NodeModel newchild = Predef.getNodeModel(child.getVariable(), new ArrayList<Variable>(newparents), Predef.getType(child));
                    additions.add(newchild);
                }
            }
        }
        for (NodeModel deleteMe : deletions) {
            nodems.remove(deleteMe.getName());
        }
        for (NodeModel addMe : additions) {
            nodems.put(addMe.getName(), addMe);
        }
    }

    public void removeParent(NodeModel child, Variable parent) {
        List<EnumVariable> parents = child.getParents();
        List<NodeModel> additions = new ArrayList<NodeModel>();
        List<NodeModel> deletions = new ArrayList<NodeModel>();
        if (parents != null) {
            if (parents.contains(parent)) {
                nodems.remove(child.getName());
                List<Variable> newparents = new ArrayList<Variable>();
                for (EnumVariable v : parents) {
                    if (!v.getName().equals(parent.getName())) {
                        newparents.add((Variable) v);
                    }
                }
                NodeModel newchild = Predef.getNodeModel(child.getVariable(), newparents, Predef.getType(child));
                nodems.put(newchild.getName(), newchild);
            }
        }
    }

    public void addParent(NodeModel child, Variable parent) {
        List<EnumVariable> parents = child.getParents();
        Set<Variable> newparents = new HashSet<Variable>();
        if (parents != null) {
            newparents.addAll(parents);
        }
        if (newparents.contains(parent)) // already has parent, so nothing needed
        {
            return;
        }
        newparents.add((Variable) parent);
        NodeModel newchild = Predef.getNodeModel(child.getVariable(), new ArrayList<Variable>(newparents), Predef.getType(child));
        nodems.remove(child.getName());
        nodems.put(newchild.getName(), newchild);
    }


    // Uses NodeModel implementation now
    public void loadnm(String filename, boolean replace) {
        if (replace) {
            nodems = new HashMap<String, NodeModel>();
        }
        BNet bnet = BNBuf.load(filename);
        for (BNode node : bnet.getNodes()) {
            System.out.println(""); // print out the variables.s.s.dflsf
            NodeModel nm = new NodeModel(node);
            vars.put(nm.getVariable().getName(), nm.getVariable());
            nodems.put(nm.getName(), nm);
            System.out.println("!!!NODEMODEL " + nm.getStateAsText());
            System.out.println("NODE !! " + node.getStateAsText());
        }
    }

    public void save(String filename) {
        BNBuf.save(getBNet(), filename);
    }
}
