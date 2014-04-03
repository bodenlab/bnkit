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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Defines a Bayesian Network. The class will manage efficient access to the
 * nodes, based on the structure.
 *
 * @author m.boden
 */
public class BNet implements Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * Indicate if the nodes are organised (prepared for inference) or not
     */
    private boolean compiled = true;

    /**
     * All nodes of the BN
     */
    private final Map<String, BNode> nodes = new HashMap<String, BNode>();
    private final Map<Variable, BNode> nodesByVar = new HashMap<Variable, BNode>();

    /**
     * Linking children-to-parents, happens when BN is compiled (@see
     * BNet.compile())
     */
    private Map<BNode, Set<BNode>> ch2par = new HashMap<BNode, Set<BNode>>();

    /**
     * Linking parent-to-children, happens when BN is compiled (@see
     * BNet.compile())
     */
    private Map<BNode, Set<BNode>> par2ch = new HashMap<BNode, Set<BNode>>();

    /**
     * Listing variables in a top-down ordered fashion, for quick query
     * processing.
     */
    private List<BNode> ordered = new ArrayList<BNode>();

    /**
     * Construct a BN.
     */
    public BNet() {
    }

    /**
     * Add a node to the BN. All nodes must implement the BNode interface.
     *
     * @param node the node (e.g. a CPT, or GDT)
     */
    public void add(BNode node) {
        if (nodes.containsKey(node.getName()) || nodesByVar.containsKey(node.getVariable())) {
            throw new BNetRuntimeException("Duplicate node names in BNet: " + node.getName());
        }
        compiled = false;
        nodes.put(node.getName(), node);
        nodesByVar.put(node.getVariable(), node);
    }

    /**
     * Add a list of nodes to the BN. All nodes must implement the BNode
     * interface.
     *
     * @param nodes the nodes (each a CPT, GDT or any other BNode instance)
     */
    public void add(BNode... nodes) {
        for (BNode node : nodes) {
            add(node);
        }
    }

    /**
     * Compile and finalise the information in the Bayesian network so that
     * inference and other computationally complex processes can be done. (For
     * example: linking ancestors and descendants.) Most (if not all) functions
     * will call this automatically so typically explicit calls are not
     * required.
     */
    public void compile() {
        if (!compiled) {
            for (BNode node : nodes.values()) {
                Set<BNode> parents = new HashSet<BNode>();
                List<EnumVariable> parvars = node.getParents();
                if (parvars != null) {
                    for (EnumVariable parent : parvars) {
                        String pname = parent.getName();
                        BNode pnode = nodes.get(pname);
                        if (pnode == null) {
                            throw new BNetRuntimeException("Invalid Bayesian network: node " + pname + " is not a member but referenced by " + node.getName());
                        }
                        parents.add(pnode);
                        Set<BNode> children = par2ch.get(pnode);
                        if (children == null) {
                            children = new HashSet<BNode>();
                            par2ch.put(pnode, children);
                        }
                        children.add(node);
                    }
                }
                ch2par.put(node, parents);
            }
            compiled = true; // note: must be declared compiled before full compilation ("ordered" below uses descendant links)

            for (BNode root : getRoots()) {
                List<BNode> ordered_from_root = new ArrayList<BNode>();
                if (ordered.contains(root)) {
                    continue;
                }
                ordered_from_root.add(root);
                for (BNode node : getDescendants(root)) {
                    if (ordered.contains(node)) {
                        ;
                    } else {
                        ordered_from_root.add(node);
                    }
                }
                ordered.addAll(0, ordered_from_root); // making sure that any new sub-trees are added "on-top"
            }

        }
    }

    /**
     * Retrieve the root nodes of the BN, i.e. the nodes with no parents. Useful
     * to know if traversing the network structure.
     *
     * @return the set of all root nodes
     */
    public Set<BNode> getRoots() {
        Set<BNode> roots = new HashSet<BNode>();
        for (BNode node : nodes.values()) {
            if (node.isRoot()) {
                roots.add(node);
            }
        }
        return roots;
    }

    /**
     * Retrieve all the nodes in the network.
     *
     * @return the nodes
     */
    public Collection<BNode> getNodes() {
        return nodes.values();
    }

    /**
     * Retrieve all parents of a specified node, ie the nodes that are
     * conditioning the specified node.
     *
     * @param nodeName the name of the node for which parents are sought
     * @return the parent nodes' names
     */
    public Set<String> getParents(String nodeName) {
        BNode node = this.getNode(nodeName);
        if (node == null) {
            throw new BNetRuntimeException("Node " + nodeName + " does not exist in network");
        }
        return getParents(node);
    }

    /**
     * Retrieve all parent names of a specified node, ie the nodes that are
     * conditioning the specified node.
     *
     * @param node the node (by reference)
     * @return the parents (by name), null if no parents
     */
    public Set<String> getParents(BNode node) {
        List<EnumVariable> parents = node.getParents();
        if (parents == null) {
            return null;
        }
        Set<String> str = new HashSet<String>();
        for (EnumVariable var : parents) {
            str.add(var.toString());
        }
        return str;
    }

    /**
     * Retrieve all ancestors of a specified node. That is, all nodes that are
     * "above" the specified node.
     *
     * @param nodeName the name of the node for which ancestors are sought
     * @return ancestor nodes (by name)
     */
    public List<BNode> getAncestors(String nodeName) {
        BNode node = this.getNode(nodeName);
        if (node == null) {
            throw new BNetRuntimeException("Node " + nodeName + " does not exist in network");
        }
        return getAncestors(node);
    }

    /**
     * Retrieve all ancestors of a specified node. That is, all nodes that are
     * "above" the specified node.
     *
     * @param node the node for which ancestors are sought (by reference)
     * @return ancestor nodes (by name), null if no ancestors
     */
    public List<BNode> getAncestors(BNode node) {
        if (!compiled) {
            this.compile(); // relies on ch2par: children-to-parent linking
        }
        Set<BNode> parents = ch2par.get(node);
        if (parents == null) {
            return null;
        }
        List<BNode> nonredundant = new ArrayList<BNode>();
        nonredundant.addAll(parents);
        for (BNode p : parents) {
            List<BNode> ancestors = getAncestors(p);
            if (ancestors != null) {
                for (BNode a : ancestors) {
                    if (!nonredundant.contains(a)) {
                        nonredundant.add(a);
                    }
                }
            }
        }
        return nonredundant;
    }

    /**
     * Retrieve all "children" nodes of the specified node. That is all nodes
     * that are "below" the specified node. Requires the network to be compiled.
     *
     * @param node the node for which all children are sought
     * @return the children nodes
     */
    public Set<String> getChildren(BNode node) {
        if (!compiled) {
            this.compile(); // relies on par2ch: parent to children linking
        }
        Set<BNode> children = par2ch.get(node);
        if (children == null) {
            return null;
        }
        Set<String> names = new HashSet<String>();
        for (BNode n : children) {
            names.add(n.getName());
        }
        return names;
    }

    /**
     * Retrieve all the descendants of a node (recursively)
     *
     * @param nodeName name of node
     * @return list of nodes
     */
    public List<BNode> getDescendants(String nodeName) {
        BNode node = this.getNode(nodeName);
        if (node == null) {
            throw new BNetRuntimeException("Node " + nodeName + " does not exist in network");
        }
        return getDescendants(node);
    }

    /**
     * Retrieve all the descendants of a node (recursively)
     *
     * @param node root node for search
     * @return list of nodes
     */
    public List<BNode> getDescendants(BNode node) {
        if (!compiled) {
            this.compile(); // relies on par2ch: parent to children linking
        }
        Set<BNode> children = par2ch.get(node);
        if (children == null) {
            return null;
        }
        List<BNode> nonredundant = new ArrayList<BNode>();
        nonredundant.addAll(children);
        for (BNode c : children) {
            List<BNode> descendants = getDescendants(c);
            if (descendants != null) {
                for (BNode a : descendants) {
                    if (!nonredundant.contains(a)) {
                        nonredundant.add(a);
                    }
                }
            }
        }
        return nonredundant;
    }

    /**
     * Retrieve all nodes in order (from root/s to leaves; with parallel paths
     * in arbitrary order).
     */
    public List<BNode> getOrdered() {
        if (!compiled) {
            this.compile();
        }
        return ordered;
    }

    /**
     * Get the names of all variables in the BN.
     *
     * @return the names of the nodes in the BN
     */
    public Set<String> getNames() {
        return nodes.keySet();
    }

    /**
     * The node structure via its name.
     *
     * @param name the name of the node
     * @return the class instance of the node
     */
    public BNode getNode(String name) {
        BNode node = nodes.get(name);
        if (node != null) {
            return node;
        } else {
            for (String nodename : nodes.keySet()) {
                String[] parts = nodename.split("\\.");
                if (parts.length >= 2) {
                    if (parts[0].equals(name)) {
                        return nodes.get(nodename);
                    }
                }
            }
            return null;
        }
    }

    /**
     * The the node structure via its variable
     *
     * @param var the variable
     * @return the class instance of the node
     */
    public BNode getNode(Variable var) {
        BNode node = nodesByVar.get(var);
        if (node != null) {
            return node;
        } else {
            return null;
        }
    }

    /**
     * Get the boolean key to retrieve entries in the conditional probability
     * GIVEN evidence ie parents' instantiations.
     *
     * @param node
     * @return the (partial) key that applies to the current instantiation of
     * nodes in the BN
     */
    public Object[] getEvidenceKey(BNode node) {
        List<EnumVariable> parents = node.getParents();
        Object[] bkey = new Object[parents.size()];
        for (int i = 0; i < bkey.length; i++) {
            BNode parent = this.getNode(parents.get(i).toString());
            bkey[i] = parent.getInstance();
        }
        return bkey;
    }

    /**
     * Determine which subset of CPTs that are relevant to a specific query and
     * evidence combination. "Every variable that is not an ancestor of a query
     * variable or evidence variable is irrelevant to the query" (Russell and
     * Norvig, 2002, p. 509).
     *
     * @param query the variables that are in the query
     * @return a new BN with the relevant CPTs only
     */
    public BNet getRelevant(Variable[] query) {
        BNet nbn = new BNet();
        Set<String> qset = new HashSet<String>();
        for (int i = 0; i < query.length; i++) {
            qset.add(query[i].toString());
        }
        // first add all instantiated/evidence variables and queried variables
        for (BNode node : nodes.values()) {
            if (node.getInstance() != null || qset.contains(node.getName())) {
                nbn.add(node);
            }
        }
        // next add the ancestors of those already added
        Set<BNode> aset = new HashSet<BNode>();
        Set<String> direct = nbn.getNames();
        for (String name : direct) {
            List<BNode> ancestors = getAncestors(name);
            if (ancestors != null) {
                aset.addAll(ancestors);
            }
        }
        for (BNode anode : aset) { // add ancestors to new BN
            if (!nbn.getNames().contains(anode.getName())) {
                nbn.add(anode);
            }
        }
        nbn.compile();
        return nbn;
    }
    
    /**
     * Get the Markov Blanket (mb) of a query node.
     * Returns the parents, children and parents of children 
     * for the query node
     * 
     * @param query variable to query
     * @return a new BN with relevant CPTs and GDTs
     */
    
    public BNet getMB(Variable query) {
    	BNet mbn = new BNet();
    	String qName = query.getName();
    	BNode qNode = getNode(query);
    	Set<String> parents = new HashSet<String>();
    	//get children of query node
    	Set<String> children = new HashSet<String>();
    	//add the parents of the query node unless it is a root node
    	if (getParents(qName) != null ) {
    		parents.addAll(getParents(qName));
    	}
    	//add children of query node unless it's a leaf node
    	Collection<String> test = getChildren(qNode);
    	if (getChildren(qNode) != null ) {
    		children.addAll(getChildren(qNode));
    	} else {
    		children.add(qNode.getName());
    	}
    	//for each child get its parents and add to the parents Set
    	for (String child : children) {
    		//child will always have at least one parent (query)
    		Set<String> newParents = getParents(child);
    		//query node will always be listed here as a parent
    		for (String parent : newParents) {
    			parents.add(parent);
    		}
    	}
    	for (BNode node : nodes.values()) {
    		//add all parents and children to the markov blanket network
    		//query node will be in parents set - unless no children!
    		if (parents.contains(node.getVariable().toString()) || children.contains(node.getName())){
    			mbn.add(node);
    		}
    	}
    	mbn.compile();
    	return mbn;
    }
    
    
    /**
     * Given the Markov Blanket network of a query, 
     * return a sample from the distribution of P(X|mb(X))
     * 
     * @param mbNet mb network of query
     * @param query 
     * @return sample from distribution
     */
    public Object getMBProb(BNet mbNet, BNode query) {  	
    	
    	query.resetInstance();
    	Set<FactorTable> fTables = new HashSet<FactorTable>();
    	//Query node not included in set, used for initial factor table
    	for (BNode node : mbNet.getNodes()){
    		if (node.getName() != query.getName()) {
    			String test = node.getVariable().getParams();
    			FactorTable fact = node.makeFactor(mbNet);
    			//only works when prior prob available?
    			//Nodes without CPT do not weigh in on result
    			if (fact != null ){
    				fTables.add(node.makeFactor(mbNet));
    			}
    		}
    	}
    	//Get the factor table for the query and make it the product start point
    	FactorTable ft = query.makeFactor(mbNet);
    	//For each factor table, add it to the product
    	for (FactorTable factor : fTables) {
    		ft = FactorTable.product(ft, factor);    		
    	}
    	Collection values = ft.map.values();
    	Object[] vals = values.toArray();
		double[] data = new double[values.size()];
    	for (int i = 0; i < values.size(); i++){
    		double res = (double)(vals[i]);
    		data[i] = res;
    	}
		Distrib dist = new EnumDistrib((Enumerable)query.getVariable().getDomain(), data);
		Object end = dist.sample(); 	
		
    	return end;
    }

    /**
     * Utility function to create an array from a set of String.
     *
     * @param set set
     * @return an array
     */
    private String[] makeArray(Set<String> set) {
        String[] s = new String[set.size()];
        Iterator i = set.iterator();
        int c = 0;
        while (i.hasNext()) {
            s[c] = (String) i.next();
            c++;
        }
        return s;
    }

    /**
     * Utility function to create a set from an array of String,
     *
     * @param array an array of String
     * @return set with String
     */
    private Set<String> makeSet(String[] array) {
        Set<String> s = new HashSet<String>();
        for (String q : array) {
            s.add(q);
        }
        return s;
    }

}

/**
 * Exceptions for the BNet class
 *
 * @author mikael
 */
class BNetRuntimeException extends RuntimeException {

    private static final long serialVersionUID = -6465152863174383970L;
    String message;

    public BNetRuntimeException(String string) {
        message = string;
    }
}
