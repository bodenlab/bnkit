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

import bn.alg.Query;
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
        Set<String> str = new HashSet<>();
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
        List<BNode> nonredundant = new ArrayList<>();
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
        Set<String> names = new HashSet<>();
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
        List<BNode> nonredundant = new ArrayList<>();
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
     * @return return nodes in a defined topological order
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
        Set<String> qset = new HashSet<>();
        for (Variable query1 : query) {
            qset.add(query1.toString());
        }
        // first add all instantiated/evidence variables and queried variables
        for (BNode node : nodes.values()) {
//            if (node.getInstance() != null || qset.contains(node.getName())) {
        	// query is list of variables and node name and variable name are different
        	if (node.getInstance() != null || qset.contains(node.getVariable().toString())){
                nbn.add(node);
            }
        }
        // next add the ancestors of those already added
        Set<BNode> aset = new HashSet<>();
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
     * 
     * Only need the factor tables of node and children
     * 
     * 
     * @param query variable to query
     * @return a list of relevant factor tables
     */
    
    //Should this method return a list of nodes in the mb instead?
    //Use of factor tables to get product could be influenced by cutting off parts of network?
    public List<BNode> getMB(Variable query) {
    	List<BNode> mbNodes = new ArrayList<>();
//    	BNet mbn = new BNet();
    	String qName = query.getName();
    	BNode qNode = getNode(query);
    	//get children of query node
    	Set<String> children = new HashSet<>();
    	//add children of query node unless it's a leaf node
    	if (getChildren(qNode) != null ) {
    		children.addAll(getChildren(qNode));
    	} else {
    		children.add(qNode.getName());
    	}
    	for (BNode node : nodes.values()) {
    		//markov blanket calculation requires factor tables of query node and children only
    		if (children.contains(node.getName()) || node.getVariable().equals(query)){
    			mbNodes.add(node);
    		}
    	}
//    	mbn.compile();
    	return mbNodes;
    }
    
    

    /**
     * Given the Markov Blanket network of a query, 
     * return a sample from the distribution of P(X|mb(X))
     * 
     * @param mbNodes
     * @param query 
     * @param cbn the current bayesian network
     * @return sample from distribution
     */
    public Object getMBProb(List<BNode> mbNodes, BNode query, BNet cbn) {  	
    	
    	//Store the instance incase the map is empty and you have to reset the node
    	Object qInstance = query.getInstance();
    	query.resetInstance();
    	System.out.println("NEW QUERY");
    	//Store all factor tables for query to iterate over to find product
    	Set<FactorTable> fTables = new HashSet<>();
    	//Check if root and GDT - special case
    	Boolean leafQuery = false;
    	//Query node not included in set, used for initial factor table in product
    	for (BNode node : mbNodes){
    		System.out.println(node.toString());
    		if (node.getName() != query.getName()){
    			FactorTable fact = node.makeFactor(cbn);
    			
    			//instantiated priors cannot be used to factorise
    			if (fact != null ){
    				fTables.add(fact);
    			}
    		}
    	}
    	
    	//Get the factor table for the query and make it the product start point
    	FactorTable ft = query.makeFactor(cbn);
    	//For each factor table, add it to the product
    	for (FactorTable factor : fTables) {
    		ft = FactorTable.product(ft, factor); 
    	} 	
    	 
    	//Distribution never altered by factor when cg is leaf node
    	if (ft.hasNonEnumVariables()) { 
    		Distrib d = ft.getDistrib(0, query.getVariable());
    		Object result = d.sample();
//    		System.out.println(result);
    		return result;
    	}
    	
    	//FIXME
    	//Where a network is properly parameterised, an atomic ft represents 
    	//a non-initialised real node as query - USING THIS ALGORITHM
    	//CAN IT REPRESENT ANYTHING ELSE??? 
    	//Is this accurate for training or only for querying the final network?
    	if (ft.isAtomic()) {
//    		System.out.println("Atomic Map getMBProbReal");
    		Distrib d = ft.getDistrib(0, query.getVariable());
    		Object result = d.sample();
//    		System.out.println(result);
    		return result;
    	}
    	
    	//How to get distribution from factor table?
    	//Two choices for method currently
    	
    	Collection<Double> values = ft.getValues();  
    	
//    	Double[] d = (Double[]) values.toArray(new Double[values.size()]);
//    	double[] data = new double[d.length];
//    	for (int i = 0; i < d.length; i++){
//    		double res = (double)(d[i]);
//    		data[i] = res;
//    	}
//		Distrib dist = new EnumDistrib((Enumerable)query.getVariable().getDomain(), data);
//    	Object end = dist.sample(); 

    	Map<Object, Double> nFt = new HashMap<Object, Double>();
    	for (Map.Entry<Integer, Double> entry : ft.getMapEntries()) {
    		nFt.put((Object)entry.getKey(), entry.getValue());
    	}
    	Distrib dist = new EnumDistrib(nFt, (Enumerable)query.getVariable().getDomain());
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
