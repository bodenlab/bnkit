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

import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;
import dat.Enumerable;
import java.io.Serializable;
import java.util.*;

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
    private final Map<String, BNode> nodesByName = new HashMap<>();
    private final Map<Variable, BNode> nodesByVar = new HashMap<>();

    /**
     * Linking children-to-parents, happens when BN is compiled (@see
     * BNet.compile())
     */
    private Map<BNode, Set<BNode>> ch2par = new HashMap<>();

    /**
     * Linking parent-to-children, happens when BN is compiled (@see
     * BNet.compile())
     */
    private Map<BNode, Set<BNode>> par2ch = new HashMap<>();

    /**
     * Listing nodes in a top-down ordered fashion, for quick query
     * processing.
     */
    private List<BNode> ordered = new ArrayList<>();

    /**
     * Listing node variables in a top-down ordered fashion
     */
    private List<Variable> orderedVars = new ArrayList<>();

    private HashMap<String, Set<BNode>> tagged = new HashMap<>();

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
        if (nodesByName.containsKey(node.getName()) || nodesByVar.containsKey(node.getVariable())) {
            throw new BNetRuntimeException("Duplicate node names in BNet: " + node.getName());
        }
        compiled = false;
        nodesByName.put(node.getName(), node);
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
     * Removes a node from the BN. 
     * Note that this does not de-associate the references to other nodes by variables.
     * @param node 
     */
    public void remove(BNode node) {
        compiled = false;
        nodesByName.remove(node.getName());
        nodesByVar.remove(node.getVariable());
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
            ch2par.clear();
            par2ch.clear();
            for (BNode node : nodesByVar.values()) {
                Set<BNode> parents = new HashSet<>();
                List<EnumVariable> parvars = node.getParents();
                if (parvars != null) {
                    for (EnumVariable parent : parvars) {
                        BNode pnode = nodesByVar.get(parent);
                        if (pnode == null) {
                            System.err.println("Invalid Bayesian network: node " + parent.getName() + " is not a member but referenced by " + node.getName());
                            throw new BNetRuntimeException("Invalid Bayesian network: node " + parent.getName() + " is not a member but referenced by " + node.getName());
                        }
                        parents.add(pnode);
                        Set<BNode> children = par2ch.get(pnode);
                        if (children == null) {
                            children = new HashSet<>();
                            par2ch.put(pnode, children);
                        } 
                        children.add(node);
                    }
                }
                ch2par.put(node, parents);
            }
            compiled = true; // note: must be declared compiled before full compilation ("ordered" below uses descendant links)

            ordered.clear();
            for (BNode root : getRoots()) {
                List<BNode> ordered_from_root = new ArrayList<>();
                if (ordered.contains(root) && ordered.size() == this.nodesByVar.size()) {
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
     * Force set the compiled status of the network
     * @param bool
     */
    public void setCompiled(Boolean bool){
        this.compiled = bool;
    }

    /**
     * Retrieve the root nodes of the BN, i.e. the nodes with no parents. Useful
     * to know if traversing the network structure.
     *
     * @return the set of all root nodes
     */
    public Set<BNode> getRoots() {
        Set<BNode> roots = new HashSet<>();
        for (BNode node : nodesByName.values()) {
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
        return nodesByName.values();
    }

    /**
     * Retrieve all parents of a specified node, ie the nodes that are
     * conditioning the specified node.
     *
     * @param nodeName the name of the node for which parents are sought
     * @return the parent nodes' names
     */
    public Set<String> getParentsNames(String nodeName) {
        BNode node = this.getNode(nodeName);
        if (node == null) {
            throw new BNetRuntimeException("Node " + nodeName + " does not exist in network");
        }
        return BNet.this.getParentsNames(node);
    }

    /**
     * Retrieve all parent names of a specified node, ie the nodes that are
     * conditioning the specified node.
     *
     * @param node the node (by reference)
     * @return the parents (by name), null if no parents
     */
    public Set<String> getParentsNames(BNode node) {
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
     * Retrieve all parent nodes of a specified node, ie the nodes that are
     * conditioning the specified node.
     *
     * @param node the node (by reference)
     * @return the parent nodes, null if no parents
     */
    public Set<BNode> getParents(BNode node) {
        if (!compiled) {
            this.compile(); // relies on ch2par: children-to-parent linking
        }
        Set<BNode> parents = ch2par.get(node);
        return parents;
    }

    /**
     * Retrieve all "sibling" nodes of a specified node, ie the nodes that are
     * conditioned by the same parent/s.
     *
     * @param node the node (by reference)
     * @return Set of sibling nodes, or null if no parents. If no siblings exist, the set is empty; excludes the node that is queried
     */
    public Set<BNode> getSiblings(BNode node) {
        if (!compiled) {
            this.compile(); // relies on ch2par: children-to-parent linking
        }
        Set<BNode> parents = ch2par.get(node);
        if (parents == null)
            return null;
        Set<BNode> children = new HashSet<>();
        for (BNode p : parents) {
            Set<BNode> ch = par2ch.get(p);
            children.addAll(ch);
        }
        children.remove(node);
        return children;
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
    public Set<String> getChildrenNames(BNode node) {
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
     * Retrieve all "children" nodes of the specified node. That is all nodes
     * that are "below" the specified node. Requires the network to be compiled.
     *
     * @param node the node for which all children are sought
     * @return the children nodes
     */
    public Set<BNode> getChildren(BNode node) {
        if (!compiled) {
            this.compile(); // relies on par2ch: parent to children linking
        }
        Set<BNode> children = par2ch.get(node);
        if (children == null) {
            return null;
        }
        Set<BNode> nodes = new HashSet<>();
        for (BNode n : children) {
            nodes.add(n);
        }
        return nodes;
    }

    /**
     * Check if there are "children" nodes of the specified node. That is nodes
     * that are "below" the specified node. Requires the network to be compiled.
     *
     * @param node the node for which children are sought
     * @return true if there is at least one child, false otherwise
     */
    public boolean hasChildren(BNode node) {
        if (!compiled) {
            this.compile(); // relies on par2ch: parent to children linking
        }
        Set<BNode> children = par2ch.get(node);
        return (children != null);
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
            return Collections.EMPTY_LIST;
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
     * Retrieve all node variables in order (from root/s to leaves; with parallel paths
     * in arbitrary order).
     * @return variables in a defined topological order
     */
    public List<Variable> getOrderedVariables() {
        List<BNode> nodes = this.getOrdered();
        for (BNode node : nodes){
            this.orderedVars.add(node.getVariable());
        }
        return this.orderedVars;
    }

    /**
     * Get an alphabetical ordering of nodes
     * @return list of alphabetically ordered nodes
     */
    public List<BNode> getAlphabetical() {
        List<BNode> nodes = new ArrayList<>();
        SortedSet<String> names = new TreeSet<>(this.getNames());
        for (String n : names) {
            nodes.add(this.getNode(n));
        }
        return nodes;
    }

    /**
     * Get the names of all variables in the BN.
     *
     * @return the names of the nodes in the BN
     */
    public Set<String> getNames(){
        return nodesByName.keySet();
    }

    /**
     * The node structure via its name.
     *
     * @param name the name of the node
     * @return the class instance of the node
     */
    public BNode getNode(String name) {
        BNode node = nodesByName.get(name);
        if (node != null) {
            return node;
        } else {
            for (String nodename : nodesByName.keySet()) {
                int s = nodename.lastIndexOf(".");
                /*
                FIXME? We get an index error if the nodename has no added index number. The below if check
                is a fix but does not address the root cause
                 */
                String nname;
                if (s < 0)
                    nname = nodename;
                else
                    nname = nodename.substring(0, s);
                if (nname.equals(name))
                    return nodesByName.get(nodename);
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
        if (node.isRoot())
            return null;
        List<EnumVariable> parents = node.getParents();
        Object[] bkey = new Object[parents.size()];
        for (int i = 0; i < bkey.length; i++) {
            BNode parent = this.getNode(parents.get(i)); // changed this from a string-search (based on .toString)
            bkey[i] = parent.getInstance();
        }
        return bkey;
    }

    /**
     * Use currently set instances to sample values for unset nodes.
     * The whole BN will be instantiated after this process.
     * Nodes are processed in a top-down manner; all parents before a child.
     * Sampling is done by using current parameters.
     */
    public void sampleInstance() {
        for (BNode node : this.getOrdered()) { // traverse nodes in top-down order
            if (node.getInstance() == null) { // not currently set, so we'll do that
                Object[] key = getEvidenceKey(node);
                Distrib d = node.getDistrib(key);
                Object sample = d.sample();
                node.setInstance(sample);
            } // else, ignore
        }
    }

    /**
     * Determine which subset of nodes that are relevant to a specific query and
     * evidence combination. "Every variable that is not an ancestor of a query
     * variable or evidence variable is irrelevant to the query" (Russell and
     * Norvig, 2002, p. 509).
     *
     * @param query the variables that are in the query
     * @return a new BN with the relevant CPTs only
     * @deprecated
     */
    public BNet getRelevant(Variable... query) {
        BNet nbn = new BNet();
        Set<String> qset = new HashSet<>();
        for (Variable query1 : query) {
            qset.add(query1.toString());
        }
        // first add all instantiated/evidence variables and queried variables
        for (BNode node : nodesByName.values()) {
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
     * Algorithm for finding nodes reachable from X(query) given Z(evidence) via active trails
     * Based on Algorithm 3.1 in Probabilistic Graphical Models - Principles and Techniques, Koller, D., Friedman, N., pg.75
     * @param query the variables that are in the query
     * @return a new BN with relevant CPTs only
     * TODO: Either rename to a "setter" since it actually modifies the BNet (see comment at end), or rewrite so that it only
     * returns the list of D-connected nodes.
     */
    public List<BNode> getDconnected(Variable... query) {
        //Get set of evidence
        Set<BNode> z = new HashSet<BNode>(); //Evidence
        for (BNode node : this.getNodes()) {
            if (node.getInstance() != null) {
                z.add(node);
            }
        }
     	//Record all ancestors of evidence
        //Phase I
        Set<BNode> a = new HashSet<BNode>(); //Ancestors of evidence
        for (BNode v : z) {

            List<BNode> anc = this.getAncestors(v);
            a.addAll(anc);
            a.add(v);
        }

        //PhaseII: traverse active trails starting from X (query)
        List<NodeDirection> l = new ArrayList<NodeDirection>(); //nodes to be visited
        //Have to add query as node to be visited
        //FIXME direction for the query?
        //I think a query node can always be added once with a single direction and that will always work
        for (Variable q : query) {
            BNode qNode = this.getNode(q);
            if (qNode == null) {
                throw new NullPointerException("Invalid query: node " + q.toString() + " does not exist in this network");
            }
            l.add(new NodeDirection(qNode, "up"));
        }

        Set<NodeDirection> v = new HashSet<NodeDirection>(); //node,direction marked as visited
        Set<BNode> r = new HashSet<BNode>(); //nodes reachable via active trail

        while (!l.isEmpty()) {
            NodeDirection cur = l.remove(0);
            if (!cur.within(v)) { //hasn't been visited
                r.add(cur.getNode()); //node is reachable
//     			if (!z.contains(cur.getNode())) { //isn't evidenced
//     				//FIXME when should an evidence node be included?
//     				//At final step check overlap between ancestors of added nodes and evidence?
//     				r.add(cur.getNode()); //node is reachable
//     			}
                v.add(cur); //mark node as visited
                if (cur.getDirection() == "up" && !z.contains(cur.getNode())) { //trail up through Y, active if Y not in Z
                    //get set of parent nodes
                    if (cur.getNode().getParents() != null) {
                        //find parents to be visited from bottom
                        for (Variable par : cur.getNode().getParents()) {
                            l.add(new NodeDirection(this.getNode(par), "up"));
                        }
                    }
                    if (this.getChildrenNames(cur.getNode()) != null) {
                        //find children to be visited from top
                        for (String c : this.getChildrenNames(cur.getNode())) {
                            l.add(new NodeDirection(this.getNode(c), "down"));
                        }
                    }
                } else if (cur.getDirection() == "down") { //trails down through Y
                    if (!z.contains(cur.getNode())) { //downward trails to Y's children are active
                        if (this.getChildrenNames(cur.getNode()) != null) {
                            //find children to be visited from top
                            for (String c : this.getChildrenNames(cur.getNode())) {
                                l.add(new NodeDirection(this.getNode(c), "down"));
                            }
                        }
                    }
                    if (a.contains(cur.getNode())) { //v-structure (converging) trails are active
                        if (cur.getNode().getParents() != null) {
                            //find parents to be visited from bottom
                            for (Variable par : cur.getNode().getParents()) {
                                l.add(new NodeDirection(this.getNode(par), "up"));
                            }
                        }
                    }
                }
            }
        }
     	//Add all reachable nodes IN ORDER and set all non-relevant nodes to false
        //FIXME If a parent of a relevant node is evidence should it be included?
        List<BNode> output = new ArrayList<BNode>();
        for (BNode n : this.getOrdered()) {
            if (r.contains(n)) {
                //n.setRelevant(true); // Modifies each globally accessible BNode so may cause issues when inference is multi-threaded
                output.add(n);
            } else {
                //n.setRelevant(false); // Modifies each globally accessible BNode so may cause issues when inference is multi-threaded
            }
        }
        return output;
    }

    public void resetNodes(){
        for (BNode node : this.getNodes()){
            node.resetInstance();
        }
    }
    
   
    /**
     * Get the Markov Blanket (mb) of a query node.
     * 
     * Only need the factor tables of node and children
     * 
     * 
     * @param qnode query node
     * @return a list of relevant factor tables
     */
    
    //Should this method return a list of nodes in the mb instead?
    //Use of factor tables to get product could be influenced by cutting off parts of network?
    public Set<BNode> getMB(BNode qnode) {
        if (!compiled) // relies on parent to children to parent linking, so BN must be compiled
            this.compile(); 
        Set<BNode> mbNodes = new HashSet<>();
    	// add children of query node unless it's a leaf node
        Set<BNode> children = par2ch.get(qnode);
        if (children != null) {
            mbNodes.addAll(children);
            // if there are children, then add parents of children
            for (BNode cnode : children) {
                Set<BNode> pnodes = ch2par.get(cnode);
                if (pnodes != null)
                    mbNodes.addAll(pnodes);
            }
            // remove query node that accidentally got added above
            mbNodes.remove(qnode); 
        }
        // finally add parents of query node
        Set<BNode> pnodes = ch2par.get(qnode); 
        if (pnodes != null)
            mbNodes.addAll(pnodes);
        return mbNodes;
    }
    
    /**
     * Given the Markov Blanket network of a query, 
     * return a sample from the distribution of P(X|mb(X))
     * Assumes that other nodes in the Markov blanket have been instantiated.
     * 
     * @param query 
     * @return sample from distribution
     */
    public Object getMBProb(BNode query) {  	
        if (!compiled) // relies on parent to children to parent linking, so BN must be compiled
            this.compile(); 
    	// Store the instance incase the map is empty and you have to reset the node
        Object qInstance = query.getInstance();
    	query.resetInstance();
        Object query_sample;
        try {
            // Enumerable query
            EnumVariable evar = (EnumVariable) query.getVariable();
            Enumerable edom = evar.getDomain();
            double[] dist = new double[edom.size()];
            Object[] values = edom.getValues();
            // FIXME: Some unnecessary re-calculation below
            for (int i = 0; i < values.length; i ++) {
                Object[] key = this.getEvidenceKey(query);
                dist[i] = query.get(key, values[i]);
                query.setInstance(values[i]);
        	// Query node not included in set, used for initial factor table in product
                Set<BNode> children = par2ch.get(query);
                if (children != null) {
                    for (BNode node : children) {
                        key = this.getEvidenceKey(node);
                        Double f = node.get(key, node.getInstance());
                        dist[i] *= f;
                    }
                }
            }
            EnumDistrib edist = new EnumDistrib(edom, dist);
            query_sample = edist.sample();
        } catch (ClassCastException e) {
            // Non-enumerable query, which means that there are no children of it
            Object[] key = this.getEvidenceKey(query);
            Distrib gdist = query.getDistrib(key);
            query_sample = gdist.sample();
        }
        query.setInstance(qInstance);
        return query_sample;
    }
    
    /**
     * Given the Markov Blanket network of a query, 
     * return a sample from the distribution of P(X|mb(X))
     * 
     * @param mbNodes
     * @param query 
     * @return sample from distribution
     */
    public Object getMBProb(Set<BNode> mbNodes, BNode query) {  	
    	// Store the instance incase the map is empty and you have to reset the node
    	Object qInstance = query.getInstance();
    	query.resetInstance();
        try {
            // Enumerable query
            EnumVariable evar = (EnumVariable) query.getVariable();
            Enumerable edom = evar.getDomain();
            double[] dist = new double[edom.size()];
            Object[] values = edom.getValues();
            // FIXME: Some unnecessary re-calculation below
            for (int i = 0; i < values.length; i ++) {
                Object[] key = this.getEvidenceKey(query);
                dist[i] = query.get(key, values[i]);
                query.setInstance(values[i]);
        	// Query node not included in set, used for initial factor table in product
        	for (BNode node : mbNodes) {
                    key = this.getEvidenceKey(node);
                    Double f = node.get(key, node.getInstance());
                    dist[i] *= f;
                }
            }
            EnumDistrib edist = new EnumDistrib(edom, dist);
            Object end = edist.sample();
            return end;   		    	
        } catch (ClassCastException e) {
            // Non-enumerable query, which means that there are no children of it
            Object[] key = this.getEvidenceKey(query);
            Distrib d = query.getDistrib(key);
            Object result = d.sample();
            return result;
        }
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
        Set<String> s = new HashSet<>();
        s.addAll(Arrays.asList(array));
        return s;
    }
    
    /**
     * Small class to store node and direction for dConnectedness
     * @author Alex
     *
     */
    public class NodeDirection {
    	private BNode node;
    	private String direction;
    	
    	public NodeDirection(BNode node, String direction) {
    		this.node = node;
    		this.direction = direction;
    	}
    	
    	public BNode getNode() {
    		return node;
    	}
    	
    	public String getDirection() {
    		return direction;
    	}
    	
    	public boolean within(Set<NodeDirection> input) {
    		for (NodeDirection n : input) {
    			if (this.node.equals(n.getNode()) && this.direction.equals(n.getDirection())){
    				return true;
    			}
    		}
    		return false;
    	}
    }
    
    public Set<String> getTagNames(){
        return this.tagged.keySet();
    }

    /**
     * Get all nodes that have a tag assigned
     * @return List of tagged nodes
     */
    public List<BNode> getTagged(){
        ArrayList taggedList = new ArrayList();
        Set<BNode> taggedSet = new HashSet();
        Iterator it = this.tagged.entrySet().iterator();
        while (it.hasNext()){
            Map.Entry<String, Set<BNode>> pairs = (Map.Entry)it.next();
            for (BNode node : pairs.getValue()){
                taggedSet.add(node);
            }
        }
        taggedList.addAll(taggedSet);
        return taggedList;
    }

    /**
     * Clear all current tags for all nodes
     */
    public void removeAllTags(){
        this.tagged.clear();
    }


    /**
     * Return nodes that are in all specified tag groups
     * @param tags tag names
     * @return list of nodes
     */
    public List<BNode> getTagged(String... tags){
        ArrayList taggedList = new ArrayList();
        Set<BNode> taggedSet = new HashSet();
//        List<BNode> taggedNodes = this.getTagged();
        for (BNode node : this.nodesByName.values()){
            boolean common = false;
            for (String tag: tags){
                if (!this.tagged.keySet().contains(tag)){
                    throw new IllegalArgumentException("Tag " + tag + " does not exist");
                }
                if (!this.tagged.get(tag).contains(node)){
                    common = false;
                    break;
                } else{
                    common = true;
                }
            }
            if (common){
                taggedSet.add(node);
            }
        }
        taggedList.addAll(taggedSet);
        return taggedList;
    }

    /**
     * Remove a node from a tag group
     * @param tag
     * @param node
     */
    public void removeTag(String tag, BNode node){
        this.tagged.get(tag).remove(node);
    }

    /**
     * Assign nodes to a specific tag group
     * @param tag
     * @param nodes
     */
    public void setTags(String tag, BNode... nodes){
        if (!this.tagged.containsKey(tag)){
            this.tagged.put(tag, new HashSet<BNode>());
        }
        for (BNode node : nodes){
            this.tagged.get(tag).add(node);
        }
    }

    public void setTags(String[] tags, BNode... nodes){
        for (String tag: tags){
            this.setTags(tag, nodes);
        }
    }

    /**
     * Get the tag groups in which a node belongs
     * @param node
     * @return List of tag groups
     */
    public List<String> getTags(BNode node){
        ArrayList<String> tags = new ArrayList();
        Iterator it = this.tagged.entrySet().iterator();
        while (it.hasNext()){
            Map.Entry<String, Set<BNode>> pairs = (Map.Entry)it.next();
            if (pairs.getValue().contains(node)){
                tags.add(pairs.getKey());
            }
        }
        return tags;
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
