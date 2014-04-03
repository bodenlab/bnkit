package bn.alg;

import bn.BNet;
import bn.BNode;
import bn.CPT;
import bn.CountTable;
import bn.Distrib;
import bn.EnumDistrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.Enumerable;
import bn.FactorTable;
import bn.JPT;
import bn.Variable;
import bn.alg.VarElim.VEQuery;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ApproxInfer implements Inference {
	
    public BNet bn;
    private double logLikelihood = 0;

    /** 
     * Approximate inference in Bayesian network by Gibbs algorithm (MCMC).
     * In accordance with the method described in Russell and Norvig, 
     * Artificial Intelligence: A Modern Approach, 3e, 2009.
     */
    @Override
    public void instantiate(BNet bn) {
        this.bn = bn;
        this.bn.compile();
    }
    
    /**
     * Construct the data structure for the query
     * X - query variable (non-evidence)
     * E - evidence variables in bn
     * Z - all other non-evidence variables
     */
    @SuppressWarnings("rawtypes")
    @Override
    public Query makeQuery(Variable... qvars) {
        List<Variable> X = new ArrayList<Variable>();
        List<Variable> E = new ArrayList<Variable>();
        List<Variable> Z = new ArrayList<Variable>();
        List<Variable> ordered = new ArrayList<Variable>();
        try {
            for (int i = 0; i < qvars.length; i++) {
                X.add((EnumVariable) qvars[i]);
            }
            for (BNode node : bn.getOrdered()) {
                Variable var = node.getVariable();
                ordered.add(var);
                if (node.getInstance() != null) {
                    E.add(var);
                //need to keep track of ALL non-evidence variables
                //this includes the query variable
                } else {
                    Z.add((Variable) var);
                }
            }
        } catch (RuntimeException e) {
            throw new RuntimeException("makeQuery, ApproxInfer didn't work");
        }
        
        return new AQuery(X, E, Z);
     
    }
    
    /**
     * Perform Approximate inference
     */
    @SuppressWarnings("rawtypes")
    @Override
    public JPT infer(Query query) {
    	JPT answer = null;
    	AQuery q = (AQuery) query;
    	//Take a copy of the current state of network
    	BNet cbn = bn;
    	//First set all non-evidenced nodes including query
    	//Should be done randomly**
    	for (Variable var : q.Z) {
    		String pre = var.getPredef();
    		if (pre.equals("Boolean")){
    			Object[] params = {true, false};
    			//Randomize the index
    			cbn.getNode(var).setInstance(params[0]);
    		} else {
    			String test = var.getParams();
        		String[] params = var.getParams().split(";");
        		//make index into params random?
        		cbn.getNode(var).setInstance(params[0]);
    		}
    	}
    	
    	//Initialise a structure for storing the results of the query
    	//Must have space for all possible parameters of query
    	
    	//For the query and associated hashMap
    	//Associated map stores parameters with final scores
    	Map<Variable, Map> store = new HashMap<Variable, Map>();
    	//Look at all variables in query
    	for (Variable param : q.X) {
    		//For the parameter and associated count
    		Map<Object, Integer> inner = new HashMap<Object, Integer>();
    		if (param.getPredef().equals("Boolean")){
    			Object[] parameters = {true, false};
        		for (Object parameter : parameters){
        			//Add each parameter with an initial count 0
        			inner.put(parameter, 0);    			
        		}
    		} else {
    			String[] parameters = param.getParams().split(";");
        		for (String parameter : parameters){
        			//Add each parameter with an initial count 0
        			inner.put(parameter, 0);    			
        		}
    		}

    		//Add the map with parameters and counts to the map with node as key
    		store.put(param, inner);
    	}

    	//Iterations of sampling needs to be flexible
    	//Have to set it in EM? Earlier?
    	int N = 100;
    	
    	//The main loop for the algorithm
    	for (int j = 0; j < N; j ++) {
    		for (Variable var : q.Z) { 
    			//Consider what order this is happening?
    			//Top down? Should it be randomised?
    			//Z is an ordered list
    			BNet mbVar = cbn.getMB(var);
    			Object result = mbVar.getMBProb(mbVar, mbVar.getNode(var));
    			cbn.getNode(var).setInstance(result);
    			// Update Query count here? Outside loop is better?
    		}
    		//Update Query count
    		for (Variable qVar : q.X) {
    			Object key = cbn.getNode(qVar).getInstance();
    			// try/catch?

    			if (key instanceof java.lang.Boolean){
    				Integer result = (Integer)store.get(qVar).get(key);
        			result++;
        			store.get(qVar).put(key, result);
    			} else {
    				Integer result = (Integer)store.get(qVar).get(key.toString());
        			result++;
        			store.get(qVar).put(key, result);
    			}

    		}   		
    	}
    	// Normalise the query counts
    	//Translate to distrib
    	//Method to sample from distrib
    	for (Variable qVar : q.X) {
    		Collection values = store.get(qVar).values();
    		Object[] vals = values.toArray();
    		double[] data = new double[values.size()];
        	for (int i = 0; i < values.size(); i++){
        		double res = (double)((int)vals[i]);
        		data[i] = res;
        	}
    		EnumDistrib dist = new EnumDistrib((Enumerable)qVar.getDomain(), data);
    		System.out.println(dist.toString());
    		System.out.println(qVar.toString());   		
    		
    		EnumTable update = cbn.getNode(qVar).getTable();
    		//Node currently has no associated CPT
    		//Does this always mean it's a root node?
    		//Should there be a check if it's a root?
    		if (update == null){
    			if (cbn.getNode(qVar).isRoot()) {
    				System.out.println("Create CPT for root");
    				CPT output = new CPT((EnumVariable)qVar);
        			output.put(dist);
    			} else {
    				//What if it isn't the root but the table is null?
    			}
    		} else {
    			//How to update a CPT with new distribution?
    			System.out.println("CPT already exists - not updated");
    			System.out.println(update.toString());
    			
    		}
//    		FactorTable fin = output.makeFactor();
//    		answer = new JPT(update);
    	}
    	
    	return answer;
    }
    
    public class AQuery implements Query {

        final List<Variable> X;
        final List<Variable> E;
        final List<Variable> Z;

        AQuery(List<Variable> X, List<Variable> E, List<Variable> Z) {
            this.X = X;
            this.E = E;
            this.Z = Z;
        }
    }

    
    @Override
    public double getLogLikelihood() {
        return logLikelihood;
    }    
    
    public double getMarkovBlanket(BNet mbNet, BNode query) {
    	query.resetInstance();
    	Set<FactorTable> fTables = new HashSet<FactorTable>();
    	//Don't include query node in Set
    	//Save it for initialising ft
    	for (BNode node : mbNet.getNodes()){
    		if (node.getName() != query.getName()) {
    			fTables.add(node.makeFactor(mbNet));
    		}
    	}
    	//Get the factor table for the query and make it the product start point
    	FactorTable ft = query.makeFactor(mbNet);
    	//For each factor table, add it to the product
    	for (FactorTable factor : fTables) {
    		ft = FactorTable.product(ft, factor);    		
    	}
    	ft.display();
    	return 0.0;
    }

}

///**
// * Get the Markov Blanket (mb) of a query node.
// * Returns the parents, children and parents of children 
// * for the query node
// * 
// * @param query variable to query
// * @return a new BN with relevant CPTs and GDTs
// */
//
//public BNet getMB(Variable query) {
//	BNet mbn = new BNet();
//	String qName = query.getName();
//	BNode qNode = getNode(query);
//	Set<String> parents = new HashSet<String>();
//	//get children of query node
//	Set<String> children = new HashSet<String>();
//	//add the parents of the query node unless it is a root node
//	if (getParents(qName) != null ) {
//		parents.addAll(getParents(qName));
//	}
//	//add children of query node unless it's a leaf node
//	Collection<String> test = getChildren(qNode);
//	if (getChildren(qNode) != null ) {
//		children.addAll(getChildren(qNode));
//	} else {
//		children.add(qNode.getName());
//	}
//	//for each child get its parents and add to the parents Set
//	for (String child : children) {
//		//child will always have at least one parent (query)
//		Set<String> newParents = getParents(child);
//		//query node will always be listed here as a parent
//		for (String parent : newParents) {
//			parents.add(parent);
//		}
//	}
//	for (BNode node : nodes.values()) {
//		//add all parents and children to the markov blanket network
//		//query node will be in parents set - unless no children!
//		if (parents.contains(node.getVariable().toString()) || children.contains(node.getName())){
//			mbn.add(node);
//		}
//	}
//	mbn.compile();
//	return mbn;
//}
//
//public Object getMBProb(BNet mbNet, BNode query) {  	
//	
//	query.resetInstance();
//	Set<FactorTable> fTables = new HashSet<FactorTable>();
//	//Don't include query node in Set
//	//Save it for initialising ft
//	for (BNode node : mbNet.getNodes()){
//		if (node.getName() != query.getName()) {
////			node.getTable().display();
//			String test = node.getVariable().getParams();
//			FactorTable fact = node.makeFactor(mbNet);
//			//only works when prior prob available
//			if (fact != null ){
//				fTables.add(node.makeFactor(mbNet));
////    			System.out.println(node.toString());
////    			node.makeFactor(mbNet).display();
//			}
//		}
//	}
//	//Get the factor table for the query and make it the product start point
//	FactorTable ft = query.makeFactor(mbNet);
//	//For each factor table, add it to the product
//	for (FactorTable factor : fTables) {
//		ft = FactorTable.product(ft, factor);    		
//	}
//	Collection values = ft.map.values();
//	Object[] vals = values.toArray();
//	double[] data = new double[values.size()];
//	for (int i = 0; i < values.size(); i++){
//		//WARNING - will this casting always work?
//		//I think yes because in a dif situation I was setting vals to Integer and that caused problems...
//		double res = (double)(vals[i]);
//		data[i] = res;
//	}
//	Distrib dist = new EnumDistrib((Enumerable)query.getVariable().getDomain(), data);
//	Object end = dist.sample(); 	
//	
////	System.out.println(query.toString());
////	System.out.println(end.toString());
//	return end;
//}


///////EXTRA
//b.setInstance(true);
//e.setInstance(false);
//s.setInstance(6.0);
//j.setInstance(true);
//m.setInstance(true);
//List<EnumVariable> parents = a.getParents();
//Set<String> parStr = new HashSet<String>();
//List<Object> parVar = new ArrayList<Object>();
//Collection<BNode> nodes = bn.getNodes();
//for (EnumVariable parent : parents) {
//	parStr.add(parent.toString());
//}
//for (BNode node : nodes) {
//	String test = node.getVariable().toString();
//	if (parStr.contains(node.getVariable().toString())){
//		parVar.add(node.getInstance());
//	}
//}
//EnumTable getTable = a.getTable();
////EnumDistrib distr = a.getDistrib();
//Double pos = a.get(parVar.toArray(), true);
//Double neg = a.get(parVar.toArray(), false);
//getTable.display();
//System.out.println("sldkf");
//
////////EXTRA

