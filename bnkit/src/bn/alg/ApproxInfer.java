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
import java.util.Random;
import java.util.Set;

public class ApproxInfer implements Inference {
	
    public BNet bn;
    private double logLikelihood = 1;
    private Random randomGenerator = new Random();
    private int iterations = 50;

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
    		//getParams() did not return suitable results for boolean nodes
    		//Set manually
    		if (pre.equals("Boolean")){
    			Object[] params = {true, false};
    			//Generate a pseudo random number to select parameter
    			int len = params.length;
    			int index = randomGenerator.nextInt(len);
    			//Randomly set the instance for the node
    			cbn.getNode(var).setInstance(params[index]);
    		} else {
        		String[] params = var.getParams().split(";");
        		//Generate a pseudo random number to select parameter
        		int index = randomGenerator.nextInt(params.length);
        		//Randomly set the instance for the node
        		cbn.getNode(var).setInstance(params[index]); 
    		}
    	}
    	
    	//Count Table for tracking state of query
    	List<EnumVariable> list = new ArrayList<EnumVariable>(q.X.size());
//    	Object[] test = q.X.toArray();
//    	System.out.println("Query");
    	for (Variable var : q.X) {
    		EnumVariable nVar = (EnumVariable)var;
//    		System.out.println(nVar.toString());
            list.add(nVar);
        }
    	CountTable storeTable = new CountTable(list);  	
//    	System.out.println("StoreTable");
//    	storeTable.display();    	
    	
    	//Iterations of sampling needs to be flexible
    	//Have to set it in EM? Earlier?
    	int N = iterations;
    	
    	//The main loop for the algorithm
    	for (int j = 0; j < N; j ++) {
    		for (Variable var : q.Z) { 
    			//Consider what order this is happening?
    			//Top down? Should it be randomised?
    			//Z is an ordered list
    			BNet mbVar = cbn.getMB(var);
    			Object result = mbVar.getMBProb(mbVar, mbVar.getNode(var));
    			if (result != null) {
    				cbn.getNode(var).setInstance(result);
    			}
    			
    			// Update Query count here? Outside loop is better?
    		}
    		//Update Query count
    		//Confused about what I'm doing below, think I've gotten a bit mixed up about what I'm looking for
    		for (Variable qVar : q.X) {
    			List<EnumVariable> parList = storeTable.table.getParents();
    			List<Object> instances = new ArrayList<Object>();
    			for (EnumVariable par : parList) {
    				instances.add(cbn.getNode(par).getInstance());
    			}
    			//How do you know parent query is in right order?
    			//Query is a parent
    			if (instances.contains(null)) {
    				System.out.println("Stop");
    			}
    			storeTable.count(instances.toArray());
//    		System.out.println("StoreTable complete");
//    	    storeTable.display();    			
    		}
    	}
//    	System.out.println("StoreTable complete");
////    	storeTable.display();
//    	Collection<Double> values = storeTable.table.getValues();
//    	Object[] points = list.toArray();
    	
//    	System.out.println("To JPT");
    	answer = new JPT(storeTable.table);
    	
    	logLikelihood += 1;
    	
    	//Reset all unevidenced nodes in network
    	for (BNode node : cbn.getNodes()) {
    		if (q.Z.contains(node.getVariable())){
    			node.resetInstance();
    		}
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
    
    public int getIterations() {
    	return iterations;
    }
    
    public void setIterations(int iter) {
    	iterations = iter;
    }
}

