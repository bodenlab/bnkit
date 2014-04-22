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
import bn.GaussianDistrib;
import bn.JPT;
import bn.Variable;
import bn.alg.AICon.ApproxInferRuntimeException;
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
     * 
     * Convergence of algorithm is incomplete
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
     * Z - all non-evidence variables
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
                X.add(qvars[i]);
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
     * Perform Approximate inference using Gibbs sampling algorithm
     */
    @SuppressWarnings("rawtypes")
    @Override
    public QueryResult infer(Query query) {
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
    			int index = randomGenerator.nextInt(params.length);
    			//Randomly set the instance for the node
    			cbn.getNode(var).setInstance(params[index]);
    		} else if (pre.equals("Real")){
    			// Ultimately want to sample from one of the possible distributions
    			// FIXME - better way to do this sampling? Something without odd casts
    			Object[] nodeDistribs = cbn.getNode(var).getTable().getValues().toArray();
    			int index = randomGenerator.nextInt(nodeDistribs.length);
    			//Get individual mu and sigma from random distribution
    			String[] values = nodeDistribs[index].toString().split(";");
    			//Create distribution to sample from
    			Distrib d = new GaussianDistrib(Double.parseDouble(values[0]), Double.parseDouble(values[1]));
    			cbn.getNode(var).setInstance(d.sample());
    		} else {
        		String parm = var.getParams();
        		String[] params;
        		if (parm != null) {
        			params = parm.split(";");
	        		//Generate a pseudo random number to select parameter
	        		int index = randomGenerator.nextInt(params.length);
	        		//Randomly set the instance for the node
	        		cbn.getNode(var).setInstance(params[index]); 
        		} else {
        			throw new ApproxInferRuntimeException("Node must contain parameters");
        		}
    		}
    	}
    	
    	//Create count Table for tracking state of query
    	List<EnumVariable> list = new ArrayList<EnumVariable>(q.X.size());
    	for (Variable var : q.X) {
    		EnumVariable nVar = (EnumVariable)var;
    		list.add(nVar);
        }
    	CountTable storeTable = new CountTable(list);  	   	
    	
    	//Iterations of sampling
    	int N = iterations;
    	
    	//The main loop for the algorithm
    	for (int j = 0; j < N; j ++) {
    		//Iterate over all non-evidenced nodes
    		for (Variable var : q.Z) { 
    			//Consider what order this is happening?
    			//Top down? Should it be randomised?
    			//Z is an ordered list
    			
    			//Get the children of node and node for factors
    			List<BNode> mbVar = cbn.getMB(var);
    			//Sample from the mb distribution
    			Object result = bn.getMBProb(mbVar, bn.getNode(var), bn);
    			if (result != null) {
    				cbn.getNode(var).setInstance(result);
    			}
//    			System.out.println("New Assignment");
//    			System.out.println(cbn.getNode(var).toString());
    			
    			// Update Query count here? Outside loop is better?
    		}
    		//Update Query count
    		//Unsure if the count table is being updated accurately 
    		//Unsure of how it works with multiple queries
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
    		}
    	}

    	//After sampling the distribution, convert countTable to JPT
    	answer = new JPT(storeTable.table);
    	
    	//Convergence of algorithm is incomplete
    	logLikelihood += 1;
    	
    	//Reset all unevidenced nodes in network
    	for (BNode node : cbn.getNodes()) {
    		if (q.Z.contains(node.getVariable())){
    			node.resetInstance();
    		}
    	}
    	
    	return new AResult(answer);
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

    public class AResult implements QueryResult {

        final private JPT jpt;
        public AResult(JPT jpt) {
            this.jpt = jpt;
        }
        
        @Override
        public JPT getJPT() {
            return this.jpt;
        }
        
        public Map<Variable, EnumTable<Distrib>> getNonEnum() {
            throw new RuntimeException("Not implemented");
        }
    }


    
    /**
     * Get the number of iterations sampling will complete
     * @return iterations
     */
    public int getIterations() {
    	return iterations;
    }
    
    /**
     * Set the number of iterations sampling will complete
     * @return iterations
     */
    public void setIterations(int iter) {
    	iterations = iter;
    }
    
    public class ApproxInferRuntimeException extends RuntimeException {
        private static final long serialVersionUID = 1L;

        public ApproxInferRuntimeException(String message) {
            super(message);
        }
    }
}

