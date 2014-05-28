package bn;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import bn.SampleTable.Sample;
import bn.alg.ApproxInferCont;

/**
 * A structure used by Approximate inference to query a hybrid network.
 * Allows for counts and samples to be stored together and then processed 
 * into appropriate structures based on the query.
 * 
 * Records each 'state' of the Markov chain as it proceeds and in this way
 * records the influence of non-evidenced variables.
 * 
 * Incomplete
 * 
 * @author Alex
 *
 */
public class DataSample {

	public Map<Variable, List<Object>> map;
	protected final int nParents;
	protected final List<Variable> parents;
	private List<EnumVariable> discrete; //ALL DISCRETE PARENTS
	public CountTable counts = null;
	public Map<Variable, Map<Integer, List<Object>>> nonEnumTables = null;
	private boolean allContinuous = false;
	private int[] step; //FOR DISCRETE PARENTS ONLY
	
	public DataSample(List<Variable> query) {
        this.parents = new ArrayList<Variable>(query.size());
        this.nParents = query.size();
        this.map = new HashMap<Variable, List<Object>>();
        List<EnumVariable> discrete = new ArrayList<>();
        for (Variable var : query) {
            this.parents.add(var);
            try{
            	EnumVariable e = (EnumVariable)var;
            	discrete.add(e);
            } catch (ClassCastException e) {
            	;
            }
            List<Object> nList = new ArrayList<Object>(ApproxInferCont.iterations);
        	map.put(var, nList);
        }
        this.discrete = discrete;
        this.step = new int[this.discrete.size()];
        int prod = 1;
        for (int i = 0; i < discrete.size(); i++) {
            int parent = discrete.size() - i - 1;
            this.step[parent] = prod;
            prod *= this.discrete.get(parent).size();
        }
    }
	
    public List<Variable> getParents() {
        return parents;
    }
    
    /**
     * To be called after all counting is completed
     * Use this to process the data and create appropriate structures
     */
    public void createData() {
    	if (discrete.size() == nParents) { //Only discrete nodes in query
	    	List<Variable> pars = this.getParents();
	    	EnumVariable[] enums = new EnumVariable[discrete.size()];
	    	for (int i = 0; i < nParents; i++) {
	    		Variable p = pars.get(i);
	    		enums[i] = (EnumVariable)pars.get(i);
	    	}
	    	CountTable counts = new CountTable(enums); //Count table for discrete nodes
	    	for(int i = 0; i < ApproxInferCont.iterations ; i++) {
				List<Object> curKey = new ArrayList<Object>(pars.size());
				for (Variable p : pars) {
					curKey.add(this.getData(p).get(i));
				}
				counts.count(curKey.toArray());
			}
	    	this.counts = counts; //Store the counts for this query
    	} else if (discrete.size() < nParents && discrete.size() != 0) { //Hybrid query - discrete and real
    		List<Variable> pars = this.getParents();
	    	List<Variable> nonEnums = new ArrayList<Variable>(nParents-discrete.size());
	    	for (Variable par: pars) {
	    		if(!discrete.contains(par)) {
	    			nonEnums.add(par);
	    		}
	    	}
	    	CountTable counts = new CountTable(discrete); //Count table for discrete nodes in query 
	    	Map<Variable, Map<Integer, List<Object>>> store = new HashMap<>();
	    	for (Variable v : nonEnums) {
	    		Map<Integer, List<Object>> samples = new HashMap<>(); // Structure to store key/index with list of samples
	    		//With 1 discrete node in query possible keys are t/f
	    		//all continuous values for 'states' where the discrete node is true are recorded
	    		//leaving raw data allows flexibility in which distribution you apply e.g. Gaussian, kernel density etc.
		    	for(int i = 0; i < ApproxInferCont.iterations ; i++) {
					List<Object> curKey = new ArrayList<Object>(pars.size());
					for (EnumVariable p : discrete) {
						curKey.add(this.getData(p).get(i));
					}
					counts.count(curKey.toArray());
					if (samples.containsKey(getIndex(curKey.toArray()))) {
						samples.get(getIndex(curKey.toArray())).add(this.getData(v).get(i));
					} else {
						samples.put(getIndex(curKey.toArray()), new ArrayList<Object>());
						samples.get(getIndex(curKey.toArray())).add(this.getData(v).get(i));
					}
				}
		    	store.put(v, samples);
	    	}
	    	this.nonEnumTables = store;
	    	this.counts = counts;
    	} else { //Continuous variables only in query
    		allContinuous = true; //Original storage data structure contains adequate info
    	}
    }
    
    /**
     * Retrieve the index for the specified key
     *
     * @param key the values by which the index is identified (order the same as
     * when constructing the factor table)
     * @return the index for the instantiated key
     */
    public int getIndex(Object[] key) {
        int sum = 0;
        if (key.length != discrete.size()) {
            throw new EnumTableRuntimeException("Invalid key: length is " + key.length + " not " + nParents);
        }
        for (int i = 0; i < discrete.size(); i++) {
            if (key[i] == null) {
                throw new EnumTableRuntimeException("Null in key");
            }
            sum += (discrete.get(i).getIndex(key[i]) * step[i]);
        }
        return sum;
    }
    
    /**
     * Normalize the counts generated by the enum/discrete variables in the query
     * Method can be used both with a mixed query and a query with only discrete nodes
     * If query is mixed - getMixedDistrib() will also have to be used
     * 
     * @return table containing normalized counts
     */
    //FIXME - more sophisticated normalization technique?
    public EnumTable<Double> getNormalizedCounts() {
    	if (counts != null) {
    		EnumTable<Double> normalized = new EnumTable(counts.table.getParents());
    		for (Entry<Integer, Double> entry : counts.table.getMapEntries()) {
    			double nobserv = entry.getValue().doubleValue();
    			Object[] cntkey = counts.table.getKey(entry.getKey().intValue());
    			double norm = nobserv/ApproxInferCont.iterations;
    			normalized.setValue(cntkey, norm);
    		}
    		return normalized;
    	} else {
    		System.out.println("No count table exists to normalize");
    		return null;
    	}
    }
    
    /**
     * Use this method when all query nodes are continuous
     * Calculates a gaussian distribution for the set of samples generated
     * @return map of variable and associated distribution
     */
    public Map<Variable, Distrib> getGaussianDistrib() {
    	if (allContinuous) {
    		Map<Variable, Distrib> output = new HashMap<Variable, Distrib>();
    		for (Entry<Variable, List<Object>> entry : map.entrySet()) {
    			List<Object> data = entry.getValue();
    			
    			double sum = 0; //sum of all counts
    			for (Object d : data) {
    				sum += (Double)d;
    			}
    			double mean = sum/data.size();
                double diff = 0;
                for (int jj = 0; jj < data.size(); jj++) {
                	//FIXME - weighted based on prob from JPT?
                    diff += ((mean - (Double)data.get(jj)) * (mean - (Double)data.get(jj)));
                }
                double variance = diff/data.size();
                GaussianDistrib d = new GaussianDistrib(mean, variance);
                output.put(entry.getKey(), d);
                System.out.println();
    		}
    		return output;
    	}
    	System.out.println("Cannot use this function when discrete variables exist in query");
    	return null;
    }
    
    
    /**
     * Use this method for processing samples from a hybrid query
     * @return map of variable and associated enum table representing the distribution
     */
    public Map<Variable, EnumTable<Distrib>> getMixedGDistrib() {
    	if (!allContinuous){
    		Map<Variable, EnumTable<Distrib>> result = new HashMap<>();
    		for(Entry<Variable, Map<Integer, List<Object>>> var : nonEnumTables.entrySet()) {
    			EnumTable<Distrib> output = new EnumTable(this.discrete);
	    		for (Entry<Integer, List<Object>> entry : var.getValue().entrySet()) {
	    			List<Object> data = entry.getValue();
	    			double sum = 0; //sum of all counts
	    			for (Object d : data) {
	    				sum += (Double)d;
	    			}
	    			double mean = sum/data.size();
	                double diff = 0;
	                for (int jj = 0; jj < data.size(); jj++) {
	                	//FIXME - weighted based on prob from JPT?
	                    diff += ((mean - (Double)data.get(jj)) * (mean - (Double)data.get(jj)));
	                }
	                double variance = diff/data.size();
	                if (variance < 0.01) {
	                	variance = 0.01;
	                }
	                GaussianDistrib d = new GaussianDistrib(mean, variance);
	                output.setValue(entry.getKey(), d);
	    		}
	    		result.put(var.getKey(), output);
    		}
    		return result;
    	}
    	System.out.println("Cannot use this function when discrete variables exist in query");
    	return null;
    }
    
    /**
     * Get the canonical names of the parent variables (names + "." + index)
     *
     * @return names of variables (in order)
     */
    public String[] getLabels() {
        String[] labels = new String[nParents];
        for (int i = 0; i < labels.length; i++) {
            labels[i] = this.parents.get(i).toString();
        }
        return labels;
    }
       
    /**
     * get the list of data/samples associated with a query variable
     * @param key query variable of interest
     * @return list of data/samples
     */
    public List<Object> getData(Variable key) {
    	List<Object> samples = map.get(key);
        return samples;
    }
    
    /** 
     * Add the current instance of the variable to the appropriate list
     * @param key - current variable
     * @param instance - variable's instance
     */
    public void addValue(Variable key, Object instance) {
    	map.get(key).add(instance);
    }
    
    public boolean allContinuous() {
    	return allContinuous;
    }
    
    public List<EnumVariable> getDiscreteNodes() {
    	return discrete;
    }
    
    public Map<Variable, Map<Integer, List<Object>>> getNonEnumTable() {
    	return nonEnumTables;
    }
	
}
