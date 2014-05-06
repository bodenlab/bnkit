package bn;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A structure which maintains the sample or count table for a query node
 * @author Alex
 *
 */

public class QueryTable {
	private BNode query = null;
	private CountTable count = null;
	private SampleTable sample = null;
	private MixtureDistrib mixDistrib = null;
	private BNet bn = null;
	
  
    public QueryTable(BNode query, BNet cbn) {
    	this.query = query;
    	this.bn = cbn;
    	List<EnumVariable> list = new ArrayList<EnumVariable>();
    	Boolean cgQuery = false;
		try {
			EnumVariable nVar = (EnumVariable)query.getVariable();
			List<EnumVariable> par = query.getParents();
			list.add(nVar);
			list.addAll(par);
		} catch (ClassCastException e){
			
			//Need the queries parents in the store table?
			//GDT cannot be currently trained as root
			List<EnumVariable> parents = query.getParents();
			if (parents != null) {
				cgQuery = true;
			}
		}
		
		if (cgQuery) {
			this.sample = new SampleTable(query.getParents());
		} else {
			this.count = new CountTable(list);
		}
    }
    
    public SampleTable getSampleTable() {
    	return sample;
    }
    
    public CountTable getCountTable() {
    	return count;
    }
    
    /**
     * count this observation
     * @param cbn the network
     */
    public void countInstance(BNet cbn) {
    	List<EnumVariable> parList = null;
    	parList = count.table.getParents();
		List<Object> instances = new ArrayList<Object>();
		for (EnumVariable par : parList) {
			instances.add(cbn.getNode(par).getInstance());
		}
		//How do you know parent query is in right order?
		//Query is a parent??
		count.count(instances.toArray());
    }
    
    /**
     * record this sample
     * @param cbn the network
     */
    public void sampleInstance(BNet cbn) {
    	List<EnumVariable> parList = null;
    	parList = sample.table.getParents();
		List<Object> instances = new ArrayList<Object>();
		for (EnumVariable par : parList) {
			instances.add(cbn.getNode(par).getInstance());
		}
		//How do you know parent query is in right order?
		//Query is a parent??
		sample.count(instances.toArray(), query.getInstance());
    }
       
    public BNode getQuery(){
    	return query;
    }
    
    public MixtureDistrib getMixtureDistrib() {
    	return mixDistrib;
    }
    
    public void setDistrib(MixtureDistrib dis) {
    	mixDistrib = dis;
    }
    
    public boolean hasSample() {
    	return !(sample == null);
    }
    
    public boolean hasCount() {
    	return !(count == null);
    }

}
