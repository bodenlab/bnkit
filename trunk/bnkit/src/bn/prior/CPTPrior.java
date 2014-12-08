package bn.prior;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import dat.EnumTable;
import dat.EnumVariable;
import bn.JPT;
import bn.node.CPT;
import bn.prob.EnumDistrib;

/**
 * This class is only used for check if the prior idea works all right
 * don't use it for constructing BN!
 * @author wangyufei
 *
 */
public class CPTPrior extends CPT {
	
	private static final long serialVersionUID = 1L;
	/**
	 * This enumtable is for prior
	 * each of them corresponding to the enumtable for EnumDistrib 
	 */
	private EnumTable<Prior> PriorTable;
	// prior for root
	private Prior rootPrior;
	
	public CPTPrior(EnumVariable var, List<EnumVariable> parents) {
		super(var, parents);
		if (parents != null) { // if this node has parent
            if (parents.size() > 0) {
            	this.PriorTable = new EnumTable<>(parents);
            	rootPrior = null;
            	return;
            }
		}
		this.PriorTable = null;
		rootPrior = new UniformPrior();
	}

	public CPTPrior(EnumVariable var, EnumVariable... parents) {
		super(var, parents);
		if (parents != null) {
            if (parents.length > 0) { // if this node has parent
            	this.PriorTable = new EnumTable<>(parents);
            	rootPrior = null;
            	return;
            }
		}
		this.PriorTable = null;
		rootPrior = new UniformPrior();
	}

	public CPTPrior(EnumVariable var) {
		super(var);
		this.PriorTable = null;
		rootPrior = new UniformPrior();
	}

	public CPTPrior(JPT jpt, EnumVariable var) {
		super(jpt, var);
		List<EnumVariable> cptParents = new ArrayList<>(jpt.getParents().size() - 1);
		for (int i = 0; i < jpt.getParents().size(); i++) {
            EnumVariable jptParent = jpt.getParents().get(i);
            if (jptParent != var) {
            	 cptParents.add(jptParent);
            }
        }
		if(jpt.getParents().size() == 0) {
			this.PriorTable = null;
			rootPrior = new UniformPrior();
		}else {
			
			this.PriorTable = new EnumTable<>(cptParents);
        	rootPrior = null;
		}
	}
	
	/**
	 * set prior for the root node in BN
	 * @param _prior
	 */
	public void setPrior(Prior _prior) {
		if(rootPrior != null && _prior != null) {
			rootPrior = _prior;
		}else {
			System.err.println("This CPT is conditioned on other nodes, need parents speicficed or prior is null");
		}
	}
	
	/**
	 * set prior for the conditioned node
	 * @param key, the array of parents' value
	 * @param _prior, the prior distribution
	 */
	public void setPrior(Object[] key, Prior _prior) {
		if(PriorTable != null && _prior != null) {
			this.PriorTable.setValue(key, _prior);
		}else {
			System.err.println("This CPT is root node, no keys are needed or prior is null");
		}
	}
	
	/**
	 * set prior for the conditioned node
	 * @param keyIndex, the key index for the combination of parents' value
	 * @param prior the prior distribution
	 */
	public void setPrior(int keyIndex, Prior prior) {
		if(PriorTable != null) {
			this.PriorTable.setValue(keyIndex, prior);
		}else {
			System.err.println("This CPT is root node, no keys are needed");
		}
	}
	
	@Override
    public void maximizeInstance() {
		if (count.table.isEmpty()) {
            return;
        }
		// The Integer is the index for parent value
		// The Double array is the count for enumvariable
		Map<Integer, Double[]> data = new HashMap<Integer, Double[]>();
		EnumVariable var = getVariable();
		
		
		if (table != null) { // there are parents in the CPT
            // Set all 'old' distributions in the CPT to valid = false, i.e.
            // we are marking entries so we can remove 'ghosts' after counting
            for (EnumDistrib d : this.table.getValues()) {
                d.setValid(false);
            }
            /**
             * This for loop is used for generate a map between the index of parent value 
             * and a count vector for EnumDistrib
             */
            for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
            	// get the count 
            	double nobserv = entry.getValue();
                Object[] cntkey = count.table.getKey(entry.getKey().intValue());
                Object[] cptkey = new Object[cntkey.length - 1];
                // get the array of parent value
                for (int i = 0; i < cptkey.length; i++) {
                    cptkey[i] = cntkey[i + 1];
                }
                // get index for the above array
                Integer index = new Integer(this.table.getIndex(cptkey));

                if(data.containsKey(index)) { // if the map contains this index
                	// set the count to the appropriate position in the "count vector"
                	data.get(index)[var.getIndex(cntkey[0])] = nobserv;
                } else { // otherwise... create a new count vector
                	Double[] subData = new Double[var.size()];
                	subData[var.getIndex(cntkey[0])] = nobserv;
                	data.put(index, subData);
                }
            }
            
            /**
             * this for loop is mainly used for setting the count vector
             * to the corresponding prior
             */
            Object[] enumObject = var.getDomain().getValues();
            for (Map.Entry<Integer, Double[]> entry : data.entrySet()) {
            	//use the index to find the prior
            	Prior prior = PriorTable.getValue(entry.getKey().intValue());
            	// count vector
            	double[] hist = new double[var.size()];
            	EnumDistrib resultDistrib = null;
            	// convert the double array 
            	for(int i = 0; i < var.size(); i++) {
            		hist[i] = entry.getValue()[i].doubleValue();
            	}
            	if(prior != null) { // baysian prior
            		prior.setLikelihoodDistrib(new EnumDistrib(var.getDomain()));
            		prior.learn(enumObject, hist);
            		resultDistrib = (EnumDistrib)prior.getMAPDistrib();
            		prior.resetParameters();
            	}else { // uniform prior, just ML
            		UniformPrior uniPrior = new UniformPrior();
            		uniPrior.setLikelihoodDistrib(new EnumDistrib(var.getDomain()));
            		uniPrior.learn(enumObject, hist);
            		resultDistrib = (EnumDistrib)uniPrior.getMAPDistrib();
            	}
            	table.setValue(entry.getKey().intValue(), resultDistrib);
            }
            
          //Remove 'old' (or 'ghost' entries from CPT (for which no counts
            for (Iterator<Entry<Integer, EnumDistrib>> it = table.getMapEntries().iterator(); it.hasNext(); ) {
                Entry<Integer, EnumDistrib> entry = it.next();
            	EnumDistrib obs = entry.getValue();
            	if (!obs.isValid()) 
                    it.remove();
            }
            
		} else { // the root node
			Object[] cntkey = new Object[1];
            double[] cnts = new double[var.size()];
            for (int i = 0; i < var.size(); i++) {
                cntkey[0] = var.getDomain().get(i);
                cnts[i] = count.get(cntkey);
            }
            rootPrior.setLikelihoodDistrib(new EnumDistrib(var.getDomain()));
            rootPrior.learn(var.getDomain().getValues(), cnts);
            put((EnumDistrib)rootPrior.getMAPDistrib());
            rootPrior.resetParameters();
		}
		
		count.table.setEmpty(); // reset counts
	}
	
	
	

}
