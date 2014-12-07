package bn.prior;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	// prior for each condition
	private EnumTable<Prior> PriorTable;
	// prior for root
	private Prior rootPrior;
	
	public CPTPrior(EnumVariable var, List<EnumVariable> parents) {
		super(var, parents);
		if (parents != null) {
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
            if (parents.length > 0) {
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
	
	public void setPrior(Prior _prior) {
		if(rootPrior != null && _prior != null) {
			rootPrior = _prior;
		}else {
			System.err.println("This CPT is conditioned on other nodes, need parents speicficed or prior is null");
		}
	}
	
	public void setPrior(Object[] key, Prior _prior) {
		if(PriorTable != null && _prior != null) {
			this.PriorTable.setValue(key, _prior);
		}else {
			System.err.println("This CPT is root node, no keys are needed or prior is null");
		}
	}
	
	public void setPrior(int keyIndex, Prior prior) {
		if(PriorTable != null) {
			this.PriorTable.setValue(keyIndex, prior);
		}else {
			System.err.println("This CPT is root node, no keys are needed");
		}
	}
	
	@Override
    public void maximizeInstance() {
		Map<Integer, Double[]> data = new HashMap<Integer, Double[]>();
		EnumVariable var = getVariable();
		
		if (count.table.isEmpty()) {
            return;
        }
		if (table != null) { // there are parents in the CPT
            // Set all 'old' distributions in the CPT to valid = false, i.e.
            // we are marking entries so we can remove 'ghosts' after counting
            for (EnumDistrib d : this.table.getValues()) {
                d.setValid(false);
            }
            
            for (Map.Entry<Integer, Double> entry : count.table.getMapEntries()) {
            	double nobserv = entry.getValue();
                Object[] cntkey = count.table.getKey(entry.getKey().intValue());
                Object[] cptkey = new Object[cntkey.length - 1];
                for (int i = 0; i < cptkey.length; i++) {
                    cptkey[i] = cntkey[i + 1];
                }
                Integer index = new Integer(this.table.getIndex(cptkey));
                if(data.containsKey(index)) {
                	data.get(index)[var.getIndex(cntkey[0])] = nobserv;
                } else {
                	Double[] subData = new Double[var.size()];
                	subData[var.getIndex(cntkey[0])] = nobserv;
                	data.put(index, subData);
                }
            }
            
            Object[] enumObject = var.getDomain().getValues();
            for (Map.Entry<Integer, Double[]> entry : data.entrySet()) {
            	Prior prior = PriorTable.getValue(entry.getKey().intValue());
            	double[] hist = new double[var.size()];
            	EnumDistrib resultDistrib = null;
            	for(int i = 0; i < var.size(); i++) {
            		hist[i] = entry.getValue()[i].doubleValue();
            	}
            	if(prior != null) {
            		prior.setLikelihoodDistrib(new EnumDistrib(var.getDomain()));
            		prior.learn(var.getDomain().getValues(), hist);
            		resultDistrib = (EnumDistrib)prior.getMAPDistrib();
            	}else {
            		UniformPrior uniPrior = new UniformPrior();
            		uniPrior.setLikelihoodDistrib(new EnumDistrib(var.getDomain()));
            		uniPrior.learn(var.getDomain().getValues(), hist);
            		resultDistrib = (EnumDistrib)uniPrior.getMAPDistrib();
            	}
            	table.setValue(entry.getKey().intValue(), resultDistrib);
            }
            
		} else {
			Object[] cntkey = new Object[1];
            double[] cnts = new double[var.size()];
            for (int i = 0; i < var.size(); i++) {
                cntkey[0] = var.getDomain().get(i);
                cnts[i] = count.get(cntkey);
            }
            rootPrior.setLikelihoodDistrib(new EnumDistrib(var.getDomain()));
            rootPrior.learn(var.getDomain().getValues(), cnts);
            put((EnumDistrib)rootPrior.getMAPDistrib());
		}
		
		count.table.setEmpty(); // reset counts
	}
	
	
	

}
