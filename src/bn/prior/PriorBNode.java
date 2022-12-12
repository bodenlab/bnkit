package bn.prior;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import dat.EnumTable;
import dat.EnumVariable;
import dat.Variable;
import bn.BNode;
import bn.Distrib;
import bn.Sample;
import bn.factor.AbstractFactor;
import bn.factor.Factor;

public class PriorBNode implements BNode {

	protected BNode bNode;
	
	/**
	 * This enumtable is for prior
	 * each of them corresponding to the enumtable for EnumDistrib 
	 */
	private EnumTable<Prior> PriorTable;
	// prior for root
	private Prior rootPrior;
	// user-defined uniform distribution
	private Prior uniformPrior;
	
	public PriorBNode(BNode node) {
		bNode = node; 
		uniformPrior = null;
		rootPrior = null;
		if(!bNode.isRoot()) {
			PriorTable = new EnumTable<Prior>(bNode.getParents());
		}
	}
	
	public PriorBNode(BNode node, Prior uniform) {
		bNode = node; 
		uniformPrior = null;
		rootPrior = null;
		uniformPrior = uniform;
		if(!bNode.isRoot()) {
			PriorTable = new EnumTable<Prior>(bNode.getParents());
		}
	}
	
	/**
	 * set prior for the root node in BN
	 * @param _prior
	 */
	public void setPrior(Prior prior) {
		if(isRoot() && prior != null) {
			rootPrior = prior;
		}else {
			System.err.println("This BNode is conditioned on other nodes, need parents speicficed or prior is null");
		}
	}
	
	/**
	 * set prior for the conditioned node
	 * @param key, the array of parents' value
	 * @param _prior, the prior distribution
	 */
	public void setPrior(Object[] key, Prior prior) {
		if((!isRoot()) && prior != null) {
			this.PriorTable.setValue(key, prior);
		}else {
			System.err.println("This BNode is root node, no keys are needed or prior is null");
		}
	}
	
	/**
	 * set prior for the conditioned node
	 * @param keyIndex, the key index for the combination of parents' value
	 * @param prior the prior distribution
	 */
	public void setPrior(int keyIndex, Prior prior) {
		if((!isRoot()) && prior != null) {
			this.PriorTable.setValue(keyIndex, prior);
		}else {
			System.err.println("This BNode is root node, no keys are needed or prior is null");
		}
	}
	
	public void setUniformPrior(Prior uniform) {
		uniformPrior = uniform;
	}
	
	@Override
	public String getName() {
		return bNode.getName();
	}

	@Override
	public Double get(Object[] key, Object value) {
		return bNode.get(key, value);
	}

	@Override
	public Double get(Object value, Object... key) {
		return bNode.get(value, key);
	}

	@Override
	public Double get(Object value) {
		return bNode.get(value);
	}

	@Override
	public Variable getVariable() {
		return bNode.getVariable();
	}

	@Override
	public List<EnumVariable> getParents() {
		return bNode.getParents();
	}

	@Override
	public EnumTable getTable() {
		return bNode.getTable();
	}

	@Override
	public Distrib getDistrib(Object[] key) {
		return bNode.getDistrib(key);
	}

	@Override
	public Distrib getDistrib() {
		return bNode.getDistrib();
	}

	@Override
	public void print() {
		bNode.print();
	}

	@Override
	public String getType() {
		return bNode.getType();
	}

	@Override
	public String getStateAsText() {
		return bNode.getStateAsText();
	}

	@Override
	public boolean setState(String dump) {
		return bNode.setState(dump);
	}

	@Override
	public boolean isRoot() {
		return bNode.isRoot();
	}

	@Override
	public void setInstance(Object value) {
		bNode.setInstance(value);
	}

	@Override
	public void resetInstance() {
		bNode.resetInstance();
	}

	@Override
	public Object getInstance() {
		return bNode.getInstance();
	}

	@Override
	public Distrib makeDistrib(Collection<Sample> samples) {
		return bNode.makeDistrib(samples);
	}

	@Override
	public Factor makeFactor(Map<Variable, Object> rel) {
		return bNode.makeFactor(rel);
	}

	@Override
	public AbstractFactor makeDenseFactor(Map<Variable, Object> rel) {
		return bNode.makeDenseFactor(rel);
	}

	@Override
	public void countInstance(Object[] key, Object value, Double prob) {
		bNode.countInstance(key, value, prob);

	}

	@Override
	public void countInstance(Object[] key, Object value) {
		bNode.countInstance(key, value);
	}

	@Override
	public void maximizeInstance() {
		// get the number of condition this Bnode has
		List<EnumVariable> parents = getParents();
		int conditionNum = 1;
		if(!isRoot()) {
			for(EnumVariable parent: parents) {
				conditionNum *= parent.getDomain().size();
			}
		}
		
		// linear go through all conditions
		for(int i = 0; i < conditionNum; i++) {
			int index = i;
			if(isRoot()) {
				index = -1;
			}
			// get condition data
			List<Sample> samples = this.getConditionDataset(index);
			// convert the format
			Object[] data = new Object[samples.size()];
			double[] prob = new double[samples.size()];
			int count = 0;
			for(Sample sample: samples) {
				data[count] = sample.instance;
				prob[count] = sample.prob;
				count ++;
			}
			// get corresponding prior
			Prior prior = null;
			if(isRoot()) {
				prior = this.rootPrior;
			} else {
				prior = PriorTable.getValue(i);
			}
			// if there is no prior for that condition, we use
			// user-defined uniform prior
			if(prior == null) {
				prior = uniformPrior;
				if(prior == null) {
					throw new RuntimeException("cannot find any unifor Prior for this node");
				}
			}
			// prior start to learn from data
			prior.setEstimatedDistrib(getlikelihoodDistrib());
			prior.learn(data, prob);
			if(!isRoot()) {
				put(i,prior.getEstimatedDistrib());
			} else {
				put(prior.getEstimatedDistrib());
			}
			// reset parameters
			prior.resetParameters();
		}
	}

	@Override
	public boolean isTrainable() {
		return bNode.isTrainable();
	}

	@Override
	public void randomize(long seed) {
		bNode.randomize(seed);
	}

	@Override
	public void setRelevant(boolean relevant) {
		bNode.setRelevant(relevant);
	}

	@Override
	public void setTrainable(boolean trainable) {
		bNode.setTrainable(trainable);
	}

	@Override
	public boolean isRelevant() {
		return bNode.isRelevant();
	}

	@Override
	public List<Sample> getConditionDataset(int conditionIndex) {
		return bNode.getConditionDataset(conditionIndex);
	}

	@Override
	public Distrib getlikelihoodDistrib() {
		return bNode.getlikelihoodDistrib();
	}

	@Override
	public void put(Object[] key, Distrib distr) {
		bNode.put(key, distr);
		
	}

	@Override
	public void put(Distrib prob) {
		bNode.put(prob);
		
	}

	@Override
	public void put(Distrib prob, Object... key) {
		bNode.put(prob, key);
		
	}

	@Override
	public void put(int index, Distrib distr) {
		bNode.put(index, distr);
		
	}

}
