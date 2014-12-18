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

public abstract class PriorBNode implements BNode {

	protected BNode bNode;
	
	public PriorBNode(BNode node) {
		bNode = node; 
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
	public abstract void maximizeInstance();

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

}
