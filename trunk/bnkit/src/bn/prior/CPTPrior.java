package bn.prior;

import java.util.ArrayList;
import java.util.List;

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
	
	public void setPrior(Prior prior) {
		if(rootPrior != null) {
			rootPrior = prior;
		}else {
			System.err.println("This CPT is conditioned on other nodes, need parents speicficed");
		}
	}
	
	public void setPrior(Object[] key, Prior prior) {
		if(PriorTable != null) {
			this.PriorTable.setValue(key, prior);
		}else {
			System.err.println("This CPT is root node, no keys are needed");
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
		
	}
	
	
	

}
