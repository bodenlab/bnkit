package bn.alg;

import bn.BNet;
import bn.BNode;
import dat.EnumTable;
import dat.Variable;

/**
 * learning algorithm for fully observed data 
 * Now EM & MAP are both available for training complete data 
 * @author wangyufei
 *
 */
public class MAP extends LearningAlg {

	public MAP(BNet bn) {
		super(bn);
	}

	@Override
	public void train(Object[][] values, Variable[] vars, long seed) {
		for (int i = 0; i < values.length; i++) {
			for (int j = 0; j < vars.length; j++) {
				Variable.Assignment[] evidence = Variable.Assignment.array(vars, values[i]);
				BNode node = bn.getNode(vars[j]);
				if(node == null) {
					System.err.println("variable is not in the BN");
					return;
				}
				
				if (values[i][j] != null) { 
					Object[] evid_key = null; // only applicable if the node has parents
                    if (!node.isRoot()) {
                        evid_key = EnumTable.getKey(node.getParents(), evidence);
                    }
                    node.countInstance(evid_key, values[i][j]);
				}
			}
		}
		
		for (int j = 0; j < vars.length; j++) {
			BNode node = bn.getNode(vars[j]);
			node.maximizeInstance();
		}

	}

}
