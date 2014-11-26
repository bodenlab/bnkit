package bn.example;
import dat.EnumVariable;
import dat.Enumerable;
import bn.prob.EnumDistrib;
import bn.node.CPT;
import bn.BNet;
import bn.Predef;
/**
 * used for secondary structure prediction
 * with algiment 
 * @author wangyufei
 */

public class SecondaryStructurePredictor {
	
	final static private int windowSize = 5;
	
	public static void main(String[] args) {
		EnumVariable classNode = Predef.SecondaryStructure("secondary structure class");
		BNet bn = new BNet();
	
		CPT head = new CPT(classNode);
		head.put(new EnumDistrib(classNode.getDomain(),0.2,0.3,0.5));
		
		bn.add(head);
		for(int i = 0; i < windowSize; i++) {
			EnumVariable feature = Predef.AminoAcid("feature" + String.valueOf(i));
			CPT featureNode = new CPT(feature,classNode);
			EnumDistrib distrib = new EnumDistrib(Enumerable.aacid);
			for(int j = 0; j < classNode.getDomain().size(); j++){
				for(int x = 0; x < feature.size(); x++) {
					distrib.set(feature.getDomain().get(x), 1 / feature.size());
				}
				featureNode.put(distrib,classNode.getDomain().get(j));
			}
			bn.add(featureNode);
		}
		
		
	}
	
}
