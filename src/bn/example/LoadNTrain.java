/**
 * 
 */
package bn.example;
import bn.alg.LearningAlg;
import bn.alg.EM;

import java.util.List;

import bn.*;
import bn.file.*;

/**
 * @author mikael
 *
 */
public class LoadNTrain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.println("Usage: LoadNTrain <bn-file> <data-file>");
			System.exit(1);
		}
		String bn_file = args[0];
		String data_file = args[1];
		BNet bn = BNBuf.load(bn_file);
		List<BNode> nodes = bn.getOrdered();
		Object[][] data = DataBuf.load(data_file, nodes);

		LearningAlg em = new EM(bn);
		em.train(data, nodes);
		BNBuf.save(bn, bn_file + "2.new");
		
//		TestNetwork tn = new TestNetwork(bn, data, "Variance");
	}

}
