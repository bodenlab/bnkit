package bn.example;

import java.util.List;

import bn.BNet;
import bn.BNode;
import bn.alg.EM;
import bn.alg.LearningAlg;
import bn.file.BNBuf;
import bn.file.DataBuf;

public class LoadNTrainParallel {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 4) {
			System.out.println("Usage: LoadNTrain <bn-file> <data-file> <EMOption> <ThreadCount>");
			System.exit(1);
		}
		String bn_file = args[0];
		String data_file = args[1];
		int EMOption = Integer.parseInt(args[2]);
		int threads = Integer.parseInt(args[3]);
		BNet bn = BNBuf.load(bn_file);
		List<BNode> nodes = bn.getOrdered();
		Object[][] data = DataBuf.load(data_file, nodes);

		EM em = new EM(bn);
		em.setEMOption(EMOption); //Select parallel implementation
		em.setThreadCount(threads); //Set number of threads to be used
		em.train(data, nodes);
		BNBuf.save(bn, bn_file + "2.new");
		
//		TestNetwork tn = new TestNetwork(bn, data, "Variance");
	}

}
