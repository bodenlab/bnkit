package bn.alg;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import dat.EnumTable;
import dat.EnumVariable;
import dat.Variable;

import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * learning algorithm for fully observed data 
 * Now EM & MAP are both available for training complete data 
 * @author wangyufei
 *
 */
public class MAP extends LearningAlg {

	public int THREAD_COUNT = 10;
	public int BATCH_SIZE = 1;

	public MAP(BNet bn) {
		super(bn);
	}

	public static int[][] indexMe(BNet bn, Variable[] vars) {
		// Create an index for node to parent variables;
		// this is so that values associated with each node can be encoded, repeated for each row
		List<Integer>[] parentIndex = new List[vars.length];
		int[][] index = new int[vars.length][];
		// first, initialise the index
		for (int i = 0; i < vars.length; i++)
			parentIndex[i] = new ArrayList<>();
		// next, link the indices up
		for (int i = 0; i < vars.length; i++) {
			BNode node = bn.getNode(vars[i]);
			if (node == null)
				throw new RuntimeException("Variable is not in the BN");
			else {
				List<EnumVariable> parents = node.getParents();
				if (parents != null) {
					// a trick to make this faster: two different strategies to map indices
					// TODO: try different parent sizes...
					if (parents.size() > 2) {
						int cnt = parents.size();
						Set<EnumVariable> parentset = new HashSet<>(parents);
						for (int k = 0; k < vars.length; k++) {
							if (parentset.contains(vars[k])) {
								parentIndex[i].add(k);
								if (--cnt <= 0)
									break;
							}
						}
					} else {
						for (int j = 0; j < parents.size(); j++) {
							for (int k = 0; k < vars.length; k++) {
								if (parents.get(j).equals(vars[k])) {
									parentIndex[i].add(k);
									break;
								}
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < vars.length; i++) {
			index[i] = new int[parentIndex[i].size()];
			for (int j = 0; j < parentIndex[i].size(); j++)
				index[i][j] = parentIndex[i].get(j);
		}
		return index;
	}

	@Override
	public void train(Object[][] values, Variable[] vars, long seed) {

		// Create an index for node to parent variables;
		// this is so that values associated with each node can be encoded, repeated for each row
		int[][] index = indexMe(bn, vars);

		// Create an ExecutorService with a fixed thread pool
		ExecutorService executorService = Executors.newFixedThreadPool(THREAD_COUNT);

		// Submit tasks to the ExecutorService
		// one row of values to be coupled to their variables for each thread
		// each thread goes through each variable and updates associated counts
		for (int i = 0; i < values.length; i++) {
			MAPTraining runnable = new MAPTraining(index, vars, values[i], (i + 1));
			executorService.submit(runnable);
			if ((i + 1) % 1000 == 0 || (i + 1) == values.length) {
				System.out.println("submitted  " + (i + 1) + " instances");
			}
		}
		// Shut down the executor service
		executorService.shutdown();
		try {
			if (executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS)) {
				System.out.println("All threads have finished, now maximising likelihood.");
				// finally update the probabilities in all the tables
				for (int j = 0; j < vars.length; j++) {
					BNode node = bn.getNode(vars[j]);
					node.maximizeInstance();
				}
				System.out.println("Maximum likelihood finished.");
			}
		} catch (InterruptedException e) {
			throw new MAPRuntimeException("Error waiting for threads to finish");
		}
	}

	public void test(Object[][] values, Variable[] vars, EnumVariable[] qvars) {
		//int[][] index = indexMe(bn, vars);
		VarElim ve = new VarElim();
		ve.instantiate(bn);
		for (int i = 0; i < values.length; i++) {
			bn.resetNodes();
			for (int j = 0; j < vars.length; j++) {
				BNode node = bn.getNode(vars[j]);
				if (node != null)
					node.setInstance(values[i][j]);
			}
			for (int k = 0; k < qvars.length; k++) {
				Query q = ve.makeQuery(qvars[k]);
				CGTable r = (CGTable) ve.infer(q);
				Distrib d = r.query(qvars[k]);
				for (int j = 0; j < qvars[k].getDomain().size(); j++) {
					Object instance = qvars[k].getDomain().get(j);
					Double p = d.get(instance);
					System.out.print("P(" + qvars[k].getName() + "=" + instance + ") = " + p + "\t");
				}
			}
			System.out.println();
		}
	}

	public void test(Object[][] values, Variable[] vars) {
		VarElim ve = new VarElim();
		ve.instantiate(bn);
		for (int i = 0; i < values.length; i ++) {
			bn.resetNodes();
			for (int j = 0; j < vars.length; j ++) {
				BNode node = bn.getNode(vars[j]);
				if (node != null)
					node.setInstance(values[i][j]);
			}
			Query q = ve.makeMPE();
			CGTable r = (CGTable) ve.infer(q);
			Variable.Assignment[] as = r.getMPE();
			for (Variable.Assignment a : as)
				System.out.println("\t" + a);
			System.out.println("------------------------------------");
		}
	}

	/**
	 * Runnable 'worker' for MAP learning*
	 */
	public class MAPTraining implements Runnable {
		private int[][] index;
		private Object[] row;
		private Variable[] vars;
		private int job;

		public MAPTraining(int[][] index, Variable[] vars, Object[] values, int job) {
			this.index = index;
			this.row = values;
			this.vars = vars;
			this.job = job;
		}

		@Override
		public void run() {
			if ((job + 1) % 1000 == 0 || (job + 1) == vars.length) System.out.println("Job " + (job + 1) + " is running");
			for (int i = 0; i < vars.length; i ++) {
				BNode node = bn.getNode(vars[i]);
				if (row[i] != null && node.isTrainable()) {
					Object[] evid_key = null; // only applicable if the node has parents
					if (index[i].length > 0) { // has parents
						Object[] parvals = new Object[index[i].length];
						Variable[] parvars = new EnumVariable[index[i].length];
						for (int j = 0; j < index[i].length; j++) {
							parvals[j] = row[index[i][j]];
							parvars[j] = vars[index[i][j]];
						}
						evid_key = EnumTable.getKey(node.getParents(), Variable.Assignment.array(parvars, parvals));
					}
					node.countInstance(evid_key, row[i]); // this update is synchronised
				}
			}
			if ((job + 1) % 1000 == 0 || (job + 1) == vars.length) System.out.println("Done with job " + (job + 1));
		}
	}

	public class MAPRuntimeException extends RuntimeException {

		private static final long serialVersionUID = 1L;

		public MAPRuntimeException(String message) {
			super(message);
		}
	}

}
