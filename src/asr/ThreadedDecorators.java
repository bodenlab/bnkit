package asr;

import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;

/**
 * this is the thread executor class
 * @param <E>
 */
public class ThreadedDecorators<E> { //
    private ExecutorService executor;
    private int nThreads;
    private Map<Integer, Job<E>> jobsToDo = new HashMap<>();
    private Map<Integer, Future<TreeDecor<E>>> jobsDone = new HashMap<>();

    /**
     * Create a batch of tree decorators.
     * Note that decors will be deployed (decorated) under the control of this thread executor class, queried and the outputs are stored
     * @param decors tree decors, each instantiated (if null, it will be ignored)
     * @param tis tree instantiations with observed values, one for each tree decorator
     * @param nThreads number of threads to use
     * @return
     */
    public ThreadedDecorators(TreeDecor[] decors, TreeInstance[] tis, int nThreads) {
        if (decors.length != tis.length)
            throw new ASRRuntimeException("Mismatch in batch jobs");
        this.nThreads = nThreads;
        executor = Executors.newFixedThreadPool(nThreads);
        // create jobs
        for (int i = 0; i < decors.length; i ++) {
            if (decors[i] != null) {
                jobsToDo.put(i, new Job<>(decors[i], tis[i]));
            }
        }
    }

    /**
     * Run all the jobs
     * @return a map indexed by job number, pointing to result instance that can be queried
     * @throws InterruptedException
     */
    public Map<Integer, TreeDecor<E>> runBatch() throws InterruptedException {
        if (GRASP.VERBOSE)
            System.out.println("Running jobs: " + jobsToDo.size() + " with " + nThreads + " threads");

        for (Map.Entry<Integer, Job<E>> entry : jobsToDo.entrySet()) {
//            if (GRASP.VERBOSE)
//                System.out.println("Considering thread " + job.toString() + " [" + !job.isCompleted() + "]");
            Callable<TreeDecor<E>> worker = entry.getValue();
            Future<TreeDecor<E>> submit = executor.submit(worker);
            jobsDone.put(entry.getKey(), submit);
        }
        // retrieve results
        Map<Integer, TreeDecor<E>> res = new HashMap<>();
        for (Map.Entry<Integer, Future<TreeDecor<E>>> entry : jobsDone.entrySet()) {
            Integer tag = entry.getKey();
            Future<TreeDecor<E>> future = entry.getValue();
            try {
                TreeDecor ret = future.get();
                res.put(tag, ret);
            } catch (InterruptedException | ExecutionException e) {
                System.err.println("Failed with thread for " + tag + " with future " + future);
                e.printStackTrace();
            }
        }
        executor.shutdown();
        return res;
    }

    public class Job<E> implements Callable<TreeDecor<E>> {

        private TreeDecor<E> decorator;
        private TreeInstance ti;

        public Job(TreeDecor<E> decorator, TreeInstance ti) {
            this.decorator = decorator;
            this.ti = ti;
        }

        @Override
        public TreeDecor<E> call() throws Exception {
//            if (GRASP.VERBOSE)
//                System.out.println("Running thread " + this.toString());
            decorator.decorate(ti);
            TreeDecor<E> ret = decorator;
            return ret;
        }
    }

}
