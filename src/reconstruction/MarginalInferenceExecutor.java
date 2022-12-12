package reconstruction;

import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;

/**
 * Class to concurrently run joint reconstruction using a specified number of threads.
 * Breaks each independent node inference into a unit that is executed in parallel with all the others.
 *
 * Created by mikael on 26/8/17.
 */
public class MarginalInferenceExecutor {

    public static final int DEFAULT_NTHREADS = 8;

    // map to hold all jobs before they are dispatched
    private Map<Integer, ASRPOGMarginalInference> batchInferences = new HashMap<>();

    private int nThreads; // number of threads
    private ExecutorService executor = null; // will be started to manage specified number of threads
    private Map<Integer, Future<EnumDistrib>> jobs = new HashMap<>(); // map to hold results, associated with node ID

    /**
     * Create the executor, with threads set to a default value
     */
    public MarginalInferenceExecutor() {
        this(DEFAULT_NTHREADS);
    }

    /**
     * Create the executor
     * @param nThreads number of threads that the executor has access to
     */
    public MarginalInferenceExecutor(int nThreads) {
        this.nThreads = nThreads;
        executor = Executors.newFixedThreadPool(nThreads);
    }

    /**
     * Add a job
     * @param tag label by which the job is referred
     * @param ve the inference requested
     * @return true if the job is lodged successfully, false if the job was not accepted
     */
    public void addMarginalInference(Integer tag, VarElim ve, EnumVariable queryNode) {
        if (nThreads < 1)
            throw new RuntimeException("No threads have been allocated");
        batchInferences.put(tag, new ASRPOGMarginalInference(ve, queryNode));
    }

    public boolean isFull() {
        return false;
    }

    /**
     * Run all jobs and wait for them to be finished
     * @return all results
     */
    public Map<Integer, EnumDistrib> run() throws InterruptedException {
        for (Map.Entry<Integer, ASRPOGMarginalInference> entry : batchInferences.entrySet()) {
            Callable<EnumDistrib> worker = entry.getValue();
            Future<EnumDistrib> submit = executor.submit(worker);
            jobs.put(entry.getKey(), submit);
        }
        // now retrieve the result
        Map<Integer, EnumDistrib> res = new HashMap<>();
        for (Map.Entry<Integer, Future<EnumDistrib>> entry : jobs.entrySet()) {
            Integer nodeId = entry.getKey();
            Future<EnumDistrib> future = entry.getValue();
            try {
                res.put(nodeId, future.get());
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        return res;
    }

    public class ASRPOGMarginalInference implements Callable<EnumDistrib> {

        private VarElim ve = null;
        private EnumVariable queryNode = null;
        public Variable.Assignment[] charAssignments = new Variable.Assignment[]{};
        public EnumDistrib d_marg = null;

        public ASRPOGMarginalInference(VarElim ve, EnumVariable queryNode) {
            this.ve = ve;
            this.queryNode = queryNode;
        }

        @Override
        public EnumDistrib call() {
            try {
                Query q_marg = ve.makeQuery(queryNode);
                CGTable r_marg = (CGTable)ve.infer(q_marg);
                d_marg = (EnumDistrib)r_marg.query(queryNode);
            } catch (NullPointerException npe) { //When node of interest has been removed from network of interest
                d_marg = null;//EnumDistrib.uniform(Enumerable.aacid);
            }
            return d_marg;
        }
    }
}
