package reconstruction;

import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Class to concurrently run joint reconstruction using a specified number of threads.
 * Breaks each independent node inference into a unit that is executed in parallel with all the others.
 *
 * Created by mikael on 26/8/17.
 */
public class JointInferenceExecutor {

    public static final int DEFAULT_NTHREADS = 8;

    // map to hold all jobs before they are dispatched
    private Map<Integer, ASRPOGJointInference> batchInferences = new HashMap<>();

    private int nThreads; // number of threads
    private ExecutorService executor = null; // will be started to manage specified number of threads
    private Map<Integer, Future<Variable.Assignment[]>> jobs = new HashMap<>(); // map to hold results, associated with node ID

    /**
     * Create the executor, with threads set to a default value
     */
    public JointInferenceExecutor() {
        this(DEFAULT_NTHREADS);
    }

    /**
     * Create the executor
     * @param nThreads number of threads that the executor has access to
     */
    public JointInferenceExecutor(int nThreads) {
        this.nThreads = nThreads;
        executor = Executors.newFixedThreadPool(nThreads);
    }

    /**
     * Add a job
     * @param tag label by which the job is referred
     * @param ve the inference requested
     * @return true if the job is lodged successfully, false if the job was not accepted
     */
    public void addJointInference(Integer tag, VarElim ve) {
        if (nThreads < 1)
            throw new RuntimeException("No threads have been allocated");
        batchInferences.put(tag, new ASRPOGJointInference(ve));
    }

    public boolean isFull() {
        return false;
    }

    /**
     * Run all jobs and wait for them to be finished
     * @return all results
     */
    public Map<Integer, Variable.Assignment[]> run() throws InterruptedException {
        for (Map.Entry<Integer, ASRPOGJointInference> entry : batchInferences.entrySet()) {
            Callable<Variable.Assignment[]> worker = entry.getValue();
            Future<Variable.Assignment[]> submit = executor.submit(worker);
            jobs.put(entry.getKey(), submit);
        }
        // now retrieve the result
        Map<Integer, Variable.Assignment[]> res = new HashMap<>();
        for (Map.Entry<Integer, Future<Variable.Assignment[]>> entry : jobs.entrySet()) {
            Integer nodeId = entry.getKey();
            Future<Variable.Assignment[]> future = entry.getValue();
            try {
                res.put(nodeId, future.get());
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        return res;
    }

    public class ASRPOGJointInference implements Callable<Variable.Assignment[]> {

        private VarElim ve = null;
        private EnumVariable queryNode = null;
        public Variable.Assignment[] charAssignments = new Variable.Assignment[]{};
        public EnumDistrib d_marg = null;

        public ASRPOGJointInference(VarElim ve) {
            this.ve = ve;
        }

        @Override
        public Variable.Assignment[] call() {
            Query q_joint = ve.makeMPE();
            CGTable r_joint = (CGTable) ve.infer(q_joint);
            charAssignments = r_joint.getMPE();
            return charAssignments;
        }
    }
}
