package reconstruction;

import api.PartialOrderGraph;
import dat.POGraph;

import java.util.*;
import java.util.concurrent.*;

/**
 * Class to concurrently assemble ancestor POGs using a specified number of threads.
 * Each ancestor can be assembled independently, in parallel with all the others.
 *
 * Created by mikael on 10/10/2019.
 */
public class POGAssemblyExecutor {

    public static final int DEFAULT_NTHREADS = 8;

    private int nThreads; // number of threads
    final private ASRPOG asr;

    private ExecutorService executor = null; // will be started to manage specified number of threads
    private Map<String, POGAssemblyWorker> batchMap = new HashMap<>();
    private Map<String, Future<PartialOrderGraph>> jobs = new HashMap<>(); // map to hold POGs, each keyed with an ancestor ID

    /**
     * Create the executor, with threads set to a default value
     */
    public POGAssemblyExecutor(ASRPOG asr) {
        this(asr, DEFAULT_NTHREADS);
    }

    /**
     * Create the executor
     * @param nThreads number of threads that the executor has access to
     */
    public POGAssemblyExecutor(ASRPOG asr, int nThreads) {
        this.asr = asr;
        this.nThreads = nThreads;
        executor = Executors.newFixedThreadPool(nThreads);
    }

    /**
     * Add a job
     * @param tag label by which the job is referred
     */
    public void addAssemblyJob(String tag) {
        if (nThreads < 1)
            throw new RuntimeException("No threads have been allocated");
        batchMap.put(tag, new POGAssemblyWorker(this.asr, tag));
    }
    /**
     * Run all jobs and wait for them to be finished
     * @return all results
     */
    public Map<String, PartialOrderGraph> run() throws InterruptedException {

        for (Map.Entry<String, POGAssemblyWorker> entry : batchMap.entrySet()) {
            Callable<PartialOrderGraph> worker = entry.getValue();
            Future<PartialOrderGraph> submit = executor.submit(worker);
            jobs.put(entry.getKey(), submit);
        }
        // now retrieve the result
        Map<String, PartialOrderGraph> res = new HashMap<>();
        for (Map.Entry<String, Future<PartialOrderGraph>> entry : jobs.entrySet()) {
            String ancId = entry.getKey();
            Future<PartialOrderGraph> future = entry.getValue();
            try {
                res.put(ancId, future.get());
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        executor.shutdown();
        return res;
    }

    public class POGAssemblyWorker implements Callable<PartialOrderGraph> {

        final private ASRPOG asr;
        final private String ancestor;

        public POGAssemblyWorker(ASRPOG asr, String ancestor) {
            this.asr = asr;
            this.ancestor = ancestor;
        }
        @Override
        public PartialOrderGraph call() {
            POGraph a = asr.getAncestor(ancestor);
            if (a == null) {
                throw new RuntimeException("Assembly job failed for ancestor \"" + ancestor + "\"");
            }
            PartialOrderGraph pog = new PartialOrderGraph(a);
            pog.getMostSupported(true);
            return pog;
        }
    }
}
