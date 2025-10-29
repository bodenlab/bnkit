package asr;

import bn.ctmc.GapSubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.Tree;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;

public class ThreadedPeeler {

    private final ExecutorService executor;
    private final Map<Integer, Job> jobsToDo = new HashMap<>();
    private final Map<Integer, Future<Peeler>> jobsDone = new HashMap<>();

    /**
     * Create a batch of jobs that will use each Tree instance to work
     * out the probability of observing a particular column.
     * @param peelers The Peeler class
     * @param nThreads number of threads to use
     */
    public ThreadedPeeler(Peeler[] peelers, int nThreads) {

        this.executor = Executors.newFixedThreadPool(nThreads);
        // create jobs
        for (int i = 0; i < peelers.length; i++) {
            jobsToDo.put(i, new Job(peelers[i]));
        }
    }

    public Map<Integer, Peeler> runBatch() throws InterruptedException, ExecutionException {

        for (Map.Entry<Integer, Job> entry: jobsToDo.entrySet()) {
            Callable<Peeler> worker = entry.getValue();
            Future<Peeler> submit = executor.submit(worker);
            jobsDone.put(entry.getKey(), submit);
        }

        // retrieve results
        Map<Integer, Peeler> results = new HashMap<>();
        for (Map.Entry<Integer, Future<Peeler>> entry: jobsDone.entrySet()) {
            Integer tag = entry.getKey();
            Future<Peeler> future = entry.getValue();
            try {
                Peeler res = future.get();
                results.put(tag, res);
            } catch (InterruptedException | ExecutionException e) {
                System.err.println("Failed with thread for " + tag + " with future " + future);
                e.printStackTrace();
            }
        }

        executor.shutdown();
        return results;
    }


    public static class Peeler {
        private final Tree tree;
        private final EnumSeq.Alignment<Enumerable> aln;
        private final Double rate;
        private final GapSubstModel model;
        private final int col_idx;
        private final double geometric_seq_len_param;
        private double col_prob = Double.NaN;

        public Peeler(Tree tree, EnumSeq.Alignment<Enumerable> aln, Double col_indel_rate,
                       GapSubstModel model, int col_idx, double geometricSeqLenParam) {

            this.tree = tree;
            this.aln = aln;
            this.rate = col_indel_rate;
            this.model = model;
            this.col_idx = col_idx;
            this.geometric_seq_len_param = geometricSeqLenParam;
        }

        public void decorate() {
            this.col_prob = this.tree.logProbColGivenRate(aln, rate, model,
                    col_idx, geometric_seq_len_param);
        }

        public double getDecoration() {
            return col_prob;
        }

    }

    public static class Job implements Callable<Peeler> {

        private final Peeler peeler;

        public Job(Peeler peeler) {
            this.peeler = peeler;
        }

        public Peeler call() throws Exception {
            this.peeler.decorate();
            return peeler;
        }
    }

}
