package asr;

import bn.ctmc.GapSubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.BranchPoint;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import smile.math.MathEx;
import java.util.*;
import java.util.concurrent.*;

public class ThreadedPeeler {

    private final ExecutorService executor;
    private final Map<Integer, Job> jobsToDo = new HashMap<>();
    private final Map<Integer, Future<IndelPeeler>> jobsDone = new HashMap<>();

    /**
     * Create a batch of jobs that will use each Tree instance to work
     * out the probability of observing a particular column.
     * @param peelers The Peeler class
     * @param nThreads number of threads to use
     */
    public ThreadedPeeler(IndelPeeler[] peelers, int nThreads) {

        this.executor = Executors.newFixedThreadPool(nThreads);
        // create jobs
        for (int i = 0; i < peelers.length; i++) {
            jobsToDo.put(i, new Job(peelers[i]));
        }
    }

    public Map<Integer, Double> runBatch() throws InterruptedException, ExecutionException {

        for (Map.Entry<Integer, Job> entry: jobsToDo.entrySet()) {
            Callable<IndelPeeler> worker = entry.getValue();
            Future<IndelPeeler> submit = executor.submit(worker);
            jobsDone.put(entry.getKey(), submit);
        }

        // retrieve results
        Map<Integer, Double> results = new HashMap<>();
        for (Map.Entry<Integer, Future<IndelPeeler>> entry: jobsDone.entrySet()) {
            Integer tag = entry.getKey();
            Future<IndelPeeler> future = entry.getValue();
            try {
                Double res = future.get().getDecoration();
                results.put(tag, res);

            } catch (InterruptedException | ExecutionException e) {
                System.err.println("Failed with thread for " + tag + " with future " + future);
                e.printStackTrace();
            }
        }

        executor.shutdown();
        return results;
    }


    /**
     *  Calculate the total probability of a given aligned
     *  column u i.e. P(col| Tree, Model, SeqLenParam). Uses Felsenstein's
     *  adapted pruning algorithm before extracting the probability of the root
     *  node. All possible residue assignments are summed over and also weighted
     *  by the geometric sequence length parameter (assumed to be the average sequence
     *  length of the alignment). Refer to <a href="https://doi.org/10.1371/journal.pcbi.1000172">Rivas & Eddy (2008)</a>
     *  for a description of Felsenstein's peeling algorithm extended to gaps.
     *
     * @return log (P(alignment col |Tree, Model, SeqLenParam))
     */
    public static double logProbColGivenRate(Tree tree, EnumSeq.Alignment<Enumerable> aln,
                                      Double colIndelRate, GapSubstModel model,
                                      int colIdx, double geometricSeqLenParam) {


        int totalNodes = tree.getNLeaves() + tree.getNParents();
        int alphabetSize = model.getDomain().size() - 1; // ignore gaps
        Object[] alphabet = model.getDomain().getValues();
        Double[][] nodeResidueProbs = new Double[totalNodes][alphabetSize]; // nodes x num_letters
        Double[]nodeGapProbs = new Double[totalNodes];

        // instantiate with negative infinity for subsequent log sum calculations
        for (int i = 0; i < totalNodes; i++) {
            Arrays.fill(nodeResidueProbs[i], Double.NEGATIVE_INFINITY);
            nodeGapProbs[i] = Double.NEGATIVE_INFINITY;
        }
        // update the nodeResidueProbs and nodeGapProbs arrays in place
        tree.felsensteinsExtendedPeeling(aln, colIdx, nodeResidueProbs, nodeGapProbs, colIndelRate, model);

        int ROOT_INDEX = 0;
        double rootGapProb = nodeGapProbs[ROOT_INDEX];

        // sum over possible residue assignments
        double[] residueTerms = new double[alphabetSize]; // number of alphabet letters
        for (int resIdx = 0; resIdx < alphabetSize; resIdx++) {
            double priorProb = Math.log(model.getProb(alphabet[resIdx]));
            double rootResidueProb = nodeResidueProbs[ROOT_INDEX][resIdx];
            residueTerms[resIdx] = priorProb + rootResidueProb; // weight by prior prob of residue
        }

        double weightedLogSumResidueProb = MathEx.logsumexp(residueTerms) + Math.log(geometricSeqLenParam);
        double[] finalColTerms = {rootGapProb, weightedLogSumResidueProb};
        return MathEx.logsumexp(finalColTerms);
    }

    public static class Job implements Callable<IndelPeeler> {

        private final IndelPeeler peeler;

        public Job(IndelPeeler peeler) {
            this.peeler = peeler;
        }

        public IndelPeeler call() throws Exception {
            this.peeler.decorate();
            return peeler;
        }
    }

}
