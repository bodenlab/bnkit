package asr;

import bn.ctmc.GapSubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.BranchPoint;
import dat.phylo.Tree;
import smile.math.MathEx;
import java.util.*;
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

    /**
     * Get the log likelihood of a particular alignment given the tree and gap augmented substitution model.
     * source: <a href="https://doi.org/10.1371/journal.pcbi.1000172">Rivas & Eddy, 2008</a>
     *
     * @param model the substitution model
     * @param aln the alignment
     * @param geometricSeqLenParam the geometric sequence length parameter
     * @param alpha the alphabet
     *
     * @return the log likelihood of the alignment given the tree and model
     */
    public static double calcProbAlnGivenTree(Tree tree, GapSubstModel model, EnumSeq.Alignment<Enumerable> aln,
                                              double geometricSeqLenParam, Enumerable alpha, int nThreads) {

        int numCols = aln.getWidth();
        Peeler[] peelers = new Peeler[numCols];
        // first compute each column likelihood
        for (int i = 0; i < numCols; i++) {
            GapSubstModel model_copy = model.deepCopy(); // copy to avoid race conditions
            peelers[i] = new Peeler(tree, aln, 1.0, model_copy, i, geometricSeqLenParam);
        }

        double LL = 0.0;
        ThreadedPeeler threadPool = new ThreadedPeeler(peelers, nThreads);
        try {
            Map<Integer, Peeler> ret = threadPool.runBatch();
            for (int i = 0; i < numCols; i++) {
                Peeler peeler = ret.get(i);
                LL += peeler.getDecoration(); // add each column likelihood
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        double probExtraCol = probExtraCol(tree, model, geometricSeqLenParam); // add the normalisation term

        // need to create a dummy aln containing only gaps
        List<EnumSeq.Gappy<Enumerable>> seqArray = new ArrayList<>();
        for (int i = 0; i < aln.getHeight(); i++) {
            EnumSeq<Enumerable> seq = aln.getEnumSeq(i);
            EnumSeq.Gappy<Enumerable> gap_copy = new EnumSeq.Gappy<>(alpha);
            gap_copy.set(new Character[1]); // add an empty column
            gap_copy.setName(seq.getName());
            seqArray.add(gap_copy);
        }

        EnumSeq.Alignment<Enumerable> gap_aln = new EnumSeq.Alignment<>(seqArray);
        //  Normalisation: log P★ - log(1 - P(col_gap))
        double gapColumnProb = logProbColGivenRate(tree, gap_aln, 1.0, model, 0, geometricSeqLenParam);
        double unobservedCols = MathEx.logm1exp(gapColumnProb);

        double normalisationTerm = probExtraCol - unobservedCols;

        return normalisationTerm + LL;
    }

    /**
     * This is used for normalisation of the overall column likelihood.
     * "the factorization in columns of the unconditional length
     * alignment distribution leaves some normalisation terms that we
     * gather together into what we think of as an “extra column” (★)
     * contribution. Thus, when calculating the total probability of a
     * multiple alignment as the product of l individual columns,
     * there is an additional term in the equation."
     * source: <a href="https://doi.org/10.1371/journal.pcbi.1000172">Rivas & Eddy, 2008</a>
     *
     * @return normalisation term for the alignment
     */
    public static double probExtraCol(Tree tree, GapSubstModel model, double geometricSeqLenParam) {


        int nNodes = tree.getNLeaves() + tree.getNParents();

        // First calculate the likelihood of an all gap column
        double[] pStar = new double[nNodes];
        Arrays.fill(pStar,  Double.NEGATIVE_INFINITY);
        for (int bpidx = nNodes - 1; bpidx >= 0; bpidx--) {
            BranchPoint node = tree.getBranchPoint(bpidx);

            if (node.isLeaf()) {
                pStar[bpidx] = 0.0; // log(1)
            } else {
                // ancestor
                int[] childrenBpindices = tree.getChildren(bpidx);
                double childrenLL = 0.0;
                for (int childBpidx: childrenBpindices) {
                    childrenLL += pStar[childBpidx];
                    // add probability of gap remaining
                    childrenLL += Math.log(1 - model.ksiT(tree.getDistance(childBpidx)));
                }
                pStar[bpidx] = childrenLL;
            }
        }

        double colLikelihood = 0.0;
        colLikelihood += Math.log(1 - geometricSeqLenParam); // penalise for length of sequence
        colLikelihood += pStar[0]; // add the probability of root node

        return colLikelihood;
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


    public static class Peeler {
        private final Tree tree;
        private final EnumSeq.Alignment<Enumerable> aln;
        private final Double colIndelRate;
        private final GapSubstModel model;
        private final int colIdx;
        private final double geometricSeqLenParam;
        private double col_prob;


        public Peeler(Tree tree, EnumSeq.Alignment<Enumerable> aln, Double colIndelRate,
                       GapSubstModel model, int col_idx, double geometricSeqLenParam) {

            this.tree = tree;
            this.aln = aln;
            this.colIndelRate = colIndelRate;
            this.model = model;
            this.colIdx = col_idx;
            this.geometricSeqLenParam = geometricSeqLenParam;
        }

        public void decorate() {
            this.col_prob = logProbColGivenRate(tree, aln, colIndelRate, model, colIdx, geometricSeqLenParam);
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
