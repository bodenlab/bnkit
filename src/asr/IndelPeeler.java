package asr;

import bn.ctmc.*;
import bn.ctmc.matrix.*;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.BranchPoint;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.Tree;
import dat.pog.POGTree;
import dat.pog.POGraph;
import dat.pog.SymNode;
import smile.math.Function;
import smile.math.MathEx;
import smile.math.special.Minimise;

import java.util.*;
import java.util.concurrent.ExecutionException;

public class IndelPeeler {

    private final PIPSubstModel model;
    final private IdxTree tree;
    double treeProb;
    private final int columnIdx;
    EnumSeq.Alignment<Enumerable> aln;
    PhyloBN pbn;
    Set nonGappedLeaves;
    Set ancestorsToRootFromMRCA;
    Map<String, Integer> alnMap;

    public IndelPeeler(IdxTree tree, PIPSubstModel model, int columnIdx,
                       EnumSeq.Alignment<Enumerable> aln, PhyloBN pbn,
                       Set nonGappedLeaves, Set ancestorsToRootFromMRCA, Map<String, Integer> alnMap) {
        this.tree = tree;
        this.model = model;
        this.columnIdx = columnIdx;
        this.aln = aln;
        this.pbn = pbn;
        this.nonGappedLeaves = nonGappedLeaves;
        this.ancestorsToRootFromMRCA = ancestorsToRootFromMRCA;
        this.alnMap = alnMap;
    }

//
//    /**
//     * Calculates the likelihood of observing each column (independently)
//     * given a particular indel rate.
//     *
//     * @param pogTree              the partial order graph tree representing the alignment and phylogeny
//     * @param model                the gap augmented substitution model
//     * @param geometricSeqLenParam the geometric sequence length parameter
//     * @param rates                the indel rates to calculate the column likelihoods for
//     * @param nThreads             the number of threads to use for parallelisation
//     * @return matrix of shape (numRates, numColumns) where each entry is the log likelihood of observing that column
//     * given the tree, model, geometric sequence length parameter and indel rate.
//     */
    public static double[][] computeColumnPriors(IdxTree tree, PIPSubstModel model,
                                                 double[] rates, int nThreads, int numCols,
                                                 EnumSeq.Alignment<Enumerable> aln, Set[] nonGappedLeaveSets,
                                                 Set[] ancestorsToRootFromMRCA, Map<String, Integer> alnMap) {

        int numRates = rates.length;
        double[][] columnPriors = new double[numCols][numRates];

        IndelPeeler[] peelers = new IndelPeeler[numRates * numCols];
        for (int rateIdx = 0; rateIdx < numRates; ++rateIdx) {
            PhyloBN pbn = PhyloBN.create(tree, model, rates[rateIdx]);
            for (int colIdx = 0; colIdx < numCols; ++colIdx) {
                //int idx = colIdx * numRates + rateIdx;
                int idx = rateIdx * numCols + colIdx;
                peelers[idx] = new IndelPeeler(tree, model, colIdx, aln, pbn, nonGappedLeaveSets[colIdx],
                        ancestorsToRootFromMRCA[colIdx], alnMap);
            }
        }

        // get back the results
        double[] results = runPeelingJobs(peelers, nThreads);
        for (int rateIdx = 0; rateIdx < numRates; ++rateIdx) {
            for (int colIdx = 0; colIdx < numCols; ++colIdx) {
                //int idx = colIdx * numRates + rateIdx;
                int idx = rateIdx * numCols + colIdx;
                columnPriors[colIdx][rateIdx] = results[idx];

            }
        }

        return columnPriors;
    }

    public static double calcProbAlnGivenTree(IdxTree tree, PIPSubstModel model, EnumSeq.Alignment<Enumerable> aln,
                                              int nThreads, PhyloBN pbn, double treeLength,
                                              EnumSeq.Alignment<Enumerable> gapCol, Set[] nonGappedLeaveSets,
                                              Set[] ancestorsToRootFromMRCA, Map<String, Integer> alnMap) {


        double columnLogLikelihoods = 0.0;
        int numCols = aln.getWidth();

        IndelPeeler[] peelers = createPeelingJobs(tree, model, numCols, aln, pbn,
                nonGappedLeaveSets, ancestorsToRootFromMRCA, alnMap);
        double[] columnProbs = runPeelingJobs(peelers, nThreads);
        for (int colIdx = 0; colIdx < numCols; colIdx++) {
            columnLogLikelihoods += columnProbs[colIdx];
        }

        double logPc0 = IdxTree.getColumnProb(model, 0, (Tree) tree,
                gapCol, pbn, new HashSet<>(), new HashSet<>(), alnMap);

        double nu = model.getLambda() * (treeLength + (1.0 / model.getMu()));
        // log(|m|!) via log-sum
        double logFactorialM = 0.0;
        for (int i = 1; i <= numCols; i++) {
            logFactorialM += Math.log(i);
        }

        double logPhi =  -logFactorialM
                + numCols * Math.log(nu)
                + (nu * (Math.exp(logPc0) - 1.0));

        return logPhi + columnLogLikelihoods;

    }

    private static IndelPeeler[] createPeelingJobs(IdxTree tree, PIPSubstModel model,
                                                   int numCols, EnumSeq.Alignment<Enumerable> aln, PhyloBN pbn,
                                                   Set[] nonGappedLeaveSets, Set[] ancestorsToRootFromMRCA,
                                                   Map<String, Integer> alnMap) {

        IndelPeeler[] peelers = new IndelPeeler[numCols];
        for (int colIdx = 0; colIdx < numCols; ++colIdx) {
            peelers[colIdx] = new IndelPeeler(tree, model, colIdx, aln, pbn,
                    nonGappedLeaveSets[colIdx], ancestorsToRootFromMRCA[colIdx], alnMap);
        }

        return peelers;
    }

    private static double[] runPeelingJobs(IndelPeeler[] peelers, int nThreads) {

        double[] results = new double[peelers.length];
        ThreadedPeeler thread_pool = new ThreadedPeeler(peelers, nThreads);
        try {
            Map<Integer, Double> ret = thread_pool.runBatch();
            for (int col_idx = 0; col_idx < peelers.length; ++col_idx) {
                results[col_idx] = ret.get(col_idx);
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to run peeling jobs");
        }

        return results;
    }

    /**
     * logProbColGivenRate
     *
     * @return
     */
    public double decorate() {

        treeProb = IdxTree.getColumnProb(model, columnIdx, (Tree) tree, aln, pbn,
                nonGappedLeaves, ancestorsToRootFromMRCA, alnMap);
        return treeProb;
    }

    public double getDecoration() {
        return treeProb;
    }

    /**
     * Optimises mu and lambda (insertion and deletion rates) assuming they are equal. Uses Brent's method to find
     * the optimal value that maximises the likelihood of the alignment given the tree. The likelihood is calculated
     * according to equation 29 in <a href="https://doi.org/10.1371/journal.pcbi.1000172"> Rivas & Eddy, 2008</a>
     *
     * @param min_val              smallest value to search
     * @param max_val              largest value to search
     * @param substModelName       index of the substitution model to use
     * @param tree                 phylogenetic tree
     * @param aln                  the alignment
     * @return the optimal mu and lambda value respectively
     * @throws IllegalArgumentException if the model is not supported
     */
    public static double[] optimiseMuLambda(double min_val, double max_val, String substModelName, IdxTree tree,
                                          EnumSeq.Alignment<Enumerable> aln) throws IllegalArgumentException {




        Map<String, Integer> alnMap = aln.getMap();
        Set[] nonGappedLeaveSets = new Set[aln.getWidth()];
        Set[] ancestorsToRootFromMRCA = new Set[aln.getWidth()];
        for (int colIdx = 0; colIdx < aln.getWidth(); colIdx++) {
            Set<Integer> S = IdxTree.findNonGappedLeaves(tree, aln, colIdx, alnMap);
            nonGappedLeaveSets[colIdx] = S;
            int mrca = IdxTree.findMRCA(tree, S);
            Set<Integer> A = IdxTree.getAllAncestorsToRoot(mrca, tree);
            ancestorsToRootFromMRCA[colIdx] = A;
        }

        LikelihoodEvaluator evaluator = new LikelihoodEvaluator(tree, aln, substModelName, nonGappedLeaveSets,
                ancestorsToRootFromMRCA, alnMap);

        double mu = 0.1;      // initial values
        double lambda = 10.0;
        evaluator.setLambda(lambda);
        evaluator.setMu(mu);
        double prevLogLik = Double.NEGATIVE_INFINITY;
        double tol = 1e-5;
        int maxIter = 100;

        System.out.println("Optimising Lambda, fixed Mu");
        double bestLambda = Minimise.brent(lambda_ -> {
            evaluator.setLambda(lambda_);
            evaluator.setMu(lambda_);
            return -evaluator.evaluate();  // negate for minimisation
        }, min_val, max_val);
        evaluator.setLambda(bestLambda);
        evaluator.setMu(bestLambda);

//        for (int iter = 0; iter < maxIter; iter++) {
//
//            // optimise mu with lambda fixed
//            System.out.println("Optimising Mu, fixed Lambda");
//            double bestMu = Minimise.brent(mu_ -> {evaluator.setMu(mu_); return -evaluator.evaluate();}, min_val, max_val);
//            evaluator.setMu(bestMu);
//
//            // optimise lambda with mu fixed
//            System.out.println("Optimising Lambda, fixed Mu");
//            double bestLambda = Minimise.brent(lambda_ -> {
//                evaluator.setLambda(lambda_);
//                return -evaluator.evaluate();  // negate for minimisation
//            }, min_val, max_val);
//            evaluator.setLambda(bestLambda);
//
//            // check convergence
//            double logLik = -evaluator.evaluate();
//            System.out.println("Iter=" + iter
//                    + " mu=" + bestMu
//                    + " lambda=" + bestLambda
//                    + " logLik=" + -logLik);
//
//            if (Math.abs(logLik - prevLogLik) < tol) {
//                System.out.println("Converged at iteration " + iter);
//                break;
//            }
//            prevLogLik = logLik;
//        }

        return new double[]{evaluator.mu, evaluator.lambda};
    }

    public static class LikelihoodEvaluator {

        final IdxTree tree;
        PIPSubstModel model;
        String substModelName;
        EnumSeq.Alignment<Enumerable> aln;
        PhyloBN pbn;
        private double mu = 1.0;
        private double lambda = 1.0;
        private final double totalTreeLength;
        EnumSeq.Alignment<Enumerable> gapCol = null;
        Set[] nonGappedLeaveSets;
        Set[] ancestorsToRootFromMRCA;
        Map<String, Integer> alnMap;


        public LikelihoodEvaluator(
                IdxTree tree,
                EnumSeq.Alignment<Enumerable> aln,
                String substModelName,
                Set[] nonGappedLeaveSets,
                Set[] ancestorsToRootFromMRCA,
                Map<String, Integer> alnMap
        ) {
            this.tree = tree;
            this.totalTreeLength = Arrays.stream(tree.getValidDistances()).sum();
            this.aln = aln;
            this.substModelName = substModelName;
            this.nonGappedLeaveSets = nonGappedLeaveSets;
            this.ancestorsToRootFromMRCA = ancestorsToRootFromMRCA;
            this.alnMap = alnMap;
        }

        public void setMu(double mu) { this.mu = mu; }
        public void setLambda(double lambda) { this.lambda = lambda; }


        public double evaluate() {

            // trying to maximise the log likelihood
            switch (substModelName) {
                case "JC" -> model = new JCPIP(mu, lambda);
                case "JTT" -> model = new JTTPIP(mu, lambda);
                default -> throw new ASRRuntimeException("Model not supported");
            }

            if (gapCol == null) {
                this.gapCol = IdxTree.createGapColumn(model.getDomain(), aln);
            }

            this.pbn = PhyloBN.create(tree, model, 1.0);
            return IndelPeeler.calcProbAlnGivenTree(tree, model, aln, GRASP.NTHREADS, pbn, totalTreeLength, gapCol,
                    nonGappedLeaveSets, ancestorsToRootFromMRCA, alnMap);

        }
    }
}