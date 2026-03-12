package asr;

import bn.ctmc.GapSubstModel;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JC;
import bn.ctmc.matrix.JTT;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.BranchPoint;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.pog.POGTree;
import dat.pog.POGraph;
import dat.pog.SymNode;
import smile.math.Function;
import smile.math.MathEx;
import smile.math.special.Minimise;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;

public class IndelPeeler {

    private final GapSubstModel model;
    private final POGTree pogTree;
    final private IdxTree tree;
    private final PhyloBN pbn;
    private  final int columnIdx;
    private final double[][] nodeResidueProbs;
    private final double[] nodeGapProbs;
    private double[] containsGap = null;
    private final int numResidues; // excluding gap
    private final Object[] alphabet;
    private final int alphabetSize; // including gap
    private final double geometricSeqLenParam;
    private double treeProb;

    public IndelPeeler(POGTree pogTree, SubstModel model, Double rate, int columnIdx, double geometricSeqLenParam) {
        this.pogTree = pogTree;
        this.tree = pogTree.getTree();
        this.model = (GapSubstModel) model;
        this.columnIdx = columnIdx;
        this.alphabet = model.getDomain().getValues();
        this.numResidues = model.getDomain().size() - 1;
        this.alphabetSize = model.getDomain().size();
        this.geometricSeqLenParam = geometricSeqLenParam;


        this.nodeResidueProbs = new double[tree.getSize()][model.getDomain().size() - 1]; // only residues, not gap
        this.nodeGapProbs = new double[tree.getSize()];
        // instantiate with negative infinity for subsequent log sum calculations
        for (int i = 0; i < tree.getSize(); i++) {
            Arrays.fill(nodeResidueProbs[i], Double.NEGATIVE_INFINITY);
            nodeGapProbs[i] = Double.NEGATIVE_INFINITY;
        }

        if (rate == null) {
            this.pbn = PhyloBN.create(tree, model);
        } else {
            this.pbn = PhyloBN.create(tree, model, rate);
        }

    }

    public IndelPeeler(POGTree pogTree, SubstModel model, int columnIdx, double geometricSeqLenParam) {
        this(pogTree, model, null, columnIdx, geometricSeqLenParam);
    }


    /**
     * Calculates the likelihood of observing each column (independently)
     * given a particular indel rate.
     *
     * @param pogTree the partial order graph tree representing the alignment and phylogeny
     * @param model the gap augmented substitution model
     * @param geometricSeqLenParam the geometric sequence length parameter
     * @param rates the indel rates to calculate the column likelihoods for
     * @param nThreads the number of threads to use for parallelisation
     *
     * @return matrix of shape (numRates, numColumns) where each entry is the log likelihood of observing that column
     * given the tree, model, geometric sequence length parameter and indel rate.
     */
    public static double[][] computeColumnPriors(POGTree pogTree, SubstModel model,
                                                 double geometricSeqLenParam, double[] rates, int nThreads) {

        int numRates = rates.length;
        int numCols = pogTree.getPositions();
        double[][] columnPriors = new double[numCols][numRates];

        IndelPeeler[] peelers = new IndelPeeler[numRates * numCols];
        for (int colIdx = 0; colIdx < numCols; ++colIdx) {
            for (int rateIdx = 0; rateIdx < numRates; ++rateIdx) {
                int idx = colIdx * numRates + rateIdx;
                peelers[idx] = new IndelPeeler(pogTree, model, rates[rateIdx], colIdx, geometricSeqLenParam);
            }
        }

        // get back the results
        double[] results = runPeelingJobs(peelers, nThreads);

        for (int colIdx = 0; colIdx < numCols; ++colIdx) {
            for (int rateIdx = 0; rateIdx < numRates; ++rateIdx) {
                int idx = colIdx * numRates + rateIdx;
                columnPriors[colIdx][rateIdx] = results[idx];

            }
        }

        return columnPriors;

    }

    public static double calcProbAlnGivenTree(POGTree pogTree, GapSubstModel model, EnumSeq.Alignment<Enumerable> aln,
                                              double geometricSeqLenParam, Enumerable alpha, int nThreads) {

        double logLikelihood = 0.0;
        int numCols = aln.getWidth();
        IndelPeeler[] peelers = createPeelingJobs(pogTree, model, geometricSeqLenParam, numCols, 1.0);
        double[] columnProbs = runPeelingJobs(peelers, nThreads);
        for (int colIdx = 0; colIdx < numCols; colIdx++) {
            logLikelihood += columnProbs[colIdx];
        }

        double probExtraCol = probExtraCol(pogTree.getTree(), model, geometricSeqLenParam); // add the normalisation term

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
        POGTree gapPogTree = new POGTree(gap_aln, pogTree.getTree());

        //  Normalisation: log P★ - log(1 - P(col_gap))
        IndelPeeler peeler = new IndelPeeler(gapPogTree, model, 0, geometricSeqLenParam);
        double gapColumnProb = peeler.decorate(); // get the probability of an all gap column
        double unobservedCols = MathEx.logm1exp(gapColumnProb);

        double normalisationTerm = probExtraCol - unobservedCols;

        return normalisationTerm + logLikelihood;

    }

    public static double probExtraCol(IdxTree tree, GapSubstModel model, double geometricSeqLenParam) {


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

    private static IndelPeeler[] createPeelingJobs(POGTree pogTree, SubstModel model,
                                                   double geometricSeqLenParam, int numCols, double rate) {

        IndelPeeler[] peelers = new IndelPeeler[numCols];
        for (int colIdx = 0; colIdx < numCols; ++colIdx) {
            peelers[colIdx] = new IndelPeeler(pogTree, model, rate, colIdx, geometricSeqLenParam);
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
     * @return
     */
    public double decorate() {

        felsensteinsExtendedPeeling();
        int ROOT_INDEX = 0;
        double rootGapProb = nodeGapProbs[ROOT_INDEX];

        // sum over possible residue assignments
        double[] residueTerms = new double[numResidues]; // number of alphabet letters
        for (int resIdx = 0; resIdx < numResidues; resIdx++) {
            double priorProb = Math.log(model.getProb(alphabet[resIdx]));
            double rootResidueProb = nodeResidueProbs[ROOT_INDEX][resIdx];
            residueTerms[resIdx] = priorProb + rootResidueProb; // weight by prior prob of residue
        }

        double weightedLogSumResidueProb = MathEx.logsumexp(residueTerms) + Math.log(geometricSeqLenParam);
        double[] finalColTerms = {rootGapProb, weightedLogSumResidueProb};

        double colProb = MathEx.logsumexp(finalColTerms);
        treeProb = colProb;
        return colProb;

    }

    public double getDecoration() {
        return treeProb;
    }


    private void felsensteinsExtendedPeeling() {

        containsGap = containsGap();

        // iterate through branch point indices backwards for postorder traversal
        for (int bpidx = tree.getSize() - 1; bpidx >= 0; bpidx--) {

            if (tree.isLeaf(bpidx)) {
                calcLeafPeelingProbabilities(bpidx);
            } else {
                calcAncestralPeelingProbabilities(bpidx);
            }
        }
    }

    private void calcLeafPeelingProbabilities(int bpidx) {
        POGraph extantPog = pogTree.getExtant(bpidx);
        SymNode n = (SymNode) extantPog.getNode(columnIdx);

        if (n == null) {
            nodeGapProbs[bpidx] = Double.NEGATIVE_INFINITY;
        } else {
            int resIdx = model.getDomain().getIndex(n.getValue());
            nodeResidueProbs[bpidx][resIdx] = 0.0;
        }

    }

    private void calcAncestralPeelingProbabilities(int bpidx) {

        nodeGapProbs[bpidx] = MathEx.logsumexp(calcLogProbChildrenGivenAncestralGap(bpidx));

        for (int parentResIdx = 0; parentResIdx < numResidues; parentResIdx++) {
            nodeResidueProbs[bpidx][parentResIdx] = calcLogProbChildrenGivenAncestralResidue(bpidx, parentResIdx);
        }

    }

    private double[] calcLogProbChildrenGivenAncestralGap(int bpidx) {

        int[] childrenBpIndices = tree.getChildren(bpidx);
        double[] allChildrenLogGapTerms = new double[childrenBpIndices.length];

        int current_child = 0;
        for (int childBpidx : childrenBpIndices) {
            // collect all possible terms we will marginalise over
            if (this.containsGap[childBpidx] == 1.0) {
                // need to enumerate over every residue sum(Pu(L_child,q).P(q|-,t) and
                // then add Pu(L_child, -) so full alphabet_size
                double[] childLogProbTerms = new double[alphabetSize];
                //1.0 sum(Pu(L_child,q).P(q|-,t)
                for (int resIdx = 0; resIdx < numResidues; resIdx++) {
                    double nodeResidueProb = nodeResidueProbs[childBpidx][resIdx];

                    Object residue = model.getDomain().get(resIdx);
                    double insertionProb = Math.log(getProbOfInsertion(bpidx, residue));

                    childLogProbTerms[resIdx] = nodeResidueProb + insertionProb;
                }

                childLogProbTerms[numResidues] = nodeGapProbs[childBpidx];
                allChildrenLogGapTerms[current_child] = MathEx.logsumexp(childLogProbTerms);

            } else {
                allChildrenLogGapTerms[current_child] = Double.NEGATIVE_INFINITY;
            }
            current_child++;
        }

        return allChildrenLogGapTerms;
    }

    private double calcLogProbChildrenGivenAncestralResidue(int bpidx, int parentResIdx) {

        int[] childrenBpindices = tree.getChildren(bpidx);

        double[] childrenLogProbTerms = new double[childrenBpindices.length];
        int currentChildResIdx = 0;
        for (int childBpidx : childrenBpindices) {

            // only include full alphabet with gap if all children have gaps
            double[] childLogProbTerms = (containsGap[childBpidx] == 1.0) ? new double[alphabetSize] : new double[numResidues];

            for (int childResIdx = 0; childResIdx < numResidues; childResIdx++) {
                double childNodeResidueProb = nodeResidueProbs[childBpidx][childResIdx];


                double gapAugmentedResidueProb = Math.log(getProbGapAugmented(childBpidx, alphabet[childResIdx],
                        alphabet[parentResIdx]));
                childLogProbTerms[childResIdx] = childNodeResidueProb + gapAugmentedResidueProb;
            }

            if (containsGap[childBpidx] == 1.0) {
                double deletionProb = Math.log(getProbOfGap(childBpidx));
                childLogProbTerms[numResidues] = deletionProb;
            }

            double childLogProb = MathEx.logsumexp(childLogProbTerms);

            childrenLogProbTerms[currentChildResIdx] = childLogProb;

            currentChildResIdx++;
        }

        return MathEx.sum(childrenLogProbTerms);
    }


    public double[] containsGap() {

        double[] containsGap = new double[tree.getSize()];

        // If we want to do a postorder traversal, can just iterate
        // through branch point indices backwards as these are labelled depth-first.
        for (int bpidx = tree.getSize() - 1; bpidx >= 0; bpidx--) {

            if (tree.isLeaf(bpidx)) {
                POGraph extantPog = pogTree.getExtant(bpidx);
                SymNode n = (SymNode) extantPog.getNode(columnIdx);
                if (n == null) {
                    containsGap[bpidx] = 1.0;
                } else {
                    containsGap[bpidx] = 0.0;
                }
            } else {
                //ancestor
                int[] children = tree.getChildren(bpidx);
                double all_gaps = 1.0;
                // if any of the children don't contain a gap, break and record
                for (int child : children) {
                    if (containsGap[child] == 0.0) {
                        all_gaps = 0.0;
                        break;
                    }
                }
                containsGap[bpidx] = all_gaps;
            }
        }

        return containsGap;
    }

    private double getProbOfInsertion(int bpidx, Object state) {

        SubstNode substNode = (SubstNode) pbn.getBNode(bpidx);

        double insertionProb = ksiT(substNode.getTime());
        double stationaryFreqResidue = model.getProb(state);//  substNode.getProb(state);

        return insertionProb * stationaryFreqResidue;
    }

    private double ksiT(double time) {

        if (model.getMu() == 0 && model.getLambda() == 0) {
            return 0.0;
        }
        double insertionRate = model.getLambda() / (model.getMu() + model.getLambda());
        double indelProp = 1 - Math.exp(-((model.getMu() + model.getLambda()) * time));

        return insertionRate * indelProp;
    }

    private double getProbGapAugmented(int bpidx, Object childState, Object parentState) {

        SubstNode substNode = (SubstNode) pbn.getBNode(bpidx);
        double probNoInsert = 1 - ksiT(substNode.getTime());
        double conditinalProb = substNode.getProb(childState, parentState);

        return conditinalProb * probNoInsert;
    }

    private double getProbOfGap(int bpidx) {

        SubstNode substNode = (SubstNode) pbn.getBNode(bpidx);
        double probNoInsertion = 1 - ksiT(substNode.getTime());
        double probDeletion = gammaT(substNode.getTime());

        return probNoInsertion * probDeletion;

    }

    public double gammaT(double time) {

        if (model.getMu() == 0 && model.getLambda() == 0) {
            return 0.0;
        }

        double indelTotal = model.getMu() + model.getLambda();
        double deletionProp = model.getMu() / (indelTotal);
        double indelProp = 1 - Math.exp(-((indelTotal) * time));

        return deletionProp * indelProp;
    }


    /**
     * Optimises mu and lambda (insertion and deletion rates) assuming they are equal. Uses Brent's method to find
     * the optimal value that maximises the likelihood of the alignment given the tree. The likelihood is calculated
     * according to equation 29 in <a href="https://doi.org/10.1371/journal.pcbi.1000172"> Rivas & Eddy, 2008</a>
     * @param min_val smallest value to search
     * @param max_val largest value to search
     * @param substModelName index of the substitution model to use
     * @param tree phylogenetic tree
     * @param geometricSeqLenParam geometric sequence length param for the alignment
     * @param aln the alignment
     * @return the optimal mu and lambda value
     * @throws IllegalArgumentException if the model is not supported
     */
    public static double optimiseMuLambda(double min_val, double max_val, String substModelName, IdxTree tree,
                                          double geometricSeqLenParam,
                                          EnumSeq.Alignment<Enumerable> aln) throws IllegalArgumentException {

        double[] F;
        double[][] IRM;
        Enumerable alpha;
        if (substModelName.equals("JTT")) {
            F = JTT.F;
            IRM = JTT.Q;
            alpha = new Enumerable(JTT.S);
        } else if (substModelName.equals("JC")) {
            F = JC.F(JC.S.length);
            IRM = JC.Q(1, JC.S.length);
            alpha = new Enumerable(JC.S);
        } else {
            throw new IllegalArgumentException(substModelName + " not implemented yet");
        }

        IndelPeeler.AlnLikelihood alnLikelihood = new IndelPeeler.AlnLikelihood(tree, aln, F, IRM,
                alpha, geometricSeqLenParam);

        return Minimise.brent(alnLikelihood, min_val, max_val);
    }

    /**
     * Function to calculate the likelihood of an alignment given a tree and gap augmented substitution model.
     */
    public static class AlnLikelihood implements Function {

        double[] F;
        double[][] IRM;
        Enumerable alpha;
        IdxTree tree;
        EnumSeq.Alignment<Enumerable> aln;
        double geometricSeqLenParam;

        public AlnLikelihood(IdxTree tree, EnumSeq.Alignment<Enumerable> aln, double[] F, double[][] IRM,
                             Enumerable alpha, double geometricSeqLenParam) {
            this.tree = tree;
            this.aln = aln;
            this.geometricSeqLenParam = geometricSeqLenParam;
            this.F = F;
            this.IRM = IRM;
            this.alpha = alpha;
        }

        /**
         * Get the log likelihood of a particular alignment given the tree and gap augmented substitution model.
         * This is used by a minimisation routine to find optimal params for the substitution model.
         * @param muLambda Assumes mu (deletion rate) and lambda (insertion rate) are equal
         * @return likelihood of the alignment given the tree + mu + lambda.
         */
        @Override
        public double f(double muLambda) {

            GapSubstModel newModel = new GapSubstModel(this.F, this.IRM, this.alpha, muLambda, muLambda);

            // trying to maximise the log likelihood
            POGTree pogTree = new POGTree(aln, tree);
            return -1.0 * IndelPeeler.calcProbAlnGivenTree(pogTree, newModel, aln, geometricSeqLenParam, alpha, GRASP.NTHREADS);
        }
    }


}
