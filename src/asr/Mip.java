package asr;

import bn.ctmc.GapSubstModel;
import bn.ctmc.matrix.*;
import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.pog.POAGraph;
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import com.google.ortools.linearsolver.MPVariable;
import com.google.ortools.linearsolver.MPSolver;
import dat.pog.POGTree;

public class Mip {

    private static final double MAX_PENALTY = 10.0;
    public static double MIN_MU_LAMBDA_VALUE = 0;
    public static double MAX_MU_LAMBDA_VALUE = 0.5;
    private static final int GAP = 0;
    private static final int NON_GAP = 1;
    private static final int VIRTUAL_START = -1;
    private static final int DEFAULT_GAP_PENALTY = 2;
    private static final double[] GAP_PENALTIES = new double[]{8.0, 6.0, 4.0, 2.0};
    private final HashMap<Integer, Integer[]> extantBinarySeqs;
    private final POAGraph alnPog;
    private final POGTree pogTree;
    private final int nPos;
    private final IdxTree tree;
    private final EnumSeq.Alignment<Enumerable> aln;
    private final int nThreads;
    private final String solverName;
    private final String substModelName;
    double[][] treeNeighbourAlphaPen;
    boolean useBranchLengths;
    private MPSolver solver;
    private final Set<Integer> nodesToSkip = new HashSet<>();
    private final Set<Integer> nodesConnectedToPreviousNode = new HashSet<>();
    private final int[] nodeWeights;
    private HashMap<Integer, MPVariable[]> ancestorPositionVars;
    private MPObjective objective;
    private final HashMap<EdgeKey, MPVariable> edges = new HashMap<>();
    private final HashMap<DiffKey, MPVariable> diff = new HashMap<>();
    private final int virtualEndIdx;

    public Mip(POGTree pogTree, EnumSeq.Alignment<Enumerable> aln, String solverName, String substModelName,
               int nThreads, boolean useBranchLengths) {

        this.pogTree = pogTree;
        this.tree = pogTree.getTree();
        this.extantBinarySeqs = createBinarySeqMap(aln, tree);
        this.alnPog = new POAGraph(aln);
        this.nPos =  aln.getWidth();
        this.aln = aln;
        this.nThreads = nThreads;
        this.solverName = solverName;
        this.substModelName = substModelName;
        this.useBranchLengths = useBranchLengths;
        this.nodeWeights = new int[nPos];
        this.virtualEndIdx = nPos;
        Arrays.fill(this.nodeWeights, 1);
        this.identifyNodesToSkip();
        this.treeNeighbourAlphaPen = createTreeNeighbourAlphaPen();

    }

    private void identifyNodesToSkip() {
        int previousNode = 0;
        for (int nodeIdx = 0; nodeIdx < nPos - 1; nodeIdx++) {


            if ((getForwardEdges(nodeIdx).length == 1)
                    && Arrays.stream(getForwardEdges(nodeIdx)).min().orElseThrow() == (nodeIdx + 1)
                    && getBackwardEdges(nodeIdx + 1).length == 1) {

                this.nodesToSkip.add(nodeIdx + 1);
                this.nodeWeights[previousNode] += 1;
            } else {
                previousNode = nodeIdx + 1;
                boolean nodeFound = false;
                for (int forwardEdge : getForwardEdges(nodeIdx)) {
                    if (nodeIdx + 1 == forwardEdge) {
                        nodeFound = true;
                        break;
                    }
                }
                if (nodeFound)  {
                    this.nodesConnectedToPreviousNode.add(nodeIdx + 1);
                }
            }
        }

        //System.out.println("Skipping " + nodesToSkip.toArray().length + " " + nodesConnectedToPreviousNode.toArray().length);
    }

    private static void outputAncestralSolutions(Tree tree,
                                                 HashMap<Integer, Integer[]> ancestorPositionVars,
                                                 String OUTPUT,
                                                 String PREFIX) {


        // Create output file path
        String fastaFilePath = OUTPUT + "/" + PREFIX + "_ancestral_indel.fasta";

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fastaFilePath))) {
            for (int ancestralIdx : tree.getAncestors()) {
                Integer[] nodePosVar = ancestorPositionVars.get(ancestralIdx);

                // Write FASTA header
                writer.write(">N" + tree.getLabel(ancestralIdx));
                writer.newLine();

                StringBuilder sequenceLine = new StringBuilder();
                for (int pos = 0; pos < nodePosVar.length; pos++) {
                    sequenceLine.append(nodePosVar[pos]);

                    // Optional: Add line breaks every 80 characters (standard FASTA format)
                    if ((pos + 1) % 80 == 0 && pos < nodePosVar.length - 1) {
                        writer.write(sequenceLine.toString());
                        writer.newLine();
                        sequenceLine = new StringBuilder();
                    }
                }

                // Write any remaining sequence
                if (!sequenceLine.isEmpty()) {
                    writer.write(sequenceLine.toString());
                    writer.newLine();
                }
            }

            System.out.println("Ancestral indel solutions saved to: " + fastaFilePath);

        } catch (IOException e) {
            System.err.println("Error writing ancestral sequences to file: " + e.getMessage());
            System.out.println("Sending to standard output instead.");
            printAncestralIndels(tree, ancestorPositionVars);


        }
    }

    public static void printAncestralIndels(IdxTree tree, HashMap<Integer, Integer[]> ancestorPositionVars) {

        for (int ancestralIdx : tree.getAncestors()) {
            Integer[] nodePosVar = ancestorPositionVars.get(ancestralIdx);
            StringBuilder sb = new StringBuilder();
            sb.append(">N").append(tree.getLabel(ancestralIdx)).append("\n");
            for (Integer integer : nodePosVar) {
                sb.append(integer);
            }
            System.out.println(sb);
        }

    }

    public HashMap<Integer, Integer[]> runMPSolverIndelInference() {

        Loader.loadNativeLibraries(); // link to Google-OR Tools
        this.solver = MPSolver.createSolver(solverName);

//        if (solverName.equalsIgnoreCase("SCIP")) {
//            this.solver = new MPSolver("SCIP", MPSolver.OptimizationProblemType.SCIP_MIXED_INTEGER_PROGRAMMING);
//        } else if (solverName.equalsIgnoreCase("Gurobi")) {
//            this.solver = MPSolver.createSolver(solverName);
//        }


        if (this.solver == null) {
            GRASP.usage(6, "Could not create MIP solver with " + solverName);
        }


        objective = solver.objective();

        createAncestralPositionVariables();
        createEdgeVariables();
        addPenaltyConstraints();
        objective.setMinimization();

        if (GRASP.VERBOSE) {
            solver.enableOutput();
        }

        int actualThreadsUsed = this.nThreads;
        if (solverName.equalsIgnoreCase("SCIP")) {
            // SCIP creates concurrent solvers - appears to be a bug where other workers
            // are not terminated when a solution is found.
            System.out.println("SCIP solver detected - using single thread for solving.");
            solver.setSolverSpecificParametersAsString("lp/initalgorithm = d");
            solver.setNumThreads(1);
            actualThreadsUsed = 1;
        } else if (solverName.equalsIgnoreCase("Gurobi")) {
            solver.setSolverSpecificParametersAsString("Presolve=-1, FeasibilityTol=1e-06, IntFeasTol=1e-05, OptimalityTol=1e-06, Method=1, DegenMoves=0, Threads=" + this.nThreads);
        }

        if (GRASP.VERBOSE) {
            System.out.println(solverName + "(" + solver.solverVersion().strip() + ") model build complete. Starting solver with " + actualThreadsUsed + (actualThreadsUsed == 1 ? " thread..." : " threads..."));
        }

        solver.setTimeLimit((long) GRASP.MIP_SOLVER_TIME_LIMIT_MINUTES * 60 * 1000); // convert to milliseconds


        MPSolver.ResultStatus resultStatus = solver.solve();

        if (resultStatus == MPSolver.ResultStatus.OPTIMAL) {
            System.out.println("Optimal solution found");
            System.out.println("Objective value: " + objective.value());
            System.out.printf("Indel solution found in %d min, %d sec%n", TimeUnit.MILLISECONDS.toMinutes(solver.wallTime()),
                    TimeUnit.MILLISECONDS.toSeconds(solver.wallTime()) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(solver.wallTime())));

        } else if (resultStatus == MPSolver.ResultStatus.FEASIBLE) {

            System.out.println("Feasible solution found");
            System.out.println("Objective value: " + objective.value());
            System.out.printf("Indel solution found in %d min, %d sec%n", TimeUnit.MILLISECONDS.toMinutes(solver.wallTime()),
                    TimeUnit.MILLISECONDS.toSeconds(solver.wallTime()) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(solver.wallTime())));

        } else if (resultStatus == MPSolver.ResultStatus.INFEASIBLE) {
            System.err.println("Indel model is infeasible given the input alignment.");
            return null;
        } else if (resultStatus == MPSolver.ResultStatus.NOT_SOLVED) {
            System.err.println("No solution found or a solution could not be identified with the given time limit.");
            return null;
        }

        return extractSolutionMPSolver(ancestorPositionVars);
    }

    private static HashMap<Integer, Integer[]> extractSolutionMPSolver(HashMap<Integer, MPVariable[]> ancestorPositionVars) {

        HashMap<Integer, Integer[]> ancestralIndels = new HashMap<>();

        for (int ancestralIdx : ancestorPositionVars.keySet()) {
            MPVariable[] nodePosVar = ancestorPositionVars.get(ancestralIdx);
            Integer[] ancestralSeq = new Integer[nodePosVar.length];
            for (int pos = 0; pos < nodePosVar.length; pos++) {
                ancestralSeq[pos] = (int) Math.round(nodePosVar[pos].solutionValue());
            }
            ancestralIndels.put(ancestralIdx, ancestralSeq);
        }

        return ancestralIndels;
    }

    private static boolean allColsOccupied(EnumSeq.Alignment<Enumerable> aln) {

        for (int i = 0; i < aln.getWidth(); i++) {
            int numSeqsWithContent = aln.getOccupancy(i);
            if (numSeqsWithContent == 0) {
                return false;
            }
        }

        return true;
    }

    public static HashMap<Integer, Integer[]> createBinarySeqMap(EnumSeq.Alignment<Enumerable> aln, IdxTree tree) {

        HashMap<Integer, Integer[]> binarySeqMap = new HashMap<>();

        for (int i = 0; i < aln.getHeight(); i++) {

            EnumSeq<Enumerable> seq = aln.getEnumSeq(i);
            Integer[] binSeqValues = new Integer[aln.getWidth()]; // + VIRTUAL_NODES];

            for (int j = 0; j < aln.getWidth(); j++) {
                if (seq.get(j) == null) {
                    binSeqValues[j] = GAP;
                } else {
                    binSeqValues[j] = NON_GAP;
                }
            }

            binarySeqMap.put(tree.getIndex(seq.getName()), binSeqValues);
        }

        return binarySeqMap;
    }

    /**
     * Each column in the alignment is first assigned to a particular indel rate category. Then, the branch
     * distances are adjusted for each column based on the assigned rate category. Finally, the adjusted distances
     * are binned into discrete gap penalties ranging from 2 to 8 in increments of 2. If the distance based option
     * is not chosen, then the default gap penalty is assigned to all positions.
     *
     * @return a matrix where the rows are alignment columns and the columns are the gap opening penalties for each
     * node in the tree. The array matches the order of nodes in the IdxTree.
     */
    private double[][] createTreeNeighbourAlphaPen() {

        double[][] treeNeighbourAlphaPen = new double[this.aln.getWidth()][this.tree.getSize()];

        if (useBranchLengths) {

            double geometric_seq_len_param = (double) 1 / aln.getAvgSeqLength();
            if (GRASP.VERBOSE) {
                System.out.println("Optimising indel parameters for distance-based MIP...");
            }
            double optimal_mu = IndelPeeler.optimiseMuLambda(MIN_MU_LAMBDA_VALUE, MAX_MU_LAMBDA_VALUE, substModelName,
                    tree, geometric_seq_len_param, aln);

            GapSubstModel gapModel = createGapSubstModel(optimal_mu);



            double[] rates = IndelSegmentation.MEAN_RATES.get(GRASP.INDEL_RATE);
            if (GRASP.VERBOSE) {
                System.out.println("Computing column priors under different indel rate categories...");
                for (int i = 0; i < rates.length; i++) {
                    System.out.println("Rate category " + i + ": " + rates[i]);
                }
            }
            double[][] columnPriors = IndelPeeler.computeColumnPriors(pogTree, gapModel, geometric_seq_len_param, rates, GRASP.NTHREADS);

            if (GRASP.VERBOSE) {
                System.out.println("Computing prefix sums for indel segment assignment...");
            }
            double[][] prefix_sums = IndelSegmentation.computePrefixSums(columnPriors);

            if (GRASP.VERBOSE) {
                System.out.println("Assigning optimal indel rate segments...");
            }

            int[][] segments = IndelSegmentation.assignSegments(columnPriors.length, IndelSegmentation.RATE_PRIORS,
                    prefix_sums);

            int[] columnRateCategories = IndelSegmentation.expandSegmentOrder(segments);
            if (GRASP.VERBOSE) {
                for (int i = 0; i < columnRateCategories.length; i++) {
                    System.out.println("Column " + i + ": indel rate category: " + columnRateCategories[i]);
                }
            }
            double[][] rateAdjustedDists = new double[rates.length][tree.getSize()];

            // adjust each length by the assigned rate category
            for (int rateIdx = 0; rateIdx < rates.length; rateIdx++) {
                for (int bpidx = 0; bpidx < tree.getSize(); bpidx++) {
                    // adjust each length by the assigned rate category
                    double adjustedDist = rates[rateIdx] * tree.getDistance(bpidx);
                    rateAdjustedDists[rateIdx][bpidx] = adjustedDist;
                }
            }

            for (int colIdx = 0; colIdx < aln.getWidth(); colIdx++) {
                int rateIdx = columnRateCategories[colIdx];
                for (int bpidx = 0; bpidx < tree.getSize(); bpidx++) {

                    // shorter distances get higher penalties
                    double penalty = Math.min(1 / rateAdjustedDists[rateIdx][bpidx], MAX_PENALTY);
                    treeNeighbourAlphaPen[colIdx][bpidx] = penalty;
                }
            }

        } else {
            for (int colIdx = 0; colIdx < aln.getWidth(); colIdx++) {
                Arrays.fill(treeNeighbourAlphaPen[colIdx], DEFAULT_GAP_PENALTY);
            }
        }

        return treeNeighbourAlphaPen;
    }

    private GapSubstModel createGapSubstModel(double optimal_mu) {
        GapSubstModel gapModel;
        switch (substModelName) {
            case "JTT" -> gapModel = new JTTGap(optimal_mu, optimal_mu);
            case "JC" -> gapModel = new JCGap(optimal_mu, optimal_mu);
            case "LG" -> gapModel = new LGGap(optimal_mu, optimal_mu);
            case "Dayhoff" -> gapModel = new DayhoffGap(optimal_mu, optimal_mu);
            case "WAG" -> gapModel = new WAGGap(optimal_mu, optimal_mu);
            case "Yang" -> gapModel = new YangGap(optimal_mu, optimal_mu);
            default -> throw new IllegalArgumentException("Unrecognized gap substitution model: " + substModelName);
        }
        return gapModel;
    }

    /**
     * Create integer position variables for ancestor sequences bounded between 0 and 1. These represent the
     * actual sequence positions for each ancestor in the tree.
     */
    private void createAncestralPositionVariables() {

        HashMap<Integer, MPVariable[]> ancestorSeqVars = new HashMap<>();
        for (int ancestorIdx : tree.getAncestors()) {
            MPVariable[] seqVars = new MPVariable[nPos];
            for (int j = 0; j < nPos; j++) {
                if (this.nodesToSkip.contains(j)) {
                    seqVars[j] = seqVars[j - 1];
                } else {
                    seqVars[j] = solver.makeBoolVar("");
                }
            }
            ancestorSeqVars.put(ancestorIdx, seqVars);
        }

        this.ancestorPositionVars = ancestorSeqVars;
    }

    public static int[] append(int[] arr, int val) {
        int[] out = Arrays.copyOf(arr, arr.length + 1);
        out[arr.length] = val;
        return out;
    }

    private void createEdgeVariables() {


        int[] startIndices = alnPog.getStarts();
        this.nodesConnectedToPreviousNode.remove(VIRTUAL_START + 1);
        this.nodesConnectedToPreviousNode.remove(virtualEndIdx);


        for (int ancestorIdx : tree.getAncestors()) {

            // VIRTUAL STARTS TO REAL STARTS //
            for (int edgeEnd : startIndices) {
                // Variable for each edge from virtual start to a "real" start node
                MPVariable edge = solver.makeBoolVar("");
                EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, edgeEnd, ancestorIdx);
                edges.put(edgeKey, edge);
            }

            for (int nodeIdx = 0; nodeIdx < alnPog.maxsize(); nodeIdx++) {
                // first go through and find nodes where we actually need edge variables
                if (!this.nodesToSkip.contains(nodeIdx + 1)) { // next node is not skipped
                    if (!this.nodesConnectedToPreviousNode.contains(nodeIdx + 1)) { // current node not connected to adjacent node
                        if (getForwardEdges(nodeIdx).length == 1) {
                            for (int posTo : getForwardEdges(nodeIdx)) {
                                // node has one outgoing edge, outgoing edge is implicit on this node being activated
                                EdgeKey edgeKey = new EdgeKey(nodeIdx, posTo, ancestorIdx);
                                this.edges.put(edgeKey, this.ancestorPositionVars.get(ancestorIdx)[nodeIdx]);
                            }
                        } else {
                            for (int posTo : getForwardEdges(nodeIdx)) {
                                EdgeKey edgeKey = new EdgeKey(nodeIdx, posTo, ancestorIdx);
                                MPVariable edge;
                                if (getBackwardEdges(posTo).length == 1) {
                                    edge = this.ancestorPositionVars.get(ancestorIdx)[posTo];
                                } else {
                                    edge = solver.makeBoolVar("");
                                }
                                this.edges.put(edgeKey, edge);
                            }

                            addEdgeConstraint(ancestorIdx, nodeIdx);
                        }
                    } else { // connected to the adjacent node
                        List<MPVariable> all_edges_from_pos = new ArrayList<>();
                        List<MPVariable> all_edges_to_pos1 = new LinkedList<>();

                        for (int posTo : getForwardEdges(nodeIdx)) {
                            if (posTo == (nodeIdx + 1)) {
                                continue;
                            }
                            MPVariable edge;
                            if (getBackwardEdges(posTo).length == 1) {
                                edge = this.ancestorPositionVars.get(ancestorIdx)[posTo];
                            } else {
                                EdgeKey edgeKey = new EdgeKey(nodeIdx, posTo, ancestorIdx);
                                edge = solver.makeBoolVar("");
                                this.edges.put(edgeKey, edge);
                            }
                            all_edges_from_pos.add(edge);
                        }

                        for (int posFrom : getBackwardEdges(nodeIdx + 1)) {
                            if (posFrom == nodeIdx) {
                                continue;
                            }
                            MPVariable edge = this.edges.get(new EdgeKey(posFrom, nodeIdx + 1, ancestorIdx));
                            all_edges_to_pos1.add(edge);
                        }

                        MPConstraint constraint = solver.makeConstraint(0, 0);
                        MPVariable nodePlusOne = this.ancestorPositionVars.get(ancestorIdx)[nodeIdx + 1];
                        constraint.setCoefficient(nodePlusOne, constraint.getCoefficient(nodePlusOne) + 1);

                        MPVariable node = this.ancestorPositionVars.get(ancestorIdx)[nodeIdx];
                        constraint.setCoefficient(node, constraint.getCoefficient(node) + -1);
                        addConstraintSum(constraint, all_edges_from_pos, 1);
                        addConstraintSum(constraint, all_edges_to_pos1, -1);
                    }
                }

                if (!this.nodesToSkip.contains(nodeIdx) && !this.nodesConnectedToPreviousNode.contains(nodeIdx)) {
                    if (getBackwardEdges(nodeIdx).length > 1 ||
                            getForwardEdges(Arrays.stream(getBackwardEdges(nodeIdx)).min().orElseThrow()).length == 1) {
                        addEdgeConstraint(ancestorIdx, nodeIdx);
                    }
                }
            }
        }
    }

    private void addEdgeConstraint(int ancestorIdx, int nodeIdx) {
        int iPrimeIndex = Arrays.stream(getBackwardEdges(nodeIdx)).max().orElseThrow();
        MPVariable aKIPrime = null;
        double constant = 0;
        if (iPrimeIndex == VIRTUAL_START) {
            constant = 1;
        } else {
            aKIPrime = this.ancestorPositionVars.get(ancestorIdx)[iPrimeIndex];
        }

        List<MPVariable> edgesBypassingI = new LinkedList<>();
        for (int posTo : getForwardEdges(iPrimeIndex)) {
            if (posTo <= nodeIdx) {
                continue;
            }
            MPVariable edge;
            if (getBackwardEdges(posTo).length == 1) {
                edge = this.ancestorPositionVars.get(ancestorIdx)[posTo];
            } else {
                edge = this.edges.get(new EdgeKey(iPrimeIndex, posTo, ancestorIdx));
            }
            edgesBypassingI.add(edge);
        }

        List<MPVariable> edgesBypassingIPrime = new LinkedList<>();
        for (int posFrom : getBackwardEdges(nodeIdx)) {
            if (posFrom == iPrimeIndex) {
                continue;
            }
            MPVariable edge = this.edges.get(new EdgeKey(posFrom, nodeIdx, ancestorIdx));
            edgesBypassingIPrime.add(edge);
        }

        MPConstraint constraint = solver.makeConstraint(constant, constant);
        // nodes can also exist as edges, need to account for this in the constraint coefficients
        double ancExisting = constraint.getCoefficient(this.ancestorPositionVars.get(ancestorIdx)[nodeIdx]);
        constraint.setCoefficient(this.ancestorPositionVars.get(ancestorIdx)[nodeIdx], ancExisting + 1);
        if (aKIPrime != null) {
            double existing = constraint.getCoefficient(aKIPrime);
            constraint.setCoefficient(aKIPrime, existing + -1);
        }
        addConstraintSum(constraint, edgesBypassingI, 1);
        addConstraintSum(constraint, edgesBypassingIPrime, -1);
    }

    private int[] getForwardEdges(int nodeIdx) {

        int[] forwardEdges;
        if (nodeIdx == VIRTUAL_START) {
            forwardEdges = alnPog.getStarts();
        } else {
            forwardEdges = alnPog.getNodeIndices(nodeIdx, true);
        }

        if (alnPog.isEndNode(nodeIdx)) {
            forwardEdges = append(forwardEdges, virtualEndIdx);
        }

        return forwardEdges;

    }

    private int[] getBackwardEdges(int nodeIdx) {

        int[] backwardEdges;
        if (nodeIdx == virtualEndIdx) {
            backwardEdges = alnPog.getEnds();
        } else {
            backwardEdges = alnPog.getNodeIndices(nodeIdx, false);
        }

        if (alnPog.isStartNode(nodeIdx)) {
            backwardEdges = append(backwardEdges, VIRTUAL_START);
        }

        return backwardEdges;

    }


    private void addConstraintSum(MPConstraint constraint, List<MPVariable> vars, double coefficient) {
        for (MPVariable var : vars) {
            double existing = constraint.getCoefficient(var);
            constraint.setCoefficient(var, existing + coefficient);
        }
    }


    private void addPenaltyConstraints() {


        for (int ancestralIdx : tree.getAncestors()) {

            MPVariable[] nodePosVar = ancestorPositionVars.get(ancestralIdx);

            for (int childIdx : tree.getChildren(ancestralIdx)) {

                MPVariable[] pen = new MPVariable[nPos];

                Integer[] nodeNeighbourPosVar = null;
                MPVariable[] nodeNeighborPosVarAncestor = null;
                MPVariable[] diffPos = null;

                boolean isExtant = false;
                if (tree.isLeaf(childIdx)) {
                    nodeNeighbourPosVar = extantBinarySeqs.get(childIdx);
                    isExtant = true;
                } else {
                    nodeNeighborPosVarAncestor = ancestorPositionVars.get(childIdx);
                    MPVariable lastVar = null;

                    diffPos = new MPVariable[nPos];
                    for (int pos = 0; pos < nPos; pos++) {
                        if (this.nodesToSkip.contains(pos)) {
                            diffPos[pos] = lastVar;
                        } else {
                            diffPos[pos] = solver.makeBoolVar("");
                            lastVar = diffPos[pos];
                        }
                    }
                }

                for (int pos = 0; pos < nPos; pos++) {
                    if (isExtant) {

                        if (this.nodesToSkip.contains(pos)) {
                            this.diff.put(new DiffKey(ancestralIdx, childIdx, pos),
                                    this.diff.get(new DiffKey(ancestralIdx, childIdx, pos - 1)));
                            continue;
                        }

                        // constraint - difference variables

                        if (nodeNeighbourPosVar[pos] == 1) {
                            // constraint: diff == 1 - nodePosVar[pos]
                            // rearrange: diff + nodePosVar[pos] == 1
                            MPVariable diffVar = solver.makeBoolVar("");
                            MPConstraint c = solver.makeConstraint(1, 1);
                            c.setCoefficient(diffVar, 1.0);
                            c.setCoefficient(nodePosVar[pos], 1.0);
                            this.diff.put(new DiffKey(ancestralIdx, childIdx, pos), diffVar);

                        } else if (nodeNeighbourPosVar[pos] == 0) {
                            this.diff.put(new DiffKey(ancestralIdx, childIdx, pos), nodePosVar[pos]);
                        }

                        // gap penalty constraints
                        if (pos == 0) {
                            // pen[0] == diff[0]
                            pen[pos] = solver.makeBoolVar("");
                            MPConstraint c = solver.makeConstraint(0, 0);
                            c.setCoefficient(pen[pos], 1.0);
                            c.setCoefficient(this.diff.get(new DiffKey(ancestralIdx, childIdx, pos)), -1.0);
                        } else {

                            if (nodeNeighbourPosVar[pos - 1] == 1 && nodeNeighbourPosVar[pos] == 0) {
                                pen[pos] = nodePosVar[pos];

                            } else if (nodeNeighbourPosVar[pos - 1] == 0 && nodeNeighbourPosVar[pos] == 1) {
                                pen[pos] = solver.makeBoolVar("");
                                MPConstraint constraint = solver.makeConstraint(1,1);
                                constraint.setCoefficient(pen[pos], 1);
                                constraint.setCoefficient(nodePosVar[pos], 1);

                            } else if (nodeNeighbourPosVar[pos - 1] == 1 && nodeNeighbourPosVar[pos] == 1) {
                                pen[pos] = solver.makeBoolVar("");
                                MPConstraint constraint = solver.makeConstraint(0,Double.POSITIVE_INFINITY);
                                constraint.setCoefficient(pen[pos], 1.0);
                                constraint.setCoefficient(nodePosVar[pos - 1], -1.0);
                                constraint.setCoefficient(nodePosVar[pos], 1.0);

                            } else if (nodeNeighbourPosVar[pos - 1] == 0 && nodeNeighbourPosVar[pos] == 0) {
                                pen[pos] = solver.makeBoolVar("");
                                MPConstraint constraint = solver.makeConstraint(0, Double.POSITIVE_INFINITY);
                                constraint.setCoefficient(pen[pos], 1.0);
                                constraint.setCoefficient(nodePosVar[pos], -1.0);
                                constraint.setCoefficient(nodePosVar[pos - 1], 1.0);
                            }
                        }


                        objective.setCoefficient(pen[pos],  objective.getCoefficient(pen[pos]) + treeNeighbourAlphaPen[pos][childIdx]);
                        DiffKey diffKey = new DiffKey(ancestralIdx, childIdx, pos);
                        objective.setCoefficient(this.diff.get(diffKey),  objective.getCoefficient(this.diff.get(diffKey)) + this.nodeWeights[pos]);

                    } else {

                        if (this.nodesToSkip.contains(pos)) {
                            continue;
                        }

                        //Absolute difference constraints
                        // diff_pos[pos] <= node_pos_var[pos] + node_neighbor_pos_var[pos]
                        MPConstraint c1 = solver.makeConstraint(Double.NEGATIVE_INFINITY, 0);
                        c1.setCoefficient(diffPos[pos], 1.0);
                        c1.setCoefficient(nodePosVar[pos], -1.0);
                        c1.setCoefficient(nodeNeighborPosVarAncestor[pos], -1.0);

                        // diff_pos[pos] >= node_pos_var[pos] - node_neighbor_pos_var[pos]
                        MPConstraint c2 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
                        c2.setCoefficient(diffPos[pos], 1.0);
                        c2.setCoefficient(nodePosVar[pos], -1.0);
                        c2.setCoefficient(nodeNeighborPosVarAncestor[pos], 1.0);

                        // diff_pos[pos] >= node_neighbor_pos_var[pos] - node_pos_var[pos]
                        MPConstraint c3 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
                        c3.setCoefficient(diffPos[pos], 1.0);
                        c3.setCoefficient(nodeNeighborPosVarAncestor[pos], -1.0);
                        c3.setCoefficient(nodePosVar[pos], 1.0);

                        // diff_pos[pos] <= 2 - node_neighbor_pos_var[pos] - node_pos_var[pos]
                        MPConstraint c4 = solver.makeConstraint(Double.NEGATIVE_INFINITY, 2.0);
                        c4.setCoefficient(diffPos[pos], 1.0);
                        c4.setCoefficient(nodeNeighborPosVarAncestor[pos], 1.0);
                        c4.setCoefficient(nodePosVar[pos], 1.0);

                        c1.setIsLazy(true);
                        c4.setIsLazy(true);

                        pen[pos] = solver.makeBoolVar("");
                        if (pos == 0) {
                            MPConstraint c5 = solver.makeConstraint(0, 0);
                            c5.setCoefficient(pen[pos], 1.0);
                            c5.setCoefficient(diffPos[pos], -1.0);
                        } else {
                            List<MPVariable> parentEdgesBypassingI = new LinkedList<>();
                            List<MPVariable> childEdgesBypassingI = new LinkedList<>();
                            // get the set of edges that leave i-1 AND BYPASS i
                            for (int posTo : getForwardEdges(pos - 1)) {

                                if (posTo > pos) {
                                    MPVariable parentEdge;
                                    MPVariable childEdge;
                                    if (getBackwardEdges(posTo).length == 1) {
                                        parentEdge = nodePosVar[posTo];
                                        childEdge = nodeNeighborPosVarAncestor[posTo];
                                    } else {
                                        parentEdge = this.edges.get(new EdgeKey(pos - 1, posTo, ancestralIdx));
                                        childEdge = this.edges.get(new EdgeKey(pos - 1, posTo, childIdx));
                                    }
                                    parentEdgesBypassingI.add(parentEdge);
                                    childEdgesBypassingI.add(childEdge);
                                }
                            }

                            // get the set of edges that enter i and bypass i-1
                            List<MPVariable> parentEdgesBypassingIMinus1 = new LinkedList<>();
                            List<MPVariable> childEdgesBypassingIMinus1 = new LinkedList<>();

                            for (int posFrom : getBackwardEdges(pos)) {

                                if (posFrom < pos - 1) {

                                    MPVariable parentEdge;
                                    MPVariable childEdge;
                                    if (getBackwardEdges(pos).length == 1) {
                                        parentEdge = nodePosVar[pos];
                                        childEdge = nodeNeighborPosVarAncestor[pos];
                                    } else {
                                        parentEdge = this.edges.get(new EdgeKey(posFrom, pos, ancestralIdx));
                                        childEdge = this.edges.get(new EdgeKey(posFrom, pos, childIdx));
                                    }

                                    parentEdgesBypassingIMinus1.add(parentEdge);
                                    childEdgesBypassingIMinus1.add(childEdge);
                                }
                            }

                            MPConstraint parentConstraint = solver.makeConstraint(-1, Double.POSITIVE_INFINITY);
                            parentConstraint.setCoefficient(pen[pos], 1.0);
                            parentConstraint.setCoefficient(diffPos[pos], -1.0);
                            addConstraintSum(parentConstraint, parentEdgesBypassingI, -1);
                            addConstraintSum(parentConstraint, parentEdgesBypassingIMinus1, -1);

                            MPConstraint childConstraint = solver.makeConstraint(-1, Double.POSITIVE_INFINITY);
                            childConstraint.setCoefficient(pen[pos], 1.0);
                            childConstraint.setCoefficient(diffPos[pos], -1.0);
                            addConstraintSum(childConstraint, childEdgesBypassingI, -1);
                            addConstraintSum(childConstraint, childEdgesBypassingIMinus1, -1);

                        }

                        double existingPen = objective.getCoefficient(pen[pos]);
                        objective.setCoefficient(pen[pos], existingPen + treeNeighbourAlphaPen[pos][childIdx]);
                        double existingDiff = objective.getCoefficient(diffPos[pos]);
                        objective.setCoefficient(diffPos[pos], existingDiff + this.nodeWeights[pos]);
                    }
                }
            }
        }
    }

    /**
     * Create a unique key for an edge in the POG associated with a specific ancestor.
     * Allows this to be used as a key in hash maps.
     */
    record EdgeKey(int from, int to, int ancestorIdx) {}
    record DiffKey(int ancestorIdx, int childIdx, int pos) {}

}
