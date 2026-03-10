package asr;

import bn.ctmc.GapSubstModel;
import bn.ctmc.matrix.JCGap;
import bn.ctmc.matrix.JTTGap;
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
import dat.pog.POGraph;

public class Mip {

    public static double MIN_MU_LAMBDA_VALUE = 0;
    public static double MAX_MU_LAMBDA_VALUE = 0.5;
    private static final int GAP = 0;
    private static final int NON_GAP = 1;
    private static final int VIRTUAL_START = -1;
    private static final int DEFAULT_GAP_PENALTY = 2;
    private static final double[] GAP_PENALTIES = new double[] {8.0, 6.0, 4.0, 2.0};
    private final HashMap<Integer, Integer[]> extantBinarySeqs;
    private final POAGraph alnPog;
    private final int nPos;
    private final IdxTree tree;
    private final EnumSeq.Alignment<Enumerable> aln;
    private final int nThreads;
    private String solverName;
    private String SubstModelName;
    double[][] treeNeighbourAlphaPen;
    boolean useBranchLengths;
    private MPSolver solver;
    private final Set<Integer> nodesToSkip = new HashSet<>();
    private Set<Integer> nodesConnectedToPreviousNode = new HashSet<>();
    private final int[] nodeWeights;
    private HashMap<Integer, MPVariable[]> ancestorPositionVars;
    private MPObjective objective;
    private HashMap<EdgeKey, MPVariable> edges = new HashMap<>();

    public Mip (IdxTree tree, EnumSeq.Alignment<Enumerable> aln, String solverName, String substModelName,
                int nThreads, boolean useBranchLengths) {

        this.extantBinarySeqs = createBinarySeqMap(aln, tree);
        this.alnPog = new POAGraph(aln);
        this.nPos = aln.getWidth();
        this.tree = tree;
        this.aln = aln;
        this.nThreads = nThreads;
        this.solverName = solverName;
        this.SubstModelName = substModelName;
        this.useBranchLengths = useBranchLengths;
        this.nodeWeights = new int[this.nPos];
        Arrays.fill(this.nodeWeights, 1);
        this.identifyNodesToSkip();

    }

    private void identifyNodesToSkip() {
        int previousNode = 0;
        for (int nodeIdx = 0; nodeIdx < nPos - 1; nodeIdx++) {
            int[] forwardEdges = alnPog.getNodeIndices(nodeIdx, true);
            int[] nextNodeBackwardEdges = alnPog.getNodeIndices(nodeIdx + 1, false);
            if ((forwardEdges.length == 1)
                    && Arrays.stream(forwardEdges).min().orElseThrow() == (nodeIdx + 1)
                    && nextNodeBackwardEdges.length == 1) {

                this.nodesToSkip.add(nodeIdx + 1);
                this.nodeWeights[previousNode] += 1;
            } else {
                previousNode = nodeIdx + 1;
                int nextNode = nodeIdx + 1;
                for (int forwardEdge : forwardEdges) {
                    if (forwardEdge == nextNode) {
                        this.nodesConnectedToPreviousNode.add(nextNode);
                        break;
                    }
                }
            }
        }
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

    private static boolean checkAncestralPogEdge(POGraph ancestorPog, int from, int to) {
        try {
            return ancestorPog.getEdge(from, to) != null;
        } catch (Exception e) { // caused if no edge exists for this ancestor solution
            return false;
        }
    }

    public HashMap<Integer, Integer[]> runMPSolverIndelInference() {

        // Uncomment for logging
        //solver.enableOutput();

        Loader.loadNativeLibraries(); // link to Google-OR Tools
        this.solver = MPSolver.createSolver(solverName);
        if (this.solver == null) {
            asr.GRASP.usage(6, "Could not create MIP solver with " + solverName);
        }

        createAncestralPositionVariables();
        createEdgeVariables();
        addEdgeConstraints();
        objective = solver.objective();

        addPenaltyConstraints();
        objective.setMinimization();

        int actualThreadsUsed = this.nThreads;
        if (solverName.equalsIgnoreCase("SCIP")) {
            // SCIP creates concurrent solvers - appears to be a bug where other workers
            // are not terminated when a solution is found.
            solver.setNumThreads(1);
            actualThreadsUsed = 1;
        } else {
            solver.setNumThreads(this.nThreads);
        }

        if (GRASP.VERBOSE) {
            System.out.println(solverName + " model build complete. Starting solver with " + actualThreadsUsed + (actualThreadsUsed == 1 ? " thread..." : " threads..."));
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
            System.err.println("A feasible solution could not be identified.");
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
     * @param tree the phylogenetic tree
     * @param aln the alignment
     * @param substModelName the substitution model name
     * @return a matrix where the rows are alignment columns and the columns are the gap opening penalties for each
     * node in the tree. The array matches the order of nodes in the IdxTree.
     */
    public double[][] createTreeNeighbourAlphaPen(IdxTree tree, EnumSeq.Alignment<Enumerable> aln,
            String substModelName) {

        double[][] treeNeighbourAlphaPen = new double[this.aln.getWidth()][this.tree.getSize()];

        if (GRASP.DISTANCE_BASED_MIP) {

            double geometric_seq_len_param = (double) 1 / aln.getAvgSeqLength();
            if (GRASP.VERBOSE) {
                System.out.println("Optimising indel parameters for distance-based MIP...");
            }
            double optimal_mu = IndelDist.optimiseMuLambda(MIN_MU_LAMBDA_VALUE, MAX_MU_LAMBDA_VALUE, substModelName,
                    tree, geometric_seq_len_param, aln);

            GapSubstModel gapModel;
            switch (substModelName) {
                case "JTT" -> gapModel = new JTTGap(optimal_mu, optimal_mu);
                case "JC" -> gapModel = new JCGap(optimal_mu, optimal_mu);
                default -> throw new IllegalArgumentException("Unrecognized gap substitution model: " + substModelName);
            }

            if (GRASP.VERBOSE) {
                System.out.println("Computing column priors under different rate categories...");
            }

            double[] rates = IndelDist.MEAN_RATES.get(GRASP.INDEL_RATE);
            Double[][] columnPriors = IndelDist.computeColumnPriors(gapModel, (Tree) tree, aln,
                    rates, geometric_seq_len_param);

            if (GRASP.VERBOSE) {
                System.out.println("Computing prefix sums for segment assignment...");
            }
            double[][] prefix_sums = IndelDist.computePrefixSums(columnPriors);

            if (GRASP.VERBOSE) {
                System.out.println("Assigning optimal rate segments...");
            }

            int[][] segments = IndelDist.assignSegments(columnPriors.length, IndelDist.RATE_PRIORS,
                    prefix_sums, IndelDist.RHO);

            int[] columnRateCategories = IndelDist.expandSegmentOrder(segments);

            double[][] rateAdjustedDists = new double[rates.length][tree.getSize()];

            double minDist = Double.POSITIVE_INFINITY;
            double maxDist = Double.NEGATIVE_INFINITY;
            // adjust each length by the assigned rate category
            for (int rateIdx = 0; rateIdx < rates.length; rateIdx++) {
                for (int bpidx = 0; bpidx < tree.getSize(); bpidx++) {
                    // adjust each length by the assigned rate category
                    double adjustedDist = rates[rateIdx] * tree.getDistance(bpidx);
                    rateAdjustedDists[rateIdx][bpidx] = adjustedDist;
                    if (adjustedDist < minDist) {
                        minDist = adjustedDist;
                    }
                    if (adjustedDist > maxDist) {
                        maxDist = adjustedDist;
                    }
                }
            }

            double binRange = (maxDist - minDist) / 2;

            for (int colIdx = 0; colIdx < aln.getWidth(); colIdx++) {
                int rateIdx = columnRateCategories[colIdx];
                for (int bpidx = 0; bpidx < tree.getSize(); bpidx++) {
                    int assignedBin = Math.min(GAP_PENALTIES.length - 1, (int) ((rateAdjustedDists[rateIdx][bpidx] + (binRange / 2) - minDist) / binRange));
                    treeNeighbourAlphaPen[colIdx][bpidx] = GAP_PENALTIES[assignedBin];
                }
            }

        } else {
            for (int colIdx = 0; colIdx < aln.getWidth(); colIdx++) {
                Arrays.fill(treeNeighbourAlphaPen[colIdx], DEFAULT_GAP_PENALTY);
            }
        }

        return treeNeighbourAlphaPen;
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

    private void createEdgeVariables() {


        int[] startIndices = alnPog.getStarts();
        int[] endIndices = alnPog.getEnds();
        for (int node: startIndices) {
            // don't need to perform check as will only remove if present
            this.nodesConnectedToPreviousNode.remove(node + 1);
        }
        for (int node: endIndices) {
            this.nodesConnectedToPreviousNode.remove(node);
        }

        for (int ancestorIdx : tree.getAncestors()) {

            // VIRTUAL STARTS TO REAL STARTS //
            for (int edgeEnd : startIndices) {
                // Variable for each edge from virtual start to a "real" start node
                MPVariable edge = solver.makeBoolVar("");// makeIntVar(0, 1, "");
                EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, edgeEnd, ancestorIdx);
                edges.put(edgeKey, edge);
            }

            for (int nodeIdx = 0; nodeIdx < alnPog.maxsize(); nodeIdx++) {

                // first go through and find nodes where we actually need edge variables
                if (!this.nodesToSkip.contains(nodeIdx + 1)) { // next node is not skipped

                    int[] forwardEdges = alnPog.getNodeIndices(nodeIdx, true);
                    if (!this.nodesConnectedToPreviousNode.contains(nodeIdx + 1)) { // current node not connected to adjacent node
                        if (forwardEdges.length == 1) {
                            // node has one outgoing edge, outgoing edge is implicit on this node being activated
                            EdgeKey edgeKey = new EdgeKey(nodeIdx, forwardEdges[0], ancestorIdx);
                            this.edges.put(edgeKey, this.ancestorPositionVars.get(ancestorIdx)[nodeIdx]);
                        } else {

                            for (int posTo : forwardEdges) {
                                EdgeKey edgeKey = new EdgeKey(nodeIdx, posTo, ancestorIdx);
                                MPVariable edge;
                                int[] backwardEdges = alnPog.getNodeIndices(posTo, false);
                                if (backwardEdges.length == 1) {
                                    edge = this.ancestorPositionVars.get(ancestorIdx)[posTo];
                                } else {
                                    edge = solver.makeBoolVar("");
                                }
                                this.edges.put(edgeKey, edge);
                            }
                        }
                    } else { // connected to the adjacent node
                        for (int posTo : forwardEdges) {
                            if (posTo == (nodeIdx + 1)) {
                                continue;
                            }
                            int[] backwardEdges = alnPog.getNodeIndices(posTo, false);
                            if (backwardEdges.length != 1) {
                                EdgeKey edgeKey = new EdgeKey(nodeIdx, posTo, ancestorIdx);
                                MPVariable edge = solver.makeBoolVar("");
                                this.edges.put(edgeKey, edge);
                            }
                        }
                    }
                }
            }
        }
    }

    private void addEdgeConstraints() {
        for (int ancestorIdx : tree.getAncestors()) {
            for (int pos = 0; pos < nPos; pos++) {

                if (this.nodesToSkip.contains(pos)) {
                    continue;
                }

                // iPrime is the immediate upstream node to this position (NOT always pos - 1)
                int iPrimeIndex = Arrays.stream(this.alnPog.getNodeIndices(pos, false)).max().orElseThrow();
                MPVariable nodeStateIPrime = null;
                double constant = 0.0;
                if (iPrimeIndex != VIRTUAL_START) {
                    nodeStateIPrime = this.ancestorPositionVars.get(ancestorIdx)[iPrimeIndex];
                } else {
                    constant = 1.0;
                }

                MPConstraint constraint = solver.makeConstraint(constant, constant);

                List<MPVariable> edgesBypassingI = new LinkedList<>();
                for (int posTo : this.alnPog.getNodeIndices(iPrimeIndex, true)) {
                    if (posTo > pos) {
                        MPVariable edge;
                        if (this.alnPog.getNodeIndices(posTo, true).length == 1) {
                            edge = this.ancestorPositionVars.get(ancestorIdx)[posTo];
                        } else {
                            edge = this.edges.get(new EdgeKey(iPrimeIndex, posTo, ancestorIdx));
                        }

                        edgesBypassingI.add(edge);
                    }
                }

                List<MPVariable> edgesBypassingIPrime = new LinkedList<>();
                for (int posFrom : alnPog.getNodeIndices(pos, false)) {
                    if (posFrom < iPrimeIndex) {
                        MPVariable edge = this.edges.get(new EdgeKey(posFrom, pos, ancestorIdx));
                        edgesBypassingIPrime.add(edge);
                    }
                }

                // Constraint: a_i = a_i` - sum(edges_bypassing_i) + sum(edges_bypassing_i`)
                // Constraint: a_i - a_i` + sum(edges_bypassing_i) - sum(edges_bypassing_i`) = 0
                constraint.setCoefficient(this.ancestorPositionVars.get(ancestorIdx)[pos], 1.0);
                if (nodeStateIPrime != null) {
                    constraint.setCoefficient(nodeStateIPrime, -1.0);
                }

                addConstraintSum(constraint, edgesBypassingI, 1.0);
                addConstraintSum(constraint, edgesBypassingIPrime, -1.0);
            }
        }
    }

    private void addConstraintSum(MPConstraint constraint, List<MPVariable> vars, double coefficient) {
        for (MPVariable var: vars) {
            constraint.setCoefficient(var, coefficient);
        }
    }

    private void addPenaltyConstraints() {

        for (int ancestralIdx : tree.getAncestors()) {

            MPVariable[] nodePosVar = ancestorPositionVars.get(ancestralIdx);


            for (int childIdx : tree.getChildren(ancestralIdx)) {

                MPVariable[] pen = new MPVariable[nPos];
                for (int i = 0; i < nPos; i++) {
                    pen[i] = solver.makeBoolVar("");
                }

                Integer[] nodeNeighbourPosVar = null;
                MPVariable[] nodeNeighborPosVarAncestor = null;
                MPVariable[] diffPos = null;

                if (tree.isLeaf(childIdx)) {
                    nodeNeighbourPosVar = extantBinarySeqs.get(childIdx);
                } else {
                    nodeNeighborPosVarAncestor = ancestorPositionVars.get(childIdx);

                    diffPos = new MPVariable[nPos];
                    for (int i = 0; i < nPos; i++) {
                        diffPos[i] = solver.makeBoolVar("");// 0, 1, "");
                    }
                }

                MPVariable diffVar = null;
                MPVariable prevDiffVar = null;
                for (int pos = 0; pos < nPos; pos++) {

                    if (pos > 0) {
                        prevDiffVar = diffVar;
                    }

                    if (tree.isLeaf(childIdx)) {

                        // constraint - difference variables
                        diffVar = solver.makeBoolVar("");
                        if (nodeNeighbourPosVar[pos] == 1) {
                            // constraint: diff == 1 - nodePosVar[pos]
                            // rearrange: diff + nodePosVar[pos] == 1
                            MPConstraint c = solver.makeConstraint(1, 1);
                            c.setCoefficient(diffVar, 1.0);
                            c.setCoefficient(nodePosVar[pos], 1.0);
                        } else {
                            //diffVar = nodePosVar[pos];
                            // diffVar - nodePosVar[pos] == 0
                            MPConstraint c = solver.makeConstraint(0, 0);
                            c.setCoefficient(diffVar, 1);
                            c.setCoefficient(nodePosVar[pos], -1);

                        }
                        objective.setCoefficient(diffVar, 1.0);

                        // penalty constraints
                        if (pos == 0) {
                            // pen[0] == diff[0]
                            // Rewrite as: pen[1] - diff[1] == 0
                            MPConstraint c = solver.makeConstraint(0, 0);
                            c.setCoefficient(pen[pos], 1.0);
                            c.setCoefficient(diffVar, -1.0);
                        } else {

                            // constraint: pen[pos]>= diff[(node, node_neighbor_item, pos)] - diff[(node, node_neighbor_item, pos - 1)])
                            // rearrange: pen[pos] -  diff[(node, node_neighbor_item, pos)] + diff[(node, node_neighbor_item, pos - 1)]) >= 0
                            MPConstraint c1 = solver.makeConstraint(0, Double.POSITIVE_INFINITY);
                            c1.setCoefficient(pen[pos], 1.0);

                            c1.setCoefficient(diffVar, -1.0);

                            c1.setCoefficient(prevDiffVar, 1.0);

                            // Constraint: pen[pos] >= node_neighbor_pos_var[pos - 1] + (1 - node_pos_var[pos - 1]) + (1 - node_neighbor_pos_var[pos]) + node_pos_var[pos] - 3
                            // because these neighbours are children, we can just evaluate everything on right hand side
                            // rearrange: pen[pos] + node_pos_var[pos - 1] - node_pos_var[pos] >= node_neighbor_pos_var[pos - 1] + 1 + (1 - node_neighbor_pos_var[pos]) - 3
                            // rearrange: pen[pos] + node_pos_var[pos - 1] - node_pos_var[pos] >= node_neighbor_pos_var[pos - 1] - node_neighbor_pos_var[pos] - 1
                            double rhs2 = nodeNeighbourPosVar[pos - 1] - nodeNeighbourPosVar[pos] - 1;
                            MPConstraint c2 = solver.makeConstraint(rhs2, Double.POSITIVE_INFINITY);
                            c2.setCoefficient(pen[pos], 1.0);
                            c2.setCoefficient(nodePosVar[pos - 1], 1.0);
                            c2.setCoefficient(nodePosVar[pos], -1.0);

                            // Constraint: pen[pos] >= node_neighbor_pos_var[pos] + (1 - node_pos_var[pos]) + (1 - node_neighbor_pos_var[pos - 1]) + node_pos_var[pos - 1] - 3
                            // Rearrange: pen[pos] + node_pos_var[pos] - node_pos_var[pos - 1]  >= node_neighbor_pos_var[pos] + 1 + (1 - node_neighbor_pos_var[pos - 1]) - 3
                            // Rearrange: pen[pos] + node_pos_var[pos] - node_pos_var[pos - 1]  >= node_neighbor_pos_var[pos] - node_neighbor_pos_var[pos - 1]) - 1
                            double rhs3 = nodeNeighbourPosVar[pos] - nodeNeighbourPosVar[pos - 1] - 1;
                            MPConstraint c3 = solver.makeConstraint(rhs3, Double.POSITIVE_INFINITY);
                            c3.setCoefficient(pen[pos], 1.0);
                            c3.setCoefficient(nodePosVar[pos], 1.0);
                            c3.setCoefficient(nodePosVar[pos - 1], -1.0);
                        }


                        //MPVariable objDiffVar = ((MPVariable[])diff.get(new Triplet(ancestralIdx, childIdx, pos)))[0];
                        //objective.setCoefficient(objDiffVar, 1.0);
                        //objective.setCoefficient(pen[pos], DEFAULT_GAP_PENALTY);
                        objective.setCoefficient(pen[pos], treeNeighbourAlphaPen[pos][childIdx]);

                    } else {
                        //diff_pos[pos] <= node_pos_var[pos] + node_neighbor_pos_var[pos]
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

                        if (pos == 0) {
                            MPConstraint c5 = solver.makeConstraint(0, 0);
                            c5.setCoefficient(pen[pos], 1.0);
                            c5.setCoefficient(diffPos[pos], -1.0);
                        } else {
                            MPConstraint c6 = solver.makeConstraint(0, Double.POSITIVE_INFINITY);
                            c6.setCoefficient(pen[pos], 1.0);
                            c6.setCoefficient(diffPos[pos], -1.0);
                            c6.setCoefficient(diffPos[pos - 1], 1.0);

                            MPConstraint c7 = solver.makeConstraint(-1.0, Double.POSITIVE_INFINITY);
                            c7.setCoefficient(pen[pos], 1.0);
                            c7.setCoefficient(nodeNeighborPosVarAncestor[pos - 1], -1.0);
                            c7.setCoefficient(nodePosVar[pos - 1], 1.0);
                            c7.setCoefficient(nodeNeighborPosVarAncestor[pos], 1.0);
                            c7.setCoefficient(nodePosVar[pos], -1.0);

                            MPConstraint c8 = solver.makeConstraint(-1.0, Double.POSITIVE_INFINITY);
                            c8.setCoefficient(pen[pos], 1.0);
                            c8.setCoefficient(nodeNeighborPosVarAncestor[pos], -1.0);
                            c8.setCoefficient(nodePosVar[pos], 1.0);
                            c8.setCoefficient(nodeNeighborPosVarAncestor[pos - 1], 1.0);
                            c8.setCoefficient(nodePosVar[pos - 1], -1.0);
                        }

                        objective.setCoefficient(diffPos[pos], 1.0);
                        //objective.setCoefficient(pen[pos], DEFAULT_GAP_PENALTY);
                        objective.setCoefficient(pen[pos], treeNeighbourAlphaPen[pos][childIdx]);
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

}
