package asr;

import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;
import dat.pog.POAGraph;
import java.io.*;
import java.util.HashMap;
import com.google.ortools.linearsolver.MPVariable;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.Loader;
import com.google.ortools.sat.CpModel;
import com.google.ortools.sat.CpSolver;
import com.google.ortools.sat.CpSolverStatus;
import com.google.ortools.sat.LinearExpr;
import com.google.ortools.sat.LinearExprBuilder;
import com.google.ortools.sat.Literal;

public class Mip {

    private static final int DEFAULT_MODEL_IDX = 0;
    private static int SOLVER_IDX = 0;
    private static final int GAP = 0;
    private static final int NON_GAP = 1;
    private static final int DEFAULT_GAP_PENALTY = 2;
    private static final int DEFAULT_THREAD_COUNT = 1;
    public static String[] SOLVERS = new String[] {"SCIP", "Gurobi", "CPSAT"};
    public static String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
    private static final Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid,
                                                                 Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};

    public static void main(String[] args) {
        HashMap<String, Object> argParser = createArgMap(args);
        String ALIGNMENT = (String) argParser.get("ALIGNMENT");
        String NEWICK = (String) argParser.get("NEWICK");
        Integer MODEL_IDX = (Integer) argParser.get("MODEL_IDX");
        String OUTPUT = (String) argParser.get("OUTPUT");
        String PREFIX = (String) argParser.get("PREFIX");
        int NUM_THREADS = argParser.get("THREADS") == null ? DEFAULT_THREAD_COUNT : (Integer) argParser.get("THREADS");

        Tree tree = null;
        EnumSeq.Alignment<Enumerable> aln = null;
        try {
            tree = Utils.loadTree(NEWICK);
            if (MODEL_IDX == null) {
                MODEL_IDX = Mip.DEFAULT_MODEL_IDX;
            }
            aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
            Utils.checkData(aln, tree);
        } catch (ASRException e) {
            usage(IndelDist.ERROR.ASR.getCode(), IndelDist.ERROR.ASR.getDescription() + e.getMessage());
        } catch (IOException e) {
            usage(IndelDist.ERROR.IO.getCode(), IndelDist.ERROR.IO.getDescription() + e.getMessage());
        }
        assert tree != null;
        assert aln != null;

        if (!allColsOccupied(aln)) {
            usage(ERROR.ALN.getCode(), "Alignment has at least one column with no sequence content. Please remove these before running indel inference.");
        }

        long programStartTime = System.currentTimeMillis();
        Loader.loadNativeLibraries();
        POAGraph alnPog = new POAGraph(aln);
        HashMap<Integer, Integer[]> extantBinarySeqs = createBinarySeqMap(aln, tree);

        // Next we create the ancestor penalty array to mirror the tree's children array
        double[][] treeNeighbourAlphaPen = createTreeNeighbourAlphaPen(tree);

        System.out.println("Constructing MIP model...");
        String solverName = SOLVERS[SOLVER_IDX];
        if (solverName.equalsIgnoreCase("CPSAT")) {
            runCPSolverIndelInference(tree, aln, alnPog, NUM_THREADS, extantBinarySeqs, treeNeighbourAlphaPen,
                    OUTPUT, PREFIX);
        } else {
            runMPSolverIndelInference(tree, alnPog, extantBinarySeqs, treeNeighbourAlphaPen, aln,
                    NUM_THREADS, OUTPUT, PREFIX, programStartTime);
        }
    }



    private static void runCPSolverIndelInference(Tree tree, EnumSeq.Alignment<Enumerable> aln ,
                                                  POAGraph alnPog, int num_threads,
                                                  HashMap<Integer, Integer[]> extantBinarySeqs,
                                                  double[][] treeNeighbourAlphaPen, String output, String prefix) {

        CpModel model = new CpModel();
        HashMap<Integer, Literal[]> ancestorPositionVars = createAncestralPositionVariablesCPModel(tree, model, aln.getWidth());

        addEdgeConstraintsAncestorsCPSolver(tree, alnPog, model, ancestorPositionVars);

        LinearExprBuilder objective = LinearExpr.newBuilder();
        addPenaltyConstraintsCPSolver(tree, extantBinarySeqs, treeNeighbourAlphaPen, ancestorPositionVars,
                                      aln.getWidth(), model, objective);

        model.minimize(objective);

        CpSolver solver = new CpSolver();
        solver.getParameters().setNumWorkers(num_threads);
        CpSolverStatus status = solver.solve(model);

        if (status == CpSolverStatus.OPTIMAL) {
            System.out.println("Total cost: " + solver.objectiveValue());
            System.out.println("Problem solved in " + solver.wallTime() + " seconds");
        } else if (status == CpSolverStatus.FEASIBLE) {
            System.out.println("Feasible objective: " + solver.objectiveValue() + "\n");
            System.out.println("Problem solved in " + solver.wallTime() + " seconds");
        } else {
            System.err.println("No solution found.");
        }

        outputAncestralSolutionsCPSolver(tree, ancestorPositionVars, aln, output, prefix, solver);
    }

    private static void outputAncestralSolutionsCPSolver(Tree tree,
                                                 HashMap<Integer, Literal[]> ancestorPositionVars,
                                                 EnumSeq.Alignment<Enumerable> aln,
                                                 String OUTPUT,
                                                 String PREFIX, CpSolver solver) {


        // Create output file path
        String fastaFilePath = OUTPUT + "/" + PREFIX + "_ancestral_indel.fasta";

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fastaFilePath))) {
            for (int ancestralIdx : tree.getAncestors()) {
                Literal[] nodePosVar = ancestorPositionVars.get(ancestralIdx);

                // Write FASTA header
                writer.write(">N" + tree.getLabel(ancestralIdx));
                writer.newLine();

                StringBuilder sequenceLine = new StringBuilder();
                for (int pos = 0; pos < aln.getWidth(); pos++) {
                    sequenceLine.append(solver.value(nodePosVar[pos]));

                    // Optional: Add line breaks every 80 characters (standard FASTA format)
                    if ((pos + 1) % 80 == 0 && pos < aln.getWidth() - 1) {
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

            for (int ancestralIdx : tree.getAncestors()) {
                Literal[] nodePosVar = ancestorPositionVars.get(ancestralIdx);
                StringBuilder sb = new StringBuilder();
                sb.append(">N").append(tree.getLabel(ancestralIdx)).append("\n");
                for (int pos = 0; pos < aln.getWidth(); pos++) {
                    sb.append(solver.value(nodePosVar[pos]));
                }
                System.out.println(sb);
            }
        }
    }

    private static HashMap<Integer, Literal[]> createAncestralPositionVariablesCPModel(Tree tree, CpModel model,
                                                                                       int seqLen)  {

        HashMap<Integer, Literal[]> ancestorSeqVars = new HashMap<>();
        for (int ancestorIdx : tree.getAncestors()) {
            Literal[] seqVars = new Literal[seqLen];
            for (int j = 0; j < seqLen; j++) {
                seqVars[j] = model.newBoolVar("");
            }

            ancestorSeqVars.put(ancestorIdx, seqVars);
        }

        return ancestorSeqVars;
    }

    private static void addEdgeConstraintsAncestorsCPSolver(Tree tree,
                                                                                    POAGraph alnPOG,
                                                                                    CpModel model,
                                                                                    HashMap<Integer, Literal[]> ancestralPositionVars) {

        int VIRTUAL_START = -1;
        int VIRTUAL_END = alnPOG.maxsize();

        HashMap<EdgeKey, Literal> allEdgeVars = new HashMap<>();

        int[] ancestors = tree.getAncestors();
        for (int ancestorIdx : ancestors) {
            // VIRTUAL STARTS TO REAL STARTS //
            int[] startIndices = alnPOG.getStarts();
            Literal[] allEdgesFromVirtualStart = new Literal[startIndices.length];
            for (int positionTo = 0; positionTo < startIndices.length; positionTo++) {

                // Variable for each edge from start node
                Literal edge = model.newBoolVar("");// makeIntVar(0, 1, "");
                allEdgesFromVirtualStart[positionTo] = edge;

                // save unique edge/ancestor combination to a map for later use
                EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, startIndices[positionTo], ancestorIdx);
                allEdgeVars.put(edgeKey, edge);
            }

            // Constraint: exactly one edge must be chosen from virtual start to actual starts
            model.addExactlyOne(allEdgesFromVirtualStart);

            // FULLY CONNECTED NODES //
            int[] endIndices = alnPOG.getEnds();
            for (int nodeIdx = 0; nodeIdx < alnPOG.maxsize(); nodeIdx++) {

                // EDGES GOING OUT //
                int[] forwardEdges = alnPOG.getNodeIndices(nodeIdx, true);
                int numEdgesToEnd = (alnPOG.isEndNode(nodeIdx) ? 1 : 0);
                int totalForwardEdges = forwardEdges.length + numEdgesToEnd;

                Literal[] forwardEdgesFromNode = new Literal[totalForwardEdges];
                for (int positionTo = 0; positionTo < forwardEdges.length; positionTo++) {

                    Literal edge = model.newBoolVar(""); //0, 1, "");
                    int edgeEnd = forwardEdges[positionTo];
                    forwardEdgesFromNode[positionTo] = edge;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, edgeEnd, ancestorIdx);
                    allEdgeVars.put(edgeKey, edge);
                }

                if (numEdgesToEnd == 1) {
                    Literal edgeToEnd = model.newBoolVar(""); //0, 1, "");
                    forwardEdgesFromNode[totalForwardEdges - 1] = edgeToEnd;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, VIRTUAL_END, ancestorIdx);
                    allEdgeVars.put(edgeKey, edgeToEnd);
                }

                // sum(edges) == positionVar
                Literal ancestralNode = ancestralPositionVars.get(ancestorIdx)[nodeIdx];
                model.addEquality(ancestralNode, LinearExpr.sum(forwardEdgesFromNode));

                // EDGES COMING IN //
                int[] backwardEdges = alnPOG.getNodeIndices(nodeIdx, false);
                int numEdgesFromStart = (alnPOG.isStartNode(nodeIdx) ? 1 : 0);
                int totalBackwardEdges = backwardEdges.length + numEdgesFromStart;

                Literal[] backwardEdgesFromNode = new Literal[totalBackwardEdges];
                for (int positionTo = 0; positionTo < backwardEdges.length; positionTo++) {

                    int edgeComingIn = backwardEdges[positionTo];

                    EdgeKey edgeKey = new EdgeKey(edgeComingIn, nodeIdx, ancestorIdx);
                    Literal edge = allEdgeVars.get(edgeKey);
                    backwardEdgesFromNode[positionTo] = edge;
                }

                if (numEdgesFromStart == 1) {
                    EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, nodeIdx, ancestorIdx);
                    Literal edge = allEdgeVars.get(edgeKey);

                    backwardEdgesFromNode[totalBackwardEdges - 1] = edge;
                }

                // Constraint: sum(edges) == position variable
                model.addEquality(ancestralNode, LinearExpr.sum(backwardEdgesFromNode));

            }

            // REAL ENDS TO VIRTUAL ENDS //
            Literal[] allEdgesToVirtualEnd = new Literal[endIndices.length];
            for (int positionFrom = 0; positionFrom < endIndices.length; positionFrom++) {
                EdgeKey edgeKey = new EdgeKey(endIndices[positionFrom], VIRTUAL_END, ancestorIdx);
                allEdgesToVirtualEnd[positionFrom] = allEdgeVars.get(edgeKey);
            }

            model.addExactlyOne(allEdgesToVirtualEnd);
        }
    }

    private static void addPenaltyConstraintsCPSolver(Tree tree,
                                                      HashMap<Integer, Integer[]> extantBinarySeqs,
                                                      double[][] treeNeighbourAlphaPen,
                                                      HashMap<Integer, Literal[]> ancestorPositionVars,
                                                      int seqLength,
                                                      CpModel model,
                                                      LinearExprBuilder objective) {

        for (int ancestralIdx : tree.getAncestors()) {

            Literal[] nodePosVar = ancestorPositionVars.get(ancestralIdx);


            for (int childIdx : tree.getChildren(ancestralIdx)) {
                Literal[] pen = new Literal[seqLength];
                for (int i = 0; i < seqLength; i++) {
                    pen[i] = model.newBoolVar("");//0, 1, "");
                }

                Integer[] nodeNeighbourPosVar = null;
                Literal[] nodeNeighborPosVarAncestor = null;
                Literal[] diffPos = null;

                if (tree.isLeaf(childIdx)) {
                    nodeNeighbourPosVar = extantBinarySeqs.get(childIdx);
                } else {
                    nodeNeighborPosVarAncestor = ancestorPositionVars.get(childIdx);

                    diffPos = new Literal[seqLength];
                    for (int i = 0; i < seqLength; i++) {
                        diffPos[i] = model.newBoolVar("");// 0, 1, "");
                    }
                }

                Literal diffVar = null;
                Literal prevDiffVar = null;
                for (int pos = 0; pos < seqLength; pos++) {

                    if (pos > 0) {
                        prevDiffVar = diffVar;
                    }

                    if (tree.isLeaf(childIdx)) {

                        // constraint - difference variables
                        diffVar = model.newBoolVar("");
                        if (nodeNeighbourPosVar[pos] == 1) {
                            // constraint: diff + nodePosVar[pos] == 1
                            model.addEquality(LinearExpr.sum(new Literal[]{diffVar, nodePosVar[pos]}), 1);

                        } else {
                            // constraint: diffVar = nodePosVar[pos];
                            model.addEquality(diffVar, nodePosVar[pos]);
                        }
                        //diff.put(diffKey, new Literal[]{diffVar});
                        objective.add(diffVar);

                        // penalty constraints
                        if (pos == 0) {
                            // pen[0] == diff[0]
                            model.addEquality(diffVar, pen[pos]);

                        } else {

                            // constraint: pen[pos]>= diff[(node, node_neighbor_item, pos)] - diff[(node, node_neighbor_item, pos - 1)])
                            model.addGreaterOrEqual(LinearExpr.weightedSum(
                                            new Literal[]{
                                                    pen[pos],
                                                    diffVar,
                                                    prevDiffVar},
                                            new long[]{
                                                    1,
                                                    -1,
                                                    1}),
                                        0);

                            // Constraint: pen[pos] >= node_neighbor_pos_var[pos - 1] + (1 - node_pos_var[pos - 1]) + (1 - node_neighbor_pos_var[pos]) + node_pos_var[pos] - 3
                            // because these neighbours are children, we can just evaluate everything on right hand side
                            // rearrange: pen[pos] + node_pos_var[pos - 1] - node_pos_var[pos] >= node_neighbor_pos_var[pos - 1] + 1 + (1 - node_neighbor_pos_var[pos]) - 3
                            // rearrange: pen[pos] + node_pos_var[pos - 1] - node_pos_var[pos] >= node_neighbor_pos_var[pos - 1] - node_neighbor_pos_var[pos] - 1
                            long rhs2 = nodeNeighbourPosVar[pos - 1] - nodeNeighbourPosVar[pos] - 1;
                            model.addGreaterOrEqual(LinearExpr.weightedSum(
                                    new Literal[]{
                                            pen[pos],
                                            nodePosVar[pos - 1],
                                            nodePosVar[pos]},
                                    new long[]{
                                            1,
                                            1,
                                            -1}),
                                    rhs2);

                            // Constraint: pen[pos] >= node_neighbor_pos_var[pos] + (1 - node_pos_var[pos]) + (1 - node_neighbor_pos_var[pos - 1]) + node_pos_var[pos - 1] - 3
                            // Rearrange: pen[pos] + node_pos_var[pos] - node_pos_var[pos - 1]  >= node_neighbor_pos_var[pos] + 1 + (1 - node_neighbor_pos_var[pos - 1]) - 3
                            // Rearrange: pen[pos] + node_pos_var[pos] - node_pos_var[pos - 1]  >= node_neighbor_pos_var[pos] - node_neighbor_pos_var[pos - 1]) - 1
                            long rhs3 = nodeNeighbourPosVar[pos] - nodeNeighbourPosVar[pos - 1] - 1;
                            model.addGreaterOrEqual(LinearExpr.weightedSum(
                                    new Literal[]{
                                            pen[pos],
                                            nodePosVar[pos],
                                            nodePosVar[pos - 1]},
                                    new long[]{
                                            1,
                                            1,
                                            -1}),
                                    rhs3);

                        }

                        objective.addTerm(pen[pos], DEFAULT_GAP_PENALTY);

                    } else { // Ancestor


                        //diff_pos[pos] <= node_pos_var[pos] + node_neighbor_pos_var[pos]
                        model.addLessOrEqual(LinearExpr.weightedSum(
                                new Literal[]{
                                        diffPos[pos],
                                        nodePosVar[pos],
                                        nodeNeighborPosVarAncestor[pos]},
                                new long[]{
                                        1,
                                        -1,
                                        -1}),
                                0);

                        // diff_pos[pos] >= node_pos_var[pos] - node_neighbor_pos_var[pos]
                        model.addGreaterOrEqual(LinearExpr.weightedSum(
                                new Literal[]{
                                        diffPos[pos],
                                        nodePosVar[pos],
                                        nodeNeighborPosVarAncestor[pos]},
                                new long[]{
                                        1,
                                        -1,
                                        1}),
                                0);

                        // diff_pos[pos] >= node_neighbor_pos_var[pos] - node_pos_var[pos]
                        model.addGreaterOrEqual(LinearExpr.weightedSum(
                                new Literal[]{
                                        diffPos[pos],
                                        nodeNeighborPosVarAncestor[pos],
                                        nodePosVar[pos]},
                                new long[]{
                                        1,
                                        -1,
                                        1}),
                                0);

                        // diff_pos[pos] <= 2 - node_neighbor_pos_var[pos] - node_pos_var[pos]
                        model.addLessOrEqual(LinearExpr.weightedSum(
                                new Literal[]{
                                        diffPos[pos],
                                        nodeNeighborPosVarAncestor[pos],
                                        nodePosVar[pos]},
                                new long[]{
                                        1,
                                        1,
                                        1}),
                                2);

                        if (pos == 0) {
                            // constraint: pen[0] == diff[0]
                            model.addEquality(pen[pos], diffPos[pos]);

                        } else {
                            // pen[pos] >= diff_pos[pos] - diff_pos[pos-1]
                            model.addGreaterOrEqual(LinearExpr.weightedSum(
                                    new Literal[]{
                                            pen[pos],
                                            diffPos[pos],
                                            diffPos[pos - 1]},
                                    new long[]{
                                            1,
                                            -1,
                                            1}),
                                    0
                            );

                            // pen[pos] >= node_neighbor_pos_var[pos-1] + (1 - node_pos_var[pos-1]) + (1 - node_neighbor_pos_var[pos]) + node_pos_var[pos] -3
                            model.addGreaterOrEqual(LinearExpr.weightedSum(
                                    new Literal[]{
                                            pen[pos],
                                            nodeNeighborPosVarAncestor[pos - 1],
                                            nodePosVar[pos - 1],
                                            nodeNeighborPosVarAncestor[pos],
                                            nodePosVar[pos]},
                                    new long[]{1, -1, 1, 1, -1}), -1);

                            // pen[pos] >= node_neighbor_pos_var[pos] + (1 - node_pos_var[pos]) + (1 - node_neighbor_pos_var[pos-1]) + node_pos_var[pos-1] - 3
                            model.addGreaterOrEqual(LinearExpr.weightedSum(
                                    new Literal[]{
                                            pen[pos],
                                            nodeNeighborPosVarAncestor[pos],
                                            nodePosVar[pos],
                                            nodeNeighborPosVarAncestor[pos - 1],
                                            nodePosVar[pos - 1]},
                                    new long[]{1, -1, 1, 1, -1}), -1);
                        }

                        objective.add(diffPos[pos]);
                        objective.addTerm(pen[pos], DEFAULT_GAP_PENALTY);
                    }
                }
            }
        }
    }

    private static void runMPSolverIndelInference(Tree tree, POAGraph alnPog,
                                                  HashMap<Integer, Integer[]> extantBinarySeqs,
                                                  double[][] treeNeighbourAlphaPen, EnumSeq.Alignment<Enumerable> aln,
                                                  int num_threads, String OUTPUT, String PREFIX, long programStartTime) {

        MPSolver solver = MPSolver.createSolver(SOLVERS[SOLVER_IDX]);
        if (solver == null) {
            usage(ERROR.MIP_ENGINE.getCode(), ERROR.MIP_ENGINE.getDescription() + SOLVERS[SOLVER_IDX]);
        }
        assert solver != null;

        solver.setNumThreads(num_threads);

        HashMap<Integer, MPVariable[]> ancestorPositionVars = createAncestralPositionVariables(tree, solver, aln.getWidth());
        addEdgeConstraintsAncestors(tree, alnPog, solver, ancestorPositionVars);

        MPObjective objective = solver.objective();
        addPenaltyConstraints(tree, extantBinarySeqs, treeNeighbourAlphaPen, ancestorPositionVars,
                aln.getWidth(), solver, objective);
        long endModelBuild = System.currentTimeMillis();
        //System.out.println("Constructed MIP model in " + (endModelBuild - startModelBuild) + " ms");

        objective.setMinimization();

        MPSolver.ResultStatus resultStatus = solver.solve();

        long progEndTime = System.currentTimeMillis();
        System.out.println("Total run time: " + (progEndTime - programStartTime) / 1000 + " seconds");
        if (resultStatus == MPSolver.ResultStatus.OPTIMAL) {
            System.out.println("Total cost: " + solver.objective().value());
            System.out.println("Problem solved in " + solver.wallTime() / 1000 + " seconds");
            //System.out.println("Problem solved in " + solver.iterations() + " iterations");
            //System.out.println(solver.nodes() + " branch-and-bound nodes");
            //System.out.println("Number of variables " + solver.numConstraints());
        } else if (resultStatus == MPSolver.ResultStatus.FEASIBLE) {
            System.out.println("Feasible objective: " + objective.value() + "\n");
        } else {
            System.err.println("No solution found.");
        }

        outputAncestralSolutions(tree, ancestorPositionVars, aln, OUTPUT, PREFIX);

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

    public static HashMap<Integer, Integer[]> createBinarySeqMap(EnumSeq.Alignment<Enumerable> aln, Tree tree) {

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

    private static double[][] createTreeNeighbourAlphaPen(Tree tree) {

        double[][] tree_neighbour_alpha_pen = new double[tree.getNParents() + tree.getNLeaves()][];
        for (int i = 0; i < tree_neighbour_alpha_pen.length; i++) {
            int[] children = tree.getChildren(i);
            tree_neighbour_alpha_pen[i] = new double[children.length];
            for (int j = 0; j < children.length; j++) {
                tree_neighbour_alpha_pen[i][j] = DEFAULT_GAP_PENALTY;
            }
        }

        return tree_neighbour_alpha_pen;
    }

    /**
     * Create integer position variables for ancestor sequences bounded between 0 and 1. These represent the
     * actual sequence positions for each ancestor in the tree.
     *
     * @param tree the phylogenetic tree
     * @param solver the MIP solver
     * @param seqLength the length of the sequences in the alignment
     * @return a map of ancestor indices to their corresponding position variables
     */
    private static HashMap<Integer, MPVariable[]> createAncestralPositionVariables(Tree tree, MPSolver solver, int seqLength) {

        HashMap<Integer, MPVariable[]> ancestorSeqVars = new HashMap<>();
        for (int ancestorIdx : tree.getAncestors()) {
            MPVariable[] seqVars = new MPVariable[seqLength];
            for (int j = 0; j < seqLength; j++) {
                // NOTE: May want to check that create labels won't create large memory usage.
                seqVars[j] = solver.makeBoolVar(""); //0, 1, "");
            }
            ancestorSeqVars.put(ancestorIdx, seqVars);
        }

        return ancestorSeqVars;

    }

    private static void addEdgeConstraintsAncestors(Tree tree, POAGraph alnPOG,
                                                                            MPSolver solver,
                                                                            HashMap<Integer, MPVariable[]> ancestralPositionVars) {

        int VIRTUAL_START = -1;
        int VIRTUAL_END = alnPOG.maxsize();

        HashMap<EdgeKey, MPVariable> allEdgeVars = new HashMap<>();

        int[] ancestors = tree.getAncestors();
        for (int ancestorIdx : ancestors) {
            // VIRTUAL STARTS TO REAL STARTS //
            int[] startIndices = alnPOG.getStarts();
            MPVariable[] allEdgesFromVirtualStart = new MPVariable[startIndices.length];
            for (int positionTo = 0; positionTo < startIndices.length; positionTo++) {

                // Variable for each edge from start node
                MPVariable edge = solver.makeBoolVar("");// makeIntVar(0, 1, "");
                allEdgesFromVirtualStart[positionTo] = edge;

                // save unique edge/ancestor combination to a map for later use
                EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, startIndices[positionTo], ancestorIdx);
                allEdgeVars.put(edgeKey, edge);
            }

            // Constraint: exactly one edge must be chosen from virtual start to actual starts
            MPConstraint virtualStartConstraint = solver.makeConstraint(1, 1);
            for (MPVariable edgeVar : allEdgesFromVirtualStart) {
                virtualStartConstraint.setCoefficient(edgeVar, 1);
            }

            // FULLY CONNECTED NODES //
            int[] endIndices = alnPOG.getEnds();
            for (int nodeIdx = 0; nodeIdx < alnPOG.maxsize(); nodeIdx++) {

                // EDGES GOING OUT //
                int[] forwardEdges = alnPOG.getNodeIndices(nodeIdx, true);
                int numEdgesToEnd = (alnPOG.isEndNode(nodeIdx) ? 1 : 0);
                int totalForwardEdges = forwardEdges.length + numEdgesToEnd;

                MPVariable[] forwardEdgesFromNode = new MPVariable[totalForwardEdges];
                for (int positionTo = 0; positionTo < forwardEdges.length; positionTo++) {

                    MPVariable edge = solver.makeBoolVar(""); //0, 1, "");
                    int edgeEnd = forwardEdges[positionTo];
                    forwardEdgesFromNode[positionTo] = edge;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, edgeEnd, ancestorIdx);
                    allEdgeVars.put(edgeKey, edge);
                }

                if (numEdgesToEnd == 1) {
                    MPVariable edgeToEnd = solver.makeBoolVar(""); //0, 1, "");
                    forwardEdgesFromNode[totalForwardEdges - 1] = edgeToEnd;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, VIRTUAL_END, ancestorIdx);
                    allEdgeVars.put(edgeKey, edgeToEnd);
                }

                // Constraint: sum(edges) - position variable == 0
                // Which is equivalent to: sum(edges) == positionVar
                MPConstraint forwardConstraint = solver.makeConstraint(0, 0);
                for (MPVariable edgeVar : forwardEdgesFromNode) {
                    forwardConstraint.setCoefficient(edgeVar, 1);
                }

                MPVariable ancestralNode = ancestralPositionVars.get(ancestorIdx)[nodeIdx];
                forwardConstraint.setCoefficient(ancestralNode, -1);

                // EDGES COMING IN //
                int[] backwardEdges = alnPOG.getNodeIndices(nodeIdx, false);
                int numEdgesFromStart = (alnPOG.isStartNode(nodeIdx) ? 1 : 0);
                int totalBackwardEdges = backwardEdges.length + numEdgesFromStart;

                MPVariable[] backwardEdgesFromNode = new MPVariable[totalBackwardEdges];
                for (int positionTo = 0; positionTo < backwardEdges.length; positionTo++) {

                    int edgeComingIn = backwardEdges[positionTo];

                    EdgeKey edgeKey = new EdgeKey(edgeComingIn, nodeIdx, ancestorIdx);
                    MPVariable edge = allEdgeVars.get(edgeKey);
                    backwardEdgesFromNode[positionTo] = edge;
                }

                if (numEdgesFromStart == 1) {
                    EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, nodeIdx, ancestorIdx);
                    MPVariable edge = allEdgeVars.get(edgeKey);

                    backwardEdgesFromNode[totalBackwardEdges - 1] = edge;
                }

                // Constraint: sum(edges) - position variable == 0
                // Which is equivalent to: sum(edges) == positionVar
                MPConstraint backwardConstraint = solver.makeConstraint(0, 0);
                for (MPVariable edgeVar : backwardEdgesFromNode) {
                    backwardConstraint.setCoefficient(edgeVar, 1);
                }
                backwardConstraint.setCoefficient(ancestralNode, -1);
            }

            // REAL ENDS TO VIRTUAL ENDS //
            MPVariable[] allEdgesToVirtualEnd = new MPVariable[endIndices.length];
            for (int positionFrom = 0; positionFrom < endIndices.length; positionFrom++) {
                EdgeKey edgeKey = new EdgeKey(endIndices[positionFrom], VIRTUAL_END, ancestorIdx);
                allEdgesToVirtualEnd[positionFrom] = allEdgeVars.get(edgeKey);
            }

            MPConstraint virtualEndConstraint = solver.makeConstraint(1, 1);
            for (MPVariable edgeVar : allEdgesToVirtualEnd) {
                virtualEndConstraint.setCoefficient(edgeVar, 1);
            }
        }
    }


    private static void addPenaltyConstraints(Tree tree,
                                              HashMap<Integer, Integer[]> extantBinarySeqs,
                                              double[][] treeNeighbourAlphaPen,
                                              HashMap<Integer, MPVariable[]> ancestorPositionVars,
                                              int seqLength,
                                              MPSolver solver,
                                              MPObjective objective) {



        for (int ancestralIdx : tree.getAncestors()) {

            MPVariable[] nodePosVar = ancestorPositionVars.get(ancestralIdx);


            for (int childIdx : tree.getChildren(ancestralIdx)) {

                MPVariable[] pen = new MPVariable[seqLength];
                for (int i = 0; i < seqLength; i++) {
                    pen[i] = solver.makeBoolVar("");//0, 1, "");
                }

                Integer[] nodeNeighbourPosVar = null;
                MPVariable[] nodeNeighborPosVarAncestor = null;
                MPVariable[] diffPos = null;

                if (tree.isLeaf(childIdx)) {
                    nodeNeighbourPosVar = extantBinarySeqs.get(childIdx);
                } else {
                    nodeNeighborPosVarAncestor = ancestorPositionVars.get(childIdx);

                    diffPos = new MPVariable[seqLength];
                    for (int i = 0; i < seqLength; i++) {
                        diffPos[i] = solver.makeBoolVar("");// 0, 1, "");
                    }
                }

                MPVariable diffVar = null;
                MPVariable prevDiffVar = null;
                for (int pos = 0; pos < seqLength; pos++) {

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
                        objective.setCoefficient(pen[pos], DEFAULT_GAP_PENALTY);

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
                        objective.setCoefficient(pen[pos], DEFAULT_GAP_PENALTY);
                    }
                }
            }
        }
    }

    private static void outputAncestralSolutions(Tree tree,
                                                 HashMap<Integer, MPVariable[]> ancestorPositionVars,
                                                 EnumSeq.Alignment<Enumerable> aln,
                                                 String OUTPUT,
                                                 String PREFIX) {


        // Create output file path
        String fastaFilePath = OUTPUT + "/" + PREFIX + "_ancestral_indel.fasta";

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fastaFilePath))) {
            for (int ancestralIdx : tree.getAncestors()) {
                MPVariable[] nodePosVar = ancestorPositionVars.get(ancestralIdx);

                // Write FASTA header
                writer.write(">N" + tree.getLabel(ancestralIdx));
                writer.newLine();

                StringBuilder sequenceLine = new StringBuilder();
                for (int pos = 0; pos < aln.getWidth(); pos++) {
                    sequenceLine.append((int) nodePosVar[pos].solutionValue());

                    // Optional: Add line breaks every 80 characters (standard FASTA format)
                    if ((pos + 1) % 80 == 0 && pos < aln.getWidth() - 1) {
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

            for (int ancestralIdx : tree.getAncestors()) {
                MPVariable[] nodePosVar = ancestorPositionVars.get(ancestralIdx);
                StringBuilder sb = new StringBuilder();
                sb.append(">N").append(tree.getLabel(ancestralIdx)).append("\n");
                for (int pos = 0; pos < aln.getWidth(); pos++) {
                    sb.append((int)nodePosVar[pos].solutionValue());
                }
                System.out.println(sb);
            }
        }
    }


    /**
     * Create a unique key for an edge in the POG associated with a specific ancestor.
     * Allows this to be used as a key in hash maps.
     */
    public record EdgeKey(int from, int to, int ancestorIdx) {}

    /**
     * Error codes for the program
     */
    public enum ERROR {
        SUCCESS(0, null),
        ALN(1, "Must specify alignment (--aln <Clustal or FASTA file>)"),
        NWK( 2, "Must specify phylogenetic tree (Newick file) or previously saved folder (--input-folder <folder>)"),
        UNKNOWN(3, "Unknown option or missing required argument: "),
        IO(4, "Failed to read or write files: "),
        SUB_MODEL(5, "Could not find model with ID "),
        MIP_ENGINE(6, "Could not create MIP solver with "),
        INVALID_SOLVER(7, " is not a valid solver name for option --solver");

        private final int code;
        private final String description;

        ERROR(int code, String description) {
            this.code = code;
            this.description = description;
        }

        public String getDescription() {
            return description;
        }

        public int getCode() {
            return code;
        }
    }

    public static void usage(int error_code, String msg) {
        PrintStream out = System.out;
        if (error_code != 0)
            out = System.err;

        out.println("""
                Usage: asr.Mip\s
                \t[-a | --aln <filename>]
                \t[-n | --nwk <filename>]
                \t{-s | --substitution-model <JTT(default)}
                \t{-t | --threads <number>}
                \t{-pre | --prefix <stub>}
                \t{-o | --output-folder <foldername>} (default is current working folder, or input folder if available)
                \t{-so | --solver <SCIP (default), Gurobi, CPSAT>}
                \t-h (or --help) will print out this screen
                """
        );

        if (msg != null) {
            out.println("\n" + msg + "(Error code " + error_code + ")");
        }
        System.exit(error_code);
    }

    public static void usage() {
        usage(ERROR.SUCCESS.getCode(), ERROR.SUCCESS.getDescription());
    }

    public static HashMap<String, Object> createArgMap(String[] args) {

        HashMap<String, Object> argMap = new HashMap<>();

        parseArgs(args, argMap);

        checkArgsValid(argMap);

        return argMap;
    }

    private static void parseArgs(String[] args, HashMap<String, Object> argMap) {
        // Read in all the arguments
        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (((arg.equalsIgnoreCase("-aln")) || (arg.equalsIgnoreCase("a"))) && args.length > a + 1) {
                    argMap.put("ALIGNMENT", args[++a]);
                } else if ((arg.equalsIgnoreCase("-output-folder") || arg.equalsIgnoreCase("o")) && args.length > a + 1) {
                    argMap.put("OUTPUT", args[++a]);
                } else if (((arg.equalsIgnoreCase("-nwk")) || (arg.equalsIgnoreCase("n"))) && args.length > a + 1) {
                    argMap.put("NEWICK", args[++a]);
                } else if ((arg.equalsIgnoreCase("-prefix") || arg.equalsIgnoreCase("pre")) && args.length > a + 1) {
                    argMap.put("PREFIX", args[++a]);
                } else if ((arg.equalsIgnoreCase("-solver") || arg.equalsIgnoreCase("so")) && args.length > a + 1) {
                    boolean foundSolver = false;
                    for (int i = 0; i < SOLVERS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(SOLVERS[i])) {
                            SOLVER_IDX = i;
                            foundSolver = true;
                        }
                    }
                    if (!foundSolver)
                        usage(ERROR.INVALID_SOLVER.getCode(), args[a + 1] + ERROR.INVALID_SOLVER.getDescription());

                } else if ((arg.equalsIgnoreCase("-threads") || arg.equalsIgnoreCase("t")) && args.length > a + 1) {
                    try {
                        argMap.put("THREADS", Integer.parseInt(args[++a]));
                    } catch (NumberFormatException e) {
                        System.err.println("Failed to set number of threads for option --threads: " + args[a] + " is not a valid integer");
                    }
                } else if ((arg.equalsIgnoreCase("-substitution-model") || arg.equalsIgnoreCase("s")) && args.length > a + 1) {
                    boolean model_found = false;
                    for (int i = 0; i < MODELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(MODELS[i])) {
                            model_found = true;
                            argMap.put("MODEL_IDX", i);
                        }
                    }

                    if (!model_found) {
                        usage(ERROR.SUB_MODEL.getCode(), ERROR.SUB_MODEL.getDescription() + args[a + 1]);
                    }
                } else if ((arg.equalsIgnoreCase("-input-folder") || arg.equalsIgnoreCase("i")) && args.length > a + 1) {
                    argMap.put("INPUT", args[++a]);
                } else {
                        usage(ERROR.UNKNOWN.getCode(), ERROR.UNKNOWN.getDescription() + "\" + args[a] + \"");
                    }
            }
        }
    }

    private static void checkArgsValid(HashMap<String, Object> argMap) {

        String ALIGNMENT = (String) argMap.get("ALIGNMENT");
        String INPUT = (String) argMap.get("INPUT");
        String OUTPUT = (String) argMap.get("OUTPUT");
        String PREFIX = (String) argMap.get("PREFIX");

        if (!argMap.containsKey("ALIGNMENT")) {
            usage(ERROR.ALN.getCode(), ERROR.ALN.getDescription());
        } else if (!argMap.containsKey("NEWICK")) {
            usage(ERROR.NWK.getCode(), ERROR.NWK.getDescription());
        } else if (OUTPUT == null) {
            OUTPUT = INPUT == null ? "." : INPUT;
        }

        if (PREFIX == null) { // default prefix is the (prefix of) alignment filename
            int idx2 = ALIGNMENT == null ? 0 : ALIGNMENT.lastIndexOf(".");
            if (idx2 == -1)
                idx2 = ALIGNMENT.length();
            int idx1 = ALIGNMENT == null ? 0 : ALIGNMENT.lastIndexOf("/") + 1;
            PREFIX = ALIGNMENT == null ? "" : ALIGNMENT.substring(idx1, idx2);
        }

        argMap.put("OUTPUT", OUTPUT);
        argMap.put("PREFIX", PREFIX);

    }
}
