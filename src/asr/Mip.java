package asr;

import com.google.ortools.linearsolver.MPConstraint;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;
import dat.pog.POAGraph;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import com.google.ortools.linearsolver.MPVariable;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.Loader;

public class Mip {

    private static final int DEFAULT_MODEL = 0;
    private static final int VIRTUAL_NODES = 2;
    private static final int GAP = 0;
    private static final int NON_GAP = 1;
    private static final int DEFAULT_GAP_PENALTY = 2;

    public static String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
    private static final Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid,
                                                                 Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};

    public static void main(String[] args) {
        HashMap<String, Object> argParser = createArgMap(args);
        String ALIGNMENT = (String) argParser.get("ALIGNMENT");
        String NEWICK = (String) argParser.get("NEWICK");
        Integer MODEL_IDX = (Integer) argParser.get("MODEL_IDX");

        Tree tree = null;
        EnumSeq.Alignment<Enumerable> aln = null;
        try {
            tree = Utils.loadTree(NEWICK);
            if (MODEL_IDX == null) {
                MODEL_IDX = DEFAULT_MODEL;
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

        // Create the linear solver with the SCIP backend.
        Loader.loadNativeLibraries();
        MPSolver solver = MPSolver.createSolver("SCIP");
        if (solver == null) {
            usage(ERROR.MIP_ENGINE.getCode(), ERROR.MIP_ENGINE.getDescription());
        }

        // have tree, pog with edges, and binary sequences for extant taxa
        POAGraph alnPog = new POAGraph(aln);

        HashMap<String, Enumerable> extantBinarySeqs = createBinarySeqMap(aln);

        // Next we create the ancestor penalty array to mirror the tree children array
        double[][] tree_neighbour_alpha_pen = createTreeNeighbourAlphaPen(tree);


        HashMap<Integer, MPVariable[]> ancestorPositionVars = createAncestralPositionVariables(tree, solver, aln.getWidth());

        HashMap<EdgeKey, MPVariable> allEdges = addEdgeConstraintsAncestors(tree, alnPog, solver, ancestorPositionVars);

        addPenaltyConstraints();


    }

    public static HashMap<String, Enumerable> createBinarySeqMap(EnumSeq.Alignment<Enumerable> aln) {

        int VIRTUAL_START = 0;
        int VIRTUAL_END = aln.getWidth() + 1;
        HashMap<String, Enumerable> binarySeqMap = new HashMap<>();

        for (int i = 0; i < aln.getHeight(); i++) {

            EnumSeq<Enumerable> seq = aln.getEnumSeq(i);
            Integer[] binSeqValues = new Integer[aln.getWidth() + VIRTUAL_NODES];
            binSeqValues[VIRTUAL_START] = NON_GAP;
            binSeqValues[VIRTUAL_END] = NON_GAP;

            for (int j = 0; j < aln.getWidth(); j++) {
                if (seq.get(j) == null) {
                    binSeqValues[j + 1] = GAP;
                } else {
                    binSeqValues[j + 1] = NON_GAP;
                }
            }

            binarySeqMap.put(seq.getName(), new Enumerable(binSeqValues));
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

        int[] ancestors = tree.getAncestors();
        HashMap<Integer, MPVariable[]> ancestorSeqVars = new HashMap<>();
        for (int i = 0; i < ancestors.length; i++) {
            int ancestorIdx = ancestors[i];
            MPVariable[] seqVars = new MPVariable[seqLength];
            for (int j = 0; j < seqLength; j++) {
                // NOTE: May want to check that create labels won't create large memory usage.
                seqVars[j] = solver.makeIntVar(0, 1, "A" + ancestorIdx + "_Pos" + j);
            }
            ancestorSeqVars.put(ancestorIdx, seqVars);
        }

        return ancestorSeqVars;

    }

    private static HashMap<EdgeKey, MPVariable> addEdgeConstraintsAncestors(Tree tree, POAGraph alnPOG,
                                                                            MPSolver solver,
                                                                            HashMap<Integer, MPVariable[]> ancestralPositionVars) {

        int VIRTUAL_START = -1;
        int VIRTUAL_END = alnPOG.maxsize();

        HashMap<EdgeKey, MPVariable> allEdgeVars = new HashMap<>();

        int[] ancestors = tree.getAncestors();
        for (int i = 0; i < ancestors.length; i++) {
            int ancestorIdx = ancestors[i];

            // VIRTUAL STARTS TO REAL STARTS //
            int[] startIndices = alnPOG.getStarts();
            MPVariable[] allEdgesFromVirtualStart = new MPVariable[startIndices.length];
            for (int positionTo = 0; positionTo < startIndices.length; positionTo++) {

                // Variable for each edge from start node
                MPVariable edge = solver.makeIntVar(0, 1, "");
                allEdgesFromVirtualStart[positionTo] = edge;

                // save unique edge/ancestor combination to a map for later use
                EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, startIndices[positionTo], ancestorIdx);
                allEdgeVars.put(edgeKey, edge);
            }

            // Constraint: exactly one edge must be chosen from virtual start to actual starts
            MPConstraint virtualStartConstraint = solver.makeConstraint(1, 1, "");
            for (MPVariable edgeVar : allEdgesFromVirtualStart) {
                virtualStartConstraint.setCoefficient(edgeVar, 1);
            }

            // FULLY CONNECTED NODES //
            //alnPOG.getNodeIndices()

            for (int nodeIdx = 0; nodeIdx < alnPOG.maxsize(); nodeIdx++) {

                // EDGES GOING OUT //
                int[] forwardEdges = alnPOG.getNodeIndices(nodeIdx, true);
                int numEdgesToEnd = (alnPOG.isEndNode(nodeIdx) ? 1 : 0);
                int totalForwardEdges = forwardEdges.length + numEdgesToEnd;

                if (totalForwardEdges == 0) {
                    continue;
                }

                MPVariable[] forwardEdgesFromNode = new MPVariable[totalForwardEdges];
                for (int positionTo = 0; positionTo < forwardEdges.length; positionTo++) {

                    MPVariable edge = solver.makeIntVar(0, 1, "");
                    int edgeEnd = forwardEdges[positionTo];
                    forwardEdgesFromNode[positionTo] = edge;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, edgeEnd, ancestorIdx);
                    allEdgeVars.put(edgeKey, edge);
                }

                if (numEdgesToEnd == 1) {
                    MPVariable edgeToEnd = solver.makeIntVar(0, 1, "");
                    forwardEdgesFromNode[totalForwardEdges - 1] = edgeToEnd;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, VIRTUAL_END, ancestorIdx);
                    allEdgeVars.put(edgeKey, edgeToEnd);
                }

                // Constraint: sum(edges) - position variable == 0
                // Which is equivalent to: sum(edges) == positionVar
                MPConstraint forwardConstraint = solver.makeConstraint(0, 0, "");
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

                    MPVariable edge = solver.makeIntVar(0, 1, "");
                    int edgeEnd = backwardEdges[positionTo];
                    backwardEdgesFromNode[positionTo] = edge;

                    EdgeKey edgeKey = new EdgeKey(nodeIdx, edgeEnd, ancestorIdx);
                    allEdgeVars.put(edgeKey, edge);
                }

                if (numEdgesFromStart == 1) {
                    MPVariable edgeFromStart = solver.makeIntVar(0, 1, "");
                    backwardEdgesFromNode[totalBackwardEdges - 1] = edgeFromStart;

                    EdgeKey edgeKey = new EdgeKey(VIRTUAL_START, nodeIdx, ancestorIdx);
                    allEdgeVars.put(edgeKey, edgeFromStart);
                }

                // Constraint: sum(edges) - position variable == 0
                // Which is equivalent to: sum(edges) == positionVar
                MPConstraint backwardConstraint = solver.makeConstraint(0, 0, "");
                for (MPVariable edgeVar : backwardEdgesFromNode) {
                    backwardConstraint.setCoefficient(edgeVar, 1);
                }
                backwardConstraint.setCoefficient(ancestralNode, -1);

            }


        }

        return allEdgeVars;
    }

    private static void constrainNumberOfConnectedEdges(int[] outgoingEdges, int isTerminal, MPSolver solver,
                                                        int nodeIdx, int ancestorIdx,
                                                        HashMap<EdgeKey, MPVariable> allEdgeVars,
                                                        HashMap<Integer, MPVariable[]> ancestralPositionVars,
                                                        int virtualEnd){

        int totalOutgoingEdges = outgoingEdges.length + isTerminal;
        MPVariable[] outgoingEdgesFromNode = new MPVariable[totalOutgoingEdges];
        for (int positionTo = 0; positionTo < outgoingEdges.length; positionTo++) {

            MPVariable edge = solver.makeIntVar(0, 1, "");
            int edgeEnd = outgoingEdges[positionTo];
            outgoingEdgesFromNode[positionTo] = edge;

            EdgeKey edgeKey = new EdgeKey(nodeIdx, edgeEnd, ancestorIdx);
            allEdgeVars.put(edgeKey, edge);
        }

        if (isTerminal == 1) {
            MPVariable edgeToEnd = solver.makeIntVar(0, 1, "");
            outgoingEdgesFromNode[totalOutgoingEdges - 1] = edgeToEnd;

            EdgeKey edgeKey = new EdgeKey(nodeIdx, virtualEnd, ancestorIdx);
            allEdgeVars.put(edgeKey, edgeToEnd);
        }

        // Constraint: sum(edges) - position variable == 0
        // Which is equivalent to: sum(edges) == positionVar
        MPConstraint forwardConstraint = solver.makeConstraint(0, 0, "");
        for (MPVariable edgeVar : outgoingEdgesFromNode) {
            forwardConstraint.setCoefficient(edgeVar, 1);
        }
        MPVariable ancestralNode = ancestralPositionVars.get(ancestorIdx)[nodeIdx];
        forwardConstraint.setCoefficient(ancestralNode, -1);
    }


    private static void addPenaltyConstraints() {

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
        MIP_ENGINE(6, "Could not create MIP solver with SCIP"),;

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
                Usage: asr.IndelDist\s
                \t[-a | --aln <filename>]
                \t[-n | --nwk <filename>]
                \t{-s | --substitution-model <JTT(default)}
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
                } else if (((arg.equalsIgnoreCase("-nwk")) || (arg.equalsIgnoreCase("n"))) && args.length > a + 1) {
                    argMap.put("NEWICK", args[++a]);
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
                } else {
                        usage(ERROR.UNKNOWN.getCode(), ERROR.UNKNOWN.getDescription() + "\" + args[a] + \"");
                    }
            }
        }
    }

    private static void checkArgsValid(HashMap<String, Object> argMap) {

        if (!argMap.containsKey("ALIGNMENT")) {
            usage(IndelDist.ERROR.ALN.getCode(), IndelDist.ERROR.ALN.getDescription());
        } else if (!argMap.containsKey("NEWICK")) {
            usage(IndelDist.ERROR.NWK.getCode(), IndelDist.ERROR.NWK.getDescription());
        }
    }
}
