package asr;

import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.*;
import dat.phylo.BranchPoint;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import dat.pog.*;
import stats.*;

import java.io.*;
import java.util.*;

import static bn.prob.GammaDistrib.calcAlpha;

/**
 * Track Ancestor Via Indels and Substitutions.
 *
 * This application is intended to simulate the evolution of biological sequences FROM a given ancestor, TO extants via
 * substitutions, insertions and deletions. These events are tracked, meaning that intermediate ancestors are recorded and
 * an alignment that includes all extants and ancestors can be created. Preceding the alignment, the program generates a
 * partial order graph that helps consolidate multiple insertions and deletions at and across ancestor branchpoints.
 *
 * First, a phylogenetic tree is generated from user-specified parameters, using a gamma-like distribution for setting
 * evolutionary distances on each branch, in a bi- or multi-furcating manner, until a user-specified number of extants
 * have been mapped as leaves.
 *
 * Mutation events are determined stochastically following and amending principles described in
 * Cartwright R. Problems and Solutions for Estimating Indel Rates and Length Distributions.
 * Mol. Biol. Evol. 26(2):473–480. 2009. <a href="https://doi.org/10.1093/molbev/msn275">...</a>
 *
 * For a given ancestor sequence (of arbitrary length) each position is looked at deciding whether to introduce
 * "indels" in the descendant, by a probability e^-(rho*t*r), where "rho*t" is the indel rate times the evolutionary distance on the branch;
 * optionally the position-specific substitution rate is also used (1 by default). There's also a "deletion" proportion that dictates whether
 * a deletion or insertion should be introduced. Then, there's the length of the indel, which is decided by another Poisson-like distribution
 * (several are provided).
 * At each position, when "matched" (not "indel"), a substitution is introduced with a probability determined by a (specified) evolutionary model.
 *
 * Ultimately, tree and alignment files are saved.
 *
 * @author Mikael Boden
 * @author Chongting Zhao
 * @author Sebastian Porras
 */
public class TrAVIS {
    public static void usage() {
        usage(0, null);
    }

    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        if (msg != null)
            out.println(msg + " (Error " + error + ")");
        out.println("Usage: asr.TrAVIS [options]\n" +
                "Options:\n" +
                "  -n0, --ancestor <sequence>                   Specify ancestor sequence (text string)\n" +
                "  -o, --output-folder <folder>                 Output directory for results (created if missing)\n" +
                "  -n, --nwk <file>                             Phylogenetic tree in Newick format\n" +
                "  --extants <number>                           Number of extant leaves (default: 5)\n" +
                "  -pre, --prefix <stub>                        file prefix for outputs\n" +
                "  --no-recon                                   Bypasses reconstruction and just produces TrAVIS parameters based on input tree and alignment\n" +
                "  -s, --substitution-model <model>             Substitution model: JTT (default), Dayhoff, LG, WAG, JC, Yang\n" +
                "  -l, --length <number>                        Length of ancestor sequence (if not provided with --ancestor). \n" +
                "                                               The sequence will be randomly generated according to the substitution model\n" +
                "  -rf, --rates-file <file>                     Tabulated file with site-specific rates (IQ-TREE format)\n" +
                "      --dist-distrib <type:params>             Branch distance (t) distribution: Gamma:<SHAPE>,<SCALE> ZeroInflatedGamma:<PI>,<SHAPE>,<SCALE> MixtureGamma:<SHAPE1>,SCALE1>,<WEIGHT1>,<SHAPE2>,SCALE2>,<WEIGHT2>\n" +
                "      --leaf2root-distrib <type:params>        Distribution for leaf-to-root distances: Gaussian, GDF\n" +
                "      --subst-rate-distrib <type:params>       Substitution rate (r) distribution: Gamma, ZeroInflatedGamma, MixtureGamma\n" +
                "      --indel-rate-distrib <type:params>       Indel rate (rho) distribution: Gamma, ZeroInflatedGamma, MixtureGamma\n" +
                "      --indel-length-distrib <type:params>     Indel length distribution: ZeroTruncatedPoisson, Poisson, Zipf, Lavalette\n" +
                "      --insertion-length-distrib <type:params> Insertion length distribution (overrides indel-length)\n" +
                "      --deletion-length-distrib <type:params>  Deletion length distribution (overrides indel-length)\n" +
                "      --delprop <fraction>                     Proportion of deletions among indels (0-1, default: 0.5)\n" +
                "      -a, --aln                                Alignment to learn indel parameters from, must be used with --learn\n" +
                "      --learn                                  Don't do simulation, only infer simulation parameters\n" +
                "      --copy-tree                              Use the tree from -nwk to perform the simulation\n" +
                "      --extants-only                           Create a separate FASTA file with simulated extants only\n" +
                "      --no-gap                                 exclude gap characters in output\n" +
                "      -sa, --save-as <type>                         Output format: FASTA (default), CLUSTAL, DOT, TREE, RATES, DIR\n" +
                "      --seed <number>                          Random seed\n" +
                "      -t, --threads <number>                   Number of threads to use. \n" +
                "      --verbose                                Print details of generated events\n" +
                "  -h, --help                                   Show this help message\n");
        out.println("Notes:\n" +
                "  * Options marked with an asterisk are not fully implemented.\n" +
                "  - Substitution models for proteins: JTT, Dayhoff, LG, WAG; for DNA: JC, Yang.\n" +
                "  - Site specific rates are set from the Gamma distribution, either by specified parameters or as estimated from rates file.\n" +
                "  - If a phylogenetic tree file is given, the simulated tree is generated from a distribution estimated from it.\n" +
                "  - If no tree file is given, a random tree is generated using the specified distribution.\n" +
                "  - Distribution parameters:\n" +
                "      Gamma: <shape,scale>\n" +
                "      ZeroInflatedGamma: <pi,shape,scale>\n" +
                "      ZeroTruncatedPoisson/Poisson: <lambda>\n" +
                "      Zipf: <s[,max]>\n" +
                "      Lavalette: <a[,max]>\n" +
                "  - Output formats: FASTA, CLUSTAL, DOT, TREE, RATES, DIR\n" +
                "  - For more details, see documentation or contact the authors.\n" +
                "  ~ This is part of GRASP-Suite version " + GRASP.VERSION + " ~");
        System.exit(error);
    }

    public static Boolean VERBOSE = false;
    public static String OUTPUT = null;
    private static final int DISTRIB_NAME = 0;
    private static final int DISTRIB_PARAMS = 1;
    private static double TREE_GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
    private static double TREE_GAMMA_SCALE = 0.2;
    private static final int DESCENDANTS_MAX = 2, DESCENDANTS_MIN = 2; // Max and min of tree branching
    private static double DELETIONPROP = 0.5; // proportion of DELETIONS v INSERTIONS
    private static final int EVOL_MODEL_IDX = 0; // default model is that above indexed
    private static final String[] EVOL_MODELS = new String[]{"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
    private static String ANCSEQ = null; // ancestor sequence as a text string, provided
    private static String INPUT_TREE = null;
    private static Double GAMMA_ALPHA = null;
    private static Double SCALEDIST = null;
    private static String SRATESFILE = null;
    private static  double[] SRATES = null;
    private static long SEED = System.currentTimeMillis();
    private static int EXTANTS_N = 5;
    private static SubstModel EVOL_MODEL = null;
    private static IndelModel INDEL_LENGTH_MODEL = null;
    private static IndelModel INSERTION_LENGTH_MODEL = null;
    private static IndelModel DELETION_LENGTH_MODEL = null;
    private static RateModel SUBST_RATE_MODEL = null;
    private static RateModel INDEL_RATE_MODEL = null;
    private static RateModel TREE_DISTANCE_MODEL = null;
    private static Distrib LEAF2ROOT_DISTANCE_MODEL = null;
    private static Integer ANCSEQ_LENGTH = null;
    private static boolean GAPPY = false;
    public static final String[] TRAVIS_FORMATS = new String[]{"FASTA", "CLUSTAL", "DOT", "ALL", "RATES"};
    private static int FORMAT_IDX = 0;
    private static final int FASTA = 0;
    private static final int CLUSTAL = 1;
    private static final int DOT = 2;
    private static final int ALL = 3;
    private static final int RATES = 4;
    private static boolean COPY_TREE = false;
    private static boolean EXTANTS_ONLY = false;
    private static boolean LEARN = false;
    private static String ALIGNMENT;
    public enum LineageState {HAS_CONTENT, DELETED, NEVER_HAD_CONTENT }
    private static int ROOT_IDX = 0;

    public static void main(String[] args) {


        parseArgs(args);
        checkArgsValid();

        if (LEARN) {
            EnumSeq.Alignment aln = null;
            Tree tree = null;
            try {
                aln = Utils.loadAlignment(ALIGNMENT, EVOL_MODEL.getDomain());
                tree = Newick.load(INPUT_TREE);
                Utils.checkData(aln, tree, true);

            } catch (IOException | ASRException e) {
                usage(10, e.getMessage());
            }

            learnTreeParams(tree, 3, SEED);

            assert aln != null;
            printRootSeq(aln, tree, null);
            System.out.println("--substitution-model " + EVOL_MODEL.getName() + " \\");

            learnIndelLengthDistributions(tree, aln, SEED, null);
            learnIndelRateDistribution(tree, aln, null, SEED);

        } else {
            // Actually perform simulation
            EnumSeq rootSeq = createRootSeq(ANCSEQ, EVOL_MODEL, ANCSEQ_LENGTH, SEED);
            IdxTree tree = setupTree(
                                    INPUT_TREE,
                                    COPY_TREE,
                                    TREE_DISTANCE_MODEL,
                                    LEAF2ROOT_DISTANCE_MODEL,
                                    EXTANTS_N,
                                    SEED,
                                    TREE_GAMMA_SHAPE,
                                    TREE_GAMMA_SCALE,
                                    DESCENDANTS_MIN,
                                    DESCENDANTS_MAX,
                                    SCALEDIST);

            TrackTree.Params params = setupParams(
                                                tree,
                                                rootSeq,
                                                INDEL_LENGTH_MODEL,
                                                INSERTION_LENGTH_MODEL,
                                                DELETION_LENGTH_MODEL,
                                                INDEL_RATE_MODEL,
                                                SUBST_RATE_MODEL,
                                                EVOL_MODEL,
                                                DELETIONPROP,
                                                null,
                                                SEED);

            TrackTree tracker = new TrackTree(params, SEED);

            EnumSeq[] seqs = tracker.getSequences();

            saveOutput(seqs, tracker, tree, FORMAT_IDX, GAPPY, EXTANTS_ONLY, OUTPUT, ALIGNMENT);
        }
    }

    public static void learnIndelRateDistribution(IdxTree tree, EnumSeq.Alignment<Enumerable> aln, Object[][] ancseqs_gappy, long seed) {
        double[] rateSampleCollection = calculateColumnIndelRates(tree, aln, ancseqs_gappy);
        RateModel indelrateDist = RateModel.bestfit(rateSampleCollection, seed);
        if (indelrateDist != null) {
            System.out.println("--indel-rate-distrib " + indelrateDist.getTrAVIS());
        }
    }


    /**
     * To simulate phylogenetic trees we learn two distributions which sufficiently capture a description of the tree
     * which allows us to generate similar trees. The first is the distribution of branch distances, which we learn as
     * a gamma mixture distribution. The second is the distribution of leaf-to-root distances, which we
     * learn as a Gaussian distribution. The branch distance distribution is used to generate trees with similar branch
     * lengths, and the leaf-to-root distribution is used to shuffle the tree to have a similar distribution of
     * leaf-to-root distances as the original tree.
     *
     * @param tree The tree to learn from
     * @param nComponents The number of components to use in the gamma mixture distribution for branch lengths.
     * @param seed Random seed for reproducibility.
     */
    public static void learnTreeParams(IdxTree tree, int nComponents, long seed) {

        RateModel ddistrib = IdxTree.getGammaMixture(tree, nComponents, seed);
        System.out.println("--dist-distrib " + ddistrib.getTrAVIS() + " \\");
        GaussianDistrib l2rdistrib = tree.getLeaf2RootDistrib();
        IdxTree newrtree = Tree.Random(tree.getNLeaves(), ddistrib, 2, 2, seed);
        newrtree.fitDistances(100, l2rdistrib, seed + 202);
        System.out.println("--leaf2root-distrib " + l2rdistrib.getTrAVIS() + " \\");

    }

    public static void printRootSeq(EnumSeq.Alignment<Enumerable> aln, IdxTree tree,
                                    Object[][] ancseqs_nogap) {

        boolean allSeqsInAln = ancseqs_nogap == null;
        StringBuilder n0 = new StringBuilder();

        if (allSeqsInAln) {

            Map<String, Integer> idToAlnIndex = aln.getMap();

            Object[] n0Gapped = aln.getEnumSeq(idToAlnIndex.get(tree.getRoot().getLabel())).get();

            for (Object o : n0Gapped) {
                if (o != null) {
                    n0.append(o);
                }
            }

        } else {
            for (int j = 0; j < ancseqs_nogap[0].length; j++)
                n0.append(ancseqs_nogap[0][j]);
        }



        System.out.println("--ancestor " + n0 + " \\");
        System.out.println("--length " + n0.length() + " (Use instead of --ancestor to generate a random sequence according to the substitution model) \\");
    }


    public static String[] parseDistribParamString(String params) {

        int colonPos = params.indexOf(':');

        String distName = colonPos >= 0 ? params.substring(0, colonPos) : params;
        String parsedParams = colonPos >= 0 ? params.substring(colonPos + 1) : "";

        return new String[]{distName, parsedParams};
    }

    public static void learnIndelLengthDistributions(IdxTree tree, EnumSeq.Alignment<Enumerable> aln,
                                                      long seed, Object[][] ancseqs_gappy) {

        Map<Integer, Object[]> seqs = getAllSeqs(tree, aln, ancseqs_gappy);
        int[] ins_total = new int[0];
        int[] del_total = new int[0];

        // Go through the tree, and look at each ancestor sequence, recording predicted indel events
        for (int idx : tree) { // go through the original, user-provided tree
            int parent = tree.getParent(idx);
            if (parent != -1) {  // Non-root node, so there is a branch with distance to catch...

                // Retrieve parent/child reconstructed sequences at a node in the user-provided tree
                Object[] pseq = seqs.get(parent);// aln.getEnumSeq(idToAlnIndex.get(tree.getBranchPoint(parent).getLabel().toString())).get();
                Object[] cseq = seqs.get(idx);

                int[] insertions = getInsertionCounts(pseq, cseq);
                int[] deletions = getDeletionCounts(pseq, cseq);

                // Accumulate insertion/deletion lengths
                ins_total = mergeCounts(ins_total, insertions);
                del_total = mergeCounts(del_total, deletions);
            }
        }

        // now, turn to indel lengths...
        int[] indel_total = TrAVIS.mergeCounts(ins_total, del_total);
        int[] ins_data = TrAVIS.unfoldCounts(ins_total);
        int[] del_data = TrAVIS.unfoldCounts(del_total);
        int[] indel_data = TrAVIS.unfoldCounts(indel_total);

        // Now we can fit the indel distribution to the lengths of insertions and deletions
        // we need to try all and pick the one with greatest log-likelihood
        IndelModel indel_length_distrib = IndelModel.bestfit(indel_data, seed);
        IndelModel insertion_length_distrib = IndelModel.bestfit(ins_data, seed);
        IndelModel deletion_length_distrib = IndelModel.bestfit(del_data, seed);

        System.out.println("--indel-length-distrib " + indel_length_distrib.getTrAVIS() + " \\");
        System.out.println("--insertion-length-distrib " + insertion_length_distrib.getTrAVIS() + " \\");
        System.out.println("--deletion-length-distrib " + deletion_length_distrib.getTrAVIS() + " \\");

        int ninsertions = Arrays.stream(ins_total).sum();
        int ndeletions = Arrays.stream(del_total).sum();
        int nindel = ninsertions + ndeletions;
        double delprop = (double) ndeletions / (double) nindel;
        System.out.printf("--delprop %.2f \\\n", delprop);
    }

    private static Map<Integer, Object[]> getAllSeqs(IdxTree tree, EnumSeq.Alignment<Enumerable> aln,
                                                     Object[][] ancseqs_gappy) {

        Iterator<Integer> dfs = tree.getDepthFirstIterator();
        Map<String, Integer> alnMap = aln.getMap();

        // if ancseqs_gappy is null, we are working with an alignment containing both ancestors and extants
        boolean allSeqsInAln = ancseqs_gappy == null;

        Map<Integer, Object[]> seqs = new HashMap<>();

        while (dfs.hasNext()) {
            int bpidx = dfs.next();

            Object[] currentSeq;

            if (allSeqsInAln) {
                currentSeq = aln.getEnumSeq(alnMap.get(tree.getBranchPoint(bpidx).getLabel().toString())).get();
            } else {
                if (tree.isLeaf(bpidx)) {
                    EnumSeq.Gappy seq = aln.getEnumSeq(alnMap.get(tree.getLabel(bpidx)));
                    currentSeq = seq.get();
                } else {
                    currentSeq = ancseqs_gappy[(Integer) tree.getLabel(bpidx)];
                }
            }

            seqs.put(bpidx, currentSeq);
        }

        return seqs;
    }

    private static void processNodeForRate(int bpidx, Map<Integer, LineageState> lineageState,
                                           Map<Integer, Integer> numNodesTraversedSinceIndel,
                                           Map<Integer, Double> distTraversedSinceIndel, IdxTree tree,
                                           List<Double> rateSampleCollection, boolean currentNodeHasContent) {


        int parentIdx = tree.getParent(bpidx);

        // need to track how many nodes since indel relative to the parent
        int nodesParentTraversed = numNodesTraversedSinceIndel.get(parentIdx);
        numNodesTraversedSinceIndel.put(bpidx, nodesParentTraversed + 1);

        // same idea for distance traversed
        double distTraversedParent = distTraversedSinceIndel.get(parentIdx);
        distTraversedSinceIndel.put(bpidx, distTraversedParent + tree.getDistance(bpidx));

        // now check the state of our parent
        LineageState parentState = lineageState.get(parentIdx);
        // two possible scenarios:
        // 1) child has content: if the parent was deleted or never had content, this is an insertion.
        // TODO - deletion in parent followed by insertion is technically a violation, potentially should stop recording indels below this node
        // 2) Child does NOT have content; if parent had content we've identified a deletion.
        boolean indelEventOccurred = ((parentState == LineageState.DELETED || parentState == LineageState.NEVER_HAD_CONTENT) && currentNodeHasContent) ||
                (parentState == LineageState.HAS_CONTENT && !currentNodeHasContent);

        if (indelEventOccurred) {
            int localNodesTraversed = numNodesTraversedSinceIndel.get(bpidx);
            double localDistTraversed = distTraversedSinceIndel.get(bpidx);
            // there is 1 indel event after we traverse a certain number of nodes
            double indelRate = -Math.log(1.0 - ((double) 1 / localNodesTraversed)) / localDistTraversed;

            rateSampleCollection.add(indelRate);
            // Except for the last node, we had no indel events, which we mark as a non-event.
            for (int x = 0; x < localNodesTraversed - 1; x++) {
                rateSampleCollection.add(0.0);
            }

            // reset all the counts
            numNodesTraversedSinceIndel.put(bpidx, 1);
            distTraversedSinceIndel.put(bpidx, 0.0);
        }

        // bookkeeping so we can identify indel events.
        LineageState currentState;
        if (currentNodeHasContent) {
            currentState = LineageState.HAS_CONTENT;
        } else if (parentState == LineageState.HAS_CONTENT) {
            currentState = LineageState.DELETED;
        } else if (parentState == LineageState.DELETED) {
            currentState = LineageState.DELETED;
        } else {
            currentState = LineageState.NEVER_HAD_CONTENT;
        }

        lineageState.put(bpidx, currentState);
    }

    public static double[] calculateColumnIndelRates(IdxTree tree, EnumSeq.Alignment<Enumerable> aln,
                                                     Object[][] ancseqs_gappy) {

        List<Double> rateSampleCollection = new ArrayList<>();
        Map<Integer, Object[]> seqs = getAllSeqs(tree, aln, ancseqs_gappy);
        for (int alnPos = 0; alnPos < aln.getWidth(); alnPos++) {

            Iterator<Integer> dfs = tree.getDepthFirstIterator();

            Map<Integer, LineageState> lineageState = new HashMap<>();
            Map<Integer, Integer> numNodesTraversedSinceIndel = new HashMap<>();
            Map<Integer, Double> distTraversedSinceIndel = new HashMap<>();

            while (dfs.hasNext()) {
                int bpidx = dfs.next();
                // grab the current node sequence
                Object[] currentSeq = seqs.get(bpidx);

                boolean currentNodeHasContent = currentSeq[alnPos] != null;
                if (bpidx == 0) {
                    // special conditions for the root
                    lineageState.put(bpidx, currentNodeHasContent ? LineageState.HAS_CONTENT : LineageState.NEVER_HAD_CONTENT);
                    numNodesTraversedSinceIndel.put(bpidx, 1); // start the count
                    distTraversedSinceIndel.put(bpidx, 0.0);
                    continue;
                }

                processNodeForRate(bpidx, lineageState, numNodesTraversedSinceIndel, distTraversedSinceIndel,
                        tree, rateSampleCollection, currentNodeHasContent);

            }
        }

        double[] colRateArray = new double[rateSampleCollection.size()];
        for (int jj = 0; jj < rateSampleCollection.size(); jj++)
            colRateArray[jj] = rateSampleCollection.get(jj);

        return colRateArray;
    }

    private static void parseArgs(String[] args) {

        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("n0") || arg.equalsIgnoreCase("-ancestor") && args.length > a + 1) {
                    ANCSEQ = args[++a];
                } else if ((arg.equalsIgnoreCase("-aln") || arg.equalsIgnoreCase("a")) && args.length > a + 1) {
                        ALIGNMENT = args[++ a];
                } else if (arg.equalsIgnoreCase("-nwk")  || arg.equalsIgnoreCase("n") && args.length > a + 1) {
                    INPUT_TREE = args[++a];
                } else if (arg.equalsIgnoreCase("o") || arg.equalsIgnoreCase("-output-folder") && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if (arg.equalsIgnoreCase("-seed") && args.length > a + 1) {
                    SEED = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-extants") && args.length > a + 1) {
                    EXTANTS_N = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-dgamma") && args.length > a + 2) {
                    TREE_GAMMA_SHAPE = Double.parseDouble(args[a+1]);
                    TREE_GAMMA_SCALE = Double.parseDouble(args[a+2]);
                    if (args.length > a + 3 && !args[a + 3].startsWith("-")) {
                        SCALEDIST = Double.parseDouble(args[a+3]); // brdist_scale
                    }
                } else if (arg.equalsIgnoreCase("-gap")) {
                    GAPPY = true;
                } else if (arg.equalsIgnoreCase("-learn")) {
                    LEARN = true;
                } else if (arg.equalsIgnoreCase("-verbose")) {
                    VERBOSE = true;
                } else if (arg.equalsIgnoreCase("-help") || arg.equalsIgnoreCase("h")) {
                    usage();
                } else if ((arg.equalsIgnoreCase("s") || arg.equalsIgnoreCase("-substitution-model")) && args.length > a + 1) {
                    EVOL_MODEL = SubstModel.createModel(args[a+1]);
                    if (EVOL_MODEL == null)
                        usage(1, args[a + 1] + " is not a valid model name");
                } else if ((arg.equalsIgnoreCase("rf") && args.length > a + 1)) {
                    SRATESFILE = args[++a];
                } else if (arg.equalsIgnoreCase("-subst-rate-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a+1]);
                    SUBST_RATE_MODEL = RateModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                } else if (arg.equalsIgnoreCase("-indel-rate-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a+1]);
                    INDEL_RATE_MODEL = RateModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                } else if (arg.equalsIgnoreCase("-indel-length-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a+1]);
                    INDEL_LENGTH_MODEL = IndelModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                } else if (arg.equalsIgnoreCase("-insertion-length-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a+1]);
                    INSERTION_LENGTH_MODEL = IndelModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                } else if (arg.equalsIgnoreCase("-deletion-length-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a+1]);
                    DELETION_LENGTH_MODEL = IndelModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                } else if (arg.equalsIgnoreCase("-delprop") && args.length > a + 1) {
                    DELETIONPROP = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("-length") || arg.equalsIgnoreCase("l") && args.length > a + 1) {
                    ANCSEQ_LENGTH = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-dist-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a+1]);
                    TREE_DISTANCE_MODEL = RateModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                } else if (arg.equalsIgnoreCase("-leaf2root-distrib") && args.length > a + 1) {
                    String[] params = parseDistribParamString(args[a + 1]);
                    LEAF2ROOT_DISTANCE_MODEL = Distrib.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS]);
                } else if (arg.equalsIgnoreCase("-copy-tree") && args.length > a + 1) {
                    COPY_TREE = true;
                } else if (arg.equalsIgnoreCase("-extants-only") && args.length > a + 1) {
                    EXTANTS_ONLY = true;
                } else if (arg.equalsIgnoreCase("-save-as") || arg.equalsIgnoreCase("sa") && args.length > a + 1) {
                    boolean found_format = false;
                    for (int i = 0; i < TRAVIS_FORMATS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(TRAVIS_FORMATS[i])) {
                            FORMAT_IDX = i;
                            found_format = true;
                        }
                    }
                    if (!found_format)
                        usage(1, args[a + 1] + " is not a valid format name");
                } else if (arg.equalsIgnoreCase("-help") || arg.equalsIgnoreCase("h")) {
                    usage();
                }
            }
        }

    }

    private static void checkArgsValid() {

        if (OUTPUT == null) {
            usage(25, "Output file or folder must be specified");
        }

        if (SRATESFILE != null) {
            try {
                SRATES = TSVFile.loadSubstitutionRatesFile(SRATESFILE);
            } catch (IOException e) {
                usage(24, e.getMessage());
            } catch (NumberFormatException e) {
                usage(23, e.getMessage());
            }

        }

        if (SRATES != null) { // position-specific rates available
            GAMMA_ALPHA = calcAlpha(SRATES);
        }

        EVOL_MODEL = SubstModel.createModel(EVOL_MODELS[EVOL_MODEL_IDX]);
        if (EVOL_MODEL == null) {
            usage(1, "Model " + EVOL_MODELS[EVOL_MODEL_IDX] + " could not be created");
        }

        if (TRAVIS_FORMATS[FORMAT_IDX].equalsIgnoreCase("CLUSTAL")) // Clustal files can only be "gappy"
            GAPPY = true;

    }

    public static IdxTree setupTree(String newickFile, boolean copyTree, RateModel treeDistanceModel,
                                     Distrib leaf2rootDistanceModel, int numExtants, long seed,
                                     double treeGammaShape, double treeGammaScale, int descendantsMin,
                                     int descendantsMax, Double scaleDist) {

        IdxTree tree = null;

        // load in a user provided tree
        if (newickFile != null) {
            try {
                tree = Newick.load(newickFile);

                if (copyTree) {
                    return tree;
                }

            } catch (IOException e) {
                usage(26, "Input tree file " + newickFile + " is invalid.");
            }
        }


        if (tree == null) {
            // no tree provided but parameters given
            if (treeDistanceModel != null)
                tree = IdxTree.generateTreeFromDistrib(treeDistanceModel, leaf2rootDistanceModel, numExtants, seed, 100);
            else {
                // otherwise create a completely random tree
                tree = Tree.Random(numExtants, seed, treeGammaShape, 1.0 / treeGammaScale, descendantsMax, descendantsMin);
                if (leaf2rootDistanceModel != null)
                    tree = IdxTree.shuffleWithLeaf2RootDistrib(tree, leaf2rootDistanceModel, seed, 100);
            }
            if (scaleDist != null)
                tree.adjustDistances(scaleDist);
        } else {
            tree = Tree.generateTreeFromMixture(tree,3, seed,100);
        }

        return tree;
    }

    /**
     * Create the root sequence for the simulation, either by parsing a user-provided sequence or by
     * generating a random sequence according to the model's equilibrium distribution.
     *
     * @return The root sequence as an EnumSeq object.
     */
    public static EnumSeq createRootSeq(String rootString, SubstModel model,
                                         Integer rootSeqLength, long seed) {

        EnumSeq ancseq = null;

        // sequence provided, try parse it according to the model
        if (rootString != null) {
            if (model.getDomain().equals(Enumerable.aacid)) {
                ancseq = EnumSeq.parseProtein(rootString);
            } else if (model.getDomain().equals(Enumerable.nacid)) {
                ancseq = EnumSeq.parseDNA(rootString);
            } else if (model.getDomain().equals(Enumerable.nacidRNA)) {
                ancseq = EnumSeq.parseRNA(rootString);
            } else {
                usage(5, "Model \"" + model + "\" alphabet is not valid");
            }
        }

        // user needs to provide a sequence length otherwise
        if (rootString == null && rootSeqLength != null) {

            Enumerable domain = model.getDomain();
            Random rand = new Random(seed);

            // sample from the model's equilibrium distribution to create a random sequence of the specified length
            String[] ancSeq = new String[rootSeqLength];
            for (int i = 0; i < rootSeqLength; i++) {
                Object nchar = null;
                double sump = 0;
                double tossAgain = rand.nextDouble();
                for (Object c : domain.getValues()) {
                    sump += model.getProb(c);
                    if (sump >= tossAgain) {
                        nchar = c;
                        break;
                    }
                }
                assert nchar != null;
                ancSeq[i] = nchar.toString();
            }

            String finalAncSeq = String.join("", ancSeq);
            if (model.getDomain().equals(Enumerable.aacid)) {
                ancseq = EnumSeq.parseProtein(finalAncSeq);
            } else if (model.getDomain().equals(Enumerable.nacid)) {
                ancseq = EnumSeq.parseDNA(finalAncSeq);
            } else if (model.getDomain().equals(Enumerable.nacidRNA)) {
                ancseq = EnumSeq.parseRNA(finalAncSeq);
            }
        }

        if (ancseq == null) {
            usage(4, "Invalid ancestor sequence \"" + rootString + "\" for model " + EVOL_MODELS[EVOL_MODEL_IDX]);
        }
        ancseq.setName("N0");

        return ancseq;
    }

    /**
     * Build the Params object to be used for the TrackTree simulation, which includes all the necessary information
     * about the tree, root sequence, and evolutionary models.
     *
     * @param tree The phylogenetic tree to simulate along.
     * @param rootSeq The root sequence from which to start the simulation.
     * @return A TrackTree.Params object containing all the parameters for the simulation.
     */
    public static TrackTree.Params setupParams(IdxTree tree, EnumSeq rootSeq,
                                               IndelModel indelLengthModel,
                                               IndelModel insertionLengthModel,
                                               IndelModel deletionLengthModel,
                                               RateModel indelRateModel,
                                               RateModel substRateModel,
                                               SubstModel model,
                                               double deletionProportion,
                                               double[] substRates,
                                               long seed) {

        // we've got an ancestor to track down the tree
        TrackTree.Params params = new TrackTree.Params(tree, rootSeq, model, seed);

        if (indelLengthModel == null) {
            params.setIndelModel(IndelModel.create("Zipf", "1.7,50", seed));
        } else {
            params.setIndelModel(indelLengthModel);
        }

        if (substRates != null) {
            params.setSubstRates(substRates);
        }

        if (insertionLengthModel != null) {
            params.setInsertmodel(insertionLengthModel);
        }
        if (deletionLengthModel != null) {
            params.setDeletemodel(deletionLengthModel);
        }

        if (indelRateModel != null)
            params.setIndelRateModel(indelRateModel);

        if (substRateModel != null)
            params.setSubstRateModel(substRateModel);

        params.PROPORTION_DELETION = deletionProportion;

        params.setSeed(seed);

        return params;
    }

    public static void saveOutput(EnumSeq[] seqs, TrackTree tracker, IdxTree tree,
                                  int formatIdx, boolean gappy, boolean extantsOnly, String outputDir,
                                  String prefix) {

        File file = new File(outputDir);
        file.mkdirs();// true if the directory was created, false otherwise

        switch (formatIdx) {
            case FASTA: // FASTA
                try {
                    FastaWriter fw = new FastaWriter(new File(outputDir, prefix +  "_travis.fa"));
                    if (!gappy) {
                        fw.save(seqs);
                    } else { // gappy
                        EnumSeq[] aln = tracker.getAlignment();
                        fw.save(aln);
                    }
                    fw.close();

                    if (extantsOnly) {
                        EnumSeq[] aln = tracker.getAlignment();
                        EnumSeq[] extants = new EnumSeq[tree.getNLeaves()];
                        int count = 0;
                        for (EnumSeq e : aln) {
                            if (e.getName().startsWith("A")) {
                                extants[count] = e;
                                count++;
                            }
                        }
                        FastaWriter fwExtants = new FastaWriter(new File(outputDir, prefix + "_extants_travis.fa"));
                        fwExtants.save(extants);
                        fwExtants.close();
                    }
                } catch (IOException e) {
                    usage(6, "FASTA file could not be saved");
                }

                try {
                    Newick.save(tree, outputDir + "/" + prefix + "_travis.nwk", Newick.MODE_DEFAULT);
                } catch (IOException e) {
                    usage(2, "Tree file could not be saved");
                }



                break;
            case DOT: // DOT
                POAGraph poag = tracker.getPOAG();
                try {
                    poag.saveToDOT(outputDir + "/" + prefix + "_travis.dot");
                } catch (IOException e) {
                    usage(6, "DOT file could not be saved");
                }
                break;
            case CLUSTAL: // CLUSTAL
                EnumSeq[] aln = tracker.getAlignment();
                try {
                    AlnWriter aw = new AlnWriter(new File(outputDir, prefix + "_travis.aln"));
                    aw.save(aln);
                    aw.close();
                } catch (IOException e) {
                    usage(6, "CLUSTAL file could not be saved");
                }
                break;
            case ALL: // ALL in a DIRECTORY
                POAGraph poaGraph = tracker.getPOAG();
                EnumSeq[] aseqs = tracker.getAlignment();
                try {

                    FastaWriter fw = new FastaWriter(new File(outputDir,  prefix + "_travis.fa"));
                    if (!gappy) {
                        fw.save(seqs);
                    } else { // gappy
                        fw.save(aseqs);
                    }
                    fw.close();
                    AlnWriter aw = new AlnWriter(new File(outputDir, prefix +"_travis.aln"));
                    aw.save(aseqs);
                    aw.close();
                    poaGraph.saveToDOT(outputDir + "/" + prefix + "_travis.dot");
                    poaGraph.saveToMatrix(outputDir + "/" + prefix + "_travis.m");
                    Newick.save(tree, outputDir + "/" + prefix + "_travis.nwk", Newick.MODE_DEFAULT);

                    EnumSeq[] alnTracker = tracker.getAlignment();
                    EnumSeq[] extants = new EnumSeq[tree.getNLeaves()];
                    int count = 0;
                    for (EnumSeq e : alnTracker) {
                        if (e.getName().startsWith("A")) {
                            extants[count] = e;
                            count++;
                        }
                    }
                    FastaWriter fwExtants = new FastaWriter(new File(outputDir, prefix + "_extants_travis.fa"));
                    fwExtants.save(extants);
                    fwExtants.close();

                    if (tracker.INDELRATES) {
                        double[] indelRates = tracker.getIndelRates();
                        Object[][] data = new Object[indelRates.length + 1][2];
                        for (int i = 0; i <= indelRates.length; i++) {
                            if (i == 0) // header
                                data[0] = new Object[]{"Site", "Rate"};
                            else
                                data[i] = new Object[]{i, indelRates[i - 1]};
                        }
                        TSVFile ratesfile = new TSVFile(data, true);
                        ratesfile.save(outputDir + "/" + prefix + "_indel_rates_travis.tsv");
                    }

                } catch (IOException e) {
                    usage(7, "Something went wrong saving files in directory");
                }
                break;
            case RATES: // RATES
                if (tracker.INDELRATES) {
                    double[] indelRates = tracker.getIndelRates();
                    Object[][] data = new Object[indelRates.length + 1][2];
                    for (int i = 0; i <= indelRates.length; i++) {
                        if (i == 0) // header
                            data[0] = new Object[]{"Site", "Rate"};
                        else
                            data[i] = new Object[]{i, indelRates[i - 1]};
                    }
                    try {
                        TSVFile ratesfile = new TSVFile(data, true);
                        ratesfile.save(outputDir + "/" + prefix + "_indel_rates_travis.tsv");
                    } catch (IOException e) {
                        usage(6, "Indel rates file could not be saved");
                    }
                }
        }
    }


    /**
     * Calculates the counts of insertions between two sequences.
     *
     * This method compares two sequences (padded to be aligned) and determines the number of insertions
     * at each position. It returns an array where the index represents the length
     * of the insertion and the value at that index represents the count of such
     * insertions.
     *
     * @param seq1 the first sequence to compare (parent sequence)
     * @param seq2 the second sequence to compare (child sequence)
     * @return an array where the index represents the length of the insertion and
     *         the value at that index represents the count of such insertions, or
     *         null if the sequences are of different lengths
     */
    static int[] getInsertionCounts(Object[] seq1, Object[] seq2) {
        if (seq1.length != seq2.length)
            return null;
        Map<Integer, Integer> ins = new HashMap<>(); // map with insert length as key and count as value
        boolean seq2ins = false;
        int seq2cnt = 0;
        for (int i = 0; i < seq1.length; i++) {
            if (seq1[i] == null && seq2[i] == null)
                continue;
            if (seq1[i] == null && seq2[i] != null) {
                seq2cnt++;
                seq2ins = true;
            } else {
                if (seq2ins) {
                    if (ins.containsKey(seq2cnt))
                        ins.put(seq2cnt, ins.get(seq2cnt) + 1);
                    else
                        ins.put(seq2cnt, 1);
                    seq2ins = false;
                    seq2cnt = 0;
                }
            }
        }
        int max = 0;
        for (int cnt : ins.keySet()) {
            if (cnt > max)
                max = cnt;
        }
        int[] ret = new int[max];
        for (int i = 0; i < max; i++) {
            if (ins.containsKey(i + 1))
                ret[i] = ins.get(i + 1);
            else
                ret[i] = 0;
        }
        return ret;
    }

    /**
     * Calculates the counts of deletions between two sequences.
     *
     * This method compares two sequences (padded to be aligned) and determines the number of deletions
     * at each position. It returns an array where the index represents the length
     * of the deletion and the value at that index represents the count of such
     * deletions.
     *
     * @param seq1 the first sequence to compare (parent sequence)
     * @param seq2 the second sequence to compare (child sequence)
     * @return an array where the index represents the length of the deletion and
     *         the value at that index represents the count of such deletions, or
     *         null if the sequences are of different lengths
     */
    static int[] getDeletionCounts(Object[] seq1, Object[] seq2) {
        if (seq1.length != seq2.length)
            return null;
        Map<Integer, Integer> del = new HashMap<>(); // map with insert length as key and count as value
        boolean seq2del = false;
        int seq2len = 0; // length of current deletion
        for (int i = 0; i < seq1.length; i++) {
            if (seq1[i] == null && seq2[i] == null) // both parent and child are gaps, so nothing changes
                continue;
            if (seq2[i] == null && seq1[i] != null) { // parent has content, but child has gap so start/continue deletion
                seq2len ++;
                seq2del = true; // start/continue current deletion
            } else { // parent is gap, child has content, or both have content; either way, we're ending deletion (if current)
                if (seq2del) {
                    if (del.containsKey(seq2len))
                        del.put(seq2len, del.get(seq2len) + 1);
                    else
                        del.put(seq2len, 1);
                    seq2del = false;
                    seq2len = 0;
                }
            }
        }
        int max = 0;
        for (int cnt : del.keySet()) {
            if (cnt > max)
                max = cnt;
        }
        int[] ret = new int[max];
        for (int i = 0; i < max; i++) {
            if (del.containsKey(i + 1))
                ret[i] = del.get(i + 1);
            else
                ret[i] = 0;
        }
        return ret;
    }

    /**
     * Merge two count arrays as produced by the getDeletionCounts, and getInsertionCounts.
     * @param cnt1
     * @param cnt2
     * @return the merged counts
     */
    public static int[] mergeCounts(int[] cnt1, int[] cnt2) {
        int[] tmp = new int[Math.max(cnt1.length, cnt2.length)];
        for (int j = 0; j < tmp.length; j++) {
            tmp[j] += j < cnt1.length ? cnt1[j] : 0;
            tmp[j] += j < cnt2.length ? cnt2[j] : 0;
        }
        return tmp;
    }

    /**
     * Take the array of counts and make the data points represented by their index
     * (as produced by the getDeletionCounts, and getInsertionCounts).
     * @param cnts
     * @return dataset
     */
    public static int[] unfoldCounts(int[] cnts) {
        int n = Arrays.stream(cnts).sum();
        int[] ret = new int[n];
        int j = 0;
        for (int i = 0; i < cnts.length; i++) {
            for (int k = 0; k < cnts[i]; k ++) {
                ret[j ++] = (i + 1);
            }
        }
        return ret;
    }

    /**
     * Calculate indel openings, matches, and mismatches between a parent and a child sequence.
     *
     * @param seq1 The sequence of the parent node.
     * @param seq2 The sequence of the child node.
     * @return An array of three integers: indel openings, matches, and mismatches.
     */
    public static int[] calculateIndelOpening(Object[] seq1, Object[] seq2) {
        int indelOpenings = 0;
        int matches = 0;
        int mismatches = 0;

        int length = Math.min(seq1.length, seq2.length);
        boolean inIndel = false;

        for (int i = 0; i < length; i++) {
            boolean isGap = (seq1[i] == null && seq2[i] != null) || (seq2[i] == null && seq1[i] != null);

            if (isGap) {
                if (!inIndel) {
                    indelOpenings++;  // New indel opening
                    inIndel = true;
                }
            } else {
                inIndel = false;
                if (seq1[i] != null && seq2[i] != null) {
                    if (seq1[i].equals(seq2[i])) {
                        matches++;
                    } else {
                        mismatches++;
                    }
                }
            }
        }

        return new int[]{indelOpenings, matches, mismatches};
    }

    /**
     * Class to track ancestor sequence to extants via intermediate ancestors.
     * Matches/substitutions are determined by a probability p=exp^-rt where rt is the rate times the evolutionary distance from the ancestor to the descendant.
     * If not a match/substitution, insertions and deletions are equally probable, i.e. (1-p)/2 each.
     * The length of an insertion or deletion is determined by a Poisson with mean (lambda) 1; note that this means that 0.37 of indels are length 0.
     * The implementation is inspired by rules extracted from.
     * Position specific rates can be supplied to the constructor.
     *  <a href="https://doi.org/10.1093/molbev/msn275">Cartwright R. Problems and Solutions for Estimating Indel Rates and Length Distributions.
     *  Mol. Biol. Evol. 26(2):473–480. 2009.</a>
     */
    static class TrackTree {

        static class Params {
            public final IdxTree tree;              // the tree that is being used for generating sequences
            public final EnumSeq ancseq;            // the ultimate ancestor sequence
            public final SubstModel substmodel;     // evolutionary model for substitution
            public double[] ancrates = null;        // rates for ancestor to override randomly set rates for the corresponding position

            public RateModel substratemodel = null; // distribution which specifies site-specific rates that modulates substitution at each site/position
            public RateModel indelratemodel = null; // distribution which specifies site-specific rates that modulates indel events at each site/position

            public IndelModel insertmodel = null;   // distribution from which insertion lengths are sampled
            public IndelModel deletemodel = null;   // distribution from which deletion lengths are sampled

            public double PROPORTION_DELETION = 0.5;   // the proportion of deletion events (as opposed to insertion) amongst all indels
            public Random rand = null;

            public Params(IdxTree tree, EnumSeq ancseq, SubstModel substmodel, long SEED) {
                this.tree = tree;
                this.ancseq = ancseq;
                this.substmodel = substmodel;
                setSeed(SEED);
            }

            public static Params seed(Params params, long SEED) {
                params.setSeed(SEED);
                return params;
            }

            public void setInsertmodel(IndelModel insertmodel) {
                this.insertmodel = insertmodel;
            }

            public void setDeletemodel(IndelModel deletemodel) {
                this.deletemodel = deletemodel;
            }

            public void setIndelModel(IndelModel indelmodel) {
                this.insertmodel = indelmodel;
                this.deletemodel = indelmodel;
            }

            /**
             * Set the distribution that defines (variable) rates that modulate substitution at each site/position.
             * By default this is the Gamma distribution, with a mean of 1.
             * @param ratesgamma alpha and shape parameter of gamma, which is also used to specify 1/beta and scale, to ensure the mean is 1
             */
            public void setSubstRateModel(double ratesgamma) {
                this.substratemodel = new GammaDistrib(ratesgamma, ratesgamma);
            }

            /**
             * Set the distribution that defines (variable) rates that modulate substitution at each site/position.
             * Standard distribution is Gamma.
             */
            public void setSubstRateModel(RateModel substratemodel) {
                this.substratemodel = substratemodel;
            }

            /**
             * Set the substitution rates for the columns that originate in the specified ancestor sequence (at the root)
             * @param rates evolutionary rates
             */
            public void setSubstRates(double[] rates) {
                if (ancseq.length() == rates.length) {
                    this.ancrates = rates;
                }
            }

            /**
             * Set the distribution that specifies node-specific rates that modulates indel events in a sequence (NOT site specific).
             * This setter assumes that the classical Gamma distribution is used.
             * @param shape parameter of Gamma
             * @param scale parameter of Gamma
             */
            public void setIndelRateModel(double shape, double scale) {
                this.indelratemodel = new GammaDistrib(shape, scale);
            }
            /**
             * Set the distribution that specifies node-specific rates that modulates indel events in a sequence (NOT site specific).
             * For example, GammaDistrib or ZeroInflatedGamma
             */
            public void setIndelRateModel(RateModel indelratemodel) {
                this.indelratemodel = indelratemodel;
            }

            public void setSeed(long SEED) {
                this.rand = new Random(SEED);
                if (insertmodel != null) insertmodel.setSeed(SEED);
                if (deletemodel != null) deletemodel.setSeed(SEED);
                if (substratemodel != null) substratemodel.setSeed(SEED);
                if (indelratemodel != null) indelratemodel.setSeed(SEED);

            }

            public Random getRandom() {
                return rand;
            }
        }

        Enumerable myType = null;

        public final Params params;
        private TreeInstance ti_seqs;
        private EnumNode[][] enumNodes = null;
        private TreeInstance ti_deletions = null;
        private TreeInstance ti_insertions = null;
        private Random rand = null;
        private POAGraph poag = null;
        private EnumSeq.Gappy[] alignedseqs = null; // alignment extracted from sequences and POAG; requires POAG to have been generated
        private int[][] alignedidxs = null;         // alignment extracted from sequences and POAG, but in the form of indices from original seq idx to alignment idx
        // private int[] order = null; // topological order of nodes in POAG; when set alignment can be extracted
        private double[] alignedrates = null;
        private double[] alignedIndelRates = null;
        private double[][] substRates; // rates in tree, reference to branchpoint specific index
        private double[][] colIndelRates;
        public boolean USERATES;
        public boolean INDELRATES;

        /**
         * @param params
         */
        public TrackTree(Params params) {
            this(params, System.currentTimeMillis());
        }

        /**
         * @param params
         * @param SEED
         */
        public TrackTree(TrackTree.Params params, long SEED) {
            this.params = params;
            this.params.setSeed(SEED);
            EnumSeq ancseq = params.ancseq;

            USERATES = (params.substratemodel != null); // check if we will generate position specific rates using the model; if not, use a constant rate
            INDELRATES = (params.indelratemodel != null);

            myType = ancseq.getType();
            IdxTree tree = params.tree;
            int[][] deletions  = new int[tree.getSize()][];
            int[][] insertions = new int[tree.getSize()][];
            substRates = new double[tree.getSize()][];
            colIndelRates = new double[tree.getSize()][];
            EnumSeq[] bpseqs = new EnumSeq[tree.getSize()];
            int length_sum = 0;
            int indel_cnt = 0;
            for (int idx : tree) {
                if (idx == ROOT_IDX) {
                    bpseqs[ROOT_IDX] = ancseq;
                    substRates[ROOT_IDX] = new double[ancseq.length()];
                    colIndelRates[ROOT_IDX] = new double[ancseq.length()];
                    for (int i = 0; i < ancseq.length(); i++) {

                        if (params.ancrates != null) {
                            substRates[ROOT_IDX][i] = params.ancrates[i];
                        } else {
                            if (USERATES) {
                                substRates[ROOT_IDX][i] = params.substratemodel.sample();
                            } else {
                                substRates[ROOT_IDX][i] = 1.0;
                            }
                        }

                        colIndelRates[ROOT_IDX][i] = INDELRATES ? params.indelratemodel.sample() : 1;
                    }
                } else { // branchpoint has parents, all of which have been instantiated (iterator order ensures this, starting with branchpoint idx 0)
                    int paridx = tree.getParent(idx);           // idx of parent
                    Object[] parseq = bpseqs[paridx].get();     // sequence of parent
                    double t = tree.getDistance(idx);           // distance from parent to child
                    insertions[idx] = new int[parseq.length+1]; // insertions at this branchpoint relative to parent indices; note that insertions can happen before or after a sequence
                    deletions[idx] = new int[parseq.length];    // deletions at this branchpoint relative to parent indices
                    substRates[idx] = new double[parseq.length];     // rates at this branchpoint relative to parent indices; note that insertion rate for before and after is shared
                    colIndelRates[idx] = new double[parseq.length];

                    // determine what indels are introduced; note: we don't yet know how many indices are required for child so we use lists before moving to array
                    List<Object> child = new ArrayList<>();     // collect character states for the resulting positions, accommodating insertions and deletions
                    List<Object> tail  = new ArrayList<>();     // collect character states for the tail of the child; intended for tailing insertions
                    List<Double> childrates = new ArrayList<>();// collect character rates for the resulting positions, accommodating insertions and deletions
                    List<Double> tailrates = new ArrayList<>(); // collect character rates for the tail of the child; intended for tailing insertions
                    List<Double> childColIndelRates = new ArrayList<>(); // collect indel rates for the resulting positions
                    List<Double> childColIndelTailRates = new ArrayList<>(); // collect indel rates for the tail of the child

                    //double rho = params.indelratemodel.sample();// node specific rate of insertions and deletions
                    //rList.add(rho);                             // save it so the whole series can be recorded
                    // ----
                    // Next, loop through each site of the parent sequence
                    // ----
                    int i = 0;                                  // idx for parent position
                    while (i < parseq.length) {
                        // the toss is different for each site
                        // move through the child by incrementing the idx in the parent
                        // three possibilities:
                        // 1. no indel so "match" parent/child and potential substitution from parent to child,
                        // 2. deletion/s in the child, and
                        // 3. insertion/s in the child
                        double toss = params.rand.nextDouble();

                        // the indel rate is specific to the position in the parent sequence, i.e. the propensity for an indel event at this position
                        double p = Math.exp(-(colIndelRates[paridx][i] * t));

                        // make decision of what happens in child for the current site i
                        if (toss < p) { // 1. no indel (so match) with prob p = e^-rt, so consider substitution (r is evolutionary rate and t is branch distance)
                            EnumDistrib d = params.substmodel.getDistrib(parseq[i], substRates[paridx][i]*t); // probability of child states GIVEN parent state
                            Object nchar = null;
                            double tossagain = params.rand.nextDouble();
                            double sump = 0;
                            for (Object c : d.getDomain().getValues()) { // different character states have different substitution probs
                                sump += d.get(c);
                                if (sump >= tossagain) {
                                    nchar = c;
                                    break;
                                }
                            }
                            if (nchar != null) {
                                child.add(nchar);
                                childrates.add(substRates[paridx][i]); // stays the same
                                childColIndelRates.add(colIndelRates[paridx][i]); // same for col indel rate
                            } else
                                throw new RuntimeException("Sampling invalid distribution");
                            i += 1; // done with site, continue to next...
                        } else {
                            indel_cnt += 1; // count the indel event for report
                            double toss2 = params.rand.nextDouble();    // decide between deletion and insertion
                            if (toss2 < params.PROPORTION_DELETION) {   // 2. deletion with prob q = (1 - p)/2, so consider length of deletion
                                int indel_length = Math.min(params.deletemodel.sample(), parseq.length - i);// length, can only delete what is left of the sequence
                                length_sum += indel_length;                    // keep track of the aggregate length of indel events
                                deletions[idx][i] = indel_length;       // deletions skip characters in the parent
                                i += indel_length;                      // jump ahead as far as the deletion took us
                            } else { // 3. insertion with prob q = (1 - p)/2, so consider length of insertion
                                // special case at i==0: insertions can happen BEFORE and AFTER the sequence,
                                // so to avoid introducing a bias for LONGER EXTANTS,
                                // we place it at either end with a uniform coin toss
                                int indel_length = params.insertmodel.sample(); // length, can only delete what is left of the sequence
                                length_sum += indel_length;
                                insertions[idx][i == 0 ? (params.rand.nextBoolean() ? 0 : parseq.length) : i] += indel_length; // insertions can be on top of another
                                for (int j = 0; j < indel_length; j ++) {
                                    Object nchar = null;
                                    double tossagain = params.rand.nextDouble();
                                    double sump = 0;
                                    for (Object c : params.substmodel.getDomain().getValues()) {
                                        sump += params.substmodel.getProb(c);
                                        if (sump >= tossagain) {
                                            nchar = c;
                                            break;
                                        }
                                    }
                                    if (nchar != null) {
                                        if (i==0 && insertions[idx][parseq.length] > 0)
                                            tail.add(nchar);
                                        else
                                            child.add(nchar);

                                        if (i == 0 && insertions[idx][parseq.length] > 0) {
                                            tailrates.add(USERATES ? params.substratemodel.sample() : 1);
                                            childColIndelTailRates.add(INDELRATES ? params.indelratemodel.sample() : 1);
                                        } else {
                                            childrates.add(USERATES ? params.substratemodel.sample() : 1); // new position means new rate
                                            childColIndelRates.add(INDELRATES ? params.indelratemodel.sample() : 1);
                                        }
                                    } else
                                        throw new RuntimeException("Sampling invalid distribution");
                                }
                                // after an insertion, what do we do with the character? Currently, we enforce a match/substitution
                                // note that for a before/after insertion, the character concerned is always that at the head of the sequence
                                EnumDistrib d = params.substmodel.getDistrib(parseq[i], t);
                                Object nchar = null;
                                double tossagain = params.rand.nextDouble();
                                double sump = 0;
                                for (Object c : d.getDomain().getValues()) {
                                    sump += d.get(c);
                                    if (sump >= tossagain) {
                                        nchar = c;
                                        break;
                                    }
                                }
                                if (nchar != null) {
                                    child.add(nchar);
                                    childrates.add(substRates[paridx][i]); // stays the same
                                    childColIndelRates.add(colIndelRates[paridx][i]);
                                } else
                                    throw new RuntimeException("Sampling invalid distribution");
                                i += 1;
                            }
                        }
                    }
                    // set character states
                    Object[] chseq = new Object[child.size() + tail.size()];
                    substRates[idx] = new double[childrates.size() + tailrates.size()];
                    colIndelRates[idx] = new double[childColIndelRates.size() + childColIndelTailRates.size()];
                    for (int j = 0; j < chseq.length; j ++) {
                        if (j >= child.size()) {
                            substRates[idx][j] = tailrates.get(j - child.size());
                            colIndelRates[idx][j] = childColIndelTailRates.get(j - child.size());
                        } else {
                            substRates[idx][j] = childrates.get(j);
                            colIndelRates[idx][j] = childColIndelRates.get(j);
                        }

                        if (j >= child.size())
                            chseq[j] = tail.get(j - child.size());
                        else
                            chseq[j] = child.get(j);
                    }
                    bpseqs[idx] = new EnumSeq(ancseq.getType());
                    bpseqs[idx].set(chseq);
                    bpseqs[idx].setName(tree.getBranchPoint(idx).getLabel().toString());
                }
            }
            // System.out.println("result: sum of indel "+ sum + " num of indel " + indel_cnt);
            ti_deletions = new TreeInstance(tree, deletions);
            ti_insertions = new TreeInstance(tree, insertions);
            ti_seqs = new TreeInstance(tree, bpseqs);

//            if (VERBOSE) {
//                String outputFile = (OUTPUT != null ? OUTPUT : "")  +"_travis_report.txt";
//
//                try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {
//
//                    System.out.println(tree);
//                    pw.println(tree);
//
//                    for (int idx : tree) {
//                        BranchPoint bp = tree.getBranchPoint(idx);
//                        BranchPoint parent = bp.getParent();
//
//                        System.out.println(bp.getLabel() + "\t" + bpseqs[idx]);
//                        pw.println(bp.getLabel() + "\t" + bpseqs[idx]);
//
//                        if (idx != 0 && parent != null) {
//                            for (int i = 0; i < deletions[idx].length; i++) {
//                                if (deletions[idx][i] > 0) {
//                                    String line = "\tDELETE " + parent.getLabel() + "->"
//                                            + bp.getLabel() + "@" + i + ":" + deletions[idx][i];
//
//                                    System.out.println(line);
//                                    pw.println(line);
//                                }
//                            }
//                            for (int i = 0; i < insertions[idx].length; i++) {
//                                if (insertions[idx][i] > 0) {
//                                    String line = "\tINSERT " + parent.getLabel() + "->"
//                                            + bp.getLabel() + "@" + i + ":" + insertions[idx][i];
//
//                                    System.out.println(line);
//                                    pw.println(line);
//                                }
//                            }
//                        }
//                    }
//
//                    pw.flush();
//
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            }
        }

        public TreeInstance getTreeWithSequences() {
            return ti_seqs;
        }

        public TreeInstance getTreeWithDeletions() {
            return ti_deletions;
        }

        public TreeInstance getTreeWithInsertions() {
            return ti_insertions;
        }



        public EnumSeq[] getSequences() {
            Object[] oseqs = ti_seqs.getInstance();
            EnumSeq[] eseqs = new EnumSeq[oseqs.length];
            for (int i = 0; i < eseqs.length; i ++)
                eseqs[i] = (EnumSeq) oseqs[i];
            return eseqs;
        }

        public EnumSeq[] getLeafSequences() {
            EnumSeq[] all = getSequences();  // 获取所有原始序列
            List<EnumSeq> leafList = new ArrayList<>();
            for (int i = 0; i < params.tree.getSize(); i++) {
                if (params.tree.isLeaf(i)) {
                    leafList.add(all[i]);
                }
            }
            return leafList.toArray(new EnumSeq[0]);
        }



        public POAGraph getPOAG() {
            if (poag != null)  // already computed
                return poag;
            Set<EnumEdge> edges = new HashSet<>();
            EnumNode start = new EnumNode(myType);
            EnumNode end = new EnumNode(myType);
            EnumSeq ancseq = (EnumSeq)ti_seqs.getInstance(0);
            if (ancseq.length() < 1)
                return null;
            EnumNode[] ancestor = new EnumNode[ancseq.length()];
            for (int i = 0; i < ancseq.length(); i ++) {
                ancestor[i] = new EnumNode(myType);
                ancestor[i].add(ancseq.get(i));
            }
            enumNodes = new EnumNode[params.tree.getSize()][];
            build(0, start, ancestor, end, edges);
            HashMap<EnumNode, Integer> nodes = new HashMap<>();
            int cnt = 0;
            for (EnumEdge edge : edges) {
                for (EnumNode node : edge.getPair()) {
                    if (node != start && node != end) {
                        if (!nodes.containsKey(node))
                            nodes.put(node, cnt++);
                    }
                }
            }
            poag = new POAGraph(myType, nodes.size());
            for (Map.Entry<EnumNode, Integer> entry : nodes.entrySet())
                poag.addNode(entry.getValue(), entry.getKey());
            for (EnumEdge edge : edges) {
                if (edge.getPair()[0] == start) {
                    if (edge.getPair()[1] != end)
                        poag.addEdge(-1, nodes.get(edge.getPair()[1]));
                    else
                        poag.addEdge(-1, poag.maxsize()); // empty sequence
                } else if (edge.getPair()[1] == end) {
                    poag.addTerminalEdge(nodes.get(edge.getPair()[0]));
                } else {
                    poag.addEdge(nodes.get(edge.getPair()[0]), nodes.get(edge.getPair()[1]));
                }
            }
            return poag;
        }

        public EnumSeq[] getAlignment() {
            if (alignedseqs != null)
                return alignedseqs;
            POAGraph poag = getPOAG();
            int[] order = poag.getTopoSortDepthFirst();
            alignedseqs = new EnumSeq.Gappy[params.tree.getSize()];
            alignedidxs = new int[params.tree.getSize()][];
            HashMap<EnumNode, Integer> nodes = new HashMap<>();
            for (int i = 0; i < order.length; i ++) {
                EnumNode node = (EnumNode) poag.getNode(order[i]);
                nodes.put(node, i);
                if (VERBOSE)
                    System.out.println(i + "\t" + order[i] + "\t" + node.getLabel());
            }
            if (USERATES)
                alignedrates = new double[order.length];
            if (INDELRATES)
                alignedIndelRates = new double[order.length];
            for (int idx : params.tree) {
                alignedseqs[idx] = new EnumSeq.Gappy(myType);
                Object[] seq = new Object[order.length];
                EnumSeq orig = (EnumSeq)ti_seqs.getInstance(idx);
                alignedidxs[idx] = new int[orig.length()];
                for (int j = 0; j < enumNodes[idx].length; j ++) {
                    int pos = nodes.get(enumNodes[idx][j]);
                    seq[pos] = orig.get(j);
                    alignedidxs[idx][j] = pos;
                    if (USERATES)
                        alignedrates[pos] = substRates[idx][j];
                    if (INDELRATES) {
                        alignedIndelRates[pos] = colIndelRates[idx][j];
                    }
                }
                alignedseqs[idx].set(seq);
                alignedseqs[idx].setName(params.tree.getBranchPoint(idx).getLabel().toString());
            }
            return alignedseqs;
        }

        public int[][] getIndels() {
            EnumSeq[] ancseqs_gappy = getAlignment();
            int[] ins_total = new int[0];
            int[] del_total = new int[0];
            // Go through the tree, and look at each ancestor sequence, recording predicted indel events
            for (int idx : params.tree) { // go through the tree
                int parent = params.tree.getParent(idx);
                if (parent != -1) {  // Non-root node, so there is a branch with distance to catch...
                    double dist = params.tree.getDistance(idx);
                    // Retrieve parent/child reconstructed sequences at a node in the user-provided tree
                    Object[] pseq = ancseqs_gappy[parent].get();
                    Object[] cseq = ancseqs_gappy[idx].get();
                    // Calculate indel rate for the reconstructed sequences in the user-provided tree
                    int[] insertions = TrAVIS.getInsertionCounts(pseq, cseq);
                    int[] deletions = TrAVIS.getDeletionCounts(pseq, cseq);
                    // Accumulate insertion/deletion lengths
                    ins_total = TrAVIS.mergeCounts(ins_total, insertions);
                    del_total = TrAVIS.mergeCounts(del_total, deletions);
                }
            }
            return new int[][] { ins_total, del_total };
        }

        public EnumSeq[] getLeafAlignments() {
            EnumSeq[] all = getAlignment();  // 获取所有对齐后的序列
            List<EnumSeq> leafList = new ArrayList<>();
            for (int i = 0; i < params.tree.getSize(); i++) {
                if (params.tree.isLeaf(i)) {
                    leafList.add(all[i]);
                }
            }
            return leafList.toArray(new EnumSeq[0]);
        }


        public double[] getSubstRates() {
            if (!USERATES)
                return null;
            if (alignedrates != null)
                return alignedrates;
            getAlignment();
            return alignedrates;
        }

        public double[] getIndelRates() {
            if (!INDELRATES) {
                return null;
            }
            if (alignedIndelRates != null) {
                return alignedIndelRates;
            }
            getAlignment();
            return alignedIndelRates;
        }


        private void build(int idx, EnumNode start, EnumNode[] nodes, EnumNode end, Set<EnumEdge> edges) {
            enumNodes[idx] = nodes;
            EnumNode prev = start;
            for (int i = 0; i < nodes.length; i ++) {
                edges.add(new EnumEdge(prev, nodes[i]));
                prev = nodes[i];
            }
            edges.add(new EnumEdge(prev, end));
            for (int chidx : params.tree.getChildren(idx)) {
                int[] deletions = (int[])ti_deletions.getInstance(chidx);
                int[] insertions = (int[])ti_insertions.getInstance(chidx);
                int childptr = 0; // index to child sequence
                int parentptr = 0; // index to parent sequence
                EnumSeq chseq = (EnumSeq) ti_seqs.getInstance(chidx);
                EnumNode[] child = new EnumNode[chseq.length()];
                for (int i = 0; i <= nodes.length; i ++) {
                    if (i < nodes.length) { // NOT special case: insertion at end
                        if (deletions[i] > 0) { // remove node/s
                            i += (deletions[i] - 1);
                        } else { // new or existing nodes
                            if (insertions[i] > 0) {
                                for (int j = 0; j < insertions[i]; j++) {
                                    child[childptr] = new EnumNode(myType);
                                    child[childptr].add(chseq.get(childptr));
                                    childptr += 1;
                                }
                            } // insertion is followed by a match... which is why no "else" clause
                            child[childptr] = nodes[i + parentptr];
                            child[childptr].add(chseq.get(childptr));
                            childptr += 1;
                        }
                    } else { // special case: insertion at end
                        if (insertions[i] > 0) {
                            for (int j = 0; j < insertions[i]; j++) {
                                child[childptr] = new EnumNode(myType);
                                child[childptr].add(chseq.get(childptr));
                                childptr += 1;
                            }
                        } // tail-end insertion is NOT followed by a match
                    }
                }
                build(chidx, start, child, end, edges);
            }
        }
    }
}