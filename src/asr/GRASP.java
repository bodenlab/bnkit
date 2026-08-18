package asr;

import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.*;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.pog.IdxGraph;
import dat.pog.POAGraph;
import dat.pog.POGTree;
import dat.pog.POGraph;
import stats.*;
import java.io.IOException;
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import asr.IndelSegmentation.RATE_CATEGORY;

/**
 * Command line version of GRASP.
 * @author mikael
 * @author ariane
 * @author gabe
 * @author Sebastian
 */
public class GRASP {

    public static String VERSION = "14-Aug-2026";

    public static void usage() {
        usage(0, null);
    }

    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        out.println("""
                Usage: asr.GRASP\s
                \t[-a | --aln <filename>]
                \t[-n | --nwk <filename>]
                \t{-o | --output-folder <foldername>} (default is current working folder, or input folder if available)
                \t{-i | --input-folder <foldername>}
                \t{-pre | --prefix <stub>}
                \t{-rf | --rates-file <filename>}
                \t{-ef | --empirical-freqs <filename>}
                \t{-s | --substitution-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}
                \t{-t | --threads <number>}
                \t{-j | --joint (default)}
                \t{-m | --marginal <branchpoint-id>}
                \t{--indel-method <methodname>} (select one from BEP(default) BEML SICP SICML PSP PSML SCIP Gurobi)
                \t{--supported-path <methodname>} (select one from DIJKSTRA(default) ASTAR)
                \t{--nogap}
                \t{--seed <seed>}
                \t{--nonibble}
                \t{--exclude-noedge}
                \t{--save-as <list-of-formats>} (select multiple from FASTA CLUSTAL TREE DISTRIB ASR DOT TREES)
                \t{--save-all} (saves reconstruction with ALL formats)
                \t{--save-tree} (bypasses inference and re-saves the tree with ancestor nodes labelled as per GRASP's
                \tdepth-first labelling scheme starting with N0)
                \t{--save-poag { <branchpoint-id> } (bypasses inference and saves the input alignment as a POAG
                \t(partial order alignment graph of extant sequences under specified ancestor [default N0])
                \t{--time}{--verbose}{--help}
                """);
        out.println("""
                Inference is a two-stage process:
                \t(1) A history of indel events is inferred by either maximum likelihood or maximum parsimony and\s
                \tmapped onto the tree to determine what positions contain actual sequence content
                \t(2) For each ancestral position, the most probable character is assigned to each phylogenetic branch\s
                \tpoint when performing a joint reconstruction. Alternatively, for each\s
                \tposition at a nominated branch point, the probability distribution over all possible\s
                \tcharacters is inferred when performing a marginal reconstruction.
                \tFinally, edges are drawn to represent all inferred combinations of indels to form an ancestor POG\s
                \twith nodes that can form a valid sequence with inferred content; a preferred path
                \tthrough the POG is then inferred, nominating a single, best supported sequence.
                """);
        out.println("""
                Mode of character inference:
                \t-j (or --joint) activates joint reconstruction (default),\s
                \t-m (or --marginal) activates marginal reconstruction (requires a branch-point to be nominated)
                \t--onlyindel disengages the stage of character state inference
                """);
        out.println("Required arguments:\n" +
                "\t-a (or --aln) must specify the name of a multiple-sequence alignment file on FASTA or CLUSTAL format\n" +
                "\t-n (or --nwk) must specify the name of a phylogenetic-tree file in Newick format\n");
        out.println("Optional arguments:\n" +
                "\t-o (or --output-folder) specifies the folder that will be used to save output files,\n\t\te.g. inferred ancestor or ancestors, tree, etc. as specified by format\n" +
                "\t-i (or --input-folder) skips indel inference, and loads a previous reconstruction from specified folder\n" +
                "\t-sa (or --save-as) lists the files and formats to be generated (see below)\n\t--save-all nominates all\n" +
                "\t-pre (or --prefix) specifies a stub that is added to result filenames (default is the prefix of the alignment file)\n" +
                "\t-indel (or --indel-method) specifies what method to use for inferring indels (see below)\n" +
                "\t-s (or --substitution-model) specifies what evolutionary model to use for inferring character states (see below)\n" +
                "\t-rf (or --rates-file) specifies a tabulated file with relative, position-specific substitution rates\n\t\tWe recommend the use of this generally, but specifically for trees with great distances, and with biologically diverse entries\n\t\tAs an example, IQ-TREE produces rates on the accepted format with the --rate option (--mlrate is NOT supported yet).\n" +
                "\t-ef (or --empirical-freqs) specifies a tabulated file with the headers Character & Proportion that contain\n\t\tstationary character frequencies for the chosen substitution model. The standard stationary character frequencies of a\n\t\tchosen substitution model are used by default when -ef is not specified\n" +
                "\t--include-extants means that extants are included in output files (when the format allows)\n" +
                "\t--nogap means that the gap-character is excluded in the resulting output (when the format allows)\n" +
                "\t--nonibble de-activates the removal of indices in partial order graphs that cannot form a path from start to end\n" +
                "\t--orphans de-activates the removal of orphaned indel trees\n" +
                "\t--exclude-noedge removes non-existing edge as an option for parsimony in BEP\n" +
                "\t--solver-time-limit the maximum time the MIP solver can run for in minutes before defaulting to BEP indel inference\n" +
                "\t--verbose prints out information about steps undertaken, and --time the time it took to finish\n" +
                "\t-h (or --help) will print out this screen\n");
        out.println("""
                Files/formats:\s
                \tFASTA: sequences (most preferred path at each ancestor, gapped or not gapped)
                \tCLUSTAL: sequences (most preferred path at each ancestor, gapped)
                \tTREE: phylogenetic tree with ancestor nodes labelled
                \tDISTRIB: character distributions for each position (indexed by POG, only available for marginal reconstruction)
                \tASR: complete reconstruction as JSON, incl. POGs of ancestors and extants, and tree (ASR.json)
                \tDOT: partial-order graphs of ancestors in DOT format
                \tTREES: position-specific trees with ancestor states labelled""");
        out.println("""
                Indel-methods:\s
                \tBEP: bi-directional edge (maximum) parsimony
                \tBEML: bi-directional edge maximum likelihood (uses uniform evolutionary model akin to JC)
                \tSICP: simple indel-coding (maximum) parsimony (based on Simmons and Ochoterena)
                \tSICML: simple indel-coding maximum likelihood (uses uniform evolutionary model)
                \tPSP: position-specific (maximum) parsimony
                \tPSML: position-specific maximum likelihood (uses uniform evolutionary model)
                \tSCIP: globally optimal distance sensitive parsimony-based indel history using the open-source
                \t\tSCIP solver (https://www.scipopt.org/). Does not support multi-threading
                \tGurobi: globally optimal distance sensitive parsimony-based indel history.
                \t\tRequires local installation of Gurobi to run (https://www.gurobi.com/downloads/)
                \tAdd '*' to method name for less conservative setting (if available) or to use a default gap opening penalty of 2 for SCIP or Gurobi\s
                """);
        out.println("""
                Substitution-models:\s
                \tJTT: Jones-Taylor-Thornton (protein; default)
                \tDayhoff: Dayhoff-Schwartz-Orcutt (protein)
                \tLG: Le-Gasquel (protein)
                \tWAG: Whelan-Goldman (protein)
                \tJC: Jukes-Cantor (DNA)
                \tYang: Yang's general reversible process model (DNA)
                """);
        out.println("Notes: \n" +
                "\tGreater number of threads may improve processing time up to a point when coordination chokes performance; default is 4 threads.\n" +
                "\tRunning GRASP requires large memory and in most cases Java needs to be run with the option -Xmx20g, \n\twhere 20g specifies that 20GB of RAM should be available.\n" +
                "\n~ This is version " + VERSION + " ~");
        if (msg != null)
            out.println("\n" + msg + " (Error " + error + ")");
        System.exit(error);
    }


    // GENERAL SETTINGS
    public static boolean VERBOSE = false;
    public static boolean TIME = false;
    public static int NTHREADS = 4;
    public static boolean NIBBLE = true;

    // INPUT VARIABLES
    private static String ALIGNMENT = null;
    private static String ASRFILE = "ASR.json";
    private static String NEWICK = null;
    static String OUTPUT = null;
    private static String INPUT = null;
    static String PREFIX = null;
    private static String RATESFILE = null;
    public static double[] RATES = null;
    private static String EMPIRICAL_FREQS_FILE = null;
    private static double[] EMPIRICAL_FREQS = null;
    private static final String[] MODELS = new String[]{"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
    private static int MODEL_IDX = 0; // default model is that above indexed 0
    public static SubstModel MODEL = null;
    // Alphabet is decided by MODEL_IDX
    private static final Enumerable[] ALPHAS = new Enumerable[]{Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};

    // COLUMN INFERENCE SETTINGS
    private static Inference MODE = Inference.JOINT; // default inference mode
    // ancestor to reconstruct if inference mode is "marginal"
    private static Integer MARG_NODE = null;
    private static int SEED;
    private static boolean BYPASS = false; // bypass inference, default is false
    private static boolean NEED_CONSENSUS = false;
    public enum Inference {
        JOINT,
        MARGINAL
    }


    // OUTPUT FORMATS
    private static boolean SAVE_AS = false;
    private static boolean INCLUDE_EXTANTS = false;
    private static final String[] FORMATS = new String[]{"FASTA", "DISTRIB", "CLUSTAL", "TREE", "ASR", "DOT", "TREES", "POAG"};
    private static boolean GAPPY = true;
    private static final boolean[] SAVE_AS_IDX = new boolean[FORMATS.length];
    // select to compute consensus path for these output formats
    private static final boolean[] CONSENSUS = new boolean[]{true, false, true, false, false, false, false, false};
    private static final int FASTA = 0;
    private static final int DISTRIB = 1;
    private static final int CLUSTAL = 2;
    private static final int TREE = 3;
    private static final int POGS = 4;
    private static final int DOT = 5;
    private static final int TREES = 6;
    private static final int POAG = 7;

    // INDEL INFERENCE SETTINGS
    private static final int BEP = 0;
    private static final int BEPML = 1;
    private static final int SICP = 2;
    private static final int SICML = 3;
    private static final int PSP = 4;
    private static final int PSML = 5;
    private static final int SCIP = 6;
    private static final int GUROBI = 7;
    public static RATE_CATEGORY INDEL_RATE = RATE_CATEGORY.HIGH;
    public static int MIP_SOLVER_TIME_LIMIT_MINUTES = 720; // 12 hours
    private static final String[] INDELS = new String[]{"BEP", "BEML", "SICP", "SICML", "PSP", "PSML", "SCIP", "Gurobi"};
    private static int INDEL_IDX = 0; // default indel approach is that above indexed 0
    private static final String[] SPATH = new String[]{"DIJKSTRA", "ASTAR"};
    public static boolean RANDOM_RATES = false;
    public static boolean SIMPLE_RATES = false;
    public static boolean SEQ_RATES = false; // in development
    public static boolean COL_RATES = false; // in development
    public static boolean INDEL_CONSERVATIVE = true;
    public static boolean DISTANCE_BASED_MIP = true;
    public static int NUM_GAMMA_CATEGORIES = 20;

    // Mode for BEP
    public static boolean RECODE_NULL = true;
    public static boolean REMOVE_INDEL_ORPHANS = true;
    public static boolean ONLYINDEL = false;

    // TRAVIS PARAMS
    private static boolean RUN_TRAVIS = false;
    private static int EXTANTS_N = 5;
    private static final int DISTRIB_NAME = 0;
    private static final int DISTRIB_PARAMS = 1;
    private static double TREE_GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
    private static double TREE_GAMMA_SCALE = 0.2;
    private static final int DESCENDANTS_MAX = 2, DESCENDANTS_MIN = 2; // Max and min of tree branching
    private static double DELETIONPROP = 0.5; // proportion of DELETIONS v INSERTIONS
    private static Double SCALEDIST = null;
    private static boolean LEARN = false;
    private static int TRAVIS_FORMAT_IDX = 0;
    private static boolean EXTANTS_ONLY = false;
    private static boolean COPY_TREE = false;
    private static IndelModel INDEL_LENGTH_MODEL = null;
    private static IndelModel INSERTION_LENGTH_MODEL = null;
    private static IndelModel DELETION_LENGTH_MODEL = null;
    private static RateModel SUBST_RATE_MODEL = null;
    private static RateModel INDEL_RATE_MODEL = null;
    private static RateModel TREE_DISTANCE_MODEL = null;
    private static Distrib LEAF2ROOT_DISTANCE_MODEL = null;
    private static Integer ANCSEQ_LENGTH = null;
    private static String ANCSEQ = null; // ancestor sequence as a text string, provided
    private static boolean PERFORM_TRAVIS_SIMUL = false;
    private static boolean TRAVIS_LEARN_NO_RECON = false;

    public static void main(String[] args) {

        // Setup GRASP/TrAVIS to run
        SEED = new Random().nextInt();
        RUN_TRAVIS = travisFlagUsed(args);
        if (RUN_TRAVIS) {
            parseTravisArgs(args);
        } else {
            parseGRASPArgs(args);
        }
        checkArgsValid();

        // main variables that will hold outputs
        EnumSeq.Alignment<Enumerable> aln = null;
        Tree tree = null;
        Object[][] ancseqs_gappy = null;
        Object[][] ancseqs_nogap = null;
        String[] ancnames = null;
        POGraph[] ancestors = null;
        Prediction indelpred = (INPUT != null) ? setupIndelPrediction() : null;
        boolean reconstructionRequired = indelpred == null && !PERFORM_TRAVIS_SIMUL && !TRAVIS_LEARN_NO_RECON;

        // Basic data required for a reconstruction
        if (reconstructionRequired) {
            try {
                aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
                tree = Utils.loadTree(NEWICK);
                Utils.checkData(aln, tree);
            } catch (ASRException e) {
                if (RUN_TRAVIS) {
                    TrAVIS.usage(10, "Invalid input for ASR: " + e.getMessage() + "\n If you're learning TrAVIS params from an alignment and tree with ancestors, use --no-recon");
                } else {
                    usage(22, "Invalid input for ASR: " + e.getMessage());
                }
            } catch (IOException e) {
                if (RUN_TRAVIS) {
                    TrAVIS.usage(2, "Failed to read or write files: " + e.getMessage() + "\n");
                } else {
                    usage(2, "Failed to read or write files: " + e.getMessage());
                }
            }
        }

        long START_TIME = System.currentTimeMillis();

        // INDEL INFERENCE
        if (!BYPASS && reconstructionRequired) {
            POGTree pogtree = new POGTree(aln, tree);
            indelpred = performIndelInference(pogtree, aln);
        }

        // COLUMN INFERENCE
        if (!BYPASS && !PERFORM_TRAVIS_SIMUL && !TRAVIS_LEARN_NO_RECON) {
            if (indelpred == null)
                usage(3, INDELS[INDEL_IDX] + " is not implemented");

            performColumnInference(indelpred);

            Map<Object, POGraph> pogs = indelpred.getAncestors(MODE);
            ancestors = storePOGsInArray(pogs);

            ancnames = new String[pogs.size()];
            if (NEED_CONSENSUS) {
                ancseqs_gappy = new Object[pogs.size()][];
                ancseqs_nogap = new Object[pogs.size()][];
                extractAncestralSequences(ancseqs_gappy, ancseqs_nogap, pogs, indelpred, ancnames);
            }
        }

        if (RUN_TRAVIS) {
            runTravis(aln, tree, indelpred, ancseqs_nogap, ancseqs_gappy);
        } else {
            saveGraspOutput(ancnames, ancseqs_nogap, ancseqs_gappy, indelpred, tree, aln, ancestors);
        }

        long ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
        if (VERBOSE || TIME) {
            System.out.printf("Done in %d min, %d sec%n", TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                    TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME)));
        }
    }

    private static boolean travisFlagUsed(String[] args) {

        boolean runTravis = false;
        for (String s : args) {
            if (s.startsWith("-")) {
                String arg = s.substring(1);
                if ((arg.equalsIgnoreCase("-travis"))) {
                    runTravis = true;
                    break;
                }
            }
        }

        return runTravis;
    }

    /**
     * Parse command line arguments for TrAVIS
     * @param args command line arguments
     */
    private static void parseTravisArgs(String[] args) {

        // Because we are creating distributions as we parse the data, need to find the seed first
        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if ((arg.equalsIgnoreCase("-seed") && args.length > a + 1)) {
                    SEED = Integer.parseInt(args[++a]);
                    break;
                }
            }
        }

        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);

                if (arg.equalsIgnoreCase("-travis")) {
                    continue;
                }

                if (arg.equalsIgnoreCase("n0") || arg.equalsIgnoreCase("-ancestor") && args.length > a + 1) {
                    ANCSEQ = args[++a];
                    PERFORM_TRAVIS_SIMUL = true;
                } else if ((arg.equalsIgnoreCase("-aln") || arg.equalsIgnoreCase("a")) && args.length > a + 1) {
                    ALIGNMENT = args[++ a];
                } else if (arg.equalsIgnoreCase("-nwk")  || arg.equalsIgnoreCase("n") && args.length > a + 1) {
                    NEWICK = args[++a];
                } else if (arg.equalsIgnoreCase("-no-recon")) {
                    TRAVIS_LEARN_NO_RECON = true; // bypasses reconstruction and just produces TrAVIS parameters based on input tree and alignment
                } else if (arg.equalsIgnoreCase("o") || arg.equalsIgnoreCase("-output-folder") && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if (arg.equalsIgnoreCase("-seed") && args.length > a + 1) {
                    ++a; // can skip over as we would've already parsed above
                } else if ((arg.equalsIgnoreCase("-prefix") || arg.equalsIgnoreCase("pre")) && args.length > a + 1) {
                    PREFIX = args[++a];
                } else if (arg.equalsIgnoreCase("-extants") && args.length > a + 1) {
                    EXTANTS_N = Integer.parseInt(args[++a]);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-dgamma") && args.length > a + 2) {
                    TREE_GAMMA_SHAPE = Double.parseDouble(args[a+1]);
                    TREE_GAMMA_SCALE = Double.parseDouble(args[a+2]);
                    if (args.length > a + 3 && !args[a + 3].startsWith("-")) {
                        SCALEDIST = Double.parseDouble(args[a+3]); // brdist_scale
                    }
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-no-gap")) {
                    GAPPY = false;
                } else if (arg.equalsIgnoreCase("-learn")) {
                    LEARN = true;
                } else if (arg.equalsIgnoreCase("-verbose")) {
                    VERBOSE = true;
                } else if (arg.equalsIgnoreCase("-help") || arg.equalsIgnoreCase("h")) {
                    TrAVIS.usage();
                } else if ((arg.equalsIgnoreCase("s") || arg.equalsIgnoreCase("-substitution-model")) && args.length > a + 1) {
                    boolean found_model = false;
                    for (int i = 0; i < MODELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(MODELS[i])) {
                            MODEL_IDX = i;
                            found_model = true;
                        }
                    }
                    if (!found_model)
                        TrAVIS.usage(1, args[a + 1] + " is not a valid model name for option --substitution-model");

                } else if ((arg.equalsIgnoreCase("-rates-file") || arg.equalsIgnoreCase("rf")) && args.length > a + 1) {
                    RATESFILE = args[++a];
                } else if (arg.equalsIgnoreCase("-indel-rate-distrib") && args.length > a + 1) {
                    String[] params = TrAVIS.parseDistribParamString(args[a + 1]);
                    INDEL_RATE_MODEL = RateModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-subst-rate-distrib")) {
                    String[] params = TrAVIS.parseDistribParamString(args[a+1]);
                    SUBST_RATE_MODEL = RateModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-indel-length-distrib") && args.length > a + 1) {
                    String[] params = TrAVIS.parseDistribParamString(args[a+1]);
                    INDEL_LENGTH_MODEL = IndelModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-insertion-length-distrib") && args.length > a + 1) {
                    String[] params = TrAVIS.parseDistribParamString(args[a+1]);
                    INSERTION_LENGTH_MODEL = IndelModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-deletion-length-distrib") && args.length > a + 1) {
                    String[] params = TrAVIS.parseDistribParamString(args[a+1]);
                    DELETION_LENGTH_MODEL = IndelModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-delprop") && args.length > a + 1) {
                    DELETIONPROP = Double.parseDouble(args[++a]);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-length") || arg.equalsIgnoreCase("l") && args.length > a + 1) {
                    ANCSEQ_LENGTH = Integer.parseInt(args[++a]);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-dist-distrib") && args.length > a + 1) {
                    String[] params = TrAVIS.parseDistribParamString(args[a+1]);
                    TREE_DISTANCE_MODEL = RateModel.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS], SEED);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-leaf2root-distrib") && args.length > a + 1) {
                    String[] params = TrAVIS.parseDistribParamString(args[a + 1]);
                    LEAF2ROOT_DISTANCE_MODEL = Distrib.create(params[DISTRIB_NAME], params[DISTRIB_PARAMS]);
                    PERFORM_TRAVIS_SIMUL = true;
                } else if (arg.equalsIgnoreCase("-copy-tree") && args.length > a + 1) {
                    COPY_TREE = true;
                } else if (arg.equalsIgnoreCase("-extants-only")) {
                    EXTANTS_ONLY = true;
                } else if (arg.equalsIgnoreCase("sa") || arg.equalsIgnoreCase("-save-as") && args.length > a + 1) {
                    boolean found_format = false;
                    for (int i = 0; i < TrAVIS.TRAVIS_FORMATS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(TrAVIS.TRAVIS_FORMATS[i])) {
                            TRAVIS_FORMAT_IDX = i;
                            found_format = true;
                        }
                    }
                    if (!found_format)
                        TrAVIS.usage(1, args[a + 1] + " is not a valid format name");
                } else if ((arg.equalsIgnoreCase("-threads") || arg.equalsIgnoreCase("t")) && args.length > a + 1) {
                        try {
                            NTHREADS = Integer.parseInt(args[++a]);
                        } catch (NumberFormatException e) {
                            TrAVIS.usage(2, "Failed to set number of threads for option --threads: " + args[a] + " is not a valid integer");
                        }
                } else {
                    TrAVIS.usage(5, "Unknown option or missing required argument: \"" + args[a] + "\"");
                }
            }
        }
    }

    /**
     * Parse command line arguments for GRASP
     * @param args command line arguments
     */
    private static void parseGRASPArgs(String[] args) {
        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if ((arg.equalsIgnoreCase("-aln") || arg.equalsIgnoreCase("a")) && args.length > a + 1) {
                    ALIGNMENT = args[++a];
                } else if ((arg.equalsIgnoreCase("-nwk") || arg.equalsIgnoreCase("n")) && args.length > a + 1) {
                    NEWICK = args[++a];
                } else if ((arg.equalsIgnoreCase("-output-folder") || arg.equalsIgnoreCase("o")) && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if ((arg.equalsIgnoreCase("-input-folder") || arg.equalsIgnoreCase("i")) && args.length > a + 1) {
                    INPUT = args[++a];
                } else if ((arg.equalsIgnoreCase("-prefix") || arg.equalsIgnoreCase("pre")) && args.length > a + 1) {
                    PREFIX = args[++a];
                } else if ((arg.equalsIgnoreCase("-rates-file") || arg.equalsIgnoreCase("rf")) && args.length > a + 1) {
                    RATESFILE = args[++a];
                } else if ((arg.equalsIgnoreCase("-seed") && args.length > a + 1)) {
                    SEED = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-joint") || arg.equalsIgnoreCase("j")) {
                    MODE = Inference.JOINT;
                } else if (arg.equalsIgnoreCase("-random-rates")) {
                    RANDOM_RATES = true;
                } else if (arg.equalsIgnoreCase("-simple-rates")) {
                    SIMPLE_RATES = true;
                } else if (arg.equalsIgnoreCase("-seq-rates")) {
                    SEQ_RATES = true;
                } else if (arg.equalsIgnoreCase("-col-rates")) {
                    COL_RATES = true;
                } else if ((arg.equalsIgnoreCase("-marginal") || arg.equalsIgnoreCase("m")) && args.length > a + 1) {
                    MODE = Inference.MARGINAL;
                    String ancid = args[++a];
                    if (ancid.startsWith("N"))
                        ancid = ancid.substring(1);
                    try {
                        MARG_NODE = Integer.parseInt(ancid);
                    } catch (NumberFormatException e) {
                        usage(2, args[a] + " is not a valid ancestor name (use <number>, or \"N<number>\", where <number> starts with 0 at root, depth-first). Tip: perform joint reconstruction first to check branch point numbering in tree.");
                    }
                } else if (arg.equalsIgnoreCase("-onlyindel")) {
                    MODE = null;
                    ONLYINDEL = true;
                } else if ((arg.equalsIgnoreCase("-substitution-model") || arg.equalsIgnoreCase("s")) && args.length > a + 1) {
                    boolean found_model = false;
                    for (int i = 0; i < MODELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(MODELS[i])) {
                            MODEL_IDX = i;
                            found_model = true;
                        }
                    }
                    if (!found_model)
                        usage(1, args[a + 1] + " is not a valid model name for option --substitution-model");
                } else if ((arg.equalsIgnoreCase("-indel-method") || arg.equalsIgnoreCase("indel")) && args.length > a + 1) {
                    boolean found_indel = false;
                    for (int i = 0; i < INDELS.length; i++) {
                        if (args[a + 1].startsWith(INDELS[i])) {
                            INDEL_IDX = i;
                            found_indel = true;
                            if (args[a + 1].endsWith("*")) {
                                INDEL_CONSERVATIVE = false;
                                DISTANCE_BASED_MIP = false;
                            }
                        }
                    }
                    if (!found_indel)
                        usage(3, args[a + 1] + " is not a valid indel approach for option --indel-method");
                } else if (arg.equalsIgnoreCase("-supported-path") && args.length > a + 1) {
                    boolean found_spath = false;
                    for (int i = 0; i < SPATH.length; i++) {
                        if (args[a + 1].startsWith(SPATH[i])) {
                            // default supported path approach is that above indexed 0
                            POGraph.SUPPORTED_PATH_DEFAULT = i;
                            found_spath = true;
                        }
                    }
                    if (!found_spath)
                        usage(6, args[a + 1] + " is not a valid method for option --supported-path");
                } else if ((arg.equalsIgnoreCase("-save-as") || arg.equalsIgnoreCase("sa")) && args.length > a + 1) {
                    String format = "<none given>";
                    for (int a1 = a + 1; a1 < args.length; a1++) {
                        if (args[a1].startsWith("-"))
                            break;
                        format = args[a1];
                        boolean found_format = false;
                        for (int i = 0; i < FORMATS.length; i++) {
                            if (format.equalsIgnoreCase(FORMATS[i])) {
                                SAVE_AS_IDX[i] = true;
                                found_format = true;
                                break;
                            }
                        }
                        if (!found_format)
                            usage(1, args[a + 1] + " is not a valid format name for option --save-as");
                    }
                    SAVE_AS = true;
                } else if (arg.equalsIgnoreCase("-save-all")) {
                    for (int i = 0; i < FORMATS.length - 2; i++)
                        SAVE_AS_IDX[i] = true;
                    SAVE_AS = true;
                } else if (arg.equalsIgnoreCase("-save-tree")) {
                    BYPASS = true;
                    SAVE_AS = true;
                    SAVE_AS_IDX[3] = true;
                } else if (arg.equalsIgnoreCase("-save-poag")) {
                    MARG_NODE = 0;
                    BYPASS = true;
                    SAVE_AS = true;
                    SAVE_AS_IDX[9] = true;
                    if (a + 1 < args.length) {
                        String ancid = args[++a];
                        if (ancid.startsWith("-")) { // another option, so no ancestor given
                            a--;
                            continue;
                        } else { // ancestor specified
                            if (ancid.startsWith("N"))
                                ancid = ancid.substring(1);
                            try {
                                MARG_NODE = Integer.parseInt(ancid);
                            } catch (NumberFormatException e) {
                                usage(2, args[a] + " is not a valid ancestor name (use <number>, or \"N<number>\", where <number> starts with 0 at root, depth-first). Tip: use option --save-tree to check branch point numbering in tree.");
                            }
                        }
                    }
                } else if (arg.equalsIgnoreCase("-exclude-noedge")) {
                    RECODE_NULL = false;
                } else if (arg.equalsIgnoreCase("-include-extants")) {
                    INCLUDE_EXTANTS = true;
                } else if ((arg.equalsIgnoreCase("-empirical-freqs") || arg.equalsIgnoreCase("ef")) && args.length > a + 1) {
                    EMPIRICAL_FREQS_FILE = args[++a];
                } else if ((arg.equalsIgnoreCase("-threads") || arg.equalsIgnoreCase("t")) && args.length > a + 1) {
                    try {
                        NTHREADS = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(2, "Failed to set number of threads for option --threads: " + args[a] + " is not a valid integer");
                    }
                } else if ((arg.equalsIgnoreCase("-gamma-cat")) && args.length > a + 1) {
                    try {
                        NUM_GAMMA_CATEGORIES = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(2, "Failed to set number of discrete Gamma categories for option --gamma-cat: " + args[a] + " is not a valid integer");
                    }

                } else if (arg.equalsIgnoreCase("-nogap")) {
                    GAPPY = false;
                } else if (arg.equalsIgnoreCase("-indel-prior")) {

                    switch (args[++a].toUpperCase()) {
                        case "LOWGAP" -> {
                            INDEL_RATE = RATE_CATEGORY.LOW;
                        }
                        case "HIGHGAP" -> {
                            INDEL_RATE = RATE_CATEGORY.HIGH;
                        }
                        default ->
                                usage(25, args[a] + " is not a valid indel prior (choose from LOWGAP (UniRef30), HIGHGAP (PFAM))");
                    }

                } else if (arg.equalsIgnoreCase("-verbose")) {
                    VERBOSE = true;
                } else if (arg.equalsIgnoreCase("-time")) {
                    TIME = true;
                } else if (arg.equalsIgnoreCase("-nonibble")) {
                    NIBBLE = false;
                } else if (arg.equalsIgnoreCase("-solver-time-limit")) {
                    try {
                        MIP_SOLVER_TIME_LIMIT_MINUTES = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(2, "Failed to set time limit for MIP solver: " + args[a] + " is not a valid integer");
                    }

                } else if (arg.equalsIgnoreCase("-orphans")) {
                    REMOVE_INDEL_ORPHANS = false;
                } else if (arg.equalsIgnoreCase("-help") || arg.equalsIgnoreCase("h")) {
                    usage();
                } else {
                    usage(5, "Unknown option or missing required argument: \"" + args[a] + "\"");
                }
            }
        }
    }

    /**
     * Check that the command line arguments are valid and set up the substitution model and output formats
     */
    private static void checkArgsValid() {

        if (OUTPUT == null)
            OUTPUT = INPUT == null ? "." : INPUT;

        if (RATESFILE != null) {
            parseRatesFile();
        }

        if (EMPIRICAL_FREQS_FILE != null) {
            checkEmpiricalFreqsFile();
        }

        setupSubstModel();
        setOutputFormats();

        if (PREFIX == null) {
            setFilePrefix();
        }

        if (TrAVIS.TRAVIS_FORMATS[TRAVIS_FORMAT_IDX].equalsIgnoreCase("CLUSTAL")) // Clustal files can only be "gappy"
            GAPPY = true;

        if (!RUN_TRAVIS) {
            if (ALIGNMENT == null && INPUT == null)
                usage(3, "Must specify alignment (--aln <CLUSTAL or FASTA file>) or previously saved folder (--input-folder <folder>");
            else if (NEWICK == null && INPUT == null)
                usage(4, "Must specify phylogenetic tree (Newick file) or previously saved folder (--input-folder <folder>");
        } else {
            if (LEARN && (ALIGNMENT == null || NEWICK == null)) {
                TrAVIS.usage(2, "Must specify alignment (--aln) and tree (--nwk) to learn parameters");
            } else if (!LEARN) {
                if (INDEL_RATE_MODEL == null) {
                    TrAVIS.usage(3, "Must specify indel rate distribution to run simulation (--indel-rate-distrib)");
                }

                if (ANCSEQ == null && ANCSEQ_LENGTH == null) {
                    TrAVIS.usage(4, "Must specify ancestor sequence or ancestor sequence length to run simulation (--ancestor or --length)");
                } else if (ANCSEQ != null && ANCSEQ_LENGTH != null) {
                    TrAVIS.usage(4, "Cannot use --length and --ancestor together; use --length for a random root sequence or --ancestor to set the root as the consensus root sequence from a reconstruction");
                }

                if (NEWICK == null && TREE_DISTANCE_MODEL == null && LEAF2ROOT_DISTANCE_MODEL == null) {
                    TrAVIS.usage(5, "Must specify tree (--nwk) or tree distance distribution (--dist-distrib or --leaf2root-distrib) to run simulation");
                } else if ((TREE_DISTANCE_MODEL != null && LEAF2ROOT_DISTANCE_MODEL == null) || (TREE_DISTANCE_MODEL == null && LEAF2ROOT_DISTANCE_MODEL != null)) {
                    TrAVIS.usage(5, "Must specify --leaf2root-distrib or --dist-distrib to run simulation if --nwk is not specified");
                }

                if (ALIGNMENT != null) {
                    TrAVIS.usage(6, "An alignment cannot be used when performing a simulation");
                }

            }
        }
    }

    private static void setFilePrefix() {
        int idx2 = ALIGNMENT == null ? 0 : ALIGNMENT.lastIndexOf(".");
        if (idx2 == -1)
            idx2 = ALIGNMENT.length();
        int idx1 = ALIGNMENT == null ? 0 : ALIGNMENT.lastIndexOf("/") + 1;
        PREFIX = ALIGNMENT == null ? "" : ALIGNMENT.substring(idx1, idx2);

    }

    private static void checkEmpiricalFreqsFile() {
        double totalFreq = 0.0;
        try {
            EMPIRICAL_FREQS = TSVFile.loadEmpiricalFreqFile(EMPIRICAL_FREQS_FILE, MODELS[MODEL_IDX]);

            for (double empiricalFreq : EMPIRICAL_FREQS) {
                totalFreq += empiricalFreq;
            }

            double tolerance = 1e-5;
            if (Math.abs(totalFreq - 1.0) >= tolerance) {
                System.out.println("WARNING: Empirical frequencies do not sum to 1.0 (sum is " + totalFreq + ")\n Renormalizing frequencies.");
                for (int i = 0; i < EMPIRICAL_FREQS.length; i++) {
                    EMPIRICAL_FREQS[i] = EMPIRICAL_FREQS[i] / totalFreq;
                }
            }

        } catch (ClassCastException e) {
            if (RUN_TRAVIS) {
                TrAVIS.usage(29, e.getMessage());
            } else {
                usage(29, e.getMessage());
            }
        } catch (NumberFormatException e) {
            if (RUN_TRAVIS) {
                TrAVIS.usage(28, e.getMessage());
            } else {
                usage(28, e.getMessage());
            }
        } catch (IOException e) {
            if (RUN_TRAVIS) {
                TrAVIS.usage(30, "Empirical frequencies file could not be opened or read: " + EMPIRICAL_FREQS_FILE);
            } else {
                usage(30, "Empirical frequencies file could not be opened or read: " + EMPIRICAL_FREQS_FILE);
            }
        } catch (RuntimeException e) {
            if (RUN_TRAVIS) {
                TrAVIS.usage(27, e.getMessage());
            } else {
                usage(27, e.getMessage());
            }
        }
    }

    /**
     * Set up the substitution model based on the command line arguments
     */
    private static void setupSubstModel() {

        if (EMPIRICAL_FREQS_FILE != null) {
            MODEL = SubstModel.createModel(MODELS[MODEL_IDX], EMPIRICAL_FREQS);
        } else {
            MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
        }

        if (MODEL == null)
            usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");

    }

    /**
     * Set the output formats based on the command line arguments
     */
    private static void setOutputFormats() {

        if (!SAVE_AS && MODE == Inference.JOINT) { // set default files to save for joint
            SAVE_AS_IDX[0] = SAVE_AS_IDX[3] = true;
        } else if (!SAVE_AS && MODE == Inference.MARGINAL) { // set default files to save for marginal
            SAVE_AS_IDX[1] = SAVE_AS_IDX[3] = true;
        }

        for (int i = 0; i < SAVE_AS_IDX.length; i++) {
            if (SAVE_AS_IDX[i] && CONSENSUS[i]) {
                NEED_CONSENSUS = true;
                break;
            }
        }

        if (RUN_TRAVIS) {
            NEED_CONSENSUS = true;
        }
    }

    private static void parseRatesFile() {
        try {
            RATES = TSVFile.loadSubstitutionRatesFile(RATESFILE);
        } catch (IOException e) {
            usage(24, e.getMessage());
        } catch (NumberFormatException e) {
            usage(23, e.getMessage());
        }
    }

    private static Prediction setupIndelPrediction() {
        Prediction indelpred = null;
            try {
                indelpred = Prediction.load(INPUT + "/" + ASRFILE);
            } catch (ASRRuntimeException e) {
                usage(7, "Prediction failed to load: " + e.getMessage());
            } catch (IOException e) {
                usage(2, "Failed to read + " + INPUT + "/" + ASRFILE + ": "  + e.getMessage());
            }

        return indelpred;
    }

    private static Prediction performIndelInference(POGTree pogtree,
                                              EnumSeq.Alignment<Enumerable> aln) {
        return switch (INDEL_IDX) {
            case BEP -> Prediction.PredictByBidirEdgeParsimony(pogtree);
            case BEPML -> Prediction.PredictByBidirEdgeMaxLhood(pogtree);
            case SICP -> Prediction.PredictBySICP(pogtree);
            case SICML -> Prediction.PredictBySICML(pogtree);
            case PSP -> Prediction.PredictByParsimony(pogtree);
            case PSML -> Prediction.PredictByMaxLhood(pogtree);
            case SCIP, GUROBI -> Prediction.PredictByMIP(pogtree, aln, INDELS[INDEL_IDX], MODELS[MODEL_IDX],
                    GRASP.NTHREADS, GRASP.DISTANCE_BASED_MIP);
            default -> null;
        };
    }

    private static void performColumnInference(Prediction indelpred) {
        if (MODE == Inference.JOINT)
            indelpred.getJoint(MODEL, RATES);
        else if (MODE == Inference.MARGINAL) {
            if (indelpred.getTree().getIndex(MARG_NODE) < 0)
                usage(2, MARG_NODE + " is not a valid ancestor number");
            indelpred.getMarginal(MARG_NODE, MODEL, RATES);
        } else if (ONLYINDEL) {
            indelpred.saveIndelSolutionAsFasta(OUTPUT, PREFIX);
        }
    }

    private static POGraph[] storePOGsInArray(Map<Object, POGraph> pogs) {
        POGraph[] ancestors = new POGraph[pogs.size()];
        try {
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                if (MODE == Inference.MARGINAL) {
                    ancestors[0] = entry.getValue();
                    break;
                }
                ancestors[(Integer) entry.getKey()] = entry.getValue();
            }
        } catch (NumberFormatException exc) {
            int ii = 0;
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
                ancestors[ii++] = entry.getValue();
        }

        return ancestors;
    }

    public static void extractAncestralSequences(Object[][] ancseqs_gappy, Object[][] ancseqs_nogap,
                                                  Map<Object, POGraph> pogs, Prediction indelpred, String[] ancnames) {

        int ii = 0;
        try {
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                if (MODE == Inference.MARGINAL) {
                    ancnames[0] = "N" + entry.getKey().toString();
                    ancseqs_gappy[0] = indelpred.getSequence(entry.getKey(), MODE, true);
                    ancseqs_nogap[0] = indelpred.getSequence(entry.getKey(), MODE, false);
                    break;
                }
                ancnames[(Integer) entry.getKey()] = "N" + entry.getKey().toString();
                ancseqs_gappy[(Integer) entry.getKey()] = indelpred.getSequence(entry.getKey(), MODE, true);
                ancseqs_nogap[(Integer) entry.getKey()] = indelpred.getSequence(entry.getKey(), MODE, false);
                ii++;
            }
        } catch (NumberFormatException exc) {
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                ancnames[ii] = "N" + entry.getKey().toString();
                ancseqs_gappy[ii] = indelpred.getSequence(entry.getKey(), MODE, true);
                ancseqs_nogap[ii++] = indelpred.getSequence(entry.getKey(), MODE, false);
            }
        }
    }

    private static void saveGraspOutput(String[] ancnames, Object[][] ancseqs_nogap,  Object[][] ancseqs_gappy,
                                        Prediction indelpred, IdxTree tree, EnumSeq.Alignment<Enumerable> aln,
                                        POGraph[] ancestors) {
        File file = new File(OUTPUT);
        file.mkdirs();// true if the directory was created, false otherwise

        try {
            for (int i = 0; i < SAVE_AS_IDX.length; i++) {
                if (!SAVE_AS_IDX[i])
                    continue;
                switch (i) {
                    case FASTA: // FASTA
                        if (!BYPASS && MODE != null) {
                            saveGraspOutputAsFasta(ancnames, ancseqs_nogap, ancseqs_gappy);
                        }
                        break;
                    case DISTRIB:
                        if (!BYPASS && MODE == Inference.MARGINAL) { // must be true for this format
                            saveGraspOutputAsDistrib(indelpred);
                        }
                        break;
                    case CLUSTAL:
                        if (!BYPASS && MODE != null) {
                            saveGraspOutputAsClustal(ancnames, ancseqs_gappy);
                        }
                        break;
                    case TREE:
                        saveGraspOutputAsTree(indelpred, tree, false);
                        break;
                    case POGS:
                        if (!BYPASS) {
                            indelpred.save(OUTPUT + "/" + ASRFILE);
                        }
                        break;
                    case DOT:
                        if (!BYPASS) {
                            saveGraspOutputAsDOT(ancestors);
                        }
                        break;
                    case TREES:
                        saveGraspOutputAsTree(indelpred, tree, true);
                        break;

                    case POAG:
                        if (BYPASS) {
                            saveGraspOutputAsPOAG(tree, aln);
                        }
                        break;
                }
            }
        } catch (ASRException e) {
            usage(22, "Invalid input for ASR: " + e.getMessage());
        } catch (IOException e) {
            usage(2, "Failed to read or write files: " + e.getMessage());
        }
    }


    private static void saveGraspOutputAsFasta(String[] ancnames, Object[][] ancseqs_nogap,
                                               Object[][] ancseqs_gappy) throws IOException {
        FastaWriter fw = null;
        if (MODE == Inference.MARGINAL) // just one sequence
            fw = new FastaWriter(new File(OUTPUT, PREFIX + "_N" + MARG_NODE + ".fa"));
        else if (MODE == Inference.JOINT)
            fw = new FastaWriter(new File(OUTPUT, PREFIX + "_ancestors.fa"));
        if (GAPPY)
            fw.save(ancnames, ancseqs_gappy);
        else
            fw.save(ancnames, ancseqs_nogap);
        fw.close();
    }

    private static void saveGraspOutputAsDistrib(Prediction indelpred) throws IOException {
        EnumDistrib[] d = indelpred.getMarginal(MARG_NODE, MODEL, RATES);
        if (d != null) {
            Object[][] m = new Object[d.length + 1][];
            for (int j = 0; j < d.length; j++) {
                if (d[j] != null) {
                    m[j + 1] = new Object[MODEL.getDomain().size() + 1];
                    m[j + 1][0] = j + 1;
                    if (m[0] == null) {
                        m[0] = new Object[MODEL.getDomain().size() + 1];
                        m[0][0] = "Index";
                    }
                    for (int jj = 0; jj < m[j + 1].length - 1; jj++) {
                        m[j + 1][jj + 1] = d[j].get(jj);
                        if (m[0][jj + 1] == null)
                            m[0][jj + 1] = MODEL.getDomain().get(jj);
                    }
                }
            }
            for (int j = 0; j < d.length; j++) {
                if (d[j] == null) {
                    m[j + 1] = new Object[m[0].length];
                    m[j + 1][0] = j + 1;
                    for (int jj = 0; jj < m[j + 1].length - 1; jj++)
                        m[j + 1][jj + 1] = null;
                }
            }
            TSVFile.saveObjects(OUTPUT + "/" + PREFIX + "_N" + MARG_NODE + ".tsv", m);
        } else {
            usage(8, "Invalid ancestor node label: " + MARG_NODE);
        }
    }

    private static void saveGraspOutputAsClustal(String[] ancnames, Object[][] ancseqs_gappy) throws IOException {
        AlnWriter aw;
        if (MODE == Inference.MARGINAL) // just one sequence
            aw = new AlnWriter(new File(OUTPUT, PREFIX + "_N" + MARG_NODE + ".aln"));
        else
            aw = new AlnWriter(new File(OUTPUT, PREFIX + "_ancestors.aln"));
        aw.save(ancnames, ancseqs_gappy);
        aw.close();
    }

    private static void saveGraspOutputAsTree(Prediction indelpred, IdxTree tree, boolean saveAllTrees) throws IOException, ASRException {

        if (saveAllTrees) {
            if (MODE == Inference.JOINT)
                indelpred.saveTreeInstances(OUTPUT);
            else if (MODE == Inference.MARGINAL)
                usage(9, "Instantiations of position specific trees not available from marginal inference");
        } else {
            if (indelpred == null)
                Newick.save(tree, OUTPUT + "/" + PREFIX + "_ancestors.nwk", Newick.MODE_ANCESTOR);
            else
                Newick.save(indelpred.getTree(), OUTPUT + "/" + PREFIX + "_ancestors.nwk", Newick.MODE_ANCESTOR);
        }
    }

    private static void saveGraspOutputAsDOT(POGraph[] ancestors) throws ASRException, IOException {
        Map<Object, IdxGraph> saveme2 = new HashMap<>();
        for (int idx = 0; idx < ancestors.length; idx++) {
            ancestors[idx].setName("N" + idx);
            saveme2.put("N" + idx, ancestors[idx]);
        }
        IdxGraph.saveToDOT(OUTPUT, saveme2);
    }

    private static void saveGraspOutputAsPOAG(IdxTree tree, EnumSeq.Alignment<Enumerable> aln) throws IOException {
        int bpidx = 0; // default root
        if (MARG_NODE != null)
            bpidx = tree.getIndex(MARG_NODE);
        if (bpidx > 0) {
            List<EnumSeq> select = new ArrayList<>();
            String[] names = aln.getNames();
            for (int idx : tree.getLeaves(bpidx)) {
                Object label = tree.getLabel(idx);
                for (int ii = 0; ii < names.length; ii++) {
                    if (names[ii].equals(label.toString())) {
                        EnumSeq.Gappy seq = aln.getEnumSeq(ii);
                        select.add(seq);
                    }
                }
            }
            aln = new EnumSeq.Alignment(select);
        }
        POAGraph poag = new POAGraph(aln);
        if (VERBOSE)
            System.out.println("Saved POAG with " + aln.getHeight() + " sequences, under ancestor N" + MARG_NODE);
        poag.saveToDOT(OUTPUT + "/" + PREFIX + "_POAGunderN" + MARG_NODE + ".dot");
    }

    private static void runTravis( EnumSeq.Alignment<Enumerable> aln, Tree tree, Prediction indelpred,
                                   Object[][] ancseqs_nogap, Object[][] ancseqs_gappy) {
        if (LEARN) {
            if (TRAVIS_LEARN_NO_RECON) {
                try {
                    aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
                    tree = Utils.loadTree(NEWICK);
                    Utils.checkData(aln, tree, true);

                } catch (ASRException e) {
                    TrAVIS.usage(22, "Invalid input for TrAVIS " + e.getMessage());
                } catch (IOException e) {
                    TrAVIS.usage(2, "Failed to read or write files: " + e.getMessage());
                }

                assert aln != null;
                assert tree != null;
                TrAVIS.learnTreeParams(tree, 3, SEED);
                TrAVIS.printRootSeq(aln, tree, null);
                TrAVIS.learnIndelLengthDistributions(tree, aln, SEED, null);
                TrAVIS.learnIndelRateDistribution(tree, aln, null, SEED);
            } else {
                IdxTree mytree = indelpred.getTree();
                TrAVIS.learnTreeParams(mytree, 3, SEED);
                TrAVIS.printRootSeq(aln, mytree, ancseqs_nogap);
                System.out.println("--substitution-model " + MODEL.getName() + " \\");
                TrAVIS.learnIndelRateDistribution(tree, aln, ancseqs_gappy, SEED);
                TrAVIS.learnIndelLengthDistributions(tree, aln, SEED, ancseqs_gappy);
            }
        }

        if (PERFORM_TRAVIS_SIMUL) {
            // actually run the simulation
            EnumSeq rootSeq = TrAVIS.createRootSeq(ANCSEQ, MODEL, ANCSEQ_LENGTH, SEED);
            IdxTree simTree = TrAVIS.setupTree(
                    NEWICK,
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

            TrAVIS.TrackTree.Params params = TrAVIS.setupParams(
                    simTree,
                    rootSeq,
                    INDEL_LENGTH_MODEL,
                    INSERTION_LENGTH_MODEL,
                    DELETION_LENGTH_MODEL,
                    INDEL_RATE_MODEL,
                    SUBST_RATE_MODEL,
                    MODEL,
                    DELETIONPROP,
                    RATES,
                    SEED);

            TrAVIS.TrackTree tracker = new TrAVIS.TrackTree(params, SEED);
            EnumSeq[] seqs = tracker.getSequences();
            TrAVIS.saveOutput(seqs, tracker, simTree, TRAVIS_FORMAT_IDX, GAPPY, EXTANTS_ONLY, OUTPUT, PREFIX);
        }
    }
}
