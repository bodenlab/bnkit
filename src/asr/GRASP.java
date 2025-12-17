package asr;

import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;
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

/**
 * Command line version of GRASP.
 * @author mikael
 * @author ariane
 * @author gabe
 */
public class GRASP {

    public static String VERSION = "17-Dec-2025";

    public static boolean VERBOSE  = false;
    public static boolean TIME     = false;
    public static int     NTHREADS = 4;
    public static boolean NIBBLE   = true;
    public static boolean INDEL_CONSERVATIVE = true;
    // Mode for BEP
    public static boolean RECODE_NULL = true;
    public static int MIP_SOLVER_TIME_LIMIT_MINUTES = 720; // 12 hours
    public static boolean REMOVE_INDEL_ORPHANS = true;
    public enum Inference {
        JOINT,
        MARGINAL
    }

    public static void usage() {
        usage(0, null);
    }
    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        out.println("Usage: asr.GRASP \n" +
                "\t[-a | --aln <filename>]\n" +
                "\t[-n | --nwk <filename>]\n" +
                "\t{-o | --output-folder <foldername>} (default is current working folder, or input folder if available)\n" +
                "\t{-i | --input-folder <foldername>}\n" +
                "\t{-pre | --prefix <stub>}\n" +
                "\t{-rf | --rates-file <filename>}\n" +
                "\t{-s | --substitution-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}\n" +
                "\t{-t | --threads <number>}\n" +
                "\t{-j | --joint (default)}\n" +
                "\t{-m | --marginal <branchpoint-id>}\n" +
                "\t{--indel-method <methodname>} (select one from BEP(default) BEML SICP SICML PSP PSML SCIP Gurobi CPSAT)\n" +
                "\t{* --indel-prior <LOWGAP|MEDGAP|HIGHGAP>}\n" +
                "\t{--indel-rate-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>}\n" +
                "\t{--copy-rates Copy substitution rates from reconstructed ancestor\n" +
                "\t{--conflate-rates Modulate indel rate (rho) by site-specific substitution rate (r): p=e^(rho*r*t)\n" +
                "\t{--indel-length-distrib <ZeroTruncatedPoisson|ZTP|Poisson|Zipf|Lavalette>}\n"+
                "\t{--supported-path <methodname>} (select one from DIJKSTRA(default) ASTAR)\n" +
                "\t{--nogap}\n" +
                "\t{--seed <seed>}\n" +
                "\t{--nonibble}\n" +
                "\t{--exclude-noedge}\n" +
                "\t{--save-as <list-of-formats>} (select multiple from FASTA CLUSTAL TREE DISTRIB ASR DOT TREES TrAVIS SIMUL)\n" +
                "\t{--save-all} (saves reconstruction with ALL formats)\n" +
                "\t{--save-tree} (bypasses inference and re-saves the tree with ancestor nodes labelled as per GRASP's\n\tdepth-first labelling scheme starting with N0)\n" +
                "\t{--save-poag { <branchpoint-id> } (bypasses inference and saves the input alignment as a POAG\n\t(partial order alignment graph of extant sequences under specified ancestor [default N0])\n" +
                "\t{--time}{--verbose}{--help}\n");
        out.println("Inference is a two-stage process:\n" +
                "\t(1) A history of indel events is inferred by either maximum likelihood or maximum parsimony and \n\tmapped onto the tree to determine what positions contain actual sequence content\n" +
                "\t(2) For each ancestral position, the most probable character is assigned to each phylogenetic branch \n\tpoint when performing a joint reconstruction. Alternatively, for each \n\tposition at a nominated branch point, the probability distribution over all possible \n\tcharacters is inferred when performing a marginal reconstruction.\n" +
                "\tFinally, edges are drawn to represent all inferred combinations of indels to form an ancestor POG \n\twith nodes that can form a valid sequence with inferred content; a preferred path\n\tthrough the POG is then inferred, nominating a single, best supported sequence.\n");
        out.println("Mode of character inference:\n" +
                "\t-j (or --joint) activates joint reconstruction (default), \n\t-m (or --marginal) activates marginal reconstruction (requires a branch-point to be nominated)\n" +
                "\t--onlyindel disengages the stage of character state inference\n");
        out.println("Required arguments:\n" +
                "\t-a (or --aln) must specify the name of a multiple-sequence alignment file on FASTA or CLUSTAL format\n" +
                "\t-n (or --nwk) must specify the name of a phylogenetic-tree file on Newick format\n");
        out.println("Optional arguments:\n" +
                "\t-o (or --output-folder) specifies the folder that will be used to save output files,\n\t\te.g. inferred ancestor or ancestors, tree, etc. as specified by format\n" +
                "\t-i (or --input-folder) skips indel inference, and loads a previous reconstruction from specified folder\n" +
                "\t-sa (or --save-as) lists the files and formats to be generated (see below)\n\t--save-all nominates all\n" +
                "\t-pre (or --prefix) specifies a stub that is added to result filenames (default is the prefix of the alignment file)\n" +
                "\t-indel (or --indel-method) specifies what method to use for inferring indels (see below)\n" +
                "\t-s (or --substitution-model) specifies what evolutionary model to use for inferring character states (see below)\n" +
                "\t-rf (or --rates-file) specifies a tabulated file with relative, position-specific substitution rates\n\t\tWe recommend the use of this generally, but specifically for trees with great distances, and with biologically diverse entries\n\t\tAs an example, IQ-TREE produces rates on the accepted format with the --rate option (--mlrate is NOT supported yet).\n" +
                "\t--indel-prior (not implemented but intended for TrAVIS) specifies Gamma priors pre-determined from Pfam alignments that have few, moderate, or large numbers of gaps.\n" +
                "\t--indel-length-distrib specifies the indel length distribution function, which serves to model both deletions and insertions in TrAVIS\n" +
                "\t--indel-rate-distrib the indel rate distribution, which serves to model both insertions and deletions in TrAVIS\n"+
                "\t--include-extants means that extants are included in output files (when the format allows)\n" +
                "\t--nogap means that the gap-character is excluded in the resulting output (when the format allows)\n" +
                "\t--nonibble de-activates the removal of indices in partial order graphs that cannot form a path from start to end\n" +
                "\t--orphans de-activates the removal of orphaned indel trees\n" +
                "\t--exclude-noedge removes non-existing edge as an option for parsimony in BEP\n" +
                "\t--solver-time-limit the maximum time the MIP solver can run for in minutes before defaulting to BEP indel inference\n" +
                "\t--verbose prints out information about steps undertaken, and --time the time it took to finish\n" +
                "\t-h (or --help) will print out this screen\n");
        out.println("Files/formats: \n" +
                "\tFASTA: sequences (most preferred path at each ancestor, gapped or not gapped)\n" +
                "\tCLUSTAL: sequences (most preferred path at each ancestor, gapped)\n" +
                "\tTREE: phylogenetic tree with ancestor nodes labelled\n" +
                "\tDISTRIB: character distributions for each position (indexed by POG, only available for marginal reconstruction)\n" +
                "\tASR: complete reconstruction as JSON, incl. POGs of ancestors and extants, and tree (ASR.json)\n" +
                "\tDOT: partial-order graphs of ancestors in DOT format\n" +
                "\tTREES: position-specific trees with ancestor states labelled\n" +
                "\tTrAVIS: Produce commandline parameters for running TrAVIS\n" +
                "\tSIMUL: Run TrAVIS based on parameters from joint reconstruction, and save tree and alignments");
        out.println("Indel-methods: \n" +
                "\tBEP: bi-directional edge (maximum) parsimony\n" +
                "\tBEML: bi-directional edge maximum likelihood (uses uniform evolutionary model akin to JC)\n" +
                "\tSICP: simple indel-coding (maximum) parsimony (based on Simmons and Ochoterena)\n" +
                "\tSICML: simple indel-coding maximum likelihood (uses uniform evolutionary model)\n" +
                "\tPSP: position-specific (maximum) parsimony\n" +
                "\tPSML: position-specific maximum likelihood (uses uniform evolutionary model)\n" +
                "\tSCIP: infer a globally optimal indel history using the open-source SCIP solver (https://www.scipopt.org/). Does not support multi-threading\n" +
                "\tGurobi: uses the Gurobi solver to infer a globally optimal indel history. Requires local installation of Gurobi to run (https://www.gurobi.com/downloads/)\n" +
                "\tCPSAT: infer a globally optimal indel history using Google's open-source CP-SAT solver (https://developers.google.com/optimization). Should use a minimum of 8 threads for reliable performance\n" +
                "\tAdd '*' to method name for less conservative setting (if available)\n");
        out.println("Substitution-models: \n" +
                "\tJTT: Jones-Taylor-Thornton (protein; default)\n" +
                "\tDayhoff: Dayhoff-Schwartz-Orcutt (protein)\n" +
                "\tLG: Le-Gasquel (protein)\n" +
                "\tWAG: Whelan-Goldman (protein)\n" +
                "\tJC: Jukes-Cantor (DNA)\n" +
                "\tYang: Yang's general reversible process model (DNA)\n");
        out.println("Notes: \n" +
                "\tGreater number of threads may improve processing time up to a point when coordination chokes performance; default is 4 threads.\n" +
                "\tRunning GRASP requires large memory and in most cases Java needs to be run with the option -Xmx20g, \n\twhere 20g specifies that 20GB of RAM should be available.\n" +
                "\n~ This is version " + VERSION + " ~");
        if (msg != null)
            out.println("\n" + msg + " (Error " + error + ")");
        System.exit(error);
    }


    public static void main(String[] args) {

        boolean BYPASS = false; // bypass inference, default is false
        String ASRFILE = "ASR.json";
        String ALIGNMENT = null;
        String NEWICK = null;
        String OUTPUT = null;
        String INPUT = null;
        String PREFIX = null;
        String RATESFILE = null;
        double[] RATES = null;
        boolean CONFLATE_RATES = false;
        boolean COPY_SUBST_RATES = false;

        String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel MODEL = null;
        // Alphabet is decided by MODEL_IDX
        Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
        // Indel approaches:
        String[] INDELS = new String[] {"BEP", "BEML", "SICP", "SICML", "PSP", "PSML", "SCIP", "Gurobi", "CPSAT"};
        int INDEL_IDX = 0; // default indel approach is that above indexed 0
        String INDEL_RATE_DISTRIB = null;
        String INDEL_LENGTH_DISTRIB = null;
        String[] SPATH = new String[] {"DIJKSTRA", "ASTAR"};
        int SPATH_IDX = 0; // default supported path approach is that above indexed 0
        boolean GAPPY = true;
        // output formats
        boolean SAVE_AS = false;
        boolean INCLUDE_EXTANTS = false;
        String[]  FORMATS    = new String[]  {"FASTA", "DISTRIB", "CLUSTAL", "TREE", "ASR", "DOT", "TREES", "MATLAB", "LATEX", "POAG", "TrAVIS", "SIMUL"};
        // select these, default for "joint reconstruction"
        boolean[] SAVE_AS_IDX = new boolean[FORMATS.length];
        // select to compute consensus path for these output formats
        boolean[] CONSENSUS = new boolean[]  {true,    false,     true,      false,  false, false, false,   false,    false,   false,  true,    true  };
        // default inference mode
        Inference MODE = Inference.JOINT;
        // ancestor to reconstruct if inference mode is "marginal"
        Integer MARG_NODE = null;
        int SEED = new Random().nextInt();

        long START_TIME = System.currentTimeMillis(), ELAPSED_TIME;

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if ((arg.equalsIgnoreCase("-aln") || arg.equalsIgnoreCase("a")) && args.length > a + 1) {
                    ALIGNMENT = args[++ a];
                } else if ((arg.equalsIgnoreCase("-nwk") || arg.equalsIgnoreCase("n")) && args.length > a + 1) {
                    NEWICK = args[++ a];
                } else if ((arg.equalsIgnoreCase("-output-folder") || arg.equalsIgnoreCase("o")) && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if ((arg.equalsIgnoreCase("-input-folder") || arg.equalsIgnoreCase("i")) && args.length > a + 1) {
                    INPUT = args[++a];
                } else if ((arg.equalsIgnoreCase("-prefix") || arg.equalsIgnoreCase("pre")) && args.length > a + 1) {
                    PREFIX = args[++ a];
                } else if ((arg.equalsIgnoreCase("-rates-file") || arg.equalsIgnoreCase("rf")) && args.length > a + 1) {
                    RATESFILE = args[++a];
                } else if ((arg.equalsIgnoreCase("-seed")  && args.length > a + 1)) {
                    SEED = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-joint") || arg.equalsIgnoreCase("j")) {
                    MODE = Inference.JOINT;
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
                            if (args[a + 1].endsWith("*"))
                                INDEL_CONSERVATIVE = false;
                        }
                    }
                    if (!found_indel)
                        usage(3, args[a + 1] + " is not a valid indel approach for option --indel-method");
                } else if (arg.equalsIgnoreCase("-supported-path") && args.length > a + 1) {
                    boolean found_spath = false;
                    for (int i = 0; i < SPATH.length; i++) {
                        if (args[a + 1].startsWith(SPATH[i])) {
                            SPATH_IDX = i;
                            found_spath = true;
                        }
                    }
                    if (!found_spath)
                        usage(6, args[a + 1] + " is not a valid method for option --supported-path");
                } else if ((arg.equalsIgnoreCase("-save-as") || arg.equalsIgnoreCase("sa")) && args.length > a + 1) {
                    String format = "<none given>";
                    for (int a1 = a + 1; a1 < args.length; a1 ++) {
                        if (args[a1].startsWith("-"))
                            break;
                        format = args[a1];
                        boolean found_format = false;
                        for (int i = 0; i < FORMATS.length; i ++) {
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
                    for (int i = 0; i < FORMATS.length - 2; i ++)
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
                            a --;
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

                    /*
                     */
                } else if (arg.equalsIgnoreCase("-indel-rate-distrib") && args.length > a + 1) {
                    //  --indel-rate-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>
                    INDEL_RATE_DISTRIB = args[++a];
                } else if (arg.equalsIgnoreCase("-indel-length-distrib") && args.length > a + 1) {
                    //  --indel-length-distrib <ZeroTruncatedPoisson|ZTP|Poisson|Zipf|Lavalette>
                    INDEL_LENGTH_DISTRIB = args[++a];
                } else if (arg.equalsIgnoreCase("-copy-rates")) {
                    COPY_SUBST_RATES = true;
                } else if (arg.equalsIgnoreCase("-conflate-rates")) {
                    CONFLATE_RATES = true;
                } else if ((arg.equalsIgnoreCase("-threads") || arg.equalsIgnoreCase("t")) && args.length > a + 1) {
                    try {
                        NTHREADS = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        usage(2, "Failed to set number of threads for option --threads: " + args[a] + " is not a valid integer");
                    }
                } else if (arg.equalsIgnoreCase("-nogap")) {
                    GAPPY = false;
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

        if (ALIGNMENT == null && INPUT == null)
            usage(3, "Must specify alignment (--aln <Clustal or FASTA file>) or previously saved folder (--input-folder <folder>");
        else if (NEWICK == null && INPUT == null)
            usage(4, "Must specify phylogenetic tree (Newick file) or previously saved folder (--input-folder <folder>");
        else if (OUTPUT == null)
            OUTPUT = INPUT == null ? "." : INPUT;

        if (PREFIX == null) { // default prefix is the (prefix of) alignment filename
            int idx2 = ALIGNMENT == null ? 0 : ALIGNMENT.lastIndexOf(".");
            if (idx2 == -1)
                idx2 = ALIGNMENT.length();
            int idx1 = ALIGNMENT == null ? 0 : ALIGNMENT.lastIndexOf("/") + 1;
            PREFIX = ALIGNMENT == null ? "" : ALIGNMENT.substring(idx1, idx2);
        }
        MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
        if (MODEL == null)
            usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
        if (!SAVE_AS && MODE == Inference.JOINT) { // set default files to save for joint
            SAVE_AS_IDX[0] = SAVE_AS_IDX[3] = true;
        } else if (!SAVE_AS && MODE == Inference.MARGINAL) { // set default files to save for marginal
            SAVE_AS_IDX[1] = SAVE_AS_IDX[3] = true;
        }
        boolean NEED_CONSENSUS = false;
        for (int i = 0; i < SAVE_AS_IDX.length; i ++) {
            if (SAVE_AS_IDX[i] && CONSENSUS[i]) {
                NEED_CONSENSUS = true;
                break;
            }
        }
        if (RATESFILE != null) {
            try {
                TSVFile ratesfile = new TSVFile(RATESFILE, true);
                int rates_col = ratesfile.getColumn("Rate");
                int index_col = ratesfile.getColumn("Site");
                if (rates_col == -1)  // not there
                    rates_col = 0;
                Object[] rateobjs = ratesfile.getCol(rates_col);
                Object[] idxobjs = null;
                if (index_col != -1)
                    idxobjs = ratesfile.getCol(index_col);
                RATES = new double[rateobjs.length];
                for (int i = 0; i < RATES.length; i++) {
                    try {
                        int index = index_col == -1 ? i : (Integer) idxobjs[i] - 1; // starts with 1, so subtract "1" to use as position index
                        RATES[index] = (Double) rateobjs[i];
                    } catch (NumberFormatException e0) {
                        usage(23, "Rates file has invalid number format:" + rateobjs[i]);
                    }
                }
            } catch (IOException e) {
                usage(24, "Rates file could not be opened or read: " + RATESFILE);
            }
        }

        Object[][] ancseqs_gappy = null;
        Object[][] ancseqs_nogap = null;
        String[] ancnames = null;
        POGraph[] ancestors = null;
        Prediction indelpred = null;
        EnumSeq.Alignment aln = null;
        Tree tree = null;
        POGTree pogtree = null;

        try {
            if (INPUT != null) {
                try {
                    indelpred = Prediction.load(INPUT + "/" + ASRFILE);
                } catch (ASRRuntimeException e) {
                    usage(7, "Prediction failed to load: " + e.getMessage());
                }
            }
            START_TIME = System.currentTimeMillis();
            if (indelpred == null) {
                aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
                tree = Utils.loadTree(NEWICK);
                Utils.checkData(aln, tree);
            }
        } catch (ASRException e) {
            usage(22, "Invalid input for ASR: " + e.getMessage());
        } catch (IOException e) {
            usage(2, "Failed to read or write files: " + e.getMessage());
        }

        if (!BYPASS && indelpred == null) {
            // if we are past the above, we can assume that the data are good to process
            pogtree = new POGTree(aln, tree);
            switch (INDEL_IDX) {
                case 0:
                    indelpred = Prediction.PredictByBidirEdgeParsimony(pogtree);
                    break;
                case 1:
                    indelpred = Prediction.PredictByBidirEdgeMaxLhood(pogtree);
                    break;
                case 2:
                    indelpred = Prediction.PredictBySICP(pogtree);
                    break;
                case 3:
                    indelpred = Prediction.PredictBySICML(pogtree);
                    break;
                case 4:
                    indelpred = Prediction.PredictByParsimony(pogtree);
                    break;
                case 5:
                    indelpred = Prediction.PredictByMaxLhood(pogtree);
                    break;
                case 6, 7, 8:
                    indelpred = Prediction.PredictByMIP(pogtree, aln, INDELS[INDEL_IDX]);
                default:
                    break;
            }
        }

        if (!BYPASS) {
            if (indelpred == null)
                usage(3, INDELS[INDEL_IDX] + " is not implemented");
            if (MODE == Inference.JOINT)
                indelpred.getJoint(MODEL, RATES);
            else if (MODE == Inference.MARGINAL) {
                if (indelpred.getTree().getIndex(MARG_NODE) < 0)
                    usage(2, MARG_NODE + " is not a valid ancestor number");
                indelpred.getMarginal(MARG_NODE, MODEL, RATES);
            }
            POGraph.SUPPORTED_PATH_DEFAULT = SPATH_IDX;
            Map<Object, POGraph> pogs = indelpred.getAncestors(MODE);
            ancestors = new POGraph[pogs.size()];
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
            ancnames = new String[pogs.size()];
            if (NEED_CONSENSUS) {
                ancseqs_gappy = new Object[pogs.size()][];
                ancseqs_nogap = new Object[pogs.size()][];
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
            File file = new File(OUTPUT);
            if (file.mkdirs()) { // true if the directory was created, false otherwise
            } else {
                // System.err.println("Directory " + OUTPUT + " already exists");
                // throw new ASRException("Directory " + directory + " already exists");
            }
        }

        try {
            for (int i = 0; i < SAVE_AS_IDX.length; i++) {
                if (!SAVE_AS_IDX[i])
                    continue;
                switch (i) { // {"FASTA", "DISTRIB", "CLUSTAL", "TREE", "POGS", "DOT", "TREES", "MATLAB", "LATEX", "POAG", "TrAVIS", "SIMUL"};
                    case 0: // FASTA
                        if (!BYPASS && MODE != null) {
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
                        break;
                    case 1: // DISTRIB
                        if (!BYPASS && MODE == Inference.MARGINAL) { // must be true for this format
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
                            } else
                                usage(8, "Invalid ancestor node label: " + MARG_NODE);
                        }
                        break;
                    case 2: // CLUSTAL
                        if (!BYPASS && MODE != null) {
                            AlnWriter aw = null;
                            if (MODE == Inference.MARGINAL) // just one sequence
                                aw = new AlnWriter(new File(OUTPUT, PREFIX + "_N" + MARG_NODE + ".aln"));
                            else
                                aw = new AlnWriter(new File(OUTPUT, PREFIX + "_ancestors.aln"));
                            aw.save(ancnames, ancseqs_gappy);
                            aw.close();
                        }
                        break;
                    case 3: // TREE
                        if (indelpred == null)
                            Newick.save(tree, OUTPUT + "/" + PREFIX + "_ancestors.nwk", Newick.MODE_ANCESTOR);
                        else
                            Newick.save(indelpred.getTree(), OUTPUT + "/" + PREFIX + "_ancestors.nwk", Newick.MODE_ANCESTOR);
                        break;
                    case 4: // POGS
                        if (!BYPASS) {
                            String filename = OUTPUT + "/" + ASRFILE;
                            indelpred.save(filename);
                        }
                        break;
                    case 5: // DOT
                        if (!BYPASS) {
                            Map<Object, IdxGraph> saveme2 = new HashMap<>();
                            for (int idx = 0; idx < ancestors.length; idx++) {
                                ancestors[idx].setName("N" + idx);
                                saveme2.put("N" + idx, ancestors[idx]);
                            }
                            IdxGraph.saveToDOT(OUTPUT, saveme2);
                        }
                        break;
                    case 6: // TREES
                        if (MODE == Inference.JOINT)
                            indelpred.saveTreeInstances(OUTPUT);
                        else if (MODE == Inference.MARGINAL)
                            usage(9, "Instantiations of position specific trees not available from marginal inference");
                        break;
/*
                    case 7: // MATLAB
                        Map<Object, IdxGraph> saveme = new HashMap<>();
                        saveme.put("Extants", new POAGraph(aln));
                        for (Object name : aln.getNames())
                            saveme.put(name, pogtree.getExtant(name));
                        for (int idx = 0; idx < ancestors.length; idx++)
                            saveme.put("N" + idx, ancestors[idx]);
                        IdxGraph.saveToMatrix(OUTPUT, saveme);
                        break;
                    case 8: // LATEX
                        Map<Object, IdxGraph> saveme3 = new HashMap<>();
                        saveme3.put("*", new POAGraph(aln));
                        for (Object name : aln.getNames())
                            saveme3.put(name, pogtree.getExtant(name));
                        for (int idx = 0; idx < ancestors.length; idx++)
                            saveme3.put("N" + idx, ancestors[idx]);
                        IdxGraph.saveToLaTeX(OUTPUT, saveme3);
                        break;
 */
                    case 9: // POAG
                        if (BYPASS) {
                            int bpidx = 0; // default root
                            if (MARG_NODE != null)
                                bpidx = tree.getIndex(MARG_NODE);
                            if (bpidx > 0) {
                                List<EnumSeq> select = new ArrayList<>();
                                String[] names = aln.getNames();
                                for (int idx : tree.getLeaves(bpidx)) {
                                    Object label = tree.getLabel(idx);
                                    for (int ii = 0; ii < names.length; ii ++) {
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
                        break;
                    case 10: // TrAVIS: figure out params to run TrAVIS, produce report
                        if (!BYPASS) {
                            Random random = new Random(SEED);

                            // the user-provided tree for which ancestor sequences have been determined
                            IdxTree mytree = indelpred.getTree();
                            RateModel ddistrib = IdxTree.getGammaMixture(mytree, 3, SEED);
                            System.out.println("--dist-distrib " + ddistrib.getTrAVIS() + " \\");
                            GaussianDistrib l2rdistrib = mytree.getLeaf2RootDistrib();
                            IdxTree newrtree = Tree.Random(mytree.getNLeaves(), ddistrib, 2, 2, SEED);
                            newrtree.fitDistances(100, l2rdistrib, SEED + 202);
                            System.out.println("--leaf2root-distrib " + l2rdistrib.getTrAVIS() + " \\");

                            // calculate some v basic stats from the alignment itself
                            double gap_prop = (double) aln.getGapCount() / (double) (aln.getWidth() * aln.getHeight());
                            double gap_open_prop = (double) aln.getGapStartCount() / (double) (aln.getWidth() * aln.getHeight());
                            double gap_length = aln.getMeanGapLength();
                            if (VERBOSE) {
                                System.out.println("Gap proportion= " + gap_prop);
                                System.out.println("Gap opening proportion= " + gap_open_prop);
                                System.out.println("Mean gap length= " + gap_length);
                            }
                            StringBuilder n0 = new StringBuilder();
                            for (int j = 0; j < ancseqs_nogap[0].length; j++)
                                n0.append(ancseqs_nogap[0][j]);
                            System.out.println("--ancestor " + n0.toString() + " \\");

                            System.out.println("--substitution-model " + MODEL.getName() + " \\");
                            RateModel subst_rate_distrib = null;
                            if (RATES != null) { // position-specific rates available
                                double rates_alpha = GammaDistrib.getAlpha(RATES);
                                subst_rate_distrib = new GammaDistrib(rates_alpha, rates_alpha);
                                System.out.println("--subst-rate-distrib " + subst_rate_distrib.getTrAVIS() + " \\");
                            }
                            // GammaDistrib.Mixture mixture = Tree.getGammaMixture4Distances(mytree, 3);
                            // retrieve sequence names to alignment index map (extants only of course)
                            Map<String, Integer> extmap = aln.getMap();
                            // collect evolutionary distances for each branch so that a distribution can be estimated
                            Map<Integer, Double> dists = new HashMap<>();
                            // collect reconstructed/observed indel rates (based on the proportion of sequence they occupy and
                            // the distance from ancestor) on each branch so that a distribution can be estimated
                            Map<Integer, Double> indelrates = new HashMap<>();
                            //
                            int[] ins_total = new int[0];
                            int[] del_total = new int[0];
                            // Go through the tree, and look at each ancestor sequence, recording predicted indel events
                            for (int idx : mytree) { // go through the original, user-provided tree
                                int parent = mytree.getParent(idx);
                                if (parent != -1) {  // Non-root node, so there is a branch with distance to catch...
                                    double dist = mytree.getDistance(idx);
                                    dists.put(idx, dist); // collect evolutionary distance on the node (to its parent)
                                    // Retrieve parent/child reconstructed sequences at a node in the user-provided tree
                                    Object[] pseq = ancseqs_gappy[(Integer) mytree.getLabel(parent)];
                                    Object[] cseq;
                                    if (mytree.isLeaf(idx)) {
                                        EnumSeq.Gappy seq = aln.getEnumSeq(extmap.get(mytree.getLabel(idx)));
                                        cseq = seq.get();
                                    } else {
                                        cseq = ancseqs_gappy[(Integer) mytree.getLabel(idx)];
                                    }
                                    // Calculate indel rate for the reconstructed sequences in the user-provided tree
                                    int[] insertions = TrAVIS.getInsertionCounts(pseq, cseq);
                                    int nInsertions = Arrays.stream(insertions).sum();
                                    int[] deletions = TrAVIS.getDeletionCounts(pseq, cseq);
                                    int nDeletions = Arrays.stream(deletions).sum();
                                    double event_prop = (double) (nInsertions + nDeletions) / ancseqs_nogap[(Integer) mytree.getLabel(parent)].length;
                                    double indelrate = -Math.log(1.0 - event_prop) / dist;
                                    indelrates.put(idx, indelrate);
                                    // Accumulate insertion/deletion lengths
                                    ins_total = TrAVIS.mergeCounts(ins_total, insertions);
                                    del_total = TrAVIS.mergeCounts(del_total, deletions);
                                }
                            }
                            //
                            List<Double> collect = new ArrayList<>();
                            for (Map.Entry<Integer, Double> entry : indelrates.entrySet()) {
                                double d = dists.get(entry.getKey());
                                int nfractions = (int) Math.floor(d * 50); // could do 100, or 10... to adjust distance resolution (greater n means higher res)
                                double rate = entry.getValue();
                                for (int jj = 0; jj < nfractions; jj++) {
                                    collect.add(rate);
                                }
                            }
                            double[] rarray = new double[collect.size()];
                            for (int jj = 0; jj < rarray.length; jj++)
                                rarray[jj] = collect.get(jj);
                            // fit a zero-inflated or standard gamma distribution to all collected rates (from the reconstruction); seems to work OK but sometime less well than the mixture below
                            RateModel indelrateDist = RateModel.bestfit(rarray, SEED);
                            // Two options:
                            if (INDEL_RATE_DISTRIB != null) {
                                // (1) an indel rate distribution has been nominated so we will use it
                                try {
                                    indelrateDist = RateModel.bestfit(INDEL_RATE_DISTRIB, rarray, SEED);
                                } catch (RuntimeException e) {
                                    throw new RuntimeException("Invalid indel length distribution " + indelrateDist.getTrAVIS());
                                }
                            } else
                                // (2) we need to try all and pick the one with greatest log-likelihood
                                indelrateDist = RateModel.bestfit(rarray, SEED);
                            // fit a 3-component mixture distribution to all collected rates (from the reconstruction)
                            // ZeroInflatedGammaMix indelrateDist = ZeroInflatedGammaMix.fit(rarray, 2);
                            if (indelrateDist != null)
                                System.out.println("--indel-rate-distrib " + indelrateDist.getTrAVIS() + " \\");
                            // now, turn to indel lengths...
                            int[] indel_total = TrAVIS.mergeCounts(ins_total, del_total);
                            int[] ins_data = TrAVIS.unfoldCounts(ins_total);
                            int[] del_data = TrAVIS.unfoldCounts(del_total);
                            int[] indel_data = TrAVIS.unfoldCounts(indel_total);

                            // Now we can fit the indel distribution to the lengths of insertions and deletions
                            IndelModel indel_length_distrib = null;
                            IndelModel insertion_length_distrib = null;
                            IndelModel deletion_length_distrib = null;
                            // Two options:
                            if (INDEL_LENGTH_DISTRIB != null) {
                                // (1) an indel distribution has been nominated so we will use it
                                try {
                                    indel_length_distrib = IndelModel.bestfit(INDEL_LENGTH_DISTRIB, indel_data, SEED);
                                    insertion_length_distrib = IndelModel.bestfit(INDEL_LENGTH_DISTRIB, ins_data, SEED);
                                    deletion_length_distrib = IndelModel.bestfit(INDEL_LENGTH_DISTRIB, del_data, SEED);
                                } catch (RuntimeException e) {
                                    throw new RuntimeException("Invalid indel length distribution " + indel_length_distrib.getTrAVIS());
                                }
                            } else {
                                // (2) we need to try all and pick the one with greatest log-likelihood
                                indel_length_distrib = IndelModel.bestfit(indel_data, SEED);
                                insertion_length_distrib = IndelModel.bestfit(ins_data, SEED);
                                deletion_length_distrib = IndelModel.bestfit(del_data, SEED);
                            }
                            System.out.println("--indel-length-distrib " + indel_length_distrib.getTrAVIS() + " \\");
                            System.out.println("--insertion-length-distrib " + insertion_length_distrib.getTrAVIS() + " \\");
                            System.out.println("--deletion-length-distrib " + deletion_length_distrib.getTrAVIS() + " \\");
                            int ninsertions = Arrays.stream(ins_total).sum();
                            int ndeletions = Arrays.stream(del_total).sum();
                            int nindel = ninsertions + ndeletions;
                            double delprop = (double) ndeletions / (double) nindel;
                            System.out.printf("--delprop %.2f \\\n", delprop);
                            /* Previously the following distributions (with params) were used to determine data likelihood
                                Zipf                    1, 2, 5, 15
                                Lavalette               1, 2, 5, 15
                                ZeroTruncatedPoisson    0.01, 0.1, 1, 2                     */
                            if (VERBOSE)
                                System.out.print("Parent\tChild\tDist \t#Insert\t#Delete\t#Indel\tLength\tRate\tP(orig)\t");
                            int NSIMUL = 10; // number of simulations
                            for (int k = 0; k < NSIMUL; k ++)
                                if (VERBOSE)
                                    System.out.printf("%s\t%5s\t%s\t%s\t", "rate/"+(k+1), "p/"+(k+1), "#ind/"+(k+1), "len/"+(k+1));
                            int[] totIndels = new int[NSIMUL + 1];
                            int[] totLengths = new int[NSIMUL + 1];
                            for (int idx : mytree) { // go through the original, user-provided tree
                                int parent = mytree.getParent(idx);
                                if (parent != -1) {  // Non-root node, so there is a branch with distance to catch...
                                    double dist = mytree.getDistance(idx);
                                    // Retrieve parent/child reconstructed sequences at a node in the user-provided tree
                                    Object[] pseq = ancseqs_gappy[(Integer) mytree.getLabel(parent)];
                                    Object[] cseq;
                                    if (mytree.isLeaf(idx)) {
                                        EnumSeq.Gappy seq = aln.getEnumSeq(extmap.get(mytree.getLabel(idx)));
                                        cseq = seq.get();
                                    } else {
                                        cseq = ancseqs_gappy[(Integer) mytree.getLabel(idx)];
                                    }
                                    // Calculate indel rate for the reconstructed sequences in the user-provided tree
                                    int[] insertions = TrAVIS.getInsertionCounts(pseq, cseq);
                                    int nInsertions = Arrays.stream(insertions).sum();
                                    int[] deletions = TrAVIS.getDeletionCounts(pseq, cseq);
                                    int nDeletions = Arrays.stream(deletions).sum();
                                    int seqlength = ancseqs_nogap[(Integer) mytree.getLabel(parent)].length;
                                    double event_prop = (double) (nInsertions + nDeletions) / seqlength;
                                    double indelrate = -Math.log(1.0 - event_prop) / dist;
                                    double p_before = Math.exp(-indelrate * dist);
                                    int[] indels = TrAVIS.mergeCounts(insertions, deletions);
                                    int lengths_total = Arrays.stream(TrAVIS.unfoldCounts(indels)).sum();
                                    if (VERBOSE)
                                        System.out.printf("%5d\t%5d\t%5.3f\t%5d\t%5d\t%5d\t%5d\t%+5.3f\t%+5.3f\t", parent, idx, dist, nInsertions, nDeletions, nInsertions + nDeletions, lengths_total, indelrate, p_before);
                                    totIndels[0] += (nInsertions + nDeletions);
                                    totLengths[0] += lengths_total;
                                    for (int k = 0; k < NSIMUL; k ++) {
                                        double sample = indelrateDist.sample();
                                        double p = Math.exp(-(sample * dist));
                                        int nIndel = 0;
                                        int lenIndel = 0;
                                        for (int pos = 0; pos < seqlength; pos++) {
                                            if (random.nextDouble() >= p) {
                                                nIndel += 1;
                                                lenIndel += indel_length_distrib.sample();
                                            }
                                        }
                                        totIndels[k + 1] += nIndel;
                                        totLengths[k + 1] += lenIndel;
                                        if (VERBOSE)
                                            System.out.printf("%5.3f\t%5.3f\t%5d\t%5d\t", sample, p, nIndel, lenIndel);
                                    }
                                }
                                if (VERBOSE)
                                    System.out.println();
                            }
                            if (VERBOSE)
                                System.out.printf("     \t     \t     \t     \t      \t%5d\t%5d\t      \t      \t", totIndels[0], totLengths[0]);
                            double mean_indel = 0, mean_length = 0;
                            for (int k = 0; k < NSIMUL; k ++) {
                                if (VERBOSE)
                                    System.out.printf("     \t     \t%5d\t%5d\t", totIndels[k + 1], totLengths[k + 1]);
                                mean_indel += totIndels[k + 1] / NSIMUL;
                                mean_length += totLengths[k + 1] / NSIMUL;
                            }
                            double sd_indel = 0, sd_length = 0;
                            for (int k = 0; k < NSIMUL; k ++) {
                                sd_indel += Math.pow(totIndels[k + 1] - mean_indel, 2) / NSIMUL;
                                sd_length += Math.pow(totLengths[k + 1] - mean_length, 2) / NSIMUL;
                            }
                            sd_indel = Math.sqrt(sd_indel);
                            sd_length = Math.sqrt(sd_length);
                            if (VERBOSE)
                                System.out.printf("\t\t\t%5.3f\t%5.3f\t(%5.3f)\t%5.3f\t%5.3f\t(%5.3f)\n", mean_indel, Math.abs(mean_indel - (double)totIndels[0]), sd_indel, mean_length, Math.abs(mean_length - (double)totLengths[0]), sd_length);
                        }
                        break;

                    case 11: // SIMUL: simulate a complete reconstruction with the same pseudo-biological proerties as the present reconstruction
                        if (!BYPASS) {
                            if (MODE == Inference.JOINT) {
                                RateModel ddistrib = IdxTree.getGammaMixture(tree, 3, SEED);
                                GaussianDistrib l2rdistrib = tree.getLeaf2RootDistrib();
                                IdxTree newrtree = Tree.Random(tree.getNLeaves(), ddistrib, 2, 2, SEED);
                                newrtree.fitDistances(100, l2rdistrib, SEED + 202);
//                                IdxTree newrtree = IdxTree.generateTreeFromMixture(tree, 3, SEED, 100);
                                IdxTree mytree = indelpred.getTree();
                                StringBuilder n0 = new StringBuilder();
                                for (int j = 0; j < ancseqs_nogap[0].length; j++)
                                    n0.append(ancseqs_nogap[0][j]);
                                if (VERBOSE)
                                    System.out.println("N0= " + n0.toString());
                                Double rates_alpha = 0.01;
                                if (RATES != null) { // position-specific rates available
                                    rates_alpha = GammaDistrib.getAlpha(RATES);
                                    if (VERBOSE)
                                        System.out.println("Position-specific rates Gamma shape= " + rates_alpha + " scale= " + 1.0/rates_alpha);
                                }
                                // retrieve sequence names to alignment index map (extants only of course)
                                Map<String, Integer> extmap = aln.getMap();
                                // collect evolutionary distances for each branch so that a distribution can be estimated
                                Map<Integer, Double> dists = new HashMap<>();
                                // collect reconstructed/observed indel rates (based on the proportion of sequence they occupy and
                                // the distance from ancestor) on each branch so that a distribution can be estimated
                                Map<Integer, Double> indelrates = new HashMap<>();
                                // collect indel lengths as well, distinguishing between insertions and deletions
                                int[] ins_total = new int[0];
                                int[] del_total = new int[0];
                                // Go through the tree, and look at each ancestor sequence, recording predicted indel events
                                for (int idx : mytree) { // go through the original, user-provided tree
                                    int parent = mytree.getParent(idx);
                                    if (parent != -1) {  // Non-root node, so there is a branch with distance to catch...
                                        double dist = mytree.getDistance(idx);
                                        dists.put(idx, dist); // collect evolutionary distance on the node (to its parent)
                                        // Retrieve parent/child reconstructed sequences at a node in the user-provided tree
                                        Object[] pseq = ancseqs_gappy[(Integer) mytree.getLabel(parent)];
                                        Object[] cseq;
                                        if (mytree.isLeaf(idx)) {
                                            EnumSeq.Gappy seq = aln.getEnumSeq(extmap.get(mytree.getLabel(idx)));
                                            cseq = seq.get();
                                        } else {
                                            cseq = ancseqs_gappy[(Integer) mytree.getLabel(idx)];
                                        }
                                        // Calculate indel rate for the reconstructed sequences in the user-provided tree
                                        int[] insertions = TrAVIS.getInsertionCounts(pseq, cseq);
                                        int nInsertions = Arrays.stream(insertions).sum();
                                        int[] deletions = TrAVIS.getDeletionCounts(pseq, cseq);
                                        int nDeletions = Arrays.stream(deletions).sum();
                                        double event_prop = (double) (nInsertions + nDeletions) / ancseqs_nogap[(Integer) mytree.getLabel(parent)].length;
                                        double indelrate = -Math.log(1.0 - event_prop) / dist;
                                        indelrates.put(idx, indelrate);
                                        // Accumulate insertion/deletion lengths
                                        ins_total = TrAVIS.mergeCounts(ins_total, insertions);
                                        del_total = TrAVIS.mergeCounts(del_total, deletions);
                                    }
                                }
                                //
                                List<Double> collect = new ArrayList<>();
                                for (Map.Entry<Integer, Double> entry : indelrates.entrySet()) {
                                    double d = dists.get(entry.getKey());
                                    int nfractions = (int) Math.floor(d * 50); // could do 100, or 10...
                                    double rate = entry.getValue();
                                    for (int jj = 0; jj < nfractions; jj++) {
                                        collect.add(rate);
                                    }
                                }
                                double[] rarray = new double[collect.size()];
                                for (int jj = 0; jj < rarray.length; jj++)
                                    rarray[jj] = collect.get(jj);
                                // fit a gamma distribution to all collected rates (from the reconstruction); seems to work OK but less well than the mixture below
                                // fit a zero-inflated or standard gamma distribution to all collected rates (from the reconstruction); seems to work OK but sometime less well than the mixture below
                                RateModel indelrateDist = RateModel.bestfit(rarray, SEED);
                                // Two options:
                                if (INDEL_RATE_DISTRIB != null) {
                                    // (1) an indel rate distribution has been nominated so we will use it
                                    try {
                                        indelrateDist = RateModel.bestfit(INDEL_RATE_DISTRIB, rarray, SEED);
                                    } catch (RuntimeException e) {
                                        throw new RuntimeException("Invalid indel length distribution " + indelrateDist.getTrAVIS());
                                    }
                                } else
                                    // (2) we need to try all and pick the one with greatest log-likelihood
                                    indelrateDist = RateModel.bestfit(rarray, SEED);
                                // now, turn to indel lengths...
                                int[] indel_total = TrAVIS.mergeCounts(ins_total, del_total);
                                int[] ins_data = TrAVIS.unfoldCounts(ins_total);
                                int[] del_data = TrAVIS.unfoldCounts(del_total);
                                int[] indel_data = TrAVIS.unfoldCounts(indel_total);
                                if (VERBOSE)
                                    System.out.println("Input tree has number indels = " + indel_data.length + " sum of lengths = " + Arrays.stream(indel_data).sum());
                                // Now we can fit the indel distribution to the lengths of insertions and deletions
                                IndelModel indel_length_distrib = null;
                                IndelModel insertion_length_distrib = null;
                                IndelModel deletion_length_distrib = null;
                                // Two options:
                                if (INDEL_LENGTH_DISTRIB != null) {
                                    // (1) an indel distribution has been nominated so we will use it
                                    try {
                                        indel_length_distrib = IndelModel.bestfit(INDEL_LENGTH_DISTRIB, indel_data, SEED);
                                        insertion_length_distrib = IndelModel.bestfit(INDEL_LENGTH_DISTRIB, ins_data, SEED);
                                        deletion_length_distrib = IndelModel.bestfit(INDEL_LENGTH_DISTRIB, del_data, SEED);
                                    } catch (RuntimeException e) {
                                        throw new RuntimeException("Invalid indel length distribution " + indel_length_distrib.getTrAVIS());
                                    }
                                } else {
                                    // (2) we need to try all and pick the one with greatest log-likelihood
                                    indel_length_distrib = IndelModel.bestfit(indel_data, SEED);
                                    insertion_length_distrib = IndelModel.bestfit(ins_data, SEED);
                                    deletion_length_distrib = IndelModel.bestfit(del_data, SEED);
                                }

                                int ninsertions = Arrays.stream(ins_total).sum();
                                int ndeletions = Arrays.stream(del_total).sum();
                                int nindel = ninsertions + ndeletions;
                                double delprop = (double) ndeletions / (double) nindel;

                                TrAVIS.TrackTree tracker = null;
                                EnumSeq[] seqs = null;
                                EnumSeq[] seqs_ex = null;
                                while (tracker == null) {
                                    TrAVIS.TrackTree.Params params = new TrAVIS.TrackTree.Params(newrtree, EnumSeq.parseProtein(n0.toString()), MODEL);
                                    if (RATES != null && COPY_SUBST_RATES) {
                                        // recover rates for N0, appropriately indexed
                                        int[] n0idxs = indelpred.getConsensus(0);
                                        double[] rates = new double[n0idxs.length];
                                        for (int k = 0; k < n0idxs.length; k ++)
                                            rates[k] = RATES[n0idxs[k]];
                                        params.setSubstRates(rates);
                                    } else if (RATES != null && !COPY_SUBST_RATES) {
                                        params.setSubstRateModel(rates_alpha);
                                    }
                                    params.setIndelRateModel(indelrateDist);
                                    params.setInsertmodel(insertion_length_distrib);
                                    params.setDeletemodel(deletion_length_distrib);
                                    params.PROPORTION_DELETION = delprop;
                                    params.SUBST_RATE_INFLUENCES_INDELS = CONFLATE_RATES;
                                    params.setSeed(SEED);

                                    tracker = new TrAVIS.TrackTree(params, SEED);
                                    seqs = tracker.getSequences();
                                    seqs_ex = tracker.getLeafSequences();
                                    for (EnumSeq seq : seqs) {
                                        if (seq.length() < 1) {
                                            tracker = null;
                                            SEED +=1;
                                            break;
                                        }
                                    }
                                    int[][] indels = tracker.getIndels();
                                    int[] ins_sim = indels[0];
                                    int[] del_sim = indels[1];
                                    int[] merged = TrAVIS.mergeCounts(ins_sim, del_sim);
                                    int[] flattened = TrAVIS.unfoldCounts(merged);
                                    if (VERBOSE)
                                        System.out.println("Simul tree has number indels = " + flattened.length + " sum of lengths = " + Arrays.stream(flattened).sum());
                                }

                                EnumSeq[] aseqs = tracker.getAlignment();
                                EnumSeq[] aseqs_ex = tracker.getLeafAlignments();
                                try {
                                    File file = new File(OUTPUT);
                                    if (!file.exists()) {
                                        file.mkdirs();
                                    }
                                    FastaWriter fw = new FastaWriter(OUTPUT + "/" + PREFIX + "_full.fa");
                                    fw.save(GAPPY ? aseqs : seqs);
                                    fw.close();
                                    FastaWriter fw1 = new FastaWriter(OUTPUT + "/" + PREFIX + ".fa");
                                    fw1.save(GAPPY ? aseqs_ex : seqs_ex);
                                    fw1.close();
                                    Newick.save(newrtree, OUTPUT + "/" + PREFIX + ".nwk", Newick.MODE_DEFAULT);

                                } catch (IOException e) {
                                    usage(7, "Something went wrong saving files in directory");
                                }

                            } else if (MODE == Inference.MARGINAL)
                                usage(23, "TrAVIS reports must be based on joint reconstructions");
                        }
                        break;
                }
                ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
                if (VERBOSE || TIME) {
                    System.out.printf("Done in %d min, %d sec%n", TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                            TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME)));
                }
            }
        } catch (ASRException e) {
            usage(22, "Invalid input for ASR: " + e.getMessage());
        } catch (IOException e) {
            usage(2, "Failed to read or write files: " + e.getMessage());
        }
    }

}
