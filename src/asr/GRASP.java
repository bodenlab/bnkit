package asr;

import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.*;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.pog.IdxGraph;
import dat.pog.POAGraph;
import dat.pog.POGTree;
import dat.pog.POGraph;
import json.JSONObject;
import json.JSONArray;
import smile.stat.distribution.GammaDistribution;
import smile.stat.distribution.GaussianDistribution;
import smile.stat.distribution.Mixture;
import stats.*;
import java.io.IOException;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.DoubleStream;

import static bn.prob.GammaDistrib.getAlpha;
import static bn.prob.GammaDistrib.getBeta;
import static stats.PowerLawCon.computeMinMax;
import static stats.PowerLawCon.estimateAlpha;
import smile.stat.distribution.ExponentialFamilyMixture;

class BranchInfo {
    String from;
    String to;
    double rate;
    double dist;
    double cumulativeLength;

    int depth;
    @Override
    public String toString() {
        return from + "→" + to + " : rate=" + rate + ", dist=" + dist;
    }
}




/**
 * Command line version of GRASP.
 * @author mikael
 * @author ariane
 * @author gabe
 */
public class GRASP {

    public static String VERSION = "9-Apr-2025";

    public static boolean VERBOSE  = false;
    public static boolean TIME     = false;
    public static int     NTHREADS = 4;
    public static boolean NIBBLE   = true;
    public static boolean INDEL_CONSERVATIVE = true;
    // Mode for BEP
    public static boolean RECODE_NULL = true;

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
                "\t{--indel-method <methodname>} (select one from BEP(default) BEML SICP SICML PSP PSML)\n" +
                "\t{--supported-path <methodname>} (select one from DIJKSTRA(default) ASTAR)\n" +
                "\t{--nogap}\n" +
                "\t{--seed <seed>}\n" +
                "\t{--nonibble}\n" +
                "\t{--exclude-noedge}\n" +
                "\t{--save-as <list-of-formats>} (select multiple from FASTA CLUSTAL TREE DISTRIB ASR DOT TREES TrAVIS)\n" +
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
                "\t-rf (or --rates-file) specifies a tabulated file with relative, position-specific substitution rates\n\t\tWe recommend the use of this generally, but specifically for trees with great distances, and with biologically diverse entries\n\t\tAs an example, IQ-TREE produces rates on the accepted format\n" +
                "\t--include-extants means that extants are included in output files (when the format allows)\n" +
                "\t--nogap means that the gap-character is excluded in the resulting output (when the format allows)\n" +
                "\t--nonibble de-activates the removal of indices in partial order graphs that cannot form a path from start to end\n" +
                "\t--orphans de-activates the removal of orphaned indel trees\n" +
                "\t--exclude-noedge removes non-existing edge as an option for parsimony in BEP\n" +
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
                "\tTrAVIS: Produce a report for reconstruction\n");
        out.println("Indel-methods: \n" +
                "\tBEP: bi-directional edge (maximum) parsimony\n" +
                "\tBEML: bi-directional edge maximum likelihood (uses uniform evolutionary model akin to JC)\n" +
                "\tSICP: simple indel-coding (maximum) parsimony (based on Simmons and Ochoterena)\n" +
                "\tSICML: simple indel-coding maximum likelihood (uses uniform evolutionary model)\n" +
                "\tPSP: position-specific (maximum) parsimony\n" +
                "\tPSML: position-specific maximum likelihood (uses uniform evolutionary model)\n" +
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

    // Compute ranks for a list of values, handling ties by assigning average ranks
    public static List<Double> getRanks(List<Double> values) {
        int n = values.size();
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < n; i++) indices.add(i);

        // Sort indices based on corresponding values
        indices.sort((i, j) -> Double.compare(values.get(i), values.get(j)));

        // Initialize all ranks with 0.0
        List<Double> ranks = new ArrayList<>(Collections.nCopies(n, 0.0));

        int i = 0;
        while (i < n) {
            int j = i;
            // Find all values that are equal (tie group)
            while (j + 1 < n && Double.compare(values.get(indices.get(i)), values.get(indices.get(j + 1))) == 0) {
                j++;
            }
            // Compute average rank for the tie group
            double rank = (i + j + 2) / 2.0; // +2 because rank starts from 1
            for (int k = i; k <= j; k++) {
                ranks.set(indices.get(k), rank);
            }
            i = j + 1;
        }
        return ranks;
    }

    // Compute Spearman rank correlation between two lists
    public static double spearman(List<Double> x, List<Double> y) {
        if (x.size() != y.size()) throw new IllegalArgumentException("Lists must be the same length");
        int n = x.size();
        List<Double> rankX = getRanks(x);
        List<Double> rankY = getRanks(y);

        double sum = 0;
        for (int i = 0; i < n; i++) {
            double d = rankX.get(i) - rankY.get(i);
            sum += d * d;
        }

        // Spearman's rank correlation coefficient formula
        return 1 - (6 * sum) / (n * (n * n - 1));
    }

    // Optimize the order of rates within each depth layer to improve Spearman correlation with ground truth
    public static void localSwapOptimize(List<BranchInfo> branchInfos, Map<String, Double> groundTruthRates, int nIter, int seed) {
        Random rand = new Random(seed);

        // Group branches by depth
        Map<Integer, List<BranchInfo>> byDepth = new HashMap<>();
        for (BranchInfo bi : branchInfos) {
            byDepth.computeIfAbsent(bi.depth, k -> new ArrayList<>()).add(bi);
        }

        // Perform nIter random swaps to try to increase Spearman correlation
        for (int iter = 0; iter < nIter; iter++) {
            for (Map.Entry<Integer, List<BranchInfo>> entry : byDepth.entrySet()) {
                List<BranchInfo> layer = entry.getValue();
                int n = layer.size();
                if (n < 2) continue;

                // Randomly select two different indices within the same depth
                int i = rand.nextInt(n), j = rand.nextInt(n);
                while (i == j) j = rand.nextInt(n);

                BranchInfo bi1 = layer.get(i);
                BranchInfo bi2 = layer.get(j);

                // Compute current Spearman correlation
                List<Double> simList = new ArrayList<>();
                List<Double> groundList = new ArrayList<>();
                for (BranchInfo b : branchInfos) {
                    simList.add(b.rate);
                    String key = b.from + "→" + b.to;
                    groundList.add(groundTruthRates.getOrDefault(key, 0.0));
                }
                double currentRho = spearman(simList, groundList);

                // Swap the rates
                double tmp = bi1.rate;
                bi1.rate = bi2.rate;
                bi2.rate = tmp;

                // Compute new Spearman correlation
                simList.clear(); groundList.clear();
                for (BranchInfo b : branchInfos) {
                    simList.add(b.rate);
                    String key = b.from + "→" + b.to;
                    groundList.add(groundTruthRates.getOrDefault(key, 0.0));
                }
                double newRho = spearman(simList, groundList);

                // If the swap made correlation worse, undo the swap
                if (newRho < currentRho) {
                    bi2.rate = bi1.rate;
                    bi1.rate = tmp;
                }
            }
        }
    }

    public static void globalSwapOptimize(List<BranchInfo> branchInfos, Map<String, Double> groundTruthRates, int nIter, int seed) {
        Random rand = new Random(seed);

        int n = branchInfos.size();
        if (n < 2) return;

        for (int iter = 0; iter < nIter; iter++) {
            // Randomly select two distinct branches from the entire list
            int i = rand.nextInt(n), j = rand.nextInt(n);
            while (i == j) j = rand.nextInt(n);

            BranchInfo bi1 = branchInfos.get(i);
            BranchInfo bi2 = branchInfos.get(j);

            // Compute current Spearman correlation between simulated and ground truth rates
            List<Double> simList = new ArrayList<>();
            List<Double> groundList = new ArrayList<>();
            for (BranchInfo b : branchInfos) {
                simList.add(b.rate);
                String key = b.from + "→" + b.to;
                groundList.add(groundTruthRates.getOrDefault(key, 0.0));
            }
            double currentRho = spearman(simList, groundList);

            // Swap the rates of the two selected branches
            double tmp = bi1.rate;
            bi1.rate = bi2.rate;
            bi2.rate = tmp;

            // Recompute the Spearman correlation after the swap
            simList.clear();
            groundList.clear();
            for (BranchInfo b : branchInfos) {
                simList.add(b.rate);
                String key = b.from + "→" + b.to;
                groundList.add(groundTruthRates.getOrDefault(key, 0.0));
            }
            double newRho = spearman(simList, groundList);

            // Revert the swap if the new correlation is worse
            if (newRho < currentRho) {
                bi2.rate = bi1.rate;
                bi1.rate = tmp;
            }
        }
    }
    


    static double logLikelihood(double[] data, ZeroInflatedGammaMix model) {
        double logL = 0.0;
        for (double val : data) {
            double p = model.p(val);
            p = Math.max(p, 1e-12);
            logL += Math.log(p);
        }
        return logL;
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

        String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel MODEL = null;
        // Alphabet is decided by MODEL_IDX
        Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
        // Indel approaches:
        String[] INDELS = new String[] {"BEP", "BEML", "SICP", "SICML", "PSP", "PSML"};
        int INDEL_IDX = 0; // default indel approach is that above indexed 0
        String[] SPATH = new String[] {"DIJKSTRA", "ASTAR"};
        int SPATH_IDX = 0; // default supported path approach is that above indexed 0
        boolean GAPPY = true;
        // output formats
        boolean SAVE_AS = false;
        boolean INCLUDE_EXTANTS = false;
        String[]  FORMATS    = new String[]  {"FASTA", "DISTRIB", "CLUSTAL", "TREE", "ASR", "DOT", "TREES", "MATLAB", "LATEX", "POAG", "TrAVIS"};
        // select these, default for "joint reconstruction"
        boolean[] SAVE_AS_IDX = new boolean[FORMATS.length];
        // select to compute consensus path for these output formats
        boolean[] CONSENSUS = new boolean[]  {true,    false,     true,      false,  false, false, false,   false,    false,   false,  true  };
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
                } else if ((arg.equalsIgnoreCase("-threads") || arg.equalsIgnoreCase("t")) && args.length > a + 1) {
                    try {
                        NTHREADS = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        System.err.println("Failed to set number of threads for option --threads: " + args[a] + " is not a valid integer");
                    }
                } else if (arg.equalsIgnoreCase("-nogap")) {
                    GAPPY = false;
                } else if (arg.equalsIgnoreCase("-verbose")) {
                    VERBOSE = true;
                } else if (arg.equalsIgnoreCase("-time")) {
                    TIME = true;
                } else if (arg.equalsIgnoreCase("-nonibble")) {
                    NIBBLE = false;
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
            int idx = ALIGNMENT == null ? 0 : ALIGNMENT.indexOf(".");
            if (idx == -1)
                idx = ALIGNMENT.length();
            PREFIX = ALIGNMENT == null ? "" : ALIGNMENT.substring(0, idx);
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
                switch (i) { // {"FASTA", "DISTRIB", "CLUSTAL", "TREE", "POGS", "DOT", "TREES", "MATLAB", "LATEX", "POAG", "TRAVIS"};
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

                    case 10: // TrAVIS
                        if (!BYPASS) {
                            if (MODE == Inference.JOINT) {

                                double gap_prop = (double) aln.getGapCount() / (double) (aln.getWidth() * aln.getHeight());
                                double gap_open_prop = (double) aln.getGapStartCount() / (double) (aln.getWidth() * aln.getHeight());
                                double gap_length = aln.getMeanGapLength();
                                //double indel_factor = 10; // FIXME: this is the default value
                                if (VERBOSE) {
                                    System.out.println("Gap proportion= " + gap_prop);
                                    System.out.println("Gap opening proportion= " + gap_open_prop);
                                    System.out.println("Mean gap length= " + gap_length);
                                }
                                Tree newrtree = Tree.generateTreeFromMixture(tree, 3, SEED, 100);
                                IdxTree mytree = indelpred.getTree();
                                Set<Double> dists = new HashSet<>();
                                Map<String, Integer> extmap = aln.getMap();
                                String[] names = aln.getNames();
                                int[] ins_total = new int[0];
                                int[] del_total = new int[0];
                                List<Double> rList = new ArrayList<>();
                                List<BranchInfo> branchInfos = new ArrayList<>();
                                Map<String, Double> branchRateMap = new LinkedHashMap<>();

                                String filenamePrefix = (PREFIX == null || PREFIX.isEmpty()) ? "result" : PREFIX;
                                String out = OUTPUT + "/" + filenamePrefix;
                                try (PrintWriter writer = new PrintWriter(new FileWriter(out + "_branch_length.txt"))) {
                                    writer.println("Branch\tLength");

                                    int n = newrtree.getSize();
                                    for (int idx = 0; idx < n; idx++) {

                                        int parent = newrtree.getParent(idx);
                                        if (parent == -1) continue;

                                        String parentLabel = newrtree.getLabel(parent).toString();
                                        String childLabel = newrtree.getLabel(idx).toString();
                                        double dist = newrtree.getDistance(idx);


                                        if (parentLabel == null || parentLabel.isEmpty()) parentLabel = String.valueOf(parent);
                                        if (childLabel == null || childLabel.isEmpty()) childLabel = String.valueOf(idx);

                                        writer.printf("%s→%s\t%.5f\n", parentLabel, childLabel, dist);
                                    }

                                    System.out.println("Tree structure saved to " + out + "_branch_length.txt");

                                } catch (IOException e) {
                                    System.err.println("Failed to write tree: " + e.getMessage());
                                }


                          
                                int rootIdx = -1;
                                for (int idx : mytree) {
                                    if (mytree.getParent(idx) == -1) {
                                        rootIdx = idx;
                                        break;
                                    }
                                }
                                if (rootIdx == -1) throw new RuntimeException("Root node not found.");
                                String rootLabel = mytree.getLabel(rootIdx).toString();
                                Map<String, Double> cumulativeDistMap = new HashMap<>();
                                cumulativeDistMap.put(rootLabel, 0.0);
                                for (int idx : mytree) {
                                    int parent = mytree.getParent(idx);
                                    if (parent != -1) { // Non-root node
                                        double dist = mytree.getDistance(idx);
                                        dists.add(dist);

                                        String from = mytree.getLabel(parent).toString();
                                        String to = mytree.getLabel(idx).toString();

                                        // Accumulate branch length from root to current node
                                        double parentCumulative = cumulativeDistMap.getOrDefault(from, 0.0);
                                        double currentCumulative = parentCumulative + dist;
                                        cumulativeDistMap.put(to, currentCumulative);

                                        // Retrieve parent/child sequences
                                        Object[] pseq = ancseqs_gappy[(Integer) mytree.getLabel(parent)];
                                        Object[] cseq;
                                        if (mytree.isLeaf(idx)) {
                                            EnumSeq.Gappy seq = aln.getEnumSeq(extmap.get(mytree.getLabel(idx)));
                                            cseq = seq.get();
                                        } else {
                                            cseq = ancseqs_gappy[(Integer) mytree.getLabel(idx)];
                                        }

                                        // Calculate indel rate
                                        int[] insertions = TrAVIS.getInsertionCounts(pseq, cseq);
                                        int[] deletions = TrAVIS.getDeletionCounts(pseq, cseq);
                                        int[] Events = TrAVIS.calculateIndelOpening(pseq, cseq);
                                        double rate = TrAVIS.calculateRForNodes(Events, pseq.length, dist);
                                        rList.add(rate);

                                        // Record branch → rate
                                        String label = from + "→" + to;
                                        branchRateMap.put(label, rate);

                                        // Save to branchInfos
                                        BranchInfo bi = new BranchInfo();
                                        bi.from = from;
                                        bi.to = to;
                                        bi.rate = rate;
                                        bi.dist = dist;
                                        bi.cumulativeLength = currentCumulative;
                                        branchInfos.add(bi);

                                        // Accumulate insertions/deletions (unchanged)
                                        int[] ins_tmp = new int[Math.max(insertions.length, ins_total.length)];
                                        for (int j = 0; j < ins_tmp.length; j++) {
                                            ins_tmp[j] += j < insertions.length ? insertions[j] : 0;
                                            ins_tmp[j] += j < ins_total.length ? ins_total[j] : 0;
                                        }
                                        ins_total = ins_tmp;

                                        int[] del_tmp = new int[Math.max(deletions.length, del_total.length)];
                                        for (int j = 0; j < del_tmp.length; j++) {
                                            del_tmp[j] += j < deletions.length ? deletions[j] : 0;
                                            del_tmp[j] += j < del_total.length ? del_total[j] : 0;
                                        }
                                        del_total = del_tmp;
                                    }
                                }
                                Map<String, Integer> nodeDepthMap = new HashMap<>();
                                for (BranchInfo bi : branchInfos) {
                                    int depth = nodeDepthMap.getOrDefault(bi.to, 0); // child 决定这条 branch 的深度
                                    System.out.println(bi.from + " -> " + bi.to + " : depth = " + depth);
                                }

                                Map<String, Double> groundTruthRates = new HashMap<>();
                                for (BranchInfo bi : branchInfos) {
                                    String label = bi.from + "→" + bi.to;
                                    groundTruthRates.put(label, bi.rate);
                                }

                                int[] indel_total = new int[Math.max(ins_total.length, del_total.length)];
                                for (int j = 0; j < indel_total.length; j++) {
                                    indel_total[j] += j < ins_total.length ? ins_total[j] : 0;
                                    indel_total[j] += j < del_total.length ? del_total[j] : 0;
                                }
                                // Now we can fit the indel distribution to the lengths of insertions and deletions
                                // ...
                                int ninsertions = 0, ndeletions = 0;
                                for (int j = 0; j < indel_total.length; j++) {
                                    ninsertions += (j < ins_total.length ? ins_total[j] : 0);
                                    ndeletions += (j < del_total.length ? del_total[j] : 0);
                                }

                                try (BufferedWriter writer = new BufferedWriter(new FileWriter(out + "_branch_length_real.txt"))) {
                                    writer.write("Branch\tLength\n");
                                    for (BranchInfo bi : branchInfos) {
                                        String branch = bi.from + "→" + bi.to;
                                        writer.write(branch + "\t" + String.format("%.5f", bi.dist) + "\n");
                                    }
                                    System.out.println("Write completed: " + out + "_branch_length_real.txt");
                                } catch (IOException e) {
                                    System.err.println("Write failed: " + e.getMessage());
                                }


                            /**
                            if (VERBOSE) {
                                System.out.println(rList);
                                try (BufferedWriter writer = new BufferedWriter(new FileWriter(out +"_sample_rlist.txt"))) {
                                    for (Double r : rList) {
                                        writer.write(r.toString());
                                        writer.newLine();
                                    }
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                            }
                             **/

                                    System.out.println("=== Branch → Indel Rate Path (Root to Tips) ===");
                                    for (Map.Entry<String, Double> entry : branchRateMap.entrySet()) {
                                        System.out.printf("%s : %.5f\n", entry.getKey(), entry.getValue());
                                    }

                                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(out + "_branch_rates.txt"))) {
                                        writer.write("Branch\tRate\n");
                                        for (Map.Entry<String, Double> entry : branchRateMap.entrySet()) {
                                            writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                                        }
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }

                                int sum = 0;
                                int count = 0;
                                if (VERBOSE) {
                                    System.out.println("Indels\tInsertions\tDeletions");
                                    System.out.println("Len\tCnt\tCnt\tCnt");
                                }
                                for (int j = 0; j < indel_total.length; j++) {
                                    if (VERBOSE) {
                                        System.out.println((j + 1) + "\t" + indel_total[j] + "\t" + (j < ins_total.length ? ins_total[j] : 0) + "\t" + (j < del_total.length ? del_total[j] : 0));
                                    }
                                    sum += indel_total[j]*(j + 1);
                                    count += indel_total[j];
                                }
                                System.out.println("sample: sum of indel "+ sum + " num of indel " + count);
                                double delprop = (double) ndeletions / (double) (ninsertions + ndeletions);

                                double[] rarray = rList.stream().mapToDouble(Double::doubleValue).toArray();
                                ZeroInflatedGammaMix zig = ZeroInflatedGammaMix.fit(rarray, 3);

                                double rhoP = zig.getP();
                                List<Double> rhoShapes = new ArrayList<>();
                                List<Double> rhoScales = new ArrayList<>();
                                List<Double> rhoWeights = new ArrayList<>();

                                for (Mixture.Component c : zig.getGammaMixture().components) {
                                    GammaDistribution gamma = (GammaDistribution) c.distribution;
                                    String[] parts = gamma.toString().split(",");
                                    double scale = Double.parseDouble(parts[0].replace("Gamma Distribution(", "").trim());
                                    double shape = Double.parseDouble(parts[1].replace(")", "").trim());
                                    double weight = c.priori;

                                    rhoShapes.add(shape);
                                    rhoScales.add(scale);
                                    rhoWeights.add(weight);
                                }


                                // Repeat simulation 5 times
                                for (int k = 1; k <= 5; k++) {
                                    // For each branch, sample a new indel rate from the ZIG distribution
                                    for (BranchInfo bi : branchInfos) {
                                        double sampled_r = zig.sample();
                                        bi.rate = sampled_r;
                                    }

                                    // Optimize the branch order locally to match ground truth using Spearman correlation
                                    //localSwapOptimize(branchInfos, groundTruthRates, 1000, k+10);
                                    globalSwapOptimize(branchInfos, groundTruthRates, 1000, k+10);
                                    // Print the simulated tree's indel rates to the console
                                    System.out.println("=== Simulated Tree #" + k + " ===");
                                    for (BranchInfo bi : branchInfos) {
                                        System.out.printf("%s→%s : %.5f\n", bi.from, bi.to, bi.rate);
                                    }

                                    // Write the simulated indel rates to a file
                                    String filepath = out + "_indel_rates_" + k + ".txt";
                                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(filepath))) {
                                        writer.write("Branch\tIndelRate\n");
                                        for (BranchInfo bi : branchInfos) {
                                            writer.write(bi.from + "→" + bi.to + "\t" + String.format("%.5f", bi.rate) + "\n");
                                        }
                                        System.out.println("Simulated indel rates #" + k + " written to: " + filepath);
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                }

                                try (BufferedWriter writer = new BufferedWriter(new FileWriter(out + "_branch_cumulative_lengths.txt"))) {
                                    writer.write("Branch\tCumulativeLength\n");
                                    for (BranchInfo bi : branchInfos) {
                                        writer.write(bi.from + "→" + bi.to + "\t" + String.format("%.5f", bi.cumulativeLength) + "\n");
                                    }
                                    System.out.println("Cumulative branch lengths written to file.");
                                } catch (IOException e) {
                                    e.printStackTrace();
                                }
                                double[] rInflated = Arrays.stream(rarray).map(x -> Math.min(x * 1.05, 100)).toArray();
                                double[] rDeflated = Arrays.stream(rarray).map(x -> x * 0.95).toArray();


                                double llOriginal = logLikelihood(rarray, zig);
                                double llInflated = logLikelihood(rInflated, zig);
                                double llDeflated = logLikelihood(rDeflated, zig);

                                System.out.printf("Log-likelihood (original): %.4f%n", llOriginal);
                                System.out.printf("Log-likelihood (inflated): %.4f%n", llInflated);
                                System.out.printf("Log-likelihood (deflated): %.4f%n", llDeflated);

                                if (llOriginal >= llInflated && llOriginal >= llDeflated) {
                                    System.out.println("Sanity check passed.");
                                } else {
                                    System.out.println("Warning: model may not reflect original data distribution.");
                                }

                                if (VERBOSE)
                                    System.out.println("Deletion proportion= " + delprop);
                                /* Trying the following distributions (with params)
                                Zipf        1, 2, 5, 15
                                Lavalette   1, 2, 5, 15
                                Poisson     0.01, 0.1, 1, 2
                                 */

                                IndelModel[] inmodels = new IndelModel[] {
                                        new Zipf(1, SEED, ins_total.length),
                                        new Zipf(2, SEED, ins_total.length),
                                        new Zipf(5, SEED, ins_total.length),
                                        new Zipf(15, SEED, ins_total.length),
                                        new Lavalette(1, SEED, ins_total.length),
                                        new Lavalette(2, SEED, ins_total.length),
                                        new Lavalette(5, SEED, ins_total.length),
                                        new Lavalette(15, SEED, ins_total.length),
                                        new ZeroTruncatedPoisson(0.01, SEED),
                                        new ZeroTruncatedPoisson(0.1, SEED),
                                        new ZeroTruncatedPoisson(1.0, SEED),
                                        new ZeroTruncatedPoisson(2.0, SEED)
                                };
                                IndelModel[] delmodels = new IndelModel[] {
                                        new Zipf(1, SEED, del_total.length),
                                        new Zipf(2, SEED, del_total.length),
                                        new Zipf(5, SEED, del_total.length),
                                        new Zipf(15, SEED, del_total.length),
                                        new Lavalette(1, SEED, del_total.length),
                                        new Lavalette(2, SEED, del_total.length),
                                        new Lavalette(5, SEED, del_total.length),
                                        new Lavalette(15, SEED, del_total.length),
                                        new ZeroTruncatedPoisson(0.01, SEED),
                                        new ZeroTruncatedPoisson(0.1, SEED),
                                        new ZeroTruncatedPoisson(1.0, SEED),
                                        new ZeroTruncatedPoisson(2.0, SEED)
                                };

                                IndelModel[] bestsofar = new IndelModel[3];
                                double[] bestscore = new double[3];
                                Arrays.fill(bestscore, Double.NEGATIVE_INFINITY);
                                for (IndelModel model : inmodels) {
                                    double sumIns = 0.0;
                                    for (int j = 0; j < ins_total.length; j++) {
                                        double y = Math.log(model.p(j + 1) + 0.0001);
                                        sumIns += y * ins_total[j];
                                    }
                                    if (sumIns > bestscore[1]) {
                                        bestscore[1] = sumIns;
                                        bestsofar[1] = model;
                                    }
                                }

                                for (IndelModel model : delmodels) {
                                    double sumDel = 0.0;
                                    for (int j = 0; j < del_total.length; j++) {
                                        double y = Math.log(model.p(j + 1) + 0.0001);
                                        sumDel += y * del_total[j];
                                    }
                                    if (sumDel > bestscore[2]) {
                                        bestscore[2] = sumDel;
                                        bestsofar[2] = model;
                                    }
                                }


                                //if (VERBOSE) {
                                    //System.out.println("Best indelmodel: " + bestsofar[0].toString() + " at p= " + bestscore[0]);
                                    //System.out.println("Best inmodel: " + bestsofar[1].toString() + " at p= " + bestscore[1]);
                                    //System.out.println("Best delmodel: " + bestsofar[2].toString() + " at p= " + bestscore[2]);
                                    //System.out.println("Longest indel= " + indel_total.length);
                                //}
                                //
                                double[] gamma = mytree.getGammaParams();
                                double alpha = gamma[0];
                                double beta = gamma[1];
                                if (VERBOSE)
                                    System.out.println("Tree distance Gamma shape= " + alpha + " scale= " + 1.0/beta);
                                StringBuilder n0 = new StringBuilder();
                                for (int j = 0; j < ancseqs_nogap[0].length; j++)
                                    n0.append(ancseqs_nogap[0][j]);
                                if (VERBOSE)
                                    System.out.println("N0= " + n0.toString());
                                Double rates_alpha = 0.01;
                                if (RATES != null) { // position-specific rates available
                                    rates_alpha = getAlpha(RATES);
                                    if (VERBOSE)
                                        System.out.println("Position-specific rates Gamma shape= " + rates_alpha + " scale= " + 1.0/rates_alpha);
                                }
                                System.out.println("To reproduce pseudo-biological properties of current reconstruction, use TrAVIS with the following parameters:");
                                //System.out.println("-d " + alpha + " " + 1.0/beta + " -indelSize " + bestsofar[0].getTrAVIS() + " -maxindel " + indel_total.length + " -n0 " + n0.toString() + " -rates " + rates_alpha + " -gap -extants " + aln.getHeight() + " -delprop " + delprop + " -indelmodel " + rhoP + " " + rhoShape + " " + rhoScale);
                                //System.out.println("Or:");

                                StringBuilder indelModelParams = new StringBuilder();
                                indelModelParams.append(rhoP);
                                for (int j = 0; j < rhoShapes.size(); j++) {
                                    indelModelParams.append(" ").append(rhoWeights.get(j));
                                    indelModelParams.append(" ").append(rhoShapes.get(j));
                                    indelModelParams.append(" ").append(rhoScales.get(j));
                                }

                                System.out.println(
                                        "-d " + alpha + " " + beta +
                                                " --inssize " + bestsofar[1].getTrAVIS() +
                                                " --delsize " + bestsofar[2].getTrAVIS() +
                                                " --maxdellen " + del_total.length +
                                                " --maxinslen " + ins_total.length +
                                                " -n0 " + n0.toString() +
                                                " -m " + MODELS[MODEL_IDX] +
                                                " " + rates_alpha +
                                                " --gap --extants " + aln.getHeight() +
                                                " --delprop " + delprop +
                                                " --indelmodel " + indelModelParams.toString() +
                                                " --seed " + SEED
                                );
                                //System.out.println("Alternatively, consider specifying:");
                                //System.out.println("Gap opening propertion= " + gap_open_prop + " Gap proportion= " + gap_prop + " Mean gap length= " + gap_length);// 如果有祖先序列且提供了输出路径
                                String[] INDELMODELS = new String[]{"ZeroTruncatedPoisson", "Poisson", "Zipf", "Lavalette"};
                                int spaceIndex = bestsofar[1].getTrAVIS().indexOf(" ");
                                String firstPart = (spaceIndex != -1) ? bestsofar[1].getTrAVIS().substring(0, spaceIndex) : bestsofar[1].getTrAVIS();
                                int IN_MODEL_IDX = 0;
                                int DEL_MODEL_IDX = 0;
                                for (int l = 0; l < INDELMODELS.length; l++) {
                                    if (firstPart.equalsIgnoreCase(INDELMODELS[l]))
                                    { IN_MODEL_IDX = l; } }
                                String numberPart = bestsofar[1].getTrAVIS().substring(spaceIndex + 1).trim();
                                double LAMBDA_OF_INMODEL = Double.parseDouble(numberPart);

                                int spaceIndex2 = bestsofar[2].getTrAVIS().indexOf(" ");
                                String firstPart2 = (spaceIndex != -1) ? bestsofar[2].getTrAVIS().substring(0, spaceIndex2) : bestsofar[2].getTrAVIS();
                                for (int l = 0; l < INDELMODELS.length; l++) {
                                    if (firstPart2.equalsIgnoreCase(INDELMODELS[l]))
                                    { DEL_MODEL_IDX = l; } }
                                String numberPart2 =  bestsofar[2].getTrAVIS().substring(spaceIndex2 + 1).trim();
                                double LAMBDA_OF_DELMODEL = Double.parseDouble(numberPart2);
                                //System.out.println(DEL_MODEL_IDX);


                                TrAVIS.TrackTree tracker = null;
                                EnumSeq[] seqs = null;
                                EnumSeq[] seqs_ex = null;
                                EnumSeq[] aseqs = null;
                                EnumSeq[] aseqs_ex = null;
                                double[] rhoShapeArray = rhoShapes.stream().mapToDouble(Double::doubleValue).toArray();
                                double[] rhoScaleArray = rhoScales.stream().mapToDouble(Double::doubleValue).toArray();
                                double[] rhoWeightArray = rhoWeights.stream().mapToDouble(Double::doubleValue).toArray();
                                while (tracker == null) {
                                    tracker = new TrAVIS.TrackTree(newrtree, EnumSeq.parseProtein(n0.toString()), MODEL, SEED,
                                            rates_alpha,
                                            DEL_MODEL_IDX, IN_MODEL_IDX,
                                            LAMBDA_OF_INMODEL, LAMBDA_OF_DELMODEL,
                                            ins_total.length, del_total.length,
                                            delprop,  rhoP , rhoWeightArray,rhoShapeArray , rhoScaleArray,VERBOSE,out);
                                    seqs = tracker.getSequences();
                                    seqs_ex = tracker.getLeafSequences();

                                    for (EnumSeq seq : seqs) {
                                        if (seq.length() < 1) {
                                            tracker = null;
                                            SEED +=1;
                                            break;
                                        }
                                    }
                                }

                                aseqs = tracker.getAlignment();
                                aseqs_ex = tracker.getLeafAlignments();
                                try {
                                    File file = new File(OUTPUT);
                                    if (!file.exists()) {
                                        file.mkdirs();
                                    }

                                    FastaWriter fw = new FastaWriter(out + "_full.fa");
                                    fw.save(GAPPY ? aseqs : seqs);
                                    fw.close();

                                    FastaWriter fw1 = new FastaWriter(out + ".fa");
                                    fw1.save(GAPPY ? aseqs_ex : seqs_ex);
                                    fw1.close();

                                    newrtree.save(out + ".nwk", "nwk");

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
