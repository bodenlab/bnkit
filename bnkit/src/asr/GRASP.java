package asr;

import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.*;
import dat.phylo.Tree;
import dat.pog.IdxGraph;
import dat.pog.POAGraph;
import dat.pog.POGTree;
import dat.pog.POGraph;
import json.JSONObject;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

/**
 * Command line version of GRASP.
 * @author mikael
 * @author ariane
 * @author gabe
 */
public class GRASP {

    public static String VERSION = "17-Oct-2022";

    public static boolean VERBOSE  = false;
    public static boolean TIME     = false;
    public static int     NTHREADS = 4;
    public static boolean NIBBLE   = true;
    public static boolean INDEL_CONSERVATIVE = true;
    // Mode for BEP
    public static boolean RECODE_NULL = true;

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
                "\t{--nonibble}\n" +
                "\t{--exclude-noedge}\n" +
                "\t{--save-as <list-of-formats>} (select multiple from FASTA CLUSTAL TREE DISTRIB ASR DOT TREES)\n" +
                "\t{--save-all} (saves reconstruction with ALL formats)\n" +
                "\t{--include-extants}\n" +
                "\t{--time}{--verbose}{--help}\n");
        out.println("Inference is a two-stage process:\n" +
                "\t(1) A history of indel events is inferred by either maximum likelihood or maximum parsimony and \n\tmapped onto the tree to determine what positions contain actual sequence content\n" +
                "\t(2) For each ancestral position, the most probable character is assigned to each phylogenetic branch \n\tpoint when performing a joint reconstruction. Alternatively, for each \n\tposition at a nominated branch point, the probability distribution over all possible \n\tcharacters is inferred when performing a marginal reconstruction.\n" +
                "\tFinally, edges are drawn to represent all inferred combinations of indels to form an ancestor POG \n\twith nodes that can form a valid sequence with inferred content; a preferred path\n\tthrough the POG is then inferred, nominating a single, best supported sequence.\n");
        out.println("Mode of character inference:\n" +
                "\t-j (or --joint) activates joint reconstruction (default), \n\t-m (or --marginal) activates marginal reconstruction (requires a branch-point to be nominated)\n");
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
                "\t-rf (or --rates-file) specifies a tabulated file with relative, position-specific rates\n\t\tWe recommend the use of this generally, but specifically for trees with great distances, and with biologically diverse entries\n\t\tAs an example, IQ-TREE produces rates on the accepted format\n" +
                "\t--include-extants means that extants are included in output files (when the format allows)\n" +
                "\t--nogap means that the gap-character is excluded in the resulting output (when the format allows)\n" +
                "\t--nonibble de-activates the removal of indices in partial order graphs that cannot form a path from start to end\n" +
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
                "\tTREES: position-specific trees with ancestor states labelled\n");
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

    public static void main(String[] args) {

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
        int SPATH_IDX = 0; // default indel approach is that above indexed 0
        boolean GAPPY = true;
        // output formats
        boolean SAVE_AS = false;
        boolean INCLUDE_EXTANTS = false;
        String[]  FORMATS    = new String[]  {"FASTA", "DISTRIB", "CLUSTAL", "TREE", "ASR", "DOT", "TREES", "MATLAB", "LATEX"};
        // select these, default for "joint reconstruction"
        boolean[] SAVE_AS_IDX = new boolean[FORMATS.length];
        // select to compute consensus path for these output formats
        boolean[] CONSENSUS = new boolean[]  {true,    false,     true,      false,  false,  false, false,   false,    false  };
        // default inference mode
        Inference MODE = Inference.JOINT;
        // ancestor to reconstruct if inference mode is "marginal"
        Integer MARG_NODE = null;

        long START_TIME, ELAPSED_TIME;

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
                        e0.printStackTrace();
                    }
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        try {
            Prediction indelpred = null;
            if (INPUT != null) {
                try {
                    indelpred = Prediction.load(INPUT + "/" + ASRFILE);
                } catch (ASRRuntimeException e) {
                    usage(7, "Prediction failed to load: " + e.getMessage());
                }
            }
            START_TIME = System.currentTimeMillis();
            if (indelpred == null) {
                EnumSeq.Alignment aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
                Tree tree = Utils.loadTree(NEWICK);
                Utils.checkData(aln, tree);
                // if we are past the above, we can assume that the data are good to process
                POGTree pogtree = new POGTree(aln, tree);
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
            Object[][] ancseqs_gappy = null;
            Object[][] ancseqs_nogap = null;
            String[] ancnames = new String[pogs.size()];
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
                        ii ++;
                    }
                } catch (NumberFormatException exc) {
                    for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                        ancnames[ii] = "N" + entry.getKey().toString();
                        ancseqs_gappy[ii] = indelpred.getSequence(entry.getKey(), MODE, true);
                        ancseqs_nogap[ii ++] = indelpred.getSequence(entry.getKey(), MODE, false);
                    }
                }
            }
            File file = new File(OUTPUT);
            if (file.mkdirs()) { // true if the directory was created, false otherwise
            } else {
                // System.err.println("Directory " + OUTPUT + " already exists");
                // throw new ASRException("Directory " + directory + " already exists");
            }
            for (int i = 0; i < SAVE_AS_IDX.length; i++) {
                if (!SAVE_AS_IDX[i])
                    continue;
                switch (i) { // {"FASTA", "DISTRIB", "CLUSTAL", "TREE", "POGS", "DOT", "TREES", "MATLAB", "LATEX"};
                    case 0: // FASTA
                        FastaWriter fw = null;
                        if (MODE == Inference.MARGINAL) // just one sequence
                            fw = new FastaWriter(new File(OUTPUT, PREFIX + "_N" + MARG_NODE + ".fa"));
                        else
                            fw = new FastaWriter(new File(OUTPUT, PREFIX + "_ancestors.fa"));
                        if (GAPPY)
                            fw.save(ancnames, ancseqs_gappy);
                        else
                            fw.save(ancnames, ancseqs_nogap);
                        fw.close();
                        break;
                    case 1: // DISTRIB
                        if (MODE == Inference.MARGINAL) { // must be true for this format
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
                        AlnWriter aw = null;
                        if (MODE == Inference.MARGINAL) // just one sequence
                            aw = new AlnWriter(new File(OUTPUT, PREFIX + "_N" + MARG_NODE + ".aln"));
                        else
                            aw = new AlnWriter(new File(OUTPUT, PREFIX + "_ancestors.aln"));
                        aw.save(ancnames, ancseqs_gappy);
                        aw.close();
                        break;
                    case 3: // TREE
                        Newick.save(indelpred.getTree(), OUTPUT + "/" + PREFIX + "_ancestors.nwk", Newick.MODE_ANCESTOR);
                        break;
                    case 4: // POGS
                        String filename = OUTPUT + "/" + ASRFILE;
                        indelpred.save(filename);
                        break;
                    case 5: // DOT
                        Map<Object, IdxGraph> saveme2 = new HashMap<>();
/*
                        if (INCLUDE_EXTANTS) {
                            IdxGraph g = new POAGraph(aln);
                            g.setName("Exts");
                            saveme2.put("*", g);
                            for (Object name : aln.getNames()) {
                                g = pogtree.getExtant(name);
                                g.setName(name.toString());
                                saveme2.put(name, g);
                            }
                        }
 */
                        for (int idx = 0; idx < ancestors.length; idx++) {
                            ancestors[idx].setName("N" + idx);
                            saveme2.put("N" + idx, ancestors[idx]);
                        }
                        IdxGraph.saveToDOT(OUTPUT, saveme2);
                        break;
                    case 6: // TREES
                        if (MODE == Inference.JOINT)
                            indelpred.saveTreeInstances(OUTPUT);
                        else
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
                }
                ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
                if (VERBOSE || TIME) {
                    System.out.println(String.format("Done in %d min, %d sec", TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                            TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME))));
                }
            }
        } catch (ASRException e) {
            usage(22, "Invalid input for ASR: " + e.getMessage());
        } catch (IOException e) {
            usage(2, "Failed to read or write files: " + e.getMessage());
/*            } catch (InterruptedException e) {
            usage(6, "Process interrupted: " + e.getMessage());
*/
        }
    }
}
