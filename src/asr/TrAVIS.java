package asr;

import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.AlnWriter;
import dat.file.FastaWriter;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.phylo.BranchPoint;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import dat.pog.*;
import stats.*;

import java.io.*;
import java.util.*;

/**
 * Track Ancestor Via Indels and Substitutions.
 *
 * This application is intended to simulate the evolution of biological sequences FROM a given ancestor, TO extants via
 * substitutions, insertions and deletions. These events are tracked, meaning that intermediate ancestors are recorded and
 * an alignment that includes all extants and ancestors can be created. Preceding the alignment, the program generates a
 * partial order graph that helps consolidate multiple insertions and deletions at and across ancestor branchpoints.
 *
 * First, a phylogenetic tree is generated from user-specified parameters, using the gamma distribution for setting
 * evolutionary distances on each branch, in a bi- or multi-furcating manner, until a user-specified number of extants
 * have been mapped as leaves.
 *
 * Mutation events are determined stochastically following stripped-down principles described in
 * Cartwright R. Problems and Solutions for Estimating Indel Rates and Length Distributions.
 * Mol. Biol. Evol. 26(2):473–480. 2009. https://doi.org/10.1093/molbev/msn275
 *
 * For a given ancestor sequence (of arbitrary length) each position is looked at deciding whether to "match" or introduce
 * "indel" in the descendant, by a probability exp^-rt, where "rt" is the rate times the evolutionary distance on the branch.
 * When "matched", a substitution is introduced with a probability determined by a (specified) evolutionary model;
 * when not matched, "indels" are determined from a Poisson with mean 1 (meaning that some will have length 0, just as many
 * length 1, and fewer with longer lengths).
 *
 * Ultimately, tree and alignment files are saved.
 *
 * @author Mikael Boden
 * @author Chongting Zhao
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
        out.println("Usage: asr.TrAVIS \n" +
                "\t[-n0 <ancestor-seq>]\n" +
                "\t[-nwk <tree-file> -out <output-file-or-dir>]\n" +
//                "\t{-rf <rates-file>}\n" +
                "\t{-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}\n" +
                "\t{-rates <a>}\n" +
                "\t{-seed <random>}\n" +
                "\t{-extants <5(default)>}\n" +
                "\t{-dist <mean-extant-to-root>\n" +
                "\t{-shape <1.1(default)>}\n" +
                "\t{-scale <0.2(default)>}\n" +
                "\t{-indel <1.0(default)>}\n" +
                "\t{-delprop <0.5(default)>}\n" +
                "\t{-indelmodel <ZeroTruncatedPoisson(default)|Poisson|Zipf|Lavalette> <model-param>}\n"+
                "\t{-delmodel <ZeroTruncatedPoisson|Poisson|Zipf|Lavalette> <model-param>}\n"+
                "\t{-inmodel <ZeroTruncatedPoisson|Poisson|Zipf|Lavalette> <model-param>}\n"+
                "\t{-maxindel <max-indel-length>}\n"+
                "\t{-gap}\n" +
                "\t{-format <FASTA(default)|CLUSTAL|DOT|TREE|RATES|DIR>}\n" +
                "\t{-verbose}\n" +
                "\t{-help}\n");
        out.println("where \n" +
                "\ttree-file is a phylogenetic tree on Newick format\n" +
                "\toutput-file-or-dir is the filename or name of directory of results\n" +
//                "\trates-file is a tabulated file with relative, position-specific rates\n\t\tAs an example, IQ-TREE produces rates on the accepted format\n" +
                "\t\"-gap\" means that the gap-character is included in the resulting output (default for CLUSTAL format)\n" +
                "\t\"-indelmodel\" specifies the indel length distribution function, which serves to model both deletions and insertions\n" +
                "\t\"-delmodel\" specifies the deletion length distribution, which overrides that specified with \"-indelmodel\"\n" +
                "\t\"-inmodel\" specifies the insertion length distribution, which overrides that specified with \"-indelmodel\"\n" +
                "\t\"-maxindel\" specifies the maximum length of an indel\n" +
                "\t\"-verbose\" means that details of evolutionary events are printed out on standard-output)\n" +
                "\t\"-help\" prints out the parameters and instructions of how to use this tool\n");
        out.println("Notes: \n" +
                "\tEvolutionary models for proteins include Jones-Taylor-Thornton (default), Dayhoff-Schwartz-Orcutt, \n\tLe-Gasquel and Whelan-Goldman; \n" +
                "\tDNA models include Jukes-Cantor and Yang (general reversible process model).\n" +
                "\tTree is set to have specified extants and gets distances from the Gamma distribution, with parameters:\n\t\tshape (aka a and K)\n\t\tscale (aka b, where Lambda=1/b)\n" +
                "\tmean distance to root is used to scale distances in the tree (noting that greater number of extants \n\tindirectly amplifies the time scope of the tree).\n" +
                "\tPosition-specific evolutionary rates from a Gamma distribution with specified parameter \"a\" and mean 1;\n\t if rate is unspecified, a uniform rate 1 is used.\n" +
                "\tRates for insertions and deletions are scaled by -indel <factor> (factor > 1 reduces, factor < 1 increases chance of indels).\n" +
                "\tThe proportion of deletions relative to insertions and deletions is given by -delprop <proportion> (proportion < 0.5 means that insertions will dominate.\n" +
                "\tThe parameter for indel/deletion/insertion length distribution models for proteins are provided as an argument.\n" +
                "\t~ This is part of GRASP-Suite version " + GRASP.VERSION + " ~");
        System.exit(error);
    }

    public static Boolean VERBOSE = false;
    static Double SCALEINDEL = 1.0;
    static Double DELETIONPROP = 0.5; // proportion of DELETIONS v INSERTIONS
    static Double LAMBDA_OF_INMODEL = 1.0;
    static Double LAMBDA_OF_DELMODEL = 1.0;
    static int DEL_MODEL_IDX = 0, IN_MODEL_IDX = 0;
    static Integer MAX_INDEL_LENGTH = 100;

    public static void main(String[] args) {
        String ANCSEQ = null; // ancestor sequence as a text string, provided
        String OUTPUTTREE = null;
        String OUTPUT = null;
        Double RATESGAMMA = null;
        Double SCALEDIST = null;
        boolean LOADTREE = false;
        long SEED = System.currentTimeMillis();
        int EXTANTS_N = 5;
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        int DESCENDANTS_MAX = 2, DESCENDANTS_MIN = 2; // Max and min of tree branching

        String[] MODELS = new String[]{"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
        int MODEL_IDX = 0; // default model is that above indexed
        String[] INDELMODELS = new String[]{"ZeroTruncatedPoisson", "Poisson", "Zipf", "Lavalette"};
        SubstModel MODEL = null;
        // Alphabet is decided by MODEL_IDX
        Enumerable[] ALPHAS = new Enumerable[]{Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
        // Indel approaches:
        boolean GAPPY = false;
        String[] FORMATS = new String[]{"FASTA", "DISTRIB", "CLUSTAL", "DOT", "TREE", "DIR", "RATES"};
        int FORMAT_IDX = 0;
        String[] INDEL_METHODS = new String[] {"BEP", "BEML", "SICP", "SICML", "PSP", "PSML"};
        int INDEL_IDX = 0; // default indel approach is that above indexed 0

        for (int a = 0; a < args.length; a++) {
            if (!args[a].startsWith("-") && ANCSEQ == null) { // ancestor sequence
                ANCSEQ = args[a];
            } else if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("n0") && args.length > a + 1) {
                    ANCSEQ = args[++a];
                } else if (arg.equalsIgnoreCase("nwk") && args.length > a + 1) {
                    OUTPUTTREE = args[++a];
                } else if (arg.equalsIgnoreCase("out") && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if (arg.equalsIgnoreCase("rates") && args.length > a + 1) {
                    RATESGAMMA = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("dist") && args.length > a + 1) {
                    SCALEDIST = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("seed") && args.length > a + 1) {
                    SEED = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("extants") && args.length > a + 1) {
                    EXTANTS_N = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("shape") && args.length > a + 1) {
                    GAMMA_SHAPE = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("scale") && args.length > a + 1) {
                    GAMMA_SCALE = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("gap")) {
                    GAPPY = true;
                } else if (arg.equalsIgnoreCase("verbose")) {
                    VERBOSE = true;
                } else if (arg.equalsIgnoreCase("help") || arg.equalsIgnoreCase("h")) {
                    usage();
                } else if (arg.equalsIgnoreCase("model") && args.length > a + 1) {
                    boolean found_model = false;
                    for (int i = 0; i < MODELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(MODELS[i])) {
                            MODEL_IDX = i;
                            found_model = true;
                        }
                    }
                    if (!found_model)
                        usage(1, args[a + 1] + " is not a valid model name");
                } else if (arg.equalsIgnoreCase("indel") && args.length > a + 1) {
                    SCALEINDEL = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("delprop") && args.length > a + 1) {
                    DELETIONPROP = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("lambda") && args.length > a + 1) {
                    LAMBDA_OF_DELMODEL = LAMBDA_OF_INMODEL = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("maxindel") && args.length > a + 1) {
                    MAX_INDEL_LENGTH = Integer.parseInt(args[++a]);
                } else if ((arg.equalsIgnoreCase("indelmodel") || arg.equalsIgnoreCase("inmodel") || arg.equalsIgnoreCase("delmodel")) && args.length > a + 2) {
                    boolean found_indelmodel = false;
                    for (int i = 0; i < INDELMODELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(INDELMODELS[i])) {
                            try {
                                if (!arg.equalsIgnoreCase("delmodel")) {
                                    IN_MODEL_IDX = i;
                                    LAMBDA_OF_INMODEL = Double.parseDouble(args[a + 2]);
                                }
                                if (!arg.equalsIgnoreCase("inmodel")) {
                                    DEL_MODEL_IDX = i;
                                    LAMBDA_OF_DELMODEL = Double.parseDouble(args[a + 2]);
                                }
                                found_indelmodel = true;
                            } catch (NumberFormatException e) {
                                usage(1, args[a + 2] + " is not a valid parameter for model " + args[a + 1]);
                            }
                        }
                    }
                    if (!found_indelmodel)
                        usage(1, args[a + 1] + " is not a valid model name");
                    else
                        a += 2;
                } else if (arg.equalsIgnoreCase("format") && args.length > a + 1) {
                    boolean found_format = false;
                    for (int i = 0; i < FORMATS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(FORMATS[i])) {
                            FORMAT_IDX = i;
                            found_format = true;
                        }
                    }
                    if (!found_format)
                        usage(1, args[a + 1] + " is not a valid format name");
                } else if (arg.equalsIgnoreCase("help") && arg.equalsIgnoreCase("h")) {
                    usage();
                }
            }
        }

        Random rand = new Random(SEED);

        if (ANCSEQ == null ) { // check standard input for sequence?
            BufferedReader br = null;
            try {
                br = new BufferedReader(new InputStreamReader(System.in));
                String input = br.readLine();
                while (input != null) {
                    ANCSEQ += input.trim();
                    input = br.readLine();
                }
            } catch (IOException e) {
                System.err.println("Error in standard input");
            } finally {
                if (br != null) {
                    try {
                        br.close();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }

        MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
        if (MODEL == null) {
            usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
        }

        EnumSeq ancseq = null;
        if (ANCSEQ != null) {
            if (MODEL.getDomain().equals(Enumerable.aacid)) {
                ancseq = EnumSeq.parseProtein(ANCSEQ);
            } else if (MODEL.getDomain().equals(Enumerable.nacid)) {
                ancseq = EnumSeq.parseDNA(ANCSEQ);
            } else if (MODEL.getDomain().equals(Enumerable.nacidRNA)) {
                ancseq = EnumSeq.parseRNA(ANCSEQ);
            } else {
                usage(5, "Model \"" + MODELS[MODEL_IDX] + "\" alphabet is not valid");
            }
            if (ancseq == null)
                usage(4, "Invalid ancestor sequence \"" + ANCSEQ + "\" for model " + MODELS[MODEL_IDX]);
            ancseq.setName("N0");
        }

        if (FORMATS[FORMAT_IDX].equalsIgnoreCase("CLUSTAL")) // Clustal files can only be "gappy"
            GAPPY = true;

        /* Load or create a tree */
        Tree tree = null;
        if (LOADTREE && OUTPUTTREE != null) {
            try {
                tree = Newick.load(OUTPUTTREE);
            } catch (IOException e) {
            }
        }
        if (tree == null) {
            tree = Tree.Random(EXTANTS_N, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, DESCENDANTS_MAX, DESCENDANTS_MIN);
            if (SCALEDIST != null)
                tree.adjustDistances(SCALEDIST);
            if (OUTPUTTREE != null) {
                try {
                    tree.save(OUTPUTTREE, "nwk");
                } catch (IOException e) {
                    usage(2, "Tree file could not be saved");
                }
            }
        }
        if (ancseq != null && OUTPUT != null) { // we've got an ancestor to track down the tree
            TrackTree tracker = new TrackTree(tree, ancseq, MODEL, SEED, RATESGAMMA==null?-1:RATESGAMMA);
            EnumSeq[] seqs = tracker.getSequences();
            switch (FORMAT_IDX) {
                case 0: // FASTA
                    try {
                        FastaWriter fw = new FastaWriter(OUTPUT);
                        if (!GAPPY) {
                            fw.save(seqs);
                        } else { // gappy
                            EnumSeq[] aln = tracker.getAlignment();
                            fw.save(aln);
                        }
                        fw.close();
                    } catch (IOException e) {
                        usage(6, "FASTA file could not be saved");
                    }
                    break;
                case 3: // DOT
                    POAGraph poag = tracker.getPOAG();
                    try {
                        poag.saveToDOT(OUTPUT);
                    } catch (IOException e) {
                        usage(6, "DOT file could not be saved");
                    }
                    break;
                case 2: // CLUSTAL
                    EnumSeq[] aln = tracker.getAlignment();
                    try {
                        AlnWriter aw = new AlnWriter(OUTPUT);
                        aw.save(aln);
                        aw.close();
                    } catch (IOException e) {
                        usage(6, "CLUSTAL file could not be saved");
                    }
                    break;
                case 5: // ALL in a DIRECTORY
                    POAGraph poaGraph = tracker.getPOAG();
                    EnumSeq[] aseqs = tracker.getAlignment();
                    try {
                        File file = new File(OUTPUT);
                        if (file.mkdirs()) { // true if the directory was created, false otherwise
                        } else {
                            // System.err.println("Directory " + OUTPUT + " already exists");
                            // throw new ASRException("Directory " + directory + " already exists");
                        }
                        FastaWriter fw = new FastaWriter(OUTPUT+"/travis.fa");
                        if (!GAPPY) {
                            fw.save(seqs);
                        } else { // gappy
                            fw.save(aseqs);
                        }
                        fw.close();
                        AlnWriter aw = new AlnWriter(OUTPUT+"/travis.aln");
                        aw.save(aseqs);
                        aw.close();
                        poaGraph.saveToDOT(OUTPUT+"/travis.dot");
                        poaGraph.saveToMatrix(OUTPUT+"/travis.m");
                        tree.save(OUTPUT+"/travis.nwk", "nwk");
                        if (RATESGAMMA != null) {
                            double[] rates = tracker.getRates();
                            Object[][] data = new Object[rates.length + 1][2];
                            for (int i = 0; i <= rates.length; i++) {
                                if (i == 0) // header
                                    data[0] = new Object[]{"Site", "Rate"};
                                else
                                    data[i] = new Object[]{i, rates[i - 1]};
                            }
                            TSVFile ratesfile = new TSVFile(data, true);
                            ratesfile.save(OUTPUT + "/travis.tsv");
                        }
                    } catch (IOException e) {
                        usage(7, "Something went wrong saving files in directory");
                    }
                    break;
                case 6: // RATES
                    if (RATESGAMMA != null) {
                        double[] rates = tracker.getRates();
                        Object[][] data = new Object[rates.length + 1][2];
                        for (int i = 0; i <= rates.length; i++) {
                            if (i == 0) // header
                                data[0] = new Object[]{"Site", "Rate"};
                            else
                                data[i] = new Object[]{i, rates[i - 1]};
                        }
                        try {
                            TSVFile ratesfile = new TSVFile(data, true);
                            ratesfile.save(OUTPUT);
                        } catch (IOException e) {
                            usage(6, "RATES file could not be saved");
                        }
                    } else {
                        usage(8, "RATES are not used hence cannot be saved");
                    }
            }
        }
    }

    /**
     * Calculates the counts of insertions between two sequences.
     *
     * This method compares two sequences and determines the number of insertions
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
     * This method compares two sequences and determines the number of deletions
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
     * Class to track ancestor sequence to extants via intermediate ancestors.
     * Matches/substitutions are determined by a probability p=exp^-rt where rt is the rate times the evolutionary distance from the ancestor to the descendant.
     * If not a match/substitution, insertions and deletions are equally probable, i.e. (1-p)/2 each.
     * The length of an insertion or deletion is determined by a Poisson with mean (lambda) 1; note that this means that 0.37 of indels are length 0.
     * The implementation is inspired by rules extracted from.
     * Position specific rates can be supplied to the constructor.
     * Cartwright R. Problems and Solutions for Estimating Indel Rates and Length Distributions.
     * Mol. Biol. Evol. 26(2):473–480. 2009. https://doi.org/10.1093/molbev/msn275
     */
    static class TrackTree {

        Enumerable myType = null;
        GammaDistrib gamma = null;

        IndelModel inmodel = null;
        IndelModel delmodel = null;

        public final Tree tree;
        private final TreeInstance ti_seqs;
        private EnumNode[][] enumNodes = null;
        private TreeInstance ti_deletions = null;
        private TreeInstance ti_insertions = null;
        private Random rand = null;
        private POAGraph poag = null;
        private EnumSeq.Gappy[] alignedseqs = null; // alignment extracted from sequences and POAG; requires POAG to have been generated
        private int[][] alignedidxs = null;         // alignment extracted from sequences and POAG, but in the form of indices from original seq idx to alignment idx
        // private int[] order = null; // topological order of nodes in POAG; when set alignment can be extracted
        private double[] alignedrates = null;
        private double[][] rates = null; // rates in tree, reference to branchpoint specific index
        public boolean USERATES;

        public TrackTree(Tree tree, EnumSeq ancseq, SubstModel MODEL, long SEED) {
            this(tree, ancseq, MODEL, SEED, -1);
        }

        public TrackTree(Tree tree, EnumSeq ancseq, SubstModel MODEL, long SEED, double ratesgamma) {
            USERATES = (ratesgamma >= 0); // check if we will generate position specific rates; if not, use a constant rate 1
            rand = new Random(SEED);
            switch (IN_MODEL_IDX) { //"Zipf","PowerLaw","Lavalette"
                case 0 -> inmodel  = new ZeroTruncatedPoisson(LAMBDA_OF_INMODEL, SEED);
                case 1 -> inmodel  = new Poisson(LAMBDA_OF_INMODEL, SEED);
                case 2 -> inmodel  = new Zipf(LAMBDA_OF_INMODEL, SEED, MAX_INDEL_LENGTH);
                case 3 -> inmodel  = new Lavalette(LAMBDA_OF_INMODEL, SEED, MAX_INDEL_LENGTH);
                default -> throw new IllegalArgumentException("Invalid model index");
            }
            switch (DEL_MODEL_IDX) { //"Zipf","PowerLaw","Lavalette"
                case 0 -> delmodel  = new ZeroTruncatedPoisson(LAMBDA_OF_DELMODEL, SEED + 1);
                case 1 -> delmodel  = new Poisson(LAMBDA_OF_DELMODEL, SEED + 1);
                case 2 -> delmodel  = new Zipf(LAMBDA_OF_DELMODEL, SEED + 1, MAX_INDEL_LENGTH);
                case 3 -> delmodel  = new Lavalette(LAMBDA_OF_DELMODEL, SEED + 1, MAX_INDEL_LENGTH);
                default -> throw new IllegalArgumentException("Invalid model index");
            }

            if (USERATES) {
                gamma = new GammaDistrib(ratesgamma, ratesgamma); // mean is a/b so setting b=a
                gamma.setSeed(SEED);
            }
            double fixed_rate = USERATES ? gamma.sample() : 1;
            myType = ancseq.getType();
            this.tree = tree;
            int[][] deletions  = new int[tree.getSize()][];
            int[][] insertions = new int[tree.getSize()][];
            rates = new double[tree.getSize()][];
            EnumSeq[] bpseqs = new EnumSeq[tree.getSize()];

            for (int idx : tree) {
                if (idx == 0) {
                    bpseqs[0] = ancseq;
                    if (USERATES) {
                        rates[0] = new double[ancseq.length()];
                        for (int i = 0; i < ancseq.length(); i++)
                            rates[0][i] = USERATES ? gamma.sample() : 1;
                            //rates[0][i] =  fixed_rate;
                    }
                } else { // branchpoint has parents, all of which have been instantiated (iterator order ensures this, starting with branchpoint idx 0)
                    int paridx = tree.getParent(idx);           // idx of parent
                    Object[] parseq = bpseqs[paridx].get();     // sequence of parent
                    double t = tree.getDistance(idx);           // distance from parent to child
                    //double t = tree.getDistance(1);
                    insertions[idx] = new int[parseq.length+1]; // insertions at this branchpoint relative to parent indices; note that insertions can happen before or after a sequence
                    deletions[idx] = new int[parseq.length];    // deletions at this branchpoint relative to parent indices
                    rates[idx] = new double[parseq.length];     // rates at this branchpoint relative to parent indices; note that insertion rate for before and after is shared
                    // determine what indels are introduced
                    List<Object> child = new ArrayList<>();     // collect character states for the resulting positions, accommodating insertions and deletions
                    List<Object> tail  = new ArrayList<>();     // collect character states for the tail of the child; intended for tailing insertions
                    List<Double> childrates = new ArrayList<>();// collect character rates for the resulting positions, accommodating insertions and deletions
                    List<Double> tailrates = new ArrayList<>(); // collect character rates for the tail of the child; intended for tailing insertions
                    // note: we don't yet know how many indices are required for child so use list before moving to array
                    int i = 0;                                  // idx for parent position
                    double toss = rand.nextDouble() / TrAVIS.SCALEINDEL;
                    while (i < parseq.length) {                 // move through the child by incrementing the idx in the parent
                        // three possibilities: 1. match and potential substitution, 2. deletion/s, and 3. insertion/s
                        double p = Math.exp(-(USERATES?rates[paridx][i]:1)*t);

                        if (toss < p) { // 1. no indel (so match) with prob p = e^-rt, so consider substitution
                            EnumDistrib d = MODEL.getDistrib(parseq[i], (USERATES?rates[paridx][i]:1)*t); // bug fix 13/3/24, prev version did not multiply with site specific rate
                            Object nchar = null;
                            double tossagain = rand.nextDouble();
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
                                if (USERATES)
                                    childrates.add(rates[paridx][i]); // stays the same
                            } else
                                throw new RuntimeException("Sampling invalid distribution");
                            i += 1;
                        } else {
                            double toss2 = rand.nextDouble();

                            if (toss2 < TrAVIS.DELETIONPROP) { // 2. deletion with prob q = (1 - p)/2, so consider length of deletion
                                int k;
                                k = Math.min(delmodel.sample(), parseq.length - i);// length, can only delete what is left of the sequence
                                deletions[idx][i] = k;  // deletions skip characters in the parent
                                i += k;
                            } else { // 3. insertion with prob q = (1 - p)/2, so consider length of insertion
                                // special case at i==0: insertions can happen BEFORE and AFTER the sequence,
                                // so to avoid introducing a bias for LONGER EXTANTS,
                                // we place it at either end with a uniform coin toss
                                int k;
                                k = inmodel.sample(); // length, can only delete what is left of the sequence
                                insertions[idx][i == 0 ? (rand.nextBoolean() ? 0 : parseq.length) : i] += k; // insertions can be on top of another
                                for (int j = 0; j < k; j ++) {
                                    Object nchar = null;
                                    double tossagain = rand.nextDouble();
                                    double sump = 0;
                                    for (Object c : MODEL.getDomain().getValues()) {
                                        sump += MODEL.getProb(c);
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
                                        if (USERATES) {
                                            if (i == 0 && insertions[idx][parseq.length] > 0)
                                                tailrates.add(gamma.sample());
                                            else
                                                childrates.add(gamma.sample()); // new position means new rate
                                        }
                                    } else
                                        throw new RuntimeException("Sampling invalid distribution");
                                }
                                // after an insertion, what do we do with the character? Currently, we enforce a match/substitution
                                // note that for a before/after insertion, the character concerned is always that at the head of the sequence
                                EnumDistrib d = MODEL.getDistrib(parseq[i], t);
                                Object nchar = null;
                                double tossagain = rand.nextDouble();
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
                                    if (USERATES)
                                        childrates.add(rates[paridx][i]); // stays the same
                                } else
                                    throw new RuntimeException("Sampling invalid distribution");
                                i += 1;
                            }
                        }
                    }
                    // set character states
                    Object[] chseq = new Object[child.size() + tail.size()];
                    if (USERATES)
                        rates[idx] = new double[childrates.size() + tailrates.size()];
                    for (int j = 0; j < chseq.length; j ++) {
                        if (USERATES) {
                            if (j >= child.size())
                                rates[idx][j] = tailrates.get(j - child.size());
                            else
                                rates[idx][j] = childrates.get(j);
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
            ti_deletions = new TreeInstance(tree, deletions);
            ti_insertions = new TreeInstance(tree, insertions);
            ti_seqs = new TreeInstance(tree, bpseqs);

            if (TrAVIS.VERBOSE) {
                System.out.println(tree);
                for (int idx : tree) {
                    BranchPoint bp = tree.getBranchPoint(idx);
                    BranchPoint parent = bp.getParent(); // is null if no parent
                    System.out.println(bp.getLabel() + "\t" + bpseqs[idx]);
                    if (idx != 0 && parent != null) {
                        for (int i = 0; i < deletions[idx].length; i++) {
                            if (deletions[idx][i] > 0)
                                System.out.println("\tDELETE " + parent.getLabel() + "->" + bp.getLabel() + "@" + i + ":" + deletions[idx][i]);
                        }
                        for (int i = 0; i < insertions[idx].length; i++) {
                            if (insertions[idx][i] > 0)
                                System.out.println("\tINSERT " + parent.getLabel() + "->" + bp.getLabel() + "@" + i + ":" + insertions[idx][i]);
                        }
                    }
                }
            }
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
            enumNodes = new EnumNode[tree.getSize()][];
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
                    poag.addEdge(-1, nodes.get(edge.getPair()[1]));
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
            alignedseqs = new EnumSeq.Gappy[tree.getSize()];
            alignedidxs = new int[tree.getSize()][];
            HashMap<EnumNode, Integer> nodes = new HashMap<>();
            for (int i = 0; i < order.length; i ++) {
                EnumNode node = (EnumNode) poag.getNode(order[i]);
                nodes.put(node, i);
                if (VERBOSE)
                    System.out.println(i + "\t" + order[i] + "\t" + node.getLabel());
            }
            if (USERATES)
                alignedrates = new double[order.length];
            for (int idx : tree) {
                alignedseqs[idx] = new EnumSeq.Gappy(myType);
                Object[] seq = new Object[order.length];
                EnumSeq orig = (EnumSeq)ti_seqs.getInstance(idx);
                alignedidxs[idx] = new int[orig.length()];
                for (int j = 0; j < enumNodes[idx].length; j ++) {
                    int pos = nodes.get(enumNodes[idx][j]);
                    seq[pos] = orig.get(j);
                    alignedidxs[idx][j] = pos;
                    if (USERATES)
                        alignedrates[pos] = rates[idx][j];
                }
                alignedseqs[idx].set(seq);
                alignedseqs[idx].setName(tree.getBranchPoint(idx).getLabel().toString());
            }
            return alignedseqs;
        }

        public double[] getRates() {
            if (!USERATES)
                return null;
            if (alignedrates != null)
                return alignedrates;
            getAlignment();
            return alignedrates;
        }

        private void build(int idx, EnumNode start, EnumNode[] nodes, EnumNode end, Set<EnumEdge> edges) {
            enumNodes[idx] = nodes;
            EnumNode prev = start;
            for (int i = 0; i < nodes.length; i ++) {
                edges.add(new EnumEdge(prev, nodes[i]));
                prev = nodes[i];
            }
            edges.add(new EnumEdge(prev, end));
            for (int chidx : tree.getChildren(idx)) {
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