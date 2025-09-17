package asr;

import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.AlnWriter;
import dat.file.FastaWriter;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.phylo.BranchPoint;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import dat.pog.*;
import stats.*;

import java.io.*;
import java.util.*;

import static bn.prob.GammaDistrib.getAlpha;

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
                "\t[-n0 | --ancestor <ancestor-seq>]\n" +
                "\t[-n | --nwk <filename>]\n" +
                "\t[-o | --output-folder <foldername>]\n" +
                "\t{--extants <5(default)>}\n" +
                "\t{-s | --substitution-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}\n" +
                "\t{-rf | --rates-file <filename>}\n" +
                "\t{--dist-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>:<model-params>}\n" +
                "\t{--leaf2root-distrib <Gaussian|GDF>:<model-params>}\n" +
                "\t{--subst-rate-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>:<model-params>}\n" +
                "\t{* --indel-prior <LOWGAP,MEDGAP,HIGHGAP>}\n" +
                "\t{--indel-rate-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>:<model-params>}\n" +
                "\t\t{* --insertion-rate-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>:<model-params>}\n" +
                "\t\t{* --deletion-rate-distrib <Gamma|ZeroInflatedGamma|ZIG|MixtureGamma>:<model-params>}\n" +
                "\t{--indel-length-distrib <ZeroTruncatedPoisson|ZTP|Poisson|Zipf|Lavalette>:<model-param>}\n"+
                "\t\t{--insertion-length-distrib <ZeroTruncatedPoisson|ZTP|Poisson|Zipf|Lavalette>:<model-param>}\n"+
                "\t\t{--deletion-length-distrib <ZeroTruncatedPoisson|ZTP|Poisson|Zipf|Lavalette>:<model-param>}\n"+
                "\t{--ratesconflate}\n" +
                "\t{--delprop <0.5(default)>}\n" +
                "\t{--max-deletion <max-length>}\n"+
                "\t{--max-insertion <max-length>}\n"+
                "\t{--gap}\n" +
                "\t{--format <FASTA(default)|CLUSTAL|DOT|TREE|RATES|DIR>}\n" +
                "\t{--seed <number>}\n" +
                "\t{--verbose}\n" +
                "\t{--help}\n");
        out.println("where \n" +
                "\t\"-n\" or \"--nwk\" requires a phylogenetic tree specified with the Newick format\n" +
                "\t\"-o\" or \"--output-folder\" requires the name of a directory, for storing results (if it does not exist, it will be created)\n" +
                "\t\"--rates-file\" requires a tabulated file with relative, position-specific rates\n\t\tAs an example, IQ-TREE produces rates on the accepted format (use the --rate option, --mlrate is not yet supported)\n" +
                "\t\"--ratesconflate\" uses the site-specific substitution rate r to modulate the branchpoint-specific indel rate rho, so that indels are introduced with probability exp(-r*rho*t)\n" +
                "\t\"--gap\" means that the gap-character is included in the resulting output (default for CLUSTAL format)\n" +
                "\t\"--indel-prior\" specifies Gamma priors \"LOWGAP\", \"MEDGAP\", \"HIGHGAP\", pre-determined from Pfam alignments that have few, moderate, or large numbers of gaps.\n" +
                "\t\"--indel-length-distrib\" specifies the indel length distribution function, which serves to model both deletions and insertions\n" +
                "\t\"--deletion-length-distrib\" specifies the deletion length distribution, which overrides that specified with \"--indel-distrib\"\n" +
                "\t\"--insertion-length-distrib\" specifies the insertion length distribution, which overrides that specified with \"--indel-distrib\"\n" +
                "\t\"--distance-distrib\" specifies the distribution from which branch distances are drawn\n" +
                "\t\"--leaf2root-distrib\" specifies the distribution that guides how distances are placed relative leaves and the root\n" +
                "\t\"--indel-rate-distrib\" the indel rate distribution, which serves to model both insertions and deletions\n"+
                "\t\"--insertion-rate-distrib\" the model for insertion rate, which overrides that specified with \"--indel-rate-distrib\"\n"+
                "\t\"--deletion-rate-distrib\" the model for deletion rate, which overrides that specified with \"--indel-rate-distrib\"\n"+

                "\t\"--max-deletion\" and \"--max-insertion\" specify the maximum length of deletions and insertions, resp.\n"+
 //               "\t\"-dist\" it can scale the distances of the tree only works when tree is not provide\n" +
                "\t\"--verbose\" means that details of evolutionary events are printed out on standard-output and generate a txt file)\n" +
                "\t\"--help\" prints out the parameters and instructions of how to use this tool\n");
        out.println("Notes: \n" +
                "\t* (asterisk) next to option above means it has not yet been implemented fully.\n" +
                "\tEvolutionary, substitution models for proteins include Jones-Taylor-Thornton (default), Dayhoff-Schwartz-Orcutt, \n\tLe-Gasquel and Whelan-Goldman; \n" +
                "\tModels for DNA include Jukes-Cantor and Yang (general reversible process model).\n" +
                "\tTree is set to have specified extants and gets distances from a distribution, either by specified parameters or estimated from specified tree.\n" +
                "\tSubstitution rates are set from the Gamma distribution, either by specified parameters or estimated from specified rates file.\n" +
                "\tIndel rates are set from the ZeroInflatedGamma or Gamma distribution, either by specified parameters or by fitting ZIG from specified rates file.\n" +
                "\tThe proportion of deletions relative to all indels is given by --delprop <proportion> (when less than 0.5, insertions are more frequent than deletions.\n" +
                "\tThe parameter for indel/deletion/insertion length distribution models for proteins are provided as an argument (see below).\n" +
                "\tGamma distribution is specified <shape,scale>\n" +
                "\t\tshape (aka K, and alpha in the \"rate\" formulation)\n\t\tscale (which is 1/lambda, or 1/beta in the \"rate\" formulation); " +
                "\t\tAlso note substitution rates from a Gamma distribution with specified parameter \"shape\" and mean 1;\n\t if no rate distribution is specified, a uniform rate is used.\n" +
                "\tZeroInflatedGamma (ZIG) distribution is specified <pi,shape,scale> (where pi is the point mass at zero)\n" +
                "\tZeroTruncatedPoisson (ZTP) distribution is specified <lambda>\n" +
                "\tPoisson distribution is specified <lambda>\n" +
                "\tZipf distribution is specified <s>\n" +
                "\tLavalette distribution is specified <a>\n" +
                "\t~ This is part of GRASP-Suite version " + GRASP.VERSION + " ~");
        System.exit(error);
    }

    public static Boolean VERBOSE = false;
    public static String OUTPUT = null;

    public static void main(String[] args) {
        String ANCSEQ = null; // ancestor sequence as a text string, provided
        String OUTPUTTREE = null;
        Double RATESGAMMA = null;
        Double SCALEDIST = null;
        boolean LOADTREE = false;
        String SRATESFILE = null;
        double[] SRATES = null;
        long SEED = System.currentTimeMillis();
        int EXTANTS_N = 5;
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        int DESCENDANTS_MAX = 2, DESCENDANTS_MIN = 2; // Max and min of tree branching
        int DEL_MODEL_IDX = 0, IN_MODEL_IDX = 0;
        double DELETIONPROP = 0.5; // proportion of DELETIONS v INSERTIONS
        double LAMBDA_OF_INMODEL = 1.0;
        double LAMBDA_OF_DELMODEL = 1.0;
        int MAX_IN_LENGTH = 10;
        int MAX_DE_LENGTH = 10;

        String[] EVOL_MODELS = new String[]{"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
        int EVOL_MODEL_IDX = 0; // default model is that above indexed
        SubstModel EVOL_MODEL = null;
        IndelModel INDEL_LENGTH_MODEL = null;
        IndelModel INSERTION_LENGTH_MODEL = null;
        IndelModel DELETION_LENGTH_MODEL = null;
        RateModel SUBST_RATE_MODEL = null;
        RateModel INDEL_RATE_MODEL = null;
        RateModel INSERTION_RATE_MODEL = null;
        RateModel DELETION_RATE_MODEL = null;
        boolean SUBST_RATE_INFLUENCE_INDELS = false;
        RateModel TREE_DISTANCE_MODEL = null;
        Distrib LEAF2ROOT_DISTANCE_MODEL = null;
        // Alphabet is decided by EVOL_MODEL_IDX
        Enumerable[] ALPHAS = new Enumerable[]{Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
        // Indel approaches:
        boolean GAPPY = false;
        String[] FORMATS = new String[]{"FASTA", "DISTRIB", "CLUSTAL", "DOT", "TREE", "DIR", "RATES"};
        int FORMAT_IDX = 0;

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
                } else if (arg.equalsIgnoreCase("-seed") && args.length > a + 1) {
                    SEED = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-extants") && args.length > a + 1) {
                    EXTANTS_N = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("d") && args.length > a + 2) {
                    GAMMA_SHAPE = Double.parseDouble(args[a+1]);
                    GAMMA_SCALE = Double.parseDouble(args[a+2]);
                    if (args.length > a + 3 && !args[a + 3].startsWith("-")) {
                        SCALEDIST = Double.parseDouble(args[a+3]);// brdist_scale
                    }
                } else if (arg.equalsIgnoreCase("-gap")) {
                    GAPPY = true;
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
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    SUBST_RATE_MODEL = RateModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-indel-rate-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    INDEL_RATE_MODEL = RateModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-insertion-rate-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    INSERTION_RATE_MODEL = RateModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-deletion-rate-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    DELETION_RATE_MODEL = RateModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-indel-length-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    INDEL_LENGTH_MODEL = IndelModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-insertion-length-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    INSERTION_LENGTH_MODEL = IndelModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-deletion-length-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    DELETION_LENGTH_MODEL = IndelModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-ratesconflate")) {
                    SUBST_RATE_INFLUENCE_INDELS = true;
                } else if (arg.equalsIgnoreCase("-delprop") && args.length > a + 1) {
                    DELETIONPROP = Double.parseDouble(args[++a]);
                } else if (arg.equalsIgnoreCase("-maxdellen") && args.length > a + 1) {
                    MAX_IN_LENGTH = Integer.parseInt(args[++a]);
                } else if (arg.equalsIgnoreCase("-maxinslen") && args.length > a + 1) {
                    MAX_DE_LENGTH = Integer.parseInt(args[++a]);


        /*
        --tree-distrib MixtureGamma:5.173,0.010,0.236,13.036,0.011,0.349,1.120,0.195,0.416
        --leaf2root-distrib GDF:2.171,1.700
         */

                } else if (arg.equalsIgnoreCase("-tree-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    TREE_DISTANCE_MODEL = RateModel.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-leaf2root-distrib") && args.length > a + 1) {
                    String params = args[a+1];
                    int colonPos = params.indexOf(':');
                    String part1 = colonPos >= 0 ? params.substring(0, colonPos) : params;
                    String part2 = colonPos >= 0 ? params.substring(colonPos + 1) : "";
                    LEAF2ROOT_DISTANCE_MODEL = Distrib.create(part1, part2);
                } else if (arg.equalsIgnoreCase("-format") && args.length > a + 1) {
                    boolean found_format = false;
                    for (int i = 0; i < FORMATS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(FORMATS[i])) {
                            FORMAT_IDX = i;
                            found_format = true;
                        }
                    }
                    if (!found_format)
                        usage(1, args[a + 1] + " is not a valid format name");
                } else if (arg.equalsIgnoreCase("-help") && arg.equalsIgnoreCase("h")) {
                    usage();
                }
            }
        }
        if (SRATESFILE != null) {
            try {
                TSVFile ratesfile = new TSVFile(SRATESFILE, true);
                int rates_col = ratesfile.getColumn("Rate");
                int index_col = ratesfile.getColumn("Site");
                if (rates_col == -1)  // not there
                    rates_col = 0;
                Object[] rateobjs = ratesfile.getCol(rates_col);
                Object[] idxobjs = null;
                if (index_col != -1)
                    idxobjs = ratesfile.getCol(index_col);
                SRATES = new double[rateobjs.length];
                for (int i = 0; i < SRATES.length; i++) {
                    try {
                        int index = index_col == -1 ? i : (Integer) idxobjs[i] - 1; // starts with 1, so subtract "1" to use as position index
                        SRATES[index] = (Double) rateobjs[i];
                    } catch (NumberFormatException e0) {
                        usage(23, "Rates file has invalid number format:" + rateobjs[i]);
                    }
                }
            } catch (IOException e) {
                usage(24, "Rates file could not be opened or read: " + SRATESFILE);
            }
        }

        if (SRATES != null) { // position-specific rates available
            RATESGAMMA = getAlpha(SRATES);
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

        EVOL_MODEL = SubstModel.createModel(EVOL_MODELS[EVOL_MODEL_IDX]);
        if (EVOL_MODEL == null) {
            usage(1, "Model " + EVOL_MODELS[EVOL_MODEL_IDX] + " could not be created");
        }

        EnumSeq ancseq = null;
        if (ANCSEQ != null) {
            if (EVOL_MODEL.getDomain().equals(Enumerable.aacid)) {
                ancseq = EnumSeq.parseProtein(ANCSEQ);
            } else if (EVOL_MODEL.getDomain().equals(Enumerable.nacid)) {
                ancseq = EnumSeq.parseDNA(ANCSEQ);
            } else if (EVOL_MODEL.getDomain().equals(Enumerable.nacidRNA)) {
                ancseq = EnumSeq.parseRNA(ANCSEQ);
            } else {
                usage(5, "Model \"" + EVOL_MODELS[EVOL_MODEL_IDX] + "\" alphabet is not valid");
            }
            if (ancseq == null)
                usage(4, "Invalid ancestor sequence \"" + ANCSEQ + "\" for model " + EVOL_MODELS[EVOL_MODEL_IDX]);
            ancseq.setName("N0");
        }

        if (FORMATS[FORMAT_IDX].equalsIgnoreCase("CLUSTAL")) // Clustal files can only be "gappy"
            GAPPY = true;

        /* Load or create a tree */
        IdxTree tree = null;
        if (LOADTREE && OUTPUTTREE != null) {
            try {
                tree = Newick.load(OUTPUTTREE);
            } catch (IOException e) {
            }
        }
        if (tree == null) {
            if (TREE_DISTANCE_MODEL != null)
                tree = IdxTree.generateTreeFromDistrib(TREE_DISTANCE_MODEL, LEAF2ROOT_DISTANCE_MODEL, EXTANTS_N, SEED, 100);
            else {
                tree = Tree.Random(EXTANTS_N, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, DESCENDANTS_MAX, DESCENDANTS_MIN);
                if (LEAF2ROOT_DISTANCE_MODEL != null)
                    tree = IdxTree.shuffleWithLeaf2RootDistrib(tree, LEAF2ROOT_DISTANCE_MODEL, SEED, 100);
            }
            if (SCALEDIST != null)
                tree.adjustDistances(SCALEDIST);
            if (OUTPUTTREE != null) {
                try {
                    Newick.save(tree, OUTPUTTREE, Newick.MODE_DEFAULT);
                } catch (IOException e) {
                    usage(2, "Tree file could not be saved");
                }
            }
        } else {
            tree = Tree.generateTreeFromMixture(tree,3,SEED,100);
        }
        if (ancseq != null && OUTPUT != null) { // we've got an ancestor to track down the tree
            TrackTree.Params params = new TrackTree.Params(tree, ancseq, EVOL_MODEL);
            if (INDEL_LENGTH_MODEL != null)
                params.setIndelModel(INDEL_LENGTH_MODEL);
            if (INSERTION_LENGTH_MODEL != null)
                params.setInsertmodel(INSERTION_LENGTH_MODEL);
            if (DELETION_LENGTH_MODEL != null)
                params.setDeletemodel(DELETION_LENGTH_MODEL);
            // still need to set rate models
            params.setSubstRateModel(RATESGAMMA);
            if (INDEL_RATE_MODEL != null)
                params.setIndelRateModel(INDEL_RATE_MODEL);

            TrackTree tracker = new TrackTree(params);

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
                        Newick.save(tree, OUTPUT+"/travis.nwk", Newick.MODE_DEFAULT);
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
     * Calculate the relative indel rate (r) for each node in a tree.
     *
     * @param Events    indel, match mismatch
     * @param branchLength Branch lengths from parent to child.
     * @return indel rates for a node.
     */
    public static Double calculateRForNodes(int[] Events, double seqlength, double branchLength) {


        double rs = 0;
        int indels = Events[0];
        int matches = Events[1];    // MB: why do we need the number of matches?
        int mismatches = Events[2]; // MB: why do we need the number of mismatches?

        // Total alignment events (B)
        double B = matches + mismatches + indels;
        // Calculate r using the formula
        if (B > 0 && branchLength > 0) {
            rs = -Math.log(1 - (double) indels / seqlength) / branchLength;
        }

        //if (rs ==0) {
            //rs =1e-4;
            //};

        return rs;
    }

    /** Method to parse a string like "[1,2,3]" into a double array
     * @deprecated ?
     * @param input
     * @return
     */
    private static double[] parseArray(String input) {
        // Remove brackets and split by comma
        String cleaned = input.replaceAll("[\\[\\]]", "");
        String[] parts = cleaned.split(",");

        // Convert string parts to double values
        return Arrays.stream(parts).mapToDouble(Double::parseDouble).toArray();
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

        static class Params {
            public final IdxTree tree;              // the tree that is being used for generating sequences
            public final EnumSeq ancseq;            // the ultimate ancestor sequence
            public final SubstModel substmodel;     // evolutionary model for substitution
            public double[] ancrates = null;        // rates for ancestor to override randomly set rates for the corresponding position

            public RateModel substratemodel = null; // distribution which specifies site-specific rates that modulates substitution at each site/position
            public RateModel indelratemodel = null; // distribution which specifies node-specific rates that modulates indel events in a sequence (NOT site specific)

            public IndelModel insertmodel = null;   // distribution from which insertion lengths are sampled
            public IndelModel deletemodel = null;   // distribution from which deletion lengths are sampled

            // if true, the site-specific evolutionary rate has an impact on the number of indels at a site (mean of rates is 1 across all sites)
            public boolean SUBST_RATE_INFLUENCES_INDELS = false;

            public double PROPORTION_DELETION = 0.5;   // the proportion of deletion events (as opposed to insertion) amongst all indels
            public Random rand = new Random(System.currentTimeMillis());

            public Params(IdxTree tree, EnumSeq ancseq, SubstModel substmodel) {
                this.tree = tree;
                this.ancseq = ancseq;
                this.substmodel = substmodel;
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
                if (insertmodel != null) insertmodel.setSeed(SEED);
                if (deletemodel != null) deletemodel.setSeed(SEED);
                if (substratemodel != null) substratemodel.setSeed(SEED);
                if (indelratemodel != null) indelratemodel.setSeed(SEED);
                this.rand = new Random(SEED);
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
        private double[][] rates = null; // rates in tree, reference to branchpoint specific index
        public boolean USERATES;

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
            EnumSeq ancseq = params.ancseq;

            List<Double> rList = new ArrayList<>();
            USERATES = (params.substratemodel != null); // check if we will generate position specific rates using the model; if not, use a constant rate
            if (USERATES)
                params.substratemodel.setSeed(SEED);
            myType = ancseq.getType();
            IdxTree tree = params.tree;
            int[][] deletions  = new int[tree.getSize()][];
            int[][] insertions = new int[tree.getSize()][];
            rates = new double[tree.getSize()][];
            EnumSeq[] bpseqs = new EnumSeq[tree.getSize()];
            int length_sum = 0;
            int indel_cnt = 0;
            Random rand = params.getRandom();

            for (int idx : tree) {
                if (idx == 0) {
                    bpseqs[0] = ancseq;
                    rates[0] = new double[ancseq.length()];
                    for (int i = 0; i < ancseq.length(); i++)
                        rates[0][i] = params.ancrates != null ? params.ancrates[i] : (USERATES ? params.substratemodel.sample() : 1);
                } else { // branchpoint has parents, all of which have been instantiated (iterator order ensures this, starting with branchpoint idx 0)
                    int paridx = tree.getParent(idx);           // idx of parent
                    Object[] parseq = bpseqs[paridx].get();     // sequence of parent
                    double t = tree.getDistance(idx);           // distance from parent to child
                    //double t = tree.getDistance(1);
                    insertions[idx] = new int[parseq.length+1]; // insertions at this branchpoint relative to parent indices; note that insertions can happen before or after a sequence
                    deletions[idx] = new int[parseq.length];    // deletions at this branchpoint relative to parent indices
                    rates[idx] = new double[parseq.length];     // rates at this branchpoint relative to parent indices; note that insertion rate for before and after is shared
                    // determine what indels are introduced; note: we don't yet know how many indices are required for child so we use lists before moving to array
                    List<Object> child = new ArrayList<>();     // collect character states for the resulting positions, accommodating insertions and deletions
                    List<Object> tail  = new ArrayList<>();     // collect character states for the tail of the child; intended for tailing insertions
                    List<Double> childrates = new ArrayList<>();// collect character rates for the resulting positions, accommodating insertions and deletions
                    List<Double> tailrates = new ArrayList<>(); // collect character rates for the tail of the child; intended for tailing insertions
                    double rho = params.indelratemodel.sample();// node specific rate of insertions and deletions
                    rList.add(rho);                             // save it so the whole series can be recorded
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
                        // indel rate (rho) is separate from site-specific substitution rate, and is node-specific
                        double p = Math.exp(-(rho*t * (params.SUBST_RATE_INFLUENCES_INDELS ? rates[paridx][i]:1)));
                        // make decision of what happens in child for the current site i
                        if (toss < p) { // 1. no indel (so match) with prob p = e^-rt, so consider substitution (r is evolutionary rate and t is branch distance)
                            EnumDistrib d = params.substmodel.getDistrib(parseq[i], rates[paridx][i]*t); // probability of child states GIVEN parent state
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
                                childrates.add(rates[paridx][i]); // stays the same
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
                                insertions[idx][i == 0 ? (rand.nextBoolean() ? 0 : parseq.length) : i] += indel_length; // insertions can be on top of another
                                for (int j = 0; j < indel_length; j ++) {
                                    Object nchar = null;
                                    double tossagain = rand.nextDouble();
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
                                        if (i == 0 && insertions[idx][parseq.length] > 0)
                                            tailrates.add(USERATES ? params.substratemodel.sample() : 1);
                                        else
                                            childrates.add(USERATES ? params.substratemodel.sample() : 1); // new position means new rate
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
                                    childrates.add(rates[paridx][i]); // stays the same
                                } else
                                    throw new RuntimeException("Sampling invalid distribution");
                                i += 1;
                            }
                        }
                    }
                    // set character states
                    Object[] chseq = new Object[child.size() + tail.size()];
                    rates[idx] = new double[childrates.size() + tailrates.size()];
                    for (int j = 0; j < chseq.length; j ++) {
                        if (j >= child.size())
                            rates[idx][j] = tailrates.get(j - child.size());
                        else
                            rates[idx][j] = childrates.get(j);
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

            if (VERBOSE) {
                System.out.println(rList);
                try (BufferedWriter writer = new BufferedWriter(new FileWriter((OUTPUT!=null? OUTPUT:"") + "_result_rlist.txt"))) {
                    for (Double r : rList) {
                        writer.write(r.toString());
                        writer.newLine();
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                String outputFile = (OUTPUT!=null? OUTPUT:"")  +"_travis_report.txt";

                try (PrintWriter pw = new PrintWriter(new FileWriter(outputFile))) {

                    System.out.println(tree);
                    pw.println(tree);

                    for (int idx : tree) {
                        BranchPoint bp = tree.getBranchPoint(idx);
                        BranchPoint parent = bp.getParent();

                        System.out.println(bp.getLabel() + "\t" + bpseqs[idx]);
                        pw.println(bp.getLabel() + "\t" + bpseqs[idx]);

                        if (idx != 0 && parent != null) {
                            for (int i = 0; i < deletions[idx].length; i++) {
                                if (deletions[idx][i] > 0) {
                                    String line = "\tDELETE " + parent.getLabel() + "->"
                                            + bp.getLabel() + "@" + i + ":" + deletions[idx][i];

                                    System.out.println(line);
                                    pw.println(line);
                                }
                            }
                            for (int i = 0; i < insertions[idx].length; i++) {
                                if (insertions[idx][i] > 0) {
                                    String line = "\tINSERT " + parent.getLabel() + "->"
                                            + bp.getLabel() + "@" + i + ":" + insertions[idx][i];

                                    System.out.println(line);
                                    pw.println(line);
                                }
                            }
                        }
                    }

                    pw.flush();

                } catch (IOException e) {
                    e.printStackTrace();
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
                        alignedrates[pos] = rates[idx][j];
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