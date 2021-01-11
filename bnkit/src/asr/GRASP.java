package asr;

import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.*;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import dat.pog.IdxGraph;
import dat.pog.POGTree;
import dat.pog.POGraph;

import java.io.*;
import java.util.Map;
import java.util.concurrent.TimeUnit;

/**
 * Command line version of GRASP.
 * @author mikael
 * @author ariane
 */
public class GRASP {

    public static String VERSION = "0110.2021";
    public static boolean VERBOSE = false;
    public static boolean TIME = false;
    public static int NTHREADS = 1;

    public enum Inference {
        JOINT,
        MARGINAL
    };

    public static void usage() {
        usage(0, null);
    }
    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        if (msg != null)
            out.println(msg + " (Error " + error + ")");
        out.println("Usage: asr.GRASP \n" +
                "\t[-aln <alignment-file> -nwk <tree-file> -out <output-file-or-directory>]\n" +
                "\t{-model <JTT(default)|Dayhoff|LG|WAG|Yang>}\n" +
                "\t{-thr <n-threads>}\n" +
                "\t{-joint (default) | -marg <branchpoint-id>} \n" +
                "\t{-indel <BEP(default)|BEML|SICP|SICML|PSP|PSML>}\n" +
                "\t{-gap}\n" +
                "\t{-savetree <tree-directory>}\n" +
                "\t{-format <FASTA(default)|CLUSTAL|DISTRIB|DOT|TREE>}\n" +
                "\t{-time}{-verbose}{-help}");
        out.println("where \n" +
                "\talignment-file is a multiple-sequence alignment on FASTA or CLUSTAL format\n" +
                "\ttree-file is a phylogenetic tree on Newick format\n" +
                "\toutput-file will be populated by inferred ancestor or ancestors; directory with files if format is DOT\n" +
                "\tInference is either joint (default) or marginal (marginal requires a branch-point to be nominated)\n" +
                "\t\"-gap\" means that the gap-character is included in the resulting output (default for CLUSTAL format, not used with DISTRIB format)\n" +
                "\t\"-savetree\" re-saves the tree on Newick format with generated ancestor labels\n" +
                "\tThe output file is written on the specified format\n" +
                "\t-verbose will print out information about steps undertaken, and -time the time it took to finish");
        out.println("Notes: \n" +
                "\tGreater number of threads may improve processing time, but implies greater memory requirement (default is 1).\n" +
                "\tEvolutionary models for proteins include Jones-Taylor-Thornton (default), Dayhoff-Schwartz-Orcutt, Le-Gasquel and Whelan-Goldman; \n" +
                "\tDNA models include Jukes-Cantor and Yang (general reversible process model).\n" +
                "\tIndel approaches include Bi-directional Edge Parsimony (default), Bi-directional Edge ML, \n" +
                "\tSimple Indel Code Parsimony, Simple Indel Code ML, Position-specific Parsimony and Position-specific ML.\n" +
                "\t~ This is version " + VERSION + " ~");
        System.exit(error);
    }

    public static void main(String[] args) {

        String ALIGNMENT = null;
        String NEWICK = null;
        String OUTPUT = null;

        String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel MODEL = null;
        // Alphabet is decided by MODEL_IDX
        Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
        // Indel approaches:
        String[] INDELS = new String[] {"BEP", "BEML", "SICP", "SICML", "PSP", "PSML"};
        int INDEL_IDX = 0; // default indel approach is that above indexed 0
        boolean GAPPY = false;
        String[] FORMATS = new String[] {"FASTA", "DISTRIB", "CLUSTAL", "DOT", "TREE"};
        int FORMAT_IDX = 0;
        // To compute consensus path is determined by output format
        boolean[] CONSENSUS = new boolean[] {true, false, true, false, false};

        Inference MODE = Inference.JOINT;
        Integer MARG_NODE = null;
        String SAVE_TREE = null;

        long START_TIME, ELAPSED_TIME;

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("aln") && args.length > a + 1) {
                    ALIGNMENT = args[++ a];
                } else if (arg.equalsIgnoreCase("nwk") && args.length > a + 1) {
                    NEWICK = args[++ a];
                } else if (arg.equalsIgnoreCase("out") && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if (arg.equalsIgnoreCase("joint")) {
                    MODE = Inference.JOINT;
                } else if (arg.equalsIgnoreCase("gap")) {
                    GAPPY = true;
                } else if (arg.equalsIgnoreCase("verbose")) {
                    VERBOSE = true;
                } else if (arg.equalsIgnoreCase("time")) {
                    TIME = true;
                } else if (arg.equalsIgnoreCase("savetree") && args.length > a + 1) {
                    SAVE_TREE = args[++a];
                } else if (arg.equalsIgnoreCase("marg") && args.length > a + 1) {
                    MODE = Inference.MARGINAL;
                    String ancid = args[++a];
                    if (ancid.startsWith("N"))
                        ancid = ancid.substring(1);
                    try {
                        MARG_NODE = Integer.parseInt(ancid);
                    } catch (NumberFormatException e) {
                        usage(2, args[a] + " is not a valid ancestor name (use <number>, or \"N<number>\", where <number> increments from 0 at root)");
                    }
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
                    boolean found_indel = false;
                    for (int i = 0; i < INDELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(INDELS[i])) {
                            INDEL_IDX = i;
                            found_indel = true;
                        }
                    }
                    if (!found_indel)
                        usage(3, args[a + 1] + " is not a valid indel approach");
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
                } else if (arg.equalsIgnoreCase("thr") && args.length > a + 1) {
                    try {
                        NTHREADS = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        System.err.println("Failed to set n-threads: " + args[a] + " is not a valid integer");
                    }
                } else if (arg.equalsIgnoreCase("help")) {
                    usage();
                }
            }
        }

        MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
        if (MODEL == null) {
            usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
        }

        if (FORMATS[FORMAT_IDX].equalsIgnoreCase("CLUSTAL")) // Clustal files can only be "gappy"
            GAPPY = true;

        if (ALIGNMENT != null && NEWICK != null && OUTPUT != null) {
            try {
                EnumSeq.Alignment aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
                Tree tree = Utils.loadTree(NEWICK);
                Utils.checkData(aln, tree);
                if (MODE == Inference.MARGINAL) {
                    if (tree.getIndex(MARG_NODE) < 0)
                        usage(2, MARG_NODE + " is not a valid ancestor number");
                } else if (MODE == Inference.JOINT) {
                    if (FORMAT_IDX == 1)
                        usage(2, FORMATS[FORMAT_IDX] + " is not a valid format for joint reconstruction");
                }
                // if we are past the above, we can assume that the data are good to process
                START_TIME = System.currentTimeMillis();
                POGTree pogtree = new POGTree(aln, tree);
                Prediction indelpred = null;
                switch (INDEL_IDX) {
                    case 0: indelpred = Prediction.PredictByBidirEdgeParsimony(pogtree); break;
                    case 1: indelpred = Prediction.PredictByBidirEdgeMaxLHood(pogtree); break;
                    case 2: indelpred = Prediction.PredictByIndelParsimony(pogtree); break;
                    case 3: indelpred = Prediction.PredictByIndelMaxLhood(pogtree); break;
                    case 4: indelpred = Prediction.PredictByParsimony(pogtree); break;
                    case 5: usage(3, "PSML is not implemented"); break;
                    default: break;
                }
                if (indelpred == null)
                    usage(3, INDELS[INDEL_IDX] + " is not implemented");
                if (MODE == Inference.JOINT)
                    indelpred.getJoint(MODEL);
                else if (MODE == Inference.MARGINAL)
                    indelpred.getMarginal(MARG_NODE, MODEL);
                Map<Object, POGraph> pogs = indelpred.getAncestors(MODE);
                POGraph[] ancestors = new POGraph[pogs.size()];
                int ii = 0;
                for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
                    ancestors[ii ++] = entry.getValue();
                EnumSeq[] ancseqs = new EnumSeq[pogs.size()];
                if (CONSENSUS[FORMAT_IDX]) {
                    ii = 0;
                    for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
                        ancseqs[ii ++] = indelpred.getSequence(entry.getKey(), MODE, GAPPY);
                }
                switch (FORMAT_IDX) {
                    case 0: // FASTA
                        FastaWriter fw = new FastaWriter(OUTPUT);
                        fw.save(ancseqs);
                        fw.close();
                        break;
                    case 2: // CLUSTAL
                        AlnWriter aw = new AlnWriter(OUTPUT);
                        aw.save(ancseqs);
                        aw.close();
                        break;
                    case 3: // DOT
                        IdxGraph.saveToDOT(OUTPUT, ancestors);
                        break;
                    case 4: // TREE
                        indelpred.saveTreeInstances(OUTPUT);
                        break;
                    case 1: // DISTRIB
                        EnumDistrib[] d = indelpred.getMarginal(MARG_NODE, MODEL);
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
                            TSVFile.saveObjects(OUTPUT, m);
                        } else
                            usage(8, "Invalid ancestor node label: " + MARG_NODE);
                        break;
                }
                if (SAVE_TREE != null)
                    Newick.save(tree, SAVE_TREE, Newick.MODE_ANCESTOR);

                ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
                if (VERBOSE || TIME) {
                    System.out.println(String.format("Done in %d min, %d sec", TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                            TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME))));
                }

            } catch (ASRException e) {
                usage(22, "Invalid input for ASR: " + e.getMessage());
            } catch (IOException e) {
                usage(2, "Failed to read or write files: " + e.getMessage());
/*            } catch (InterruptedException e) {
                usage(6, "Process interrupted: " + e.getMessage());
 */
            }

        } else if (OUTPUT == null && NEWICK != null && SAVE_TREE != null) {
        } else if (ALIGNMENT == null)
                usage(3, "Need to specify alignment (Clustal or FASTA file)");
        else if (NEWICK == null)
                usage(4, "Need to specify phylogenetic tree (Newick file)");
        else if (OUTPUT == null)
                usage(5, "Need to specify output file");
    }
}
