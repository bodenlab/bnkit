package asr;

import bn.BNet;
import bn.BNode;
import bn.Distrib;
import bn.alg.EM;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.file.BNBuf;
import bn.prob.EnumDistrib;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.phylo.PhyloBN;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import json.JSONException;
import json.JSONObject;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class BPState {

    public static void usage() {
        usage(0, null);
    }
    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        if (msg != null)
            out.println(msg + " (Error " + error + ")");
        out.println("Usage: asr.BPState \n" +
                "\t[-nwk <tree-file> -in <input-file> -out <output-file>]\n" +
                "\t{-model <uniform(default)>}\n" +
                "\t{-gamma <gamma-value(default 1.0)>}\n" +
                "\t{-params <JSON-string-no-spaces>}\n" +
                "\t{-latent <#states>}\n" +
                "\t{-internal }\n" +
                "\t{-learn}\n" +
                "\t{-seed <seed>} \n" +
                "\t{-joint (default) | -marg <branchpoint-id>} \n" +
                "\t{-format <TSV(default), TREE, STDOUT>}\n");
        out.println("where \n" +
                "\ttree-file is a phylogenetic tree on Newick format\n" +
                "\tinput-file is a table headed with sequence or ancestor names, and corresponding values (empty implies not assigned) on TSV format\n" +
                "\toutput-file will contain:\n" +
                "\t- inferred branch point states on specified format (TSV by default, TREE is a labelled tree on Newick format), or\n" +
                "\tgamma-value is parameter to the uniform model (n-state generalisation of Jukes-Cantor)\n" +
                "\tJSON-string specifies distribution on leaf nodes (if latent mode is used)\n" +
                "\tlatent indicates that the tree consists of latent values only (latent mode), with specified values as extensions to the leaves.\n" +
                    "\t- #states is the number of latent states to learn (should not exceed 25, labelled A-Z).\n" +
                "\tinternal indicates that internal nodes are also extended with user-specified or learned distributions (default leaves-only).\n" +
                "\tlearn excludes inference and instead prompts EM learning of parameters, using input data as training data.\n" +
                "\tInference is either joint (default) or marginal (marginal requires a branch-point to be nominated)\n");
        out.println("Notes:\n" +
                "\tEvolutionary models of substitution are currently limited to uniform, which is an adaptation of Jukes-Cantor for arbitrary number of states.\n" +
                    "\t- gamma-value is used by this model\n" +
                "\tIf specified values are real, a conditional Gaussian mixture distribution conditioned on latent state is learned.\n" +
                "\tIf specified values are discrete, a multinomial distribution conditioned on latent state is learned.\n");
        System.exit(error);
    }

    public enum MODEL_MODE {DIRECT, LATENT};

    public static void main(String[] args) {
        String NEWICK = null;
        String OUTPUT = null;
        String INPUT = null;

        String[] MODELS = new String[] {"uniform"};
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel SUBST_MODEL = null;
        String[] FORMATS = new String[] {"TSV", "TREE", "STDOUT"};
        int FORMAT_IDX = 0;
        asr.GRASP.Inference INFERENCE_MODE = asr.GRASP.Inference.JOINT;
        MODEL_MODE mode = MODEL_MODE.DIRECT; // not latent, discrete labels
        boolean LEARN = false; // train (when true) or infer (when false, default)
        boolean LEAVES_ONLY = true; // leaves-only are equipped with distributions (when true, default), or all nodes including internal (when false)
        Object MARG_LABEL = null;
        Double GAMMA = 1.0;
        Long SEED = 1L;
        Integer NSTATES = null; // number of states if continuous; if discrete, the same as the number of used values in the input file
        JSONObject PARAMS = null; // optional parameters as JSON

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("nwk") && args.length > a + 1) {
                    NEWICK = args[++ a];
                } else if (arg.equalsIgnoreCase("out") && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if (arg.equalsIgnoreCase("in") && args.length > a + 1) {
                    INPUT = args[++a];
                } else if (arg.equalsIgnoreCase("joint")) {
                    INFERENCE_MODE = asr.GRASP.Inference.JOINT;
                } else if (arg.equalsIgnoreCase("marg") && args.length > a + 1) {
                    INFERENCE_MODE = asr.GRASP.Inference.MARGINAL;
                    String NODE_LABEL = args[++a];
                    if (NODE_LABEL.startsWith("N")) {
                        try {
                            MARG_LABEL = Integer.parseInt(NODE_LABEL.substring(1));
                        } catch (NumberFormatException e) { // failed to assign an internal ancestor node number
                            MARG_LABEL = NODE_LABEL;        // so this is probably just an extant
                        }
                    } else {
                        try {
                            MARG_LABEL = Integer.parseInt(NODE_LABEL);
                        } catch (NumberFormatException e) { // failed to assign an internal ancestor node number
                            MARG_LABEL = NODE_LABEL;        // so this is probably just an extant
                        }
                    }
                } else if (arg.equalsIgnoreCase("gamma") && args.length > a + 1) {
                    try {
                        GAMMA = Double.parseDouble(args[a + 1]);
                        if (GAMMA <= 0)
                            throw new NumberFormatException("Cannot be zero or negative");
                    } catch (NumberFormatException e) {
                        usage(3, args[a + 1] + " is not a valid gamma value (must be a non-zero, positive real value)");
                    }
                } else if (arg.equalsIgnoreCase("latent") && args.length > a + 1) {
                    try {
                        NSTATES = Integer.parseInt(args[a + 1]);
                        if (NSTATES <= 1)
                            throw new NumberFormatException("Cannot be zero, one or negative");
                        mode = MODEL_MODE.LATENT;
                    } catch (NumberFormatException e) {
                        usage(7, args[a + 1] + " is not a valid column value (must be a non-zero, positive integer)");
                    }
                } else if (arg.equalsIgnoreCase("seed") && args.length > a + 1) {
                    try {
                        SEED = Long.parseLong(args[a + 1]);
                    } catch (NumberFormatException e) {
                        usage(10, args[a + 1] + " is not a valid seed value (must be a valid long integer)");
                    }
                } else if (arg.equalsIgnoreCase("params") && args.length > a + 1) {
                    try {
                        PARAMS = new JSONObject(args[a + 1]);
                    } catch (JSONException e) {
                        usage(9, "Invalid JSON specification of distribution for external nodes: " + PARAMS);
                    }
                    if (PARAMS == null)
                        usage(9, "Invalid JSON specification of distribution for external nodes: " + PARAMS);
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
                } else if (arg.equalsIgnoreCase("learn")) {
                    LEARN = true;
                } else if (arg.equalsIgnoreCase("internal")) {
                    LEAVES_ONLY = false;
                } else if (arg.equalsIgnoreCase("help")) {
                    usage();
                }
            }
        }

        if (INPUT != null && NEWICK != null) {
            Object[][] inputs = null;
            Tree tree = null;
            try {
                inputs = TSVFile.loadObjects(INPUT);
                tree = Tree.load(NEWICK, "newick");
            } catch (IOException e) {
                usage(5, inputs == null ? "Failed to load the input-file " + INPUT : "Failed to load the tree-file " + NEWICK);
            }
            TSVFile tsv = new TSVFile(inputs, true); // put data in TSV instance, allowing for header
            Set observed;
            String[] xalpha = null;
            if (mode == MODEL_MODE.LATENT) {
                Object[] alpha = new Object[NSTATES];
                for (int i = 0; i < alpha.length; i ++)
                    alpha[i] = Character.valueOf((char)('A' + i));
                SUBST_MODEL = new JC(GAMMA, alpha);
                if (TSVFile.isDouble(tsv.getRow(0))) { // the first row of data is all Double
                    observed = null;
                    System.out.println("Detected real values in data--using Gaussian mixtures with " + NSTATES + " components");
                } else { // the column has discrete values (not real)
                    observed = tsv.getValues();
                    Object[] extalpha = new Object[observed.size()];
                    observed.toArray(extalpha);
                    xalpha = new String[extalpha.length];
                    for (int i = 0; i < extalpha.length; i ++)
                        xalpha[i] = extalpha[i].toString();
                    System.out.println("Detected " + xalpha.length + " discrete values in data--using multi-nomial distributions conditioned on " + NSTATES + " states");
                }
            } else if (mode == MODEL_MODE.DIRECT) {
                if (TSVFile.isDouble(tsv.getRow(0))) { // the column is all Double
                    usage(8, "Real-value data requires (discrete) latent states; please specify number of latent states.");
                } else { // the column has discrete values (not real)
                    observed = tsv.getValues();
                    NSTATES = observed.size();
                    if (MODELS[MODEL_IDX].equals("uniform")) {
                        Object[] alpha = new Object[observed.size()];
                        observed.toArray(alpha);
                        SUBST_MODEL = new JC(GAMMA, alpha);
                    } else { // some other model
                        SUBST_MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
                        for (Object v : observed) {
                            if (!SUBST_MODEL.getDomain().isValid(v))
                                usage(1, "Model " + MODELS[MODEL_IDX] + " does not accept value/state \"" + v.toString() + "\"");
                        }
                    }
                }
            }
            if (SUBST_MODEL == null) {
                usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
            }
            String[] headers = tsv.getHeaders();
            Object[][] rows = tsv.getRows();
            TreeInstance ti = tree.getInstance(headers, rows[0]);
            if (mode == MODEL_MODE.DIRECT && INFERENCE_MODE == GRASP.Inference.JOINT) {
                MaxLhoodJoint inf = new MaxLhoodJoint(tree, SUBST_MODEL);
                inf.decorate(ti);
                Object[][] save = new Object[tree.getSize()][2];
                Object[] vals = new Object[tree.getSize()];
                for (int bpidx : tree) {
                    Object state = inf.getDecoration(bpidx);
                    save[bpidx][0] = tree.getLabel(bpidx);
                    save[bpidx][1] = state;
                    vals[bpidx] = state;
                    if (FORMAT_IDX == 2)
                        System.out.println(bpidx + "\t" + tree.getLabel(bpidx) + "\t" + state);
                }
                if (FORMAT_IDX == 0) {
                    try {
                        TSVFile.saveObjects(OUTPUT, save);
                    } catch (IOException e) {
                        usage(6, "Failed to save output to " + OUTPUT);
                    }
                } else if (FORMAT_IDX == 1) { // TREE
                    try {
                        Newick.save(tree, OUTPUT, vals);
                    } catch (IOException e) {
                        usage(6, "Failed to save output to " + OUTPUT);
                    }
                }
            } else if (mode == MODEL_MODE.DIRECT && INFERENCE_MODE == GRASP.Inference.MARGINAL) {
                int bpidx = tree.getIndex(MARG_LABEL);
                if (bpidx < 0) {
                    usage(5, "Invalid branch point name " + MARG_LABEL);
                }
                MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(bpidx, tree, SUBST_MODEL);
                inf.decorate(ti);
                EnumDistrib distrib = inf.getDecoration(bpidx);
                Object[][] save = new Object[2][1 + distrib.getDomain().size()];
                save[0][0] = "";
                save[1][0] = tree.getLabel(bpidx);
                Object[] vals = distrib.getDomain().getValues();
                for (int i = 0; i < vals.length; i ++)
                    save[0][1 + i] = vals[i];
                for (int i = 0; i < vals.length; i ++)
                    save[1][1 + i] = distrib.get(i);
                if (FORMAT_IDX == 2)
                    System.out.println(bpidx + "\t" + tree.getLabel(bpidx) + "\t" + distrib);
                else if (FORMAT_IDX == 0) {
                    try {
                        TSVFile.saveObjects(OUTPUT, save);
                    } catch (IOException e) {
                        usage(6, "Failed to save output to " + OUTPUT);
                    }
                }
            } else if (mode == MODEL_MODE.LATENT) {
                // System.out.println(tree.toJSON());
                // create a Bayesian network from a phylogenetic tree
                PhyloBN pbn;
                if (xalpha == null)  // specified values are "real", so equip "ext" nodes as GDTs
                    pbn = PhyloBN.withGDTs(tree, SUBST_MODEL, 1, LEAVES_ONLY, SEED);
                else  // specified values are discrete and nominated as strings in "xalpha", so add multi-nomial CPTs
                    pbn = PhyloBN.withCPTs(tree, SUBST_MODEL, xalpha, 1, LEAVES_ONLY, SEED);
                if (PARAMS != null) {
                    pbn.overrideMasterJSON(PARAMS);
                    System.out.println("Using pre-set distribution: " + pbn.toString());
                }
                if (LEARN) { // learn, do not infer
                    // train the BN; headers correspond to labels of leaves or any other node that has been nominated in the TSV file
                    pbn.trainEM(headers, rows, 1L);
                    // training done
                    System.out.println("Learned parameters specified as:");
                    System.out.println("\t-params " + pbn.getMasterJSON().toString());
                } else if (INFERENCE_MODE == GRASP.Inference.MARGINAL) { // inference
                    // retrieve the branchpoint index for the nominated ancestors
                    int bpidx = tree.getIndex(MARG_LABEL);
                    if (bpidx < 0) // did not find it...
                        usage(5, "Invalid branch point name " + MARG_LABEL);
                    // inference below; first create the inference instance
                    MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(bpidx, pbn);
                    // perform marginal inference
                    inf.decorate(ti);
                    // retrieve the distribution at the node previously nominated
                    Distrib anydistrib = inf.getDecoration(bpidx);
                    Object[][] save = null;
                    try {
                        EnumDistrib distrib = (EnumDistrib) anydistrib; // if this cast succeeds, it is a discrete distribution
                        if (FORMAT_IDX == 2)
                            System.out.println(bpidx + "\t" + tree.getLabel(bpidx) + "\t" + distrib);
                        else if (FORMAT_IDX == 0) {
                            save = new Object[2][1 + distrib.getDomain().size()];
                            save[0][0] = "";
                            save[1][0] = tree.getLabel(bpidx);
                            Object[] vals = distrib.getDomain().getValues();
                            for (int i = 0; i < vals.length; i++)
                                save[0][1 + i] = vals[i];
                            for (int i = 0; i < vals.length; i++)
                                save[1][1 + i] = distrib.get(i);
                        }
                    } catch (ClassCastException e) { // GDT probably
                        if (FORMAT_IDX == 2)
                            System.out.println(bpidx + "\t" + tree.getLabel(bpidx) + "\t" +anydistrib);
                        else if (FORMAT_IDX == 0) {
                            save = new Object[1][2];
                            save[0][0] = tree.getLabel(bpidx);
                            save[0][1] = anydistrib;
                        }
                    }
                    try {
                        TSVFile.saveObjects(OUTPUT, save);
                    } catch (IOException e) {
                        usage(6, "Failed to save output to " + OUTPUT);
                    }
                }
            }
        } else {
            usage(4, "Both an input-file and a tree-file are required");
        }
    }

}
