package asr;

import bn.BNode;
import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.node.GDT;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;
import bn.prob.MixtureDistrib;
import dat.Enumerable;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.file.Utils;
import dat.phylo.PhyloBN;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import json.JSONException;
import json.JSONObject;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class TreeGazer {

    public static void usage() {
        usage(0, null);
    }
    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        if (msg != null)
            out.println(msg + " (Error " + error + ")");
        out.println("""
                Usage: asr.TreeGazer\s
                \t[-nwk <tree-file> -in {<label>{:<parser>}@}<input-file> -out <output-file>]
                \t{-model <uniform(default)>}
                \t{-gamma <value (default 1.0)>}
                \t{-params <JSON-file>}
                \t{-latent <#states>}
                \t{-internal}
                \t{-learn}
                \t{-tied}\s
                \t{-seed <seed>}\s
                \t{-joint (default) | -marg {<branchpoint-id>} }\s
                \t{-format <TSV(default), TREE, STDOUT, ITOL>}
                \t{-lambda <value (default 5.0)>}
                \t{-help|-h}
                \t{-verbose|-v}
                """);
        out.println("""
                where\s
                \ttree-file is a phylogenetic tree on Newick format
                \tinput-file is a table with sequence or ancestor names in the first column, and corresponding values
                \t(empty or None or null implies not assigned) on TSV format
                \tlabel flags that a header is used in the input-file and identifies the column with values to be modelled;
                \tif no label is given, headers are assumed absent and values from the second column will be modelled
                \tparser identifies a parser to use on the column with values (e.g. BRENDA).
                \toutput-file will contain:
                \t- inferred branch point states on specified format (TSV by default, TREE is a labelled tree on Newick format, ITOL is a dataset to decorate trees in iTOL.embl.de), or
                \tgamma-value is parameter to the uniform model (n-state generalisation of Jukes-Cantor)
                \tlambda is the multiplier for the upper confidence bound of predicted values (used only when latent mode with real values is applied)
                \tJSON-file contains a JSON string specifying the distribution for latent nodes (if latent mode is used)
                \tlatent indicates that the tree consists of latent values only (latent mode), with specified values as extensions to the leaves.
                \t- #states is the number of latent states to learn (should not exceed 25, labelled A-Z).
                \tinternal indicates that internal nodes are also extended with user-specified or learned distributions (default leaves-only).
                \tlearn excludes inference and instead prompts EM learning of parameters, using input data as training data.
                \ttied implies that the variance learned is the same across the latent states (only applicable when EM-learning GDTs; default is off).
                \thelp prints out commandline arguments (this screen).
                \tverbose completes the requested steps while printing out messages about the process.
                """);
        out.println("""
                Notes:
                \tEvolutionary models of substitution are currently limited to uniform, which is an adaptation of Jukes-Cantor for arbitrary number of states.
                \t- gamma-value is used by this model
                \tIf specified values are real, a conditional Gaussian mixture distribution conditioned on latent state is learned.
                \tIf specified values are discrete, a multinomial distribution conditioned on latent state is learned.
                \tInference is either joint (default) or marginal (marginal allows a branch-point to be nominated;\s
                \tif one is not given all uninstantiated nodes are inferred)
                """);
        System.exit(error);
    }

    public enum MODEL_MODE {DIRECT, LATENT}
    private static final int TSV = 0;
    private static final int TREE = 1;
    private static final int STDOUT = 2;
    private static final int ITOL = 3;
    private static final int NODE = 0;
    private static final int VALUE = 1;
    private static final int MEAN = 1;
    private static final int SD = 2;
    private static final int IWD_VAL = 3;
    private static final int UCB_VAL = 4;
    private static final int DEFAULT_VALUES_IDX = 1;
    private static final int DEFAULT_ENTRIES_IDX = 0;
    private static final String[] FORMATS = new String[] {"TSV", "TREE", "STDOUT", "ITOL"};
    private static final String[] MODELS = new String[] {"uniform"};


    public static void main(String[] args) {
        String NEWICK = null;
        String OUTPUT = null;
        String INPUT = null;
        String LABEL = null; // the label to model in the datafile, if null/not specified assume first column
        String COLPARSER = null; // the parser for the column (in turn identified by LABEL)
        String PARAMS_FILE = null; // the name of the file from which parameters are read and/or written
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel SUBST_MODEL = null;
        int FORMAT_IDX = 0;
        asr.GRASP.Inference INFERENCE_MODE = GRASP.Inference.JOINT; // If LEARN is not specified, this will default to Joint
        MODEL_MODE MODE = MODEL_MODE.DIRECT; // not latent, discrete labels
        boolean VERBOSE = false; // print out various messages during processing to inform the user
        boolean LEARN = false; // train (when true) or infer (when false, default)
        boolean LEAVES_ONLY = true; // leaves-only are equipped with distributions (when true, default), or all nodes including internal (when false)
        boolean TIED_VARIANCE = false;
        Object MARG_LABEL = null;
        double GAMMA = 1.0;
        double LAMBDA = 5.0; // for upper confidence bound of predicted values
        Long SEED = 1L;
        Integer NSTATES = null; // number of states if continuous; if discrete, the same as the number of used values in the input file
        JSONObject PARAMS = null; // optional parameters as JSON
        int NBINS = 10; // ITOL save double vals binned in this many
        int NSAMPLES = 500; // when estimating mean and variance of inferred mixture of Gaussians

        for (int a = 0; a < args.length; a ++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (arg.equalsIgnoreCase("nwk") && args.length > a + 1) {
                    NEWICK = args[++a];
                } else if (arg.equalsIgnoreCase("out") && args.length > a + 1) {
                    OUTPUT = args[++a];
                } else if (arg.equalsIgnoreCase("in") && args.length > a + 1) {
                    INPUT = args[++a];
                    int atsign = INPUT.indexOf('@');
                    if (atsign >= 0) { // there is a label to extract
                        LABEL = INPUT.substring(0, Math.max(0, atsign));
                        int colon = LABEL.indexOf(':');
                        if (colon >= 0) {
                            COLPARSER = LABEL.substring(colon + 1);
                            LABEL = LABEL.substring(0, Math.max(0, colon));
                        }
                        INPUT = INPUT.substring(atsign + 1);
                    }
                } else if (arg.equalsIgnoreCase("tied")) {
                    TIED_VARIANCE = true;
                } else if (arg.equalsIgnoreCase("joint")) {
                    INFERENCE_MODE = asr.GRASP.Inference.JOINT;
                } else if (arg.equalsIgnoreCase("marg")) {
                    INFERENCE_MODE = asr.GRASP.Inference.MARGINAL;
                    if (args.length > a + 1) { // possibly there is a specification of node
                        if (!args[a + 1].startsWith("-")) {
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
                } else if (arg.equalsIgnoreCase("lambda") && args.length > a + 1) {
                    try {
                        LAMBDA = Double.parseDouble(args[a + 1]);
                        if (GAMMA <= 0)
                            throw new NumberFormatException("Cannot be zero or negative");
                    } catch (NumberFormatException e) {
                        usage(3, args[a + 1] + " is not a valid lambda value (must be a non-zero, positive real value)");
                    }
                } else if (arg.equalsIgnoreCase("latent") && args.length > a + 1) {
                    try {
                        NSTATES = Integer.parseInt(args[a + 1]);
                        if (NSTATES <= 1)
                            throw new NumberFormatException("Cannot be zero, one or negative");
                        MODE = MODEL_MODE.LATENT;
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
                    PARAMS_FILE = args[a + 1];
                    String contents = null;
                    try {
                        Path path = Paths.get(PARAMS_FILE);
                        if (Files.exists(path)) {
                            contents = Utils.LoadStringsFrom(PARAMS_FILE);
                            PARAMS = new JSONObject(contents);
                        }
                    } catch (IOException e) {
                        usage(11, "Invalid JSON file for external nodes: " + args[a + 1]);
                    } catch (JSONException e) {
                        usage(9, "Invalid JSON specification of distribution for external nodes: " + contents);
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
                } else if ((arg.equalsIgnoreCase("help"))  || (arg.equalsIgnoreCase("h"))) {
                    usage();
                } else if ((arg.equalsIgnoreCase("verbose")) || (arg.equalsIgnoreCase("v"))) {
                    VERBOSE = true;
                }
            }
        }

        //*************************************************************
        //Commandline arguments are set, now check if we're good to go
        //*************************************************************/
        if (INPUT == null || NEWICK == null) {
            usage(4, "Both an input-file and a tree-file are required");
        }

        // should be OK to load tree and input at least
        Object[][] inputs = null;
        Tree tree = null;
        try {
            inputs = TSVFile.loadObjects(INPUT); // load table with known annotations (for testing or training)
            tree = Tree.load(NEWICK, "newick"); // load tree
        } catch (IOException e) {
            usage(5, inputs == null ? "Failed to load the input-file " + INPUT : "Failed to load the tree-file " + NEWICK);
        }

        if (OUTPUT == null) { // default prefix is the (prefix of) newick filename
            int idx2 = NEWICK == null ? 0 : NEWICK.lastIndexOf(".");
            if (idx2 == -1)
                idx2 = NEWICK.length();
            int idx1 = NEWICK == null ? 0 : NEWICK.lastIndexOf("/") + 1;
            OUTPUT = NEWICK == null ? "" : NEWICK.substring(idx1, idx2);
        }

        // put data in TSV instance, header always has to be applied
        TSVFile tsv = new TSVFile(inputs, true);

        //********************************************************************************
        // Two modes: direct and latent; set models accordingly (and do appropriate checks)
        // ********************************************************************************/
        int valcol = DEFAULT_VALUES_IDX; // column with values, 1 being the default if no label is specified
        String[] discreteVals = null; // if applicable, nominated discrete values conditioned on latent states

        if (LABEL != null) {
            valcol = tsv.getColumnIndex(LABEL); // label is given so make sure to adjust column
            if (valcol == -1)
                usage(13, "Invalid name of column: " + LABEL);
        } else {
            LABEL = tsv.getHeader(DEFAULT_VALUES_IDX);
        }
        Object[] ENTRIES_OBJ = tsv.getColData(DEFAULT_ENTRIES_IDX); // entry-names as an array of Object (perhaps a mix of String and Integer, since GRASP uses numeric ancestor IDs)
        String[] ENTRIES = new String[ENTRIES_OBJ.length]; // entry-names as an array of String
        Set<Object> ENTRIES_SET = new HashSet<>();
        unpackEntryNames(ENTRIES_OBJ, ENTRIES, ENTRIES_SET, tree);

        Object[] ENTRY_VALUES = tsv.getColData(valcol, COLPARSER); // values associated with entries

        // mode is "LATENT",
        // so invent latent variable with NSTATES values
        if (MODE == MODEL_MODE.LATENT) {

            if (NSTATES == null) {
                usage(14, "-latent <NSTATES> must be specified");
            }

            Object[] alpha = createLatentStateLabels(NSTATES, LABEL);
            SUBST_MODEL = new JC(GAMMA, alpha); // set the evolutionary model (based on user specified params and the alphabet established above)

            /* values are either real or discrete; we establish this from the user-provided table */
            if (TSVFile.isDoubleOrInt(ENTRY_VALUES)) { // the value column of data is all Double (real), hence must use latent states for nodes in tree
                ENTRY_VALUES = TSVFile.getDouble(ENTRY_VALUES);
                if (VERBOSE) {
                    System.out.println("Detected real values in data--using Gaussian mixture with " + NSTATES + " components");
                }
            } else { // the column has discrete values (not real)
                discreteVals = checkDiscreteDataValidForLatentModel(ENTRY_VALUES, NSTATES, VERBOSE);
            }

        // mode is "DIRECT", so reject real values (if nominated), and
        // consolidate discrete values to nominate what states are possible
        } else if (MODE == MODEL_MODE.DIRECT) {

            if (TSVFile.isDoubleOrInt(ENTRY_VALUES)) { // the col of data is all Double (real), hence must use latent states for nodes in tree
                usage(8, "Real-value (continuous) data requires the number of latent states to be specified with -latent <NSTATES>");
            }

            Set<Object> observedStates  = getDiscreteStates(ENTRY_VALUES);
            NSTATES = observedStates.size();
            SUBST_MODEL = setupDirectModel(observedStates, GAMMA, MODEL_IDX, NSTATES);
        }


        // We should now have a working model, so check this is the case...
        if (SUBST_MODEL == null)
            usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
        else if (VERBOSE) {
            System.out.println(SUBST_MODEL.toJSON());
        }

        // Other checks...
        if (LEARN && PARAMS_FILE == null) {
            usage(12, "Learnt distribution must be saved to a parameter file, use -params to specify parameter file");
        } else if (LEARN && MODE == MODEL_MODE.DIRECT) {
            usage(16, "-learn must be used with -latent <NSTATES> to learn parameters of the latent model.");
        }

        /*
         * Processes (options):
         * -- Direct mode, joint recon inference (only discrete)
         * -- Direct mode, marginal recon inference (only discrete)
         * -- Latent mode,
         *      --learning (both discrete and real)
         *      --marginal recon inference (both discrete and real)
         *      --joint recon inference (both discrete and real)
         */

        TreeInstance ti = tree.getInstance(ENTRIES_OBJ, ENTRY_VALUES);  // create a tree with values attached to name-matched nodes

        // now decide type of operation...
        if (MODE == MODEL_MODE.DIRECT && INFERENCE_MODE == GRASP.Inference.JOINT) {
            performDirectJointInference(tree, SUBST_MODEL, ti, ENTRIES_SET, FORMAT_IDX, tsv, OUTPUT, valcol);

        } else if (MODE == MODEL_MODE.DIRECT && INFERENCE_MODE == GRASP.Inference.MARGINAL) {
            performDirectMarginalInference(MARG_LABEL, tree, ENTRIES_SET, tsv, FORMAT_IDX, ti, SUBST_MODEL, OUTPUT, valcol);

        } else if (MODE == MODEL_MODE.LATENT) {

            PhyloBN pbn = setupPhyloBN(discreteVals, tree, SUBST_MODEL, LEAVES_ONLY,
                    TIED_VARIANCE, SEED, ENTRY_VALUES, PARAMS, LEARN, VERBOSE);

            if (LEARN) { // learn, do not infer
                trainPBN(pbn, ENTRIES, ENTRY_VALUES, SEED, VERBOSE,
                        FORMAT_IDX, OUTPUT, tree, tsv, valcol, PARAMS_FILE);

            } else if (INFERENCE_MODE == GRASP.Inference.MARGINAL) {

                performLatentMarginalInference(tree, MARG_LABEL, ENTRY_VALUES,  tsv,  FORMAT_IDX, ti,  pbn,
                        OUTPUT,  valcol, NSAMPLES,  LAMBDA,  NBINS);

            } else if (INFERENCE_MODE == GRASP.Inference.JOINT) {
                performLatentJointInference(tree, pbn, tsv, ti, ENTRIES_SET, FORMAT_IDX, OUTPUT, valcol);
            }
        }
    }

    /**
     * Unpack entry names from mixed Object[] (String and Integer) to String[],
     * checking against tree labels and ancestor IDs. Also fills a Set with the actual labels used.
     * @param entryObjects
     * @param entryStrings
     * @param entrySet
     * @param tree
     */
    private static void unpackEntryNames(Object[] entryObjects, String[] entryStrings,
                                         Set<Object> entrySet, Tree tree) {
        for (int i = 0; i < entryStrings.length; i ++) {
            int bpidx = tree.getIndex(entryObjects[i]); // this search for branchpoint, does not care if ancestor is specified as string with "N" prefix or as ancestor number (Integer)
            if (bpidx < 0) { // not found
                System.err.println("No node in tree with label or ancestor ID \"" + entryObjects[i] + "\". Ignored.");
                entryStrings[i] = entryObjects[i].toString();
                entrySet.add(entryObjects[i]);
            } else {
                entryStrings[i] = tree.getLabel(bpidx).toString(); //
                entrySet.add(tree.getLabel(bpidx));
            }
        }
    }

    private static Object[] createLatentStateLabels(Integer numStates, String label) {
        Object[] alpha = new Object[numStates]; // alphabet (symbols that are used to represent states)
        // label the latent states; use A, B, ... if no LABEL has been given
        if (label == null) { // no LABEL specified so values/latent states will be A1, A2, ...
            for (int i = 0; i < alpha.length; i++)
                alpha[i] = (char) ('A' + i);
        } else { // use the given LABEL plus underscore then a number, e.g. kcat_1, kcat_2, ...
            for (int i = 0; i < alpha.length; i++)
                alpha[i] = label + "_" + (i + 1);
        }

        return alpha;
    }

    private static String[] checkDiscreteDataValidForLatentModel(Object[] ENTRY_VALUES, Integer NSTATES, boolean VERBOSE) {

        Set<Object> observed = new HashSet<>();
        for (int i = 0; i < ENTRY_VALUES.length; i ++)
            if (ENTRY_VALUES[i] != null)
                observed.add(ENTRY_VALUES[i]);
        String[] discreteVals = new String[observed.size()];
        try {
            observed.toArray(discreteVals);

        } catch (ArrayStoreException e) {
            usage(15, "Discrete values cannot contain integers");
        }

        if (VERBOSE) {
            System.out.println("Detected " + discreteVals.length + " discrete values in data--using multi-nomial distributions conditioned on " + NSTATES + " states");
        }

        return discreteVals;
    }

    private static SubstModel setupDirectModel(Set<Object> observedStates, double GAMMA, int MODEL_IDX, Integer numStates) {

        if (MODELS[MODEL_IDX].equals("uniform")) {
            Object[] xxalpha = new Object[numStates];
            observedStates.toArray(xxalpha);
            return new JC(GAMMA, xxalpha);
        } else { // some other model
            SubstModel substModel = SubstModel.createModel(MODELS[MODEL_IDX]);
            for (Object v : observedStates) {
                if (!substModel.getDomain().isValid(v))
                    usage(1, "Model " + MODELS[MODEL_IDX] + " does not accept value/state \"" + v.toString() + "\"");
            }

            return substModel;
        }
    }

    private static Set<Object> getDiscreteStates(Object[] ENTRY_VALUES) {

        Set<Object> observed = new HashSet<>();
        for (Object entryValue : ENTRY_VALUES)
            if (entryValue != null)
                observed.add(entryValue);

        return observed;

    }

    private static void saveDirectOutput(int FORMAT_IDX, String OUTPUT,
                                   Tree tree, TSVFile tsv, int valCol, TSVFile tempTSV, Integer NBINS, Object[] ENTRY_VALUES) {
        try {
            switch (FORMAT_IDX) {
                case TSV:
                    tempTSV.save(OUTPUT + ".tsv"); break;
                case TREE:
                    Newick.save(tree, OUTPUT + ".nwk", tempTSV.getColData(1)); break;
                case ITOL:
                    Object[][] matrix = TSVFile.Transpose(tempTSV.getRows());
                    if (ENTRY_VALUES == null) {
                        TSVFile.save2iTOL(OUTPUT + ".itol", matrix[NODE], matrix[MEAN], tsv.getHeader(valCol));
                    } else {
                        TSVFile.save2iTOL(OUTPUT + ".itol", matrix[NODE], matrix[MEAN], matrix[SD], tsv.getHeader(valCol), NBINS, TSVFile.getMin(ENTRY_VALUES), TSVFile.getMax(ENTRY_VALUES));
                    }

                    break;
            }
        } catch (IOException e) {
            usage(6, "Failed to save output to " + OUTPUT + " using format " + FORMATS[FORMAT_IDX]);
        }
    }

    /**
     * Joint inference of all nodes in the tree under a latent model.
     */
    private static void performLatentJointInference(Tree tree, PhyloBN pbn, TSVFile tsv, TreeInstance ti,
                                                    Set<Object> ENTRIES_SET, int FORMAT_IDX, String OUTPUT, int valcol) {

        MaxLhoodJoint mlj = new MaxLhoodJoint(pbn);
        mlj.decorate(ti);
        String[] headers = new String[] {tsv.getHeader(0), tsv.getHeader(1)};
        Object[][] save = new Object[tree.getSize()][headers.length];
        int bpcnt = 0;
        for (int bpidx : tree) {
            if (!ENTRIES_SET.contains(tree.getLabel(bpidx))) {
                Object d = mlj.getDecoration(bpidx);
                if (d != null) {
                    save[bpcnt][0] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                    save[bpcnt][1] = d;
                    bpcnt += 1;
                }
            }
        }
        TSVFile tempTsv = new TSVFile(headers, save);
        saveDirectOutput(FORMAT_IDX, OUTPUT, tree, tsv, valcol, tempTsv, null, null);
    }


    /* ---------- Direct, joint recon inference ----------
     * Create an inference instance based on tree (with values attached) and model.
     * (Checks that values are valid and discrete already done.)
     */
    private static void performDirectJointInference(Tree tree, SubstModel SUBST_MODEL, TreeInstance ti,
                                                    Set<Object> ENTRIES_SET, int FORMAT_IDX, TSVFile tsv,
                                                    String OUTPUT, int valcol) {

        MaxLhoodJoint inf = new MaxLhoodJoint(tree, SUBST_MODEL);
        inf.decorate(ti); // attach states to instantiated nodes and perform (joint) inference (based on model and tree)
        Object[][] save = new Object[tree.getSize()][2]; // matrix to store all entry-names and their states as values
        int bpcnt = 0;
        for (int bpidx : tree) {
            if (!ENTRIES_SET.contains(tree.getLabel(bpidx))) {
                Object state = inf.getDecoration(bpidx); // retrieve the already inferred state for the specified node
                save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx); // if ancestor attach N-prefix
                save[bpcnt][VALUE] = state;
                if (FORMAT_IDX == STDOUT)
                    System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + state);
                bpcnt += 1;
            }
        }
        TSVFile tempTsv = new TSVFile(new String[] {tsv.getHeader(0), tsv.getHeader(valcol)}, save);
        saveDirectOutput(FORMAT_IDX, OUTPUT, tree, tsv, valcol, tempTsv, null, null);
    }

    /**
     *
     *  ---------- Direct, marginal recon inference ----------
     * Create an inference instance based on tree (with values attached) and model.
     * Need to check nominated ancestor ID or extant name from the commandline (that has been stripped from an optional N-prefix).
     * (Checks that values are valid and discrete already done.)
     */
    private static void performDirectMarginalInference(Object MARG_LABEL, Tree tree, Set<Object> ENTRIES_SET,
                                                       TSVFile tsv, int FORMAT_IDX, TreeInstance ti, SubstModel SUBST_MODEL,
                                                       String OUTPUT, int valcol) {

        Integer[] bpidxs = collectNodesForDirectInference(MARG_LABEL, tree, ENTRIES_SET);

        Object[][] save = new Object[tree.getSize()][]; // matrix to store all entry-names and their distributions
        String[] headers = null;
        int bpcnt = 0;
        for (int bpidx : bpidxs) {
            if (ti.getInstance(bpidx) != null)
                continue;
            // inference below; first create the inference instance
            MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(bpidx, tree, SUBST_MODEL);
            inf.decorate(ti); // attach states to instantiated nodes and perform (marginal) inference for ONE ancestor or extant (based on model and tree)
            EnumDistrib distrib = inf.getDecoration(bpidx); // retrieve the distribution for the nominated ancestor or extant
            if (distrib != null) {
                if (headers == null) {
                    headers = new String[1 + distrib.getDomain().size()];
                    headers[0] = tsv.getHeader(0);
                    for (int i = 0; i < distrib.getDomain().size(); i++)
                        headers[1 + i] = distrib.getDomain().get(i).toString();
                }
                if (FORMAT_IDX == STDOUT)
                    System.out.println(tree.getLabel(bpidx) + "\t" + tree.getLabel(bpidx) + "\t" + distrib);
                save[bpcnt] = new Object[1 + distrib.getDomain().size()];
                save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                for (int i = 0; i < distrib.getDomain().size(); i++)
                    save[bpcnt][1 + i] = distrib.get(i);
                bpcnt += 1;
            }
        }

        TSVFile tempTsv = new TSVFile(headers, save);
        saveDirectOutput(FORMAT_IDX, OUTPUT, tree, tsv, valcol, tempTsv, null, null);

    }

    private static Integer[] collectNodesForDirectInference(Object MARG_LABEL, Tree tree, Set<Object> ENTRIES_SET) {
        Integer[] bpidxs;
        if (MARG_LABEL != null) { // node label provided, so perform ONE round of inference
            // retrieve the branchpoint index for the nominated ancestors
            int bpidx = tree.getIndex(MARG_LABEL);
            if (bpidx < 0) // did not find it...
                usage(5, "Invalid branch point name " + MARG_LABEL);
            bpidxs = new Integer[]{bpidx};
        } else { // perform inference for ALL nodes, EXCEPT those instantiated
            List<Integer> arrlst = new ArrayList<>();
            for (int bpidx : tree) {
                Object label = tree.getLabel(bpidx);
                if (!ENTRIES_SET.contains(label))
                    arrlst.add(bpidx);
            }
            bpidxs = new Integer[arrlst.size()];
            arrlst.toArray(bpidxs);
        }

        return bpidxs;
    }

    private static Integer[] collectNodesForLatentInference(Tree tree, Object MARG_LABEL) {

        Integer[] bpidxs; // put node indices to be inferred in an array...
        if (MARG_LABEL != null) { // node label provided, so perform ONE round of inference
            // retrieve the branchpoint index for the nominated ancestors
            int bpidx = tree.getIndex(MARG_LABEL);
            if (bpidx < 0) // did not find it...
                usage(5, "Invalid branch point name " + MARG_LABEL);
            bpidxs = new Integer[]{bpidx};
        } else { // perform inference for ALL nodes,
            List<Integer> arrlst = new ArrayList<>(); // unknown nodes
            for (int bpidx : tree) {
                arrlst.add(bpidx);
            }

            bpidxs = new Integer[arrlst.size()];
            arrlst.toArray(bpidxs);
        }

        return bpidxs;
    }

    /**
     * Set up a PhyloBN instance based on whether discrete values are provided or not.
     * Also applies any pre-specified parameters if provided and if learning is not being performed.
     *
     * @param discreteVals
     * @param tree
     * @param SUBST_MODEL
     * @param LEAVES_ONLY
     * @param TIED_VARIANCE
     * @param SEED
     * @param ENTRY_VALUES
     * @param PARAMS
     * @param LEARN
     * @param VERBOSE
     * @return
     */
    private static PhyloBN setupPhyloBN(String[] discreteVals, Tree tree, SubstModel SUBST_MODEL,
                                        boolean LEAVES_ONLY, boolean TIED_VARIANCE, Long SEED, Object[] ENTRY_VALUES,
                                        JSONObject PARAMS, boolean LEARN, boolean VERBOSE) {

        PhyloBN pbn;
        if (discreteVals == null) { // specified values are "real", so equip "ext" nodes as GDTs
            pbn = PhyloBN.withGDTs(tree, SUBST_MODEL, 1, LEAVES_ONLY, SEED);
            GDT gdt = pbn.getMasterGDT();
            gdt.setTieVariances(TIED_VARIANCE ? GDT.VARIANCE_TIED_POOLED : GDT.VARIANCE_UNTIED);
            gdt.randomize(ENTRY_VALUES, SEED.intValue());
        } else  // specified values are discrete and nominated as strings in "xalpha", so add multi-nomial CPTs
            pbn = PhyloBN.withCPTs(tree, SUBST_MODEL, discreteVals, 1, LEAVES_ONLY, SEED);

        if (PARAMS != null) {
            if (!LEARN) {
                // only override if inference is being performed
                pbn.overrideMasterJSON(PARAMS);
            }

            if (VERBOSE) {
                System.out.println("Using " + (!LEARN ? "pre-set" : "initial") + " distribution: " + pbn);
            }
        }

        return pbn;
    }

    private static void trainPBN(PhyloBN pbn, String[] ENTRIES, Object[] ENTRY_VALUES,
                                 Long SEED, boolean VERBOSE, int FORMAT_IDX, String OUTPUT,
                                 Tree tree, TSVFile tsv, int valcol, String PARAMS_FILE) {

        // train the BN; headers correspond to labels of leaves or any other node that has been nominated in the input TSV file
        pbn.trainEM(ENTRIES, new Object[][]{ENTRY_VALUES}, SEED);
        // training done
        if (VERBOSE) {
            System.out.println("Learned parameters specified as:");
            System.out.println("\t-params " + pbn.getMasterJSON().toString());
        }
        if (FORMAT_IDX == ITOL && OUTPUT != null) {
            Object[][] tsave = new Object[2][ENTRIES.length];
            for (int i = 0; i < ENTRIES.length; i++) {
                int bpidx = tree.getIndex(ENTRIES[i]);
                if (bpidx == -1)
                    continue;
                tsave[0][i] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                if (ENTRY_VALUES[i] != null)
                    tsave[1][i] = ENTRY_VALUES[i];
            }
            try {
                TSVFile.save2iTOL(OUTPUT, tsave[0], tsave[1], tsv.getHeader(valcol));
            } catch (IOException e) {
                usage(6, "Failed to save output to " + OUTPUT + " using format " + FORMATS[FORMAT_IDX]);
            }
        }

        if (PARAMS_FILE != null) { // save parameters
            try {
                Utils.SaveStringTo(PARAMS_FILE, pbn.getMasterJSON().toString());
            } catch (IOException e) {
                usage(12, "Failed to save parameters to " + PARAMS_FILE);
            }
        }
    }

    private static BNode findExampleNode(Integer[] bpidxs, PhyloBN pbn) {
        BNode example = null;
        for (int bpidx : bpidxs) { // find one node which has a distribution to check what kind it is
            example = pbn.getExtNode(bpidx);
            if (example != null) // found one
                break;
        }
        if (example == null)
            throw new RuntimeException("Invalid setting with latent nodes for marginal inference (external nodes not set)");

        return example;
    }

    private static String[] buildTsvHeader(PhyloBN pbn, BNode exampleNode, TSVFile tsv) {
        String[] headers;
        if (pbn.getMasterCPT() == null && pbn.getMasterGDT() != null) { // Gaussian, so real value
            headers = new String[] {tsv.getHeader(0), tsv.getHeader(1)+" (Mean)", tsv.getHeader(1)+" (SD)", tsv.getHeader(1)+" (IWD)", tsv.getHeader(1)+" (UCB)"};

        } else if (pbn.getMasterCPT() != null && pbn.getMasterGDT() == null) { // Discrete

            // extract possible values so that we work out what probs/states that are available
            Enumerable ed = (Enumerable) exampleNode.getVariable().getDomain();
            headers = new String[1 + ed.size()];
            headers[0] = tsv.getHeader(0);
            for (int i = 0; i < ed.size(); i++)
                headers[1 + i] = ed.get(i).toString(); // header to use for column
        } else {
            throw new RuntimeException("Invalid setting with latent nodes for marginal inference (type not defined)");
        }

        return headers;
    }

    private static void processUninstantiatedNode(Distrib anydistrib, int FORMAT_IDX, Tree tree, Object[][] save,
                                                  int bpidx, int bpcnt, int NSAMPLES, double LAMBDA) {
        try {
            EnumDistrib distrib = (EnumDistrib) anydistrib; // if this cast succeeds, it is a discrete distribution
            if (FORMAT_IDX == STDOUT)
                System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + distrib);
            save[bpcnt][0] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
            for (int i = 0; i < distrib.getDomain().size(); i++)
                save[bpcnt][1 + i] = distrib.get(i);

        } catch (ClassCastException e) { // Mixture of Gaussians, probably
            try {

                double[] samples = new double[NSAMPLES];
                for (int i = 0; i < NSAMPLES; i++) {
                    samples[i] = (Double) anydistrib.sample();
                }

                GaussianDistrib gd = GaussianDistrib.estimate(samples);
                if (FORMAT_IDX == STDOUT)
                    System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + gd.getMean() + "\t" + Math.sqrt(gd.getVariance()) + "\t" + anydistrib);
                save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                save[bpcnt][MEAN] = gd.getMean();
                save[bpcnt][SD] = Math.sqrt(gd.getVariance()); // standard deviation
                save[bpcnt][UCB_VAL] = gd.getMean() + LAMBDA * Math.sqrt(gd.getVariance());


            } catch (ClassCastException ee) {

                if (FORMAT_IDX == STDOUT)
                    System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + anydistrib);
                save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                save[bpcnt][MEAN] = anydistrib;
            }
        }
    }

    private static void processInstantiatedNode(PhyloBN pbn, TreeInstance ti, int FORMAT_IDX, Tree tree,
                                                Object[][] save, int bpidx, int bpcnt, int NSAMPLES) {


        Object instance = pbn.getExtNode(bpidx).getInstance();

        double instance_dub = (Double) instance;
        MaxLhoodMarginal<EnumDistrib> instan_inf = new MaxLhoodMarginal<>(bpidx, pbn);
        ti.setInstance(bpidx, null); // remove evidence, treat as uninstantiated

        // perform marginal inference
        instan_inf.decorate(ti);
        // retrieve the distribution at the node previously nominated
        Distrib instan_anydistrib = instan_inf.getDecoration(bpidx);

        ti.setInstance(bpidx, instance); // replace the value in the tree

        // cast under assumption using a Gaussian mixture
        MixtureDistrib mixtureDistrib = (MixtureDistrib) instan_anydistrib;

        double[] samples = new double[NSAMPLES];
        for (int i = 0; i < NSAMPLES; i++) {
            samples[i] = (Double) instan_anydistrib.sample();
        }

        GaussianDistrib gd = GaussianDistrib.estimate(samples);

        double abs_error = Math.abs(instance_dub - gd.getMean());
        double iwd = mixtureDistrib.cdf(gd.getMean() + abs_error) - mixtureDistrib.cdf(gd.getMean() - abs_error);
        save[bpcnt][IWD_VAL] = iwd;
        save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
        save[bpcnt][MEAN] = instance; // here this is not actually the mean, but the instantiated value

        if (FORMAT_IDX == STDOUT)
            System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + instance);
    }

    private static void performLatentMarginalInference(Tree tree, Object MARG_LABEL,
                                                       Object[] ENTRY_VALUES, TSVFile tsv, int FORMAT_IDX,
                                                       TreeInstance ti, PhyloBN pbn, String OUTPUT, int valcol,
                                                       int NSAMPLES, double LAMBDA, Integer NBINS) {

        Object[][] save; // the matrix in which all results are stored, incl the entry-name then the probability distribution (or sample thereof)
        // marginal inference of ONE or ALL nodes in the tree
        Integer[] bpidxs = collectNodesForLatentInference(tree, MARG_LABEL);

        // find one node which has a distribution to check what kind it is
        BNode exampleNode = findExampleNode(bpidxs, pbn);
        String[] headers = buildTsvHeader(pbn, exampleNode, tsv);

        save = new Object[bpidxs.length][headers.length];
        int bpcnt = 0; // count nodes that are inferred (excl those that are null)
        for (int bpidx : bpidxs) { // go through all nodes to be inferred
            // inference below; first create the inference instance
            MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal<>(bpidx, pbn);

            // perform marginal inference
            inf.decorate(ti);
            // retrieve the distribution at the node previously nominated
            Distrib anydistrib = inf.getDecoration(bpidx);

            if (anydistrib != null) {
                processUninstantiatedNode(anydistrib, FORMAT_IDX, tree, save, bpidx,
                        bpcnt, NSAMPLES, LAMBDA);
            } else {
                // anydistrib == null; this happens when a node is instantiated, so retrieve the value accordingly
                processInstantiatedNode(pbn, ti, FORMAT_IDX, tree, save, bpidx, bpcnt, NSAMPLES);
            }

            bpcnt += 1;
        }

        TSVFile tempTsv = new TSVFile(headers, save);
        saveDirectOutput(FORMAT_IDX, OUTPUT, tree, tsv, valcol, tempTsv, NBINS, ENTRY_VALUES);

    }
}
