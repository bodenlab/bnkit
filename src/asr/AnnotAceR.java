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
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import json.JSONException;
import json.JSONObject;
import stats.ZScore;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;



public class AnnotAceR {

    public static void usage() {
        usage(0, null);
    }
    public static void usage(int error, String msg) {
        PrintStream out = System.out;
        if (error != 0)
            out = System.err;
        if (msg != null)
            out.println(msg + " (Error " + error + ")");
        out.println("Usage: asr.AnnotAceR \n" +
                "\t[-nwk <tree-file> -in {<label>{:<parser>}@}<input-file> -out <output-file>]\n" +
                "\t{-model <uniform(default)>}\n" +
                "\t{-gamma <gamma-value(default 1.0)>}\n" +
                "\t{-params <JSON-file>}\n" +
                "\t{-latent <#states>}\n" +
                "\t{-internal}\n" +
                "\t{-learn}\n" +
                "\t{-tied} \n" +
                "\t{-seed <seed>} \n" +
                "\t{-joint (default) | -marg {<branchpoint-id>} } \n" +
                "\t{-format <TSV(default), TREE, STDOUT, ITOL>}\n" +
                "\t{-help|-h}\n" +
                "\t{-verbose|-v}\n");
        out.println("where \n" +
                "\ttree-file is a phylogenetic tree on Newick format\n" +
                "\tinput-file is a table with sequence or ancestor names in the first column, and corresponding values\n\t(empty or None or null implies not assigned) on TSV format\n" +
                "\tlabel flags that a header is used in the input-file and identifies the column with values to be modelled;\n\tif no label is given, headers are assumed absent and values from the second column will be modelled\n" +
                "\tparser identifies a parser to use on the column with values (e.g. BRENDA).\n" +
                "\toutput-file will contain:\n" +
                "\t- inferred branch point states on specified format (TSV by default, TREE is a labelled tree on Newick format, ITOL is a dataset to decorate trees in iTOL.embl.de), or\n" +
                "\tgamma-value is parameter to the uniform model (n-state generalisation of Jukes-Cantor)\n" +
                "\tJSON-file contains a JSON string specifying the distribution for latent nodes (if latent mode is used)\n" +
                "\tlatent indicates that the tree consists of latent values only (latent mode), with specified values as extensions to the leaves.\n" +
                    "\t- #states is the number of latent states to learn (should not exceed 25, labelled A-Z).\n" +
                "\tinternal indicates that internal nodes are also extended with user-specified or learned distributions (default leaves-only).\n" +
                "\tlearn excludes inference and instead prompts EM learning of parameters, using input data as training data.\n" +
                "\tiwd Calculates the integral weighted distance between instantiated values and their estimate. Only applies to marginal inference for real values. \n" +
                "\ttied implies that the variance learned is the same across the latent states (only applicable when EM-learning GDTs; default is off).\n" +
                "\thelp prints out commandline arguments (this screen).\n" +
                "\tverbose completes the requested steps while printing out messages about the process.\n");
        out.println("Notes:\n" +
                "\tEvolutionary models of substitution are currently limited to uniform, which is an adaptation of Jukes-Cantor for arbitrary number of states.\n" +
                    "\t- gamma-value is used by this model\n" +
                "\tIf specified values are real, a conditional Gaussian mixture distribution conditioned on latent state is learned.\n" +
                "\tIf specified values are discrete, a multinomial distribution conditioned on latent state is learned.\n" +
                "\tInference is either joint (default) or marginal (marginal allows a branch-point to be nominated; \n\tif one is not given all uninstantiated nodes are inferred)\n");
        System.exit(error);
    }

    public enum MODEL_MODE {DIRECT, LATENT};
    private static final int TSV = 0;
    private static final int TREE = 1;
    private static final int STDOUT = 2;
    private static final int ITOL = 3;
    private static final int NODE = 0;
    private static final int VALUE = 1;
    private static final int MEAN = 1;
    private static final int SD = 2;
    private static final int IWD_VAL = 3;


    public static void main(String[] args) {
        String NEWICK = null;
        String OUTPUT = null;
        String INPUT = null;
        String LABEL = null; // the label to model in the datafile, if null/not specified assume first column
        String COLPARSER = null; // the parser for the column (in turn identified by LABEL)
        String PARAMS_FILE = null; // the name of the file from which parameters are read and/or written

        String[] MODELS = new String[] {"uniform"};
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel SUBST_MODEL = null;
        String[] FORMATS = new String[] {"TSV", "TREE", "STDOUT", "ITOL"};
        int FORMAT_IDX = 0;
        asr.GRASP.Inference INFERENCE_MODE = null; // If LEARN is not specified, this will default to Joint
        MODEL_MODE mode = MODEL_MODE.DIRECT; // not latent, discrete labels
        boolean VERBOSE = false; // print out various messages during processing to inform the user
        boolean LEARN = false; // train (when true) or infer (when false, default)
        boolean LEAVES_ONLY = true; // leaves-only are equipped with distributions (when true, default), or all nodes including internal (when false)
        boolean TIED_VARIANCE = false;
        boolean IWD = false; // if true, calculate the integral weighted distance
        Object MARG_LABEL = null;
        Double GAMMA = 1.0;
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
                } else if (arg.equalsIgnoreCase("iwd")){
                    IWD = true;
                } else if ((arg.equalsIgnoreCase("help"))  || (arg.equalsIgnoreCase("h"))) {
                    usage();
                } else if ((arg.equalsIgnoreCase("verbose")) || (arg.equalsIgnoreCase("v"))) {
                    VERBOSE = true;
                }
            }
        }

        /*************************************************************
        Commandline arguments are set, now check if we're good to go
         *************************************************************/

        if (INPUT != null && NEWICK != null) { // should be OK to load tree and input at least
            Object[][] inputs = null;
            Tree tree = null;
            try {
                inputs = TSVFile.loadObjects(INPUT); // load table with known annotations (for testing or training)
                tree = Tree.load(NEWICK, "newick"); // load tree
            } catch (IOException e) {
                usage(5, inputs == null ? "Failed to load the input-file " + INPUT : "Failed to load the tree-file " + NEWICK);
            }
            TSVFile tsv = new TSVFile(inputs, LABEL != null ? true : false); // put data in TSV instance, allowing for header if a LABEL has been specified

            /********************************************************************************
             Two modes: direct and latent; set models accordingly (and do appropriate checks)
             ********************************************************************************/
            Set observed; // this is for (the set of) observed values (as established through the complete list of values across all entries)
            int valcol = 1; // column with values, 1 being the default if no label is specified
            String[] xalpha = null; // if applicable, nominated discrete values conditioned on latent states

            if (LABEL != null) {
                valcol = tsv.getColumn(LABEL); // label is given so make sure to adjust column
                if (valcol == -1)
                    usage(13, "Invalid name of column: " + LABEL);
            }
            Object[] ENTRIES_OBJ = tsv.getCol(0);       // entry-names as an array of Object (perhaps a mix of String and Integer, since GRASP uses numeric ancestor IDs)
            String[] ENTRIES = new String[ENTRIES_OBJ.length];      // entry-names as an array of String
            Set<Object> ENTRIES_SET = new HashSet<>();
            for (int i = 0; i < ENTRIES.length; i ++) {
//                int bpidx = tree.getIndex(ENTRIES_OBJ[i]); // this search for branchpoint, does not care if ancestor is specified as string with "N" prefix or as ancestor number (Integer)
//                if (bpidx < 0) { // not found
//                    System.err.println("No node in tree with label or ancestor ID \"" + ENTRIES_OBJ[i] + "\". Ignored.");
                    ENTRIES[i] = ENTRIES_OBJ[i].toString();
                    ENTRIES_SET.add(ENTRIES_OBJ[i]);
//                } else {
//                    ENTRIES[i] = tree.getLabel(bpidx).toString(); //
//                    ENTRIES_SET.add(tree.getLabel(bpidx));
//                }
            }
            Object[] ENTRY_VALUES = tsv.getCol(valcol, COLPARSER); // values associated with entries
            HashMap<String, Integer> ENTRY_TO_VAL = new HashMap<>();
            for (int i = 0; i < ENTRY_VALUES.length; i ++) {
                ENTRY_TO_VAL.put(ENTRIES_OBJ[i].toString(), i);
            }

            // mode is "LATENT",
            // so invent latent variable with NSTATES values
            if (mode == MODEL_MODE.LATENT) {
                if (NSTATES == null) {
                    usage(14, "-latent <NSTATES> must be specified");
                }
                Object[] alpha = new Object[NSTATES]; // alphabet (symbols that are used to represent states)
                // label the latent states; use A, B, ... if no LABEL has been given
                if (LABEL == null) { // no LABEL specified so values/latent states will be A1, A2, ...
                    for (int i = 0; i < alpha.length; i++)
                        alpha[i] = Character.valueOf((char) ('A' + i));
                } else { // use the given LABEL plus underscore then a number, e.g. kcat_1, kcat_2, ...
                    for (int i = 0; i < alpha.length; i++)
                        alpha[i] = LABEL + "_" + (i + 1);
                }
                SUBST_MODEL = new JC(GAMMA, alpha); // set the evolutionary model (based on user specified params and the alphabet established above)
                /* values are either real or discrete; we establish this from the user-provided table */
                if (TSVFile.isDoubleOrInt(new Object[][] {ENTRY_VALUES})) { // the first col of data is all Double (real), hence must use latent states for nodes in tree
                    observed = null;
                    ENTRY_VALUES = TSVFile.getDouble(ENTRY_VALUES);
                    if (VERBOSE)
                        System.out.println("Detected real values in data--using Gaussian mixture with " + NSTATES + " components");
                } else { // the column has discrete values (not real)
                    observed = new HashSet();
                    for (int i = 0; i < ENTRY_VALUES.length; i ++)
                        if (ENTRY_VALUES[i] != null)
                            observed.add(ENTRY_VALUES[i]);
                    xalpha = new String[observed.size()];
                    try {
                        observed.toArray(xalpha);

                    } catch (ArrayStoreException e) {
                        usage(15, "Discrete values cannot contain integers");
                    }

                    if (VERBOSE)
                        System.out.println("Detected " + xalpha.length + " discrete values in data--using multi-nomial distributions conditioned on " + NSTATES + " states");
                }
            //
            // mode is "DIRECT",
            // so reject real values (if nominated), and
            // consolidate discrete values to nominate what states are possible
            } else if (mode == MODEL_MODE.DIRECT) {
                if (TSVFile.isDouble(ENTRY_VALUES)) { // the col of data is all Double (real), hence must use latent states for nodes in tree
                    usage(8, "Real-value data requires (discrete) latent states; please specify number of latent states.");
                } else { // the column has discrete values (not real)
                    observed = new HashSet();
                    for (int i = 0; i < ENTRY_VALUES.length; i ++)
                        if (ENTRY_VALUES[i] != null)
                            observed.add(ENTRY_VALUES[i]);
                    NSTATES = observed.size();
                    if (MODELS[MODEL_IDX].equals("uniform")) {
                        Object[] xxalpha = new Object[NSTATES];
                        observed.toArray(xxalpha);
                        SUBST_MODEL = new JC(GAMMA, xxalpha);
                    } else { // some other model
                        SUBST_MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
                        for (Object v : observed) {
                            if (!SUBST_MODEL.getDomain().isValid(v))
                                usage(1, "Model " + MODELS[MODEL_IDX] + " does not accept value/state \"" + v.toString() + "\"");
                        }
                    }
                }
            }
            // We should now have a working model, so check this is the case...
            if (SUBST_MODEL == null)
                usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
            else if (VERBOSE) {
                System.out.println(SUBST_MODEL.toJSON());
            }
            // Other checks...
            if (LEARN && PARAMS_FILE == null)
                usage(12, "Learnt distribution must be saved to a parameter file, use -params to specify parameter file");
            else if (INFERENCE_MODE != null && OUTPUT == null)
                usage(13, "Inference needs to save to an output file, use -out to specify output file.");
            else if (LEARN && INFERENCE_MODE != null)
                usage(16, "-learn and -marg/-joint cannot be used together. Learning and inference must be performed in separate steps");
            /**********************************************************************************
             * Processes (options):
             * -- Direct mode, joint recon inference (only discrete)
             * -- Direct mode, marginal recon inference (only discrete)
             * -- Latent mode,
             *      --learning (both discrete and real)
             *      --marginal recon inference (both discrete and real)
             *      --joint recon inference (both discrete and real)
             **********************************************************************************/
            if (INFERENCE_MODE == null) INFERENCE_MODE = GRASP.Inference.JOINT; // assign default inference mode if not specified
            TreeInstance ti = tree.getInstance(ENTRIES_OBJ, ENTRY_VALUES);  // create a tree with values attached to name-matched nodes

            // now decide type of operation...
            if (mode == MODEL_MODE.DIRECT && INFERENCE_MODE == GRASP.Inference.JOINT) {
                /* ---------- Direct, joint recon inference ----------
                 * Create an inference instance based on tree (with values attached) and model.
                 * (Checks that values are valid and discrete already done.)
                 */
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
                TSVFile tmptsv = new TSVFile(new String[] {tsv.getHeader(0), tsv.getHeader(valcol)}, save);
                try {
                    switch (FORMAT_IDX) {
                        case TSV:
                            tmptsv.save(OUTPUT); break;
                        case TREE:
                            Newick.save(tree, OUTPUT, tmptsv.getCol(1)); break;
                        case ITOL:
                            Object[][] matrix = TSVFile.Transpose(tmptsv.getRows());
                            TSVFile.save2iTOL(OUTPUT, matrix[0], matrix[1], tsv.getHeader(valcol)); break;
                    }
                } catch (IOException e) {
                    usage(6, "Failed to save output to " + OUTPUT + " using format " + FORMATS[FORMAT_IDX]);
                }

            } else if (mode == MODEL_MODE.DIRECT && INFERENCE_MODE == GRASP.Inference.MARGINAL) {
                /* ---------- Direct, marginal recon inference ----------
                 * Create an inference instance based on tree (with values attached) and model.
                 * Need to check nominated ancestor ID or extant name from the commandline (that has been stripped from an optional N-prefix).
                 * (Checks that values are valid and discrete already done.)
                 */
                Integer[] bpidxs = null;
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
                Object[][] save = new Object[tree.getSize()][]; // matrix to store all entry-names and their distributions
                String[] myheaders = null;
                int bpcnt = 0;
                for (int bpidx : bpidxs) {
                    if (ti.getInstance(bpidx) != null)
                        continue;
                    // inference below; first create the inference instance
                    MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(bpidx, tree, SUBST_MODEL);
                    inf.decorate(ti); // attach states to instantiated nodes and perform (marginal) inference for ONE ancestor or extant (based on model and tree)
                    EnumDistrib distrib = inf.getDecoration(bpidx); // retrieve the distribution for the nominated ancestor or extant
                    if (distrib != null) {
                        if (myheaders == null) {
                            myheaders = new String[1 + distrib.getDomain().size()];
                            myheaders[0] = tsv.getHeader(0);
                            for (int i = 0; i < distrib.getDomain().size(); i++)
                                myheaders[1 + i] = distrib.getDomain().get(i).toString();
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
                TSVFile tmptsv = new TSVFile(myheaders, save);
                try {
                    switch (FORMAT_IDX) {
                        case TSV:
                            tmptsv.save(OUTPUT); break;
                        case TREE:
                            Newick.save(tree, OUTPUT, tmptsv.getCol(1)); break;
                        case ITOL:
                            Object[][] matrix = TSVFile.Transpose(tmptsv.getRows());
                            TSVFile.save2iTOL(OUTPUT, matrix[0], matrix[1], tsv.getHeader(valcol)); break;
                    }
                } catch (IOException e) {
                    usage(6, "Failed to save output to " + OUTPUT + " using format " + FORMATS[FORMAT_IDX]);
                }

            } else if (mode == MODEL_MODE.LATENT) {
                /* ---------- Latent states (all options) ---------- */
                PhyloBN pbn;
                if (xalpha == null) { // specified values are "real", so equip "ext" nodes as GDTs
                    pbn = PhyloBN.withGDTs(tree, SUBST_MODEL, 1, LEAVES_ONLY, SEED);
                    GDT gdt = pbn.getMasterGDT();
                    gdt.setTieVariances(TIED_VARIANCE ? GDT.VARIANCE_TIED_POOLED : GDT.VARIANCE_UNTIED);
                    gdt.randomize(ENTRY_VALUES, SEED.intValue());
                } else  // specified values are discrete and nominated as strings in "xalpha", so add multi-nomial CPTs
                    pbn = PhyloBN.withCPTs(tree, SUBST_MODEL, xalpha, 1, LEAVES_ONLY, SEED);

                if (PARAMS != null) {
                    if (!LEARN) {
                        // only override if inference is being performed
                        pbn.overrideMasterJSON(PARAMS);
                        if (VERBOSE)
                            System.out.println("Using pre-set distribution: " + pbn);
                    }

                    if (VERBOSE)
                        System.out.println("Using initial distribution: " + pbn);
                }

                if (LEARN) { // learn, do not infer
                    // train the BN; headers correspond to labels of leaves or any other node that has been nominated in the input TSV file
                    pbn.trainEM(ENTRIES, new Object[][] {ENTRY_VALUES}, SEED);
                    // training done
                    if (VERBOSE) {
                        System.out.println("Learned parameters specified as:");
                        System.out.println("\t-params " + pbn.getMasterJSON().toString());
                    }
                    if (FORMAT_IDX == ITOL && OUTPUT != null) {
                        Object[][] tsave = new Object[2][ENTRIES.length];
                        for (int i = 0; i < ENTRIES.length; i ++) {
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
                } else if (INFERENCE_MODE == GRASP.Inference.MARGINAL) {
                    // marginal inference of ONE or ALL nodes in the tree
                    Integer[] bpidxs = null; // put node indices to be inferred in an array...
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
                    Object[][] save = null; // the matrix in which all results are stored, incl the entry-name then the probability distribution (or sample thereof)
                    BNode example = null;
                    for (int bpidx : bpidxs) { // find one node which has a distribution to check what kind it is
                        example = pbn.getExtNode(bpidx);
                        if (example != null) // found one
                            break;
                    }
                    if (example == null)
                        throw new RuntimeException("Invalid setting with latent nodes for marginal inference (external nodes not set)");
                    String[] myheaders = null;
                    if (pbn.getMasterCPT() == null && pbn.getMasterGDT() != null) { // Gaussian, so real value
                        if (IWD) {
                            myheaders = new String[] {tsv.getHeader(0), tsv.getHeader(1)+" (Mean)", tsv.getHeader(1)+" (SD)", tsv.getHeader(1)+" (IWD)"}; // which means only one value is available (that is the average of multiple samples)
                        } else {
                            myheaders = new String[] {tsv.getHeader(0), tsv.getHeader(1)+" (Mean)", tsv.getHeader(1)+" (SD)"};
                        }

                    } else if (pbn.getMasterCPT() != null && pbn.getMasterGDT() == null) { // Discrete
                        Enumerable ed = (Enumerable) example.getVariable().getDomain(); // extract possible values so that we work out what probs/states that are available
                        myheaders = new String[1 + ed.size()];
                        myheaders[0] = tsv.getHeader(0);
                        for (int i = 0; i < ed.size(); i++)
                            myheaders[1 + i] = ed.get(i).toString(); // header to use for column
                    } else {
                        throw new RuntimeException("Invalid setting with latent nodes for marginal inference (type not defined)");
                    }

                    save = new Object[bpidxs.length][myheaders.length];
                    int bpcnt = 0; // count nodes that are inferred (excl those that are null)
                    for (int bpidx : bpidxs) { // go through all nodes to be inferred
                        // inference below; first create the inference instance
                        MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(bpidx, pbn);

                        // perform marginal inference
                        inf.decorate(ti);
                        // retrieve the distribution at the node previously nominated
                        Distrib anydistrib = inf.getDecoration(bpidx);

                        if (anydistrib != null) {
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
                                    double sum = 0;
                                    for (int i = 0; i < NSAMPLES; i++) {
                                        samples[i] = (Double) anydistrib.sample();
                                        sum += samples[i];
                                    }

                                    GaussianDistrib gd = GaussianDistrib.estimate(samples);
                                    if (FORMAT_IDX == STDOUT)
                                        System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + gd.getMean() + "\t" + Math.sqrt(gd.getVariance()) + "\t" + anydistrib);
                                    save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                                    save[bpcnt][MEAN] = gd.getMean();
                                    save[bpcnt][SD] = Math.sqrt(gd.getVariance()); // standard deviation


                                } catch (ClassCastException ee) {
                                    if (FORMAT_IDX == 2)
                                        System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + anydistrib);
                                    save[bpcnt][0] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                                    save[bpcnt][1] = anydistrib;
                                }
                            }
                            bpcnt += 1;
                        } else { // anydistrib == null; this happens when a node is instantiated, so retrieve the value accordingly
                            Object instance = pbn.getExtNode(bpidx).getInstance();

                            if (IWD) {
                                double instance_dub = (Double) instance;
                                MaxLhoodMarginal<EnumDistrib> instan_inf = new MaxLhoodMarginal(bpidx, pbn);
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

                            }

                            save[bpcnt][NODE] = (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx);
                            save[bpcnt][MEAN] = instance; // here this is not actually the mean, but the instantiated value

                            if (FORMAT_IDX == STDOUT)
                                System.out.println(bpidx + "\t" + (tree.isLeaf(bpidx) ? "" : "N") + tree.getLabel(bpidx) + "\t" + instance);

                            bpcnt += 1;
                        }
                    }
                    TSVFile tmptsv = new TSVFile(myheaders, save);
                    try {
                        switch (FORMAT_IDX) {
                            case TSV:
                                tmptsv.save(OUTPUT); break;
                            case TREE:
                                Newick.save(tree, OUTPUT, tmptsv.getCol(0), tmptsv.getCol(1)); break;
                            case ITOL:
                                Object[][] matrix = TSVFile.Transpose(tmptsv.getRows());

                                if (IWD) {
                                    TSVFile.save2iTOL("IWD_" + OUTPUT, matrix[NODE], matrix[IWD_VAL],"IWD_" + tsv.getHeader(valcol), NBINS, TSVFile.getMin(matrix[IWD_VAL]), TSVFile.getMax(matrix[IWD_VAL]));
                                }

                                TSVFile.save2iTOL(OUTPUT, matrix[NODE], matrix[MEAN], matrix[SD], tsv.getHeader(valcol), NBINS, TSVFile.getMin(ENTRY_VALUES), TSVFile.getMax(ENTRY_VALUES)); break;
                        }
                    } catch (IOException e) {
                        usage(6, "Failed to save output to " + OUTPUT + " using format " + FORMATS[FORMAT_IDX]);
                    }

                } else if (INFERENCE_MODE == GRASP.Inference.JOINT) {
                    // joint inference of ALL nodes in the tree
                    MaxLhoodJoint mlj = new MaxLhoodJoint(pbn);
                    mlj.decorate(ti);
                    String[] myheaders = null;
                    myheaders = new String[] {tsv.getHeader(0), tsv.getHeader(1)};
                    Object[][] save = new Object[tree.getSize()][myheaders.length];
                    Map<String, Object> mlvalues = new HashMap<>();
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
                    TSVFile tmptsv = new TSVFile(myheaders, save);
                    try {
                        switch (FORMAT_IDX) {
                            case 0: // TSV
                                tmptsv.save(OUTPUT); break;
                            case 1: // TREE
                                Newick.save(tree, OUTPUT, tmptsv.getCol(1)); break;
                            case 3: // ITOL
                                Object[][] matrix = TSVFile.Transpose(tmptsv.getRows());
                                TSVFile.save2iTOL(OUTPUT, matrix[0], matrix[1], tsv.getHeader(valcol)); break;
                        }
                    } catch (IOException e) {
                        usage(6, "Failed to save output to " + OUTPUT + " using format " + FORMATS[FORMAT_IDX]);
                    }

                }
            }
        } else { // we can't do anything without an input and a tree
            usage(4, "Both an input-file and a tree-file are required");
        }
    }

}
