package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import dat.file.Newick;
import dat.file.TSVFile;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;

import java.io.IOException;
import java.io.PrintStream;
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
                "\t{-joint (default) | -marg <branchpoint-id>} \n" +
                "\t{-format <TSV(default), TREE, STDOUT>}\n");
        out.println("where \n" +
                "\ttree-file is a phylogenetic tree on Newick format\n" +
                "\tinput-file is a table with sequence or ancestor names, and corresponding values (empty implies not assigned) on TSV format\n" +
                "\toutput-file will be populated by inferred branch point states on specified format (TSV by default, TREE is a labelled tree on Newick format)\n" +
                "\tInference is either joint (default) or marginal (marginal requires a branch-point to be nominated)\n");
        out.println("Notes: \n" +
                "\tEvolutionary models are currently limited to uniform, which is an adaptation of Jukes-Cantor for arbitrary number of states.\n\tgamma-value is used by this model\n");
        System.exit(error);
    }

    public static void main(String[] args) {
        String NEWICK = null;
        String OUTPUT = null;
        String INPUT = null;

        String[] MODELS = new String[] {"uniform"};
        int MODEL_IDX = 0; // default model is that above indexed 0
        SubstModel MODEL = null;
        String[] FORMATS = new String[] {"TSV", "TREE", "STDOUT"};
        int FORMAT_IDX = 0;
        asr.GRASP.Inference MODE = asr.GRASP.Inference.JOINT;
        Object MARG_LABEL = null;
        Double GAMMA = 1.0;

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
                    MODE = asr.GRASP.Inference.JOINT;
                } else if (arg.equalsIgnoreCase("marg") && args.length > a + 1) {
                    MODE = asr.GRASP.Inference.MARGINAL;
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
            TSVFile tsv = new TSVFile(inputs);
            Set values = tsv.getValues(1);
            if (MODELS[MODEL_IDX].equals("uniform")) {
                Object[] alpha = new Object[values.size()];
                values.toArray(alpha);
                MODEL = new JC(GAMMA, alpha);
            } else {
                MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
            }
            if (MODEL == null) {
                usage(1, "Model " + MODELS[MODEL_IDX] + " could not be created");
            }
            Object[][] all = tsv.getCols(new int[] {0, 1});
            TreeInstance ti = tree.getInstance(all[0], all[1]);
            if (MODE == GRASP.Inference.JOINT) {
                MaxLhoodJoint inf = new MaxLhoodJoint(tree, MODEL);
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
            } else if (MODE == GRASP.Inference.MARGINAL) {
                int bpidx = tree.getIndex(MARG_LABEL);
                if (bpidx < 0) {
                    usage(5, "Invalid branch point name " + MARG_LABEL);
                }
                MaxLhoodMarginal inf = new MaxLhoodMarginal(bpidx, tree, MODEL);
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
            }
        } else {
            usage(4, "Both an input-file and a tree-file are required");
        }


    }

}
