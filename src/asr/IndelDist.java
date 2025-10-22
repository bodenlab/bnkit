package asr;

import bn.ctmc.SubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;
import smile.math.MathEx;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Random;


public class IndelDist {

    public enum RATE_CATEGORY {
        LOW,
        MEDIUM,
        HIGH
    }
    public static EnumMap<RATE_CATEGORY, Double[]> MEAN_RATES = new EnumMap<>(RATE_CATEGORY.class);
    public static Double RATE_PRIOR = Math.log(0.25);
    public static String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};

    public static void usage(int error_code, String msg) {
        PrintStream out = System.out;
        if (error_code != 0)
            out = System.err;

        out.println("""
                Usage: asr.IndelDist\s
                \t[-a | --aln <filename>]
                \t[-n | --nwk <filename>]
                \t{-s | --substitution-model <JTT(default)}
                \t{--seed}
                \t{-i | --input-folder <folder>}
                \t{-o | --output-folder <folder>} (default is current working folder, or input folder if available)
                """
        );

        if (msg != null) {
            out.println("\n" + msg + "(Error code " + error_code + ")");
        }
        System.exit(error_code);
    }

    private static <T> T getArg(HashMap<String, Object> args, String key, Class<T> type) {
        Object value = args.get(key);
        if (value == null) {
            return null;
        }

        if (!type.isInstance(value)) {
            throw new ClassCastException("Value for " +
                    key + " is not of type " + type.getName());
        }

        return type.cast(value);

    }

    public static HashMap<String, Object> parseArgs(String[] args) {

        HashMap<String, Object> argMap = new HashMap<>();

        argMap.put("SEED", new Random().nextInt());

        // Read in all the arguments
        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (((arg.equalsIgnoreCase("-aln")) || (arg.equalsIgnoreCase("a"))) && args.length > a+1) {
                    argMap.put("ALIGNMENT", args[++a]);
                } else if (((arg.equalsIgnoreCase("-nwk")) || (arg.equalsIgnoreCase("n"))) && args.length > a+1) {
                    argMap.put("NEWICK", args[++a]);
                } else if ((arg.equalsIgnoreCase("-substitution-model") || arg.equalsIgnoreCase("s")) && args.length > a+1) {
                    boolean model_found = false;
                    for (int i = 0; i < MODELS.length; i++) {
                        if (args[a + 1].equalsIgnoreCase(MODELS[i])) {
                            model_found = true;
                            argMap.put("MODEL_IDX", i);
                        }
                    }

                    if (!model_found) {
                        usage(1, "Could not find model with ID " + args[a + 1]);
                    }
                } else if ((arg.equalsIgnoreCase("-seed") && args.length > a+1)) {
                    argMap.put("SEED", Integer.parseInt(args[++a]));
                } else if ((arg.equalsIgnoreCase("-output-folder") || arg.equalsIgnoreCase("o")) && args.length > a + 1) {
                    argMap.put("OUTPUT", args[++a]);
                } else if ((arg.equalsIgnoreCase("-input-folder") || arg.equalsIgnoreCase("i")) && args.length > a + 1) {
                    argMap.put("INPUT", args[++a]);
                }
            }
        }

        // now perform checks
        if (argMap.get("ALIGNMENT") == null) {
            usage(3, "Must specify alignment (--aln <Clustal or FASTA file>)");
        } else if (argMap.get("NEWICK") == null) {
            usage(4, "Must specify phylogenetic tree (Newick file) or previously saved folder (--input-folder <folder>)");
        } else if (argMap.get("OUTPUT") == null) {
            String INPUT = null;
            try {
                INPUT = (String) argMap.get("INPUT");
            } catch (ClassCastException e) {
                usage(5, "Input must be a string");
            }

            argMap.put("OUTPUT", INPUT == null ? "." : INPUT);
        }

        return argMap;
    }


    public static void main(String[] args) {
        MEAN_RATES.put(RATE_CATEGORY.LOW, new Double[]{0.05109101808171236, 0.13936648207205818, 0.2548798168737909, 0.5313293496383602});
        MEAN_RATES.put(RATE_CATEGORY.MEDIUM, new Double[]{0.10822622872656855, 0.3236414936793391, 0.6232831770054337, 1.368582433909853});
        MEAN_RATES.put(RATE_CATEGORY.HIGH, new Double[]{0.11819301147978692, 0.3823224943520077, 0.7695135589733401, 1.7645916247946944} );


        HashMap<String, Object> argParser = parseArgs(args);

        Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
        String ALIGNMENT = getArg(argParser, "ALIGNMENT", String.class);
        String NEWICK = getArg(argParser, "NEWICK", String.class);
        Integer SEED = getArg(argParser, "SEED", Integer.class);
        Integer MODEL_IDX = getArg(argParser, "MODEL_IDX", Integer.class);
        String OUTPUT = getArg(argParser, "OUTPUT", String.class);
        String INPUT = getArg(argParser, "INPUT", String.class);

        System.out.println(ALPHAS);
        System.out.println(ALIGNMENT);
        System.out.println(NEWICK);
        System.out.println(SEED);
        System.out.println(MODEL_IDX);
        System.out.println(OUTPUT);
        System.out.println(INPUT);

        Tree tree = null;
        EnumSeq.Alignment<Enumerable> aln = null;
        try {
            tree = Utils.loadTree(NEWICK);
            aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
            Utils.checkData(aln, tree);
        } catch (ASRException e) {
            usage(6, "Invalid input for ASR: " + e.getMessage());
        } catch (IOException e) {
            usage(2, "Failed to read or write files: " + e.getMessage());
        }
        assert tree != null;
        assert aln != null;

        // This will be updated to a gap model later
        SubstModel MODEL = SubstModel.createModel(MODELS[MODEL_IDX]);
        Double[][] columnPriors = computeColumnPriors(MODEL, tree, aln,
                MEAN_RATES.get(RATE_CATEGORY.HIGH));

        double[] test = {1.0, 0.5, 0.0003};
        double result = MathEx.logsumexp(test);
        System.out.println(result);

    }



    /**
     * Calculates the likelihood of observing each column (independently)
     * given a particular indel rate.
     *
     * @param model Gap augmented substitution model
     * @param tree A Tree object
     * @param aln M columns with N sequences in each col.
     * @return matrix of shape (M, K): Each column, M, and likelihood with rate K
     */
    public static Double[][] computeColumnPriors(SubstModel model, Tree tree, EnumSeq.Alignment<Enumerable> aln,
                                                 Double[] mean_rates) {


        // TODO: make this multithreaded
        Double[][] column_priors = new Double[aln.getWidth()][mean_rates.length];
        System.out.println("Computing column priors");
        for (int col_idx = 0; col_idx < aln.getWidth(); ++col_idx) {
            for (int rate_idx = 0; rate_idx < mean_rates.length; ++rate_idx) {
                double rate = mean_rates[rate_idx];
                column_priors[col_idx][rate_idx] = log_prob_col_given_rate(tree, aln, rate, model, col_idx);
            }
        }

        return column_priors;
    }

    public static double log_prob_col_given_rate(Tree tree, EnumSeq.Alignment<Enumerable> aln, Double rate, SubstModel model,
                                                 int col_idx) {


        int total_nodes = tree.getNLeaves() + tree.getNParents();
        Double[][] Pu_Lk_residue = new Double[total_nodes][model.getDomain().size()]; // nodes x num_letters
        Double[]Pu_Lk_gap = new Double[total_nodes];

        // instantiate with negative infinity for subsequent log sum calculations
        for (int i = 0; i < total_nodes; i++) {
            Arrays.fill(Pu_Lk_residue[i], Double.NEGATIVE_INFINITY);
            Pu_Lk_gap[i] = Double.NEGATIVE_INFINITY;
        }
        // update the Pu_Lk_residue and Pu_Lk_gap arrays in place
        tree.felsensteins_extended_peeling(aln, col_idx, Pu_Lk_residue, Pu_Lk_gap, rate, model);

        return 0.0;
    }

}
