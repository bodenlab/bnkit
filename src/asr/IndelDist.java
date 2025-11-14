package asr;

import bn.ctmc.GapSubstModel;
import bn.ctmc.matrix.JTT;
import bn.ctmc.matrix.JTTGap;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ExecutionException;
import asr.ThreadedPeeler.Peeler;
import smile.math.phylo.AlnLikelihood;
import smile.math.special.Minimise;

public class IndelDist {

    public enum RATE_CATEGORY {
        LOW,
        MEDIUM,
        HIGH
    }

    public enum ERROR {
        SUCCESS(0, null),
        SUB_MODEL(1, "Could not find model with ID "),
        INPUT(2, "Input must be a string"),
        ALN(3, "Must specify alignment (--aln <Clustal or FASTA file>)"),
        NWK( 4, "Must specify phylogenetic tree (Newick file) or previously saved folder (--input-folder <folder>)"),
        UNKNOWN(5, "Unknown option or missing required argument: "),
        MAX_MIN_MU(6, "min-rate and max-rate must be >= 0 and min-rate must be < max-rate"),
        ASR(7, "Invalid input for ASR: "),
        IO(8, "Failed to read or write files: "),
        MODEL_AVAIL(9, " not supported for gap augmentation");



        private final int code;
        private final String description;

        ERROR(int code, String description) {
            this.code = code;
            this.description = description;
        }

        public String getDescription() {
            return description;
        }

        public int getCode() {
            return code;
        }

    }

    public static EnumMap<RATE_CATEGORY, Double[]> MEAN_RATES = new EnumMap<>(RATE_CATEGORY.class);
    public static String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
    public static int NTHREADS = 4;
    public static double MIN_MU_LAMBDA_VALUE = 0;
    public static double MAX_MU_LAMBDA_VALUE = 0;

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
                \t{-t | --threads <number>}
                \t{-mnr | --min-rate <double>} (The lower search bound when determining optimal insertions/deletion rate in the substitution model)
                \t{-mxr | --max-rate <double>} (The upper search bound when determining optimal insertions/deletion rate in the substitution model)
                \t-h (or --help) will print out this screen
                """
        );

        if (msg != null) {
            out.println("\n" + msg + "(Error code " + error_code + ")");
        }
        System.exit(error_code);
    }

    public static void usage() {
        usage(ERROR.SUCCESS.getCode(), ERROR.SUCCESS.getDescription());
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
                        usage(ERROR.SUB_MODEL.getCode(), ERROR.SUB_MODEL.getDescription() + args[a + 1]);
                    }
                } else if ((arg.equalsIgnoreCase("-seed") && args.length > a+1)) {
                    argMap.put("SEED", Integer.parseInt(args[++a]));
                } else if ((arg.equalsIgnoreCase("-output-folder") || arg.equalsIgnoreCase("o")) && args.length > a + 1) {
                    argMap.put("OUTPUT", args[++a]);
                } else if ((arg.equalsIgnoreCase("-input-folder") || arg.equalsIgnoreCase("i")) && args.length > a + 1) {
                    argMap.put("INPUT", args[++a]);
                } else if ((arg.equalsIgnoreCase("-threads") || arg.equalsIgnoreCase("t")) && args.length > a + 1) {
                    try {
                        NTHREADS = Integer.parseInt(args[++a]);
                    } catch (NumberFormatException e) {
                        System.err.println("Failed to set number of threads for option --threads: " + args[a] + " is not a valid integer");
                    }
                } else if (arg.equalsIgnoreCase("-help") || arg.equalsIgnoreCase("h")) {
                    usage();
                } else if ((arg.equalsIgnoreCase("-min-rate") || arg.equalsIgnoreCase("mnr") && args.length > a + 1)) {
                    try {
                        MIN_MU_LAMBDA_VALUE = Double.parseDouble(args[++a]);
                    } catch (NumberFormatException e) {
                        System.err.println("Failed to set minimum Mu/Lambda value for option --min-rate: " + args[a] + " is not a valid double");
                    }
                } else if ((arg.equalsIgnoreCase("-max-rate") || arg.equalsIgnoreCase("mxr") && args.length > a + 1)) {
                    try {
                        MAX_MU_LAMBDA_VALUE = Double.parseDouble(args[++a]);
                    } catch (NumberFormatException e) {
                        System.err.println("Failed to set maximum Mu/Lambda value for option --min-rate: " + args[a] + " is not a valid double");
                    }
                } else {
                    usage(ERROR.UNKNOWN.getCode(), ERROR.UNKNOWN.getDescription() + "\" + args[a] + \"");
            }

            }
        }

        // now perform checks
        if (argMap.get("ALIGNMENT") == null) {
            usage(ERROR.ALN.getCode(), ERROR.ALN.getDescription());
        } else if (argMap.get("NEWICK") == null) {
            usage(ERROR.NWK.getCode(), ERROR.NWK.getDescription());
        } else if (argMap.get("OUTPUT") == null) {
            String INPUT = null;
            try {
                INPUT = (String) argMap.get("INPUT");
            } catch (ClassCastException e) {
                usage(ERROR.INPUT.getCode(), ERROR.INPUT.getDescription());
            }

            argMap.put("OUTPUT", INPUT == null ? "." : INPUT);
        } else if (MAX_MU_LAMBDA_VALUE <= MIN_MU_LAMBDA_VALUE || MIN_MU_LAMBDA_VALUE < 0.0) {
            usage(ERROR.MAX_MIN_MU.getCode(), ERROR.MAX_MIN_MU.getDescription());
        }

        return argMap;
    }


    public static void main(String[] args) {
        MEAN_RATES.put(RATE_CATEGORY.LOW, new Double[]{0.05109101808171236, 0.13936648207205818, 0.2548798168737909, 0.5313293496383602});
        MEAN_RATES.put(RATE_CATEGORY.MEDIUM, new Double[]{0.10822622872656855, 0.3236414936793391, 0.6232831770054337, 1.368582433909853});
        MEAN_RATES.put(RATE_CATEGORY.HIGH, new Double[]{0.11819301147978692, 0.3823224943520077, 0.7695135589733401, 1.7645916247946944} );
        double[] RATE_PRIORS = {Math.log(0.25), Math.log(0.25), Math.log(0.25),Math.log(0.25)};

        HashMap<String, Object> argParser = parseArgs(args);

        Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid, Enumerable.aacid,
                                                Enumerable.nacid, Enumerable.nacid};
        String ALIGNMENT = getArg(argParser, "ALIGNMENT", String.class);
        String NEWICK = getArg(argParser, "NEWICK", String.class);
        Integer SEED = getArg(argParser, "SEED", Integer.class);
        Integer MODEL_IDX = getArg(argParser, "MODEL_IDX", Integer.class);
        String OUTPUT = getArg(argParser, "OUTPUT", String.class);
        String INPUT = getArg(argParser, "INPUT", String.class);


        long progStart = System.currentTimeMillis();
        Tree tree = null;
        EnumSeq.Alignment<Enumerable> aln = null;
        try {
            tree = Utils.loadTree(NEWICK);
            aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
            Utils.checkData(aln, tree);
        } catch (ASRException e) {
            usage(ERROR.ASR.getCode(), ERROR.ASR.getDescription() + e.getMessage());
        } catch (IOException e) {
            usage(ERROR.IO.getCode(), ERROR.IO.getDescription() + e.getMessage());
        }
        assert tree != null;
        assert aln != null;

        double geometric_seq_len_param = (double) 1 / aln.getAvgSeqLength();

        // This will be updated to a gap model later
        if (!MODELS[MODEL_IDX].equals("JTT")) {
            throw new IllegalArgumentException(MODELS[MODEL_IDX] + " not implemented yet");
        }

        Double optimal_mu = null;
        try {
            long startOptCalc = System.currentTimeMillis();
            optimal_mu = optimiseMuLambda(MIN_MU_LAMBDA_VALUE, MAX_MU_LAMBDA_VALUE, MODEL_IDX, tree, geometric_seq_len_param, aln);
            long endOptCalc = System.currentTimeMillis();
            long OptDuration = endOptCalc - startOptCalc;
            System.out.println("Model optimisation time: " + OptDuration / 1000 + " s");

        } catch (IllegalArgumentException e) {
            usage(ERROR.MODEL_AVAIL.getCode(), MODELS[MODEL_IDX] + ERROR.MODEL_AVAIL.getDescription());
        }

        GapSubstModel MODEL;
        if (MODELS[MODEL_IDX].equals("JTT")) {
            MODEL = new JTTGap(optimal_mu, optimal_mu);
        } else {
            throw new IllegalArgumentException(MODELS[MODEL_IDX] + "not implemented yet");
        }

        long startColCalc = System.currentTimeMillis();
        Double[][] columnPriors = computeColumnPriors(MODEL, tree, aln, MEAN_RATES.get(RATE_CATEGORY.HIGH), geometric_seq_len_param);
        long endColCalc = System.currentTimeMillis();
        long colCalcDuration = endColCalc - startColCalc;
        System.out.println("Column prior calculation time: " + colCalcDuration / 1000 + " s");

        double[][] prefix_sums = compute_prefix_sums(columnPriors);

        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;
        int[][] segments = assign_segments(columnPriors.length, RATE_PRIORS, prefix_sums, rho);

        long probEnd = System.currentTimeMillis();
        long totalRunTime = progStart - probEnd;

        System.out.println(Arrays.deepToString(segments));
        System.out.println("Total program run time: " + totalRunTime / 1000 + " s");

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
    public static Double[][] computeColumnPriors(GapSubstModel model, Tree tree, EnumSeq.Alignment<Enumerable> aln,
                                                 Double[] mean_rates, double geometric_seq_len_param) {


        int numCols = aln.getWidth();
        int numRates = mean_rates.length;
        int nThreads = IndelDist.NTHREADS;
        Peeler[] peelers = new Peeler[numCols * numRates];

        System.out.println("Computing column priors");
        for (int col_idx = 0; col_idx < numCols; ++col_idx) {
            for (int rate_idx = 0; rate_idx < numRates; ++rate_idx) {
                double indel_rate = mean_rates[rate_idx];
                int idx = col_idx * numRates + rate_idx;
                GapSubstModel model_copy = model.deepCopy();
                peelers[idx] = new Peeler(tree, aln, indel_rate, model_copy, col_idx, geometric_seq_len_param);
            }
        }

        Double[][] column_priors = new Double[numCols][numRates];
        ThreadedPeeler thread_pool = new ThreadedPeeler(peelers, IndelDist.NTHREADS);
        try {
            Map<Integer, Peeler> ret = thread_pool.runBatch();
            for (int col_idx = 0; col_idx < numCols; ++col_idx) {
                for (int rate_idx = 0; rate_idx < numRates; ++rate_idx) {
                    int idx = col_idx * numRates + rate_idx;
                    Peeler peeler = ret.get(idx);
                    column_priors[col_idx][rate_idx] = peeler.getDecoration();
                }
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        return column_priors;
    }

    /**
     *
     * @param column_likelihoods array of num_cols x rate
     * @return an array of (num_cols + 1) x k with the likelihood
     * of observing up to that row with rate k
     */
    public static double[][] compute_prefix_sums(Double[][] column_likelihoods) {

        int numberOfCols = column_likelihoods.length + 1;
        int numRates = column_likelihoods[0].length;

        double[][] partial_sums = new double[numberOfCols][numRates];
        for (int k = 0; k < numRates; k++) {
            double total = 0.0;
            for (int i = 1; i < numberOfCols; i++)  {
                total += column_likelihoods[i - 1][k];
                partial_sums[i][k] = total;
            }
        }

        return partial_sums;

    }

    /**
     * Uses dynamic programming to assign segments to the MSA columns. This adapts
     * the algorithm from Zhai & Alexandre (2017) <a href="https://doi.org/10.1093/sysbio/syx033">...</a>
     * but column likelihoods are calculated using a gap augmented substitution
     * model.
     *
     * @param num_cols number of alignment columns
     * @param rate_priors the prior probabilities of each rate
     * @param prefix_sums an array of (num_cols + 1) x k with the likelihood of observing up to that row with rate k
     *
     * @return A list of tuples (start, end, rate_category) representing the optimal segmentation of the MSA columns.
     */
    public static int[][] assign_segments(int num_cols, double[] rate_priors,
                                  double[][] prefix_sums, double rho) {


        int START = 0;
        int RATE_ASSIGNED = 1;
        Double[] dp_path = new Double[num_cols + 1];
        Arrays.fill(dp_path, Double.NEGATIVE_INFINITY);
        dp_path[0] = 0.0;
        // each col will record the most optimal start and end position
        int[][] back_path = new int[num_cols + 1][2];
        int K = rate_priors.length;

        int max_seg_len = 50;
        for (int j = 1; j < num_cols + 1; j++) {
            // only look back as far as the maximum segment length
            int i_min = Math.max(1, j - max_seg_len + 1);
            double best_score = Double.NEGATIVE_INFINITY;
            int[] best_entry = new int[2];
            for (int i = i_min; i < (j + 1); i++) {
                int L = j - i + 1; // segment length
                for (int k = 0; k < K; k++) {
                    double LL = prefix_sums[j][k] - prefix_sums[i - 1][k]; // prob of this segment at this indel rate
                    double len_weight = length_prior(L, rho); // longer segments penalised more heavily
                    double prior = rate_priors[k]; // prior prob of that particular weight category
                    double score = dp_path[i - 1] + LL + prior + len_weight; // include score from last most likely pos
                    if (score > best_score) {
                        best_score = score;
                        best_entry[START] = i;
                        best_entry[RATE_ASSIGNED] = k;
                    }
                }
            }

            dp_path[j] = best_score;
            back_path[j][START] = best_entry[START];
            back_path[j][RATE_ASSIGNED] = best_entry[RATE_ASSIGNED];
        }

        // perform the backtrace through segments
        ArrayList<int[]> segments = new ArrayList<>();
        int j = num_cols;
        while (j > 0) {
            int i = back_path[j][START];
            int k = back_path[j][RATE_ASSIGNED];
            segments.add(new int[] {i, j, k});
            j = i - 1;
        }

        // need to reverse the order
        int[][] optimal_segments = new int[segments.size()][3];
        int SEG_START = 0;
        int SEG_END = 1;
        int SEG_RATE = 2;
        int pos = 0;
        for (int i = segments.size() - 1; i >= 0; i--) {
            int[] segment = segments.get(i);
            // note the start and ends are inclusive and zero-indexed
            optimal_segments[pos][SEG_START] = segment[SEG_START] -1;
            optimal_segments[pos][SEG_END] = segment[SEG_END] - 1;
            optimal_segments[pos][SEG_RATE] = segment[SEG_RATE];
            pos++;
        }

        return optimal_segments;
    }

    /**
     * Want to penalise the model for having segements that are too long.
     * Can use a Geometric distribution. Expected value is 1/rho. i.e.
     * seg_len = 1/rho therefore for a segment length of 20 rho = 0.05.
     *
     * @param segment_length the length of the indel rate segment
     * @param rho the geometric sequence length param for a segment.
     * @return
     */
    public static double length_prior(int segment_length, double rho) {
        return (segment_length - 1) * Math.log(1 - rho) + Math.log(rho);
    }

    public static double optimiseMuLambda(double min_val, double max_val, int model_idx, Tree tree,
                                          double geometric_seq_len_param,
                                          EnumSeq.Alignment<Enumerable> aln) throws IllegalArgumentException {

        double[] F;
        double[][] IRM;
        Enumerable alpha;
        if (MODELS[model_idx].equals("JTT")) {
            F = JTT.F;
            IRM = JTT.Q;
            alpha = new Enumerable(JTT.S);
        } else {
            throw new IllegalArgumentException(MODELS[model_idx] + "not implemented yet");
        }

        AlnLikelihood alnLikelihood = new AlnLikelihood(tree, aln, geometric_seq_len_param, alpha, F, IRM);

        return Minimise.brent(alnLikelihood, min_val, max_val);
    }
}
