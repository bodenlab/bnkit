package asr;

import bn.ctmc.GapSubstModel;
import bn.ctmc.matrix.JC;
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
import smile.math.Function;
import smile.math.special.Minimise;

public class IndelDist {

    private static final int SEG_START = 0;
    private static final int SEG_END = 1;
    private static final int SEG_RATE = 2;
    private static final int START = 0;
    private static final int RATE_ASSIGNED = 1;
    private static final int DEFAULT_MODEL = 0;
    public static int NTHREADS = 4;
    public static double MIN_MU_LAMBDA_VALUE = 0;
    public static double MAX_MU_LAMBDA_VALUE = 0.5;
    public static String[] MODELS = new String[] {"JTT", "Dayhoff", "LG", "WAG", "Yang", "JC"};
    private static final Enumerable[] ALPHAS = new Enumerable[] {Enumerable.aacid, Enumerable.aacid, Enumerable.aacid,
                                                                 Enumerable.aacid, Enumerable.nacid, Enumerable.nacid};
    private static final double[] RATE_PRIORS = {Math.log(0.25), Math.log(0.25), Math.log(0.25),Math.log(0.25)};


    public enum RATE_CATEGORY {
        LOW,
        MEDIUM,
        HIGH
    }

    private static final Map<RATE_CATEGORY, Double[]> MEAN_RATES =
            Map.of(
                    RATE_CATEGORY.LOW, new Double[]{0.05109101808171236, 0.13936648207205818, 0.2548798168737909, 0.5313293496383602},
                    RATE_CATEGORY.MEDIUM, new Double[]{0.10822622872656855, 0.3236414936793391, 0.6232831770054337, 1.368582433909853},
                    RATE_CATEGORY.HIGH, new Double[]{0.11819301147978692, 0.3823224943520077, 0.7695135589733401, 1.7645916247946944}
            );


    /**
     * Error codes for the program
     */
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

    public static HashMap<String, Object> createArgMap(String[] args) {

        HashMap<String, Object> argMap = new HashMap<>();

        argMap.put("SEED", new Random().nextInt());

        parseArgs(args, argMap);

        checkArgsValid(argMap);

        return argMap;
    }

    private static void parseArgs(String[] args, HashMap<String, Object> argMap) {
        // Read in all the arguments
        for (int a = 0; a < args.length; a++) {
            if (args[a].startsWith("-")) {
                String arg = args[a].substring(1);
                if (((arg.equalsIgnoreCase("-aln")) || (arg.equalsIgnoreCase("a"))) && args.length > a + 1) {
                    argMap.put("ALIGNMENT", args[++a]);
                } else if (((arg.equalsIgnoreCase("-nwk")) || (arg.equalsIgnoreCase("n"))) && args.length > a + 1) {
                    argMap.put("NEWICK", args[++a]);
                } else if ((arg.equalsIgnoreCase("-substitution-model") || arg.equalsIgnoreCase("s")) && args.length > a + 1) {
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
                } else if ((arg.equalsIgnoreCase("-seed") && args.length > a + 1)) {
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
    }

    private static void checkArgsValid(HashMap<String, Object> argMap) {

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
    }


    public static void main(String[] args) throws ExecutionException, InterruptedException {

        HashMap<String, Object> argParser = createArgMap(args);


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
            if (MODEL_IDX == null) {
                MODEL_IDX = DEFAULT_MODEL;
            }

            aln = Utils.loadAlignment(ALIGNMENT, ALPHAS[MODEL_IDX]);
            Utils.checkData(aln, tree);
        } catch (ASRException e) {
            usage(ERROR.ASR.getCode(), ERROR.ASR.getDescription() + e.getMessage());
        } catch (IOException e) {
            usage(ERROR.IO.getCode(), ERROR.IO.getDescription() + e.getMessage());
        }
        assert tree != null;
        assert aln != null;


        if (!MODELS[MODEL_IDX].equals("JTT") || !MODELS[MODEL_IDX].equals("JC")) {
            throw new IllegalArgumentException(MODELS[MODEL_IDX] + " not implemented yet");
        }

        double geometric_seq_len_param = (double) 1 / aln.getAvgSeqLength();
        double optimal_mu = 0.0;
        try {
            long startOptCalc = System.currentTimeMillis();
            optimal_mu = optimiseMuLambda(MIN_MU_LAMBDA_VALUE, MAX_MU_LAMBDA_VALUE, MODEL_IDX, tree, geometric_seq_len_param, aln);
            long endOptCalc = System.currentTimeMillis();
            long OptDuration = endOptCalc - startOptCalc;
            System.out.println("Gap augmented substitution model optimisation time: " + OptDuration / 1000 + " s");

        } catch (IllegalArgumentException e) {
            usage(ERROR.MODEL_AVAIL.getCode(), MODELS[MODEL_IDX] + ERROR.MODEL_AVAIL.getDescription());
        }

        GapSubstModel MODEL;
        if (MODELS[MODEL_IDX].equals("JTT")) {
            MODEL = new JTTGap(optimal_mu, optimal_mu);
        } else {
            throw new IllegalArgumentException(MODELS[MODEL_IDX] + "not implemented yet");
        }

        System.out.println("Computing column priors under different rate categories...");
        long startColCalc = System.currentTimeMillis();
        Double[][] columnPriors = computeColumnPriors(MODEL, tree, aln, MEAN_RATES.get(RATE_CATEGORY.HIGH), geometric_seq_len_param);
        long colCalcDuration = System.currentTimeMillis() - startColCalc;

        System.out.println("Column prior calculation time: " + colCalcDuration / 1000 + " s");


        System.out.println("Computing prefix sums for segment assignment...");
        long startPrefixTime = System.currentTimeMillis();
        double[][] prefix_sums = computePrefixSums(columnPriors);
        long prefixDuration = System.currentTimeMillis() - startPrefixTime;
        System.out.println("Prefix sum calculation time: " + prefixDuration / 1000 + " s");


        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;
        System.out.println("Assigning optimal rate segments...");
        long startSegAssign = System.currentTimeMillis();
        int[][] segments = assignSegments(columnPriors.length, RATE_PRIORS, prefix_sums, rho);
        long segAssignDuration = System.currentTimeMillis() - startSegAssign;
        System.out.println("Segment assignment time: " + segAssignDuration / 1000 + " s");

        long probEnd = System.currentTimeMillis();
        long totalRunTime = progStart - probEnd;

        System.out.println("Optimal indel rate segments (start, end, rate_category):");
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
        Peeler[] peelers = createPeelingJobs(aln, numCols, numRates, mean_rates, model, tree, geometric_seq_len_param);

        Double[][] column_priors = runPeelingJobs(peelers, numCols, numRates);

        return column_priors;
    }

    private static Peeler[] createPeelingJobs(EnumSeq.Alignment<Enumerable> aln, int numCols, int numRates,
                                                Double[] mean_rates, GapSubstModel model, Tree tree,
                                                double geometric_seq_len_param) {

        Peeler[] peelers = new Peeler[numCols * numRates];

        for (int col_idx = 0; col_idx < numCols; ++col_idx) {
            for (int rate_idx = 0; rate_idx < numRates; ++rate_idx) {
                double indel_rate = mean_rates[rate_idx];
                int idx = col_idx * numRates + rate_idx;
                GapSubstModel model_copy = model.deepCopy();
                peelers[idx] = new Peeler(tree, aln, indel_rate, model_copy, col_idx, geometric_seq_len_param);
            }
        }

        return peelers;
    }


    /**
     * Runs the peeling jobs in parallel using a thread pool.
     *
     * @param peelers array of peeling jobs to run
     * @param numCols number of columns
     * @param numRates number of rates
     * @return matrix of shape (numCols, numRates) with the likelihoods of each column under each rate
     */
    private static Double[][] runPeelingJobs(Peeler[] peelers, int numCols, int numRates) {

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
            throw new RuntimeException("Failed to run peeling jobs");
        }

        return column_priors;
    }

    /**
     *
     * @param column_likelihoods array of num_cols x rate
     * @return an array of (num_cols + 1) x k with the likelihood
     * of observing up to that row with rate k
     */
    public static double[][] computePrefixSums(Double[][] column_likelihoods) {

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
     * the algorithm from <a href="https://doi.org/10.1093/sysbio/syx033">Zhai & Alexandre (2017)</a>
     * but column likelihoods are calculated using a gap augmented substitution
     * model.
     *
     * @param num_cols number of alignment columns
     * @param rate_priors the prior probabilities of each rate
     * @param prefix_sums an array of (num_cols + 1) x k with the likelihood of observing up to that row with rate k
     *
     * @return A list of tuples (start, end, rate_category) representing the optimal segmentation of the MSA columns.
     */
    public static int[][] assignSegments(int num_cols, double[] rate_priors,
                                         double[][] prefix_sums, double rho) {

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
                    // longer segments penalised more heavily
                    double segment_length_penalty = calcLogSegmentLengthPenalty(L, rho);
                    double prior_prob_rate = rate_priors[k];
                    // include score from last most likely pos
                    double score = dp_path[i - 1] + LL + prior_prob_rate + segment_length_penalty;
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
        ArrayList<int[]> segments = performBacktrace(num_cols, back_path);

        return reverseSegmentOrder(segments);
    }

    /**
     * Perform backtrace through the dynamic programming table to get the segments.
     * @param num_cols number of columns
     * @param back_path the back path table
     * @return list of segments from backtrace with data organised as (start, end, rate_assigned) zero-indexed
     */
    private static ArrayList<int[]> performBacktrace(int num_cols, int[][] back_path) {

        ArrayList<int[]> segments = new ArrayList<>();
        int j = num_cols;
        while (j > 0) {
            int i = back_path[j][START];
            int k = back_path[j][RATE_ASSIGNED];
            segments.add(new int[] {i, j, k});
            j = i - 1;
        }

        return segments;
    }

    /**
     * Because the backtrace goes from the end to the start, need to
     * reverse the order of the segments.
     * @param segments list of segments from backtrace
     * @return 2D array of segments in correct order
     */
    private static int[][] reverseSegmentOrder(ArrayList<int[]> segments) {

        int[][] optimal_segments = new int[segments.size()][3];

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
     * Want to penalise the model for having segments that are too long.
     * Can use a Geometric distribution. Expected value is 1/rho. i.e.
     * seg_len = 1/rho therefore for a segment length of 20 rho = 0.05.
     *
     * @param segment_length the length of the indel rate segment
     * @param rho the geometric sequence length param for a segment.
     * @return the log penalty for a segment of this length
     */
    public static double calcLogSegmentLengthPenalty(int segment_length, double rho) {
        return (segment_length - 1) * Math.log(1 - rho) + Math.log(rho);
    }

    /**
     * Optimises mu and lambda (insertion and deletion rates) assuming they are equal. Uses Brent's method to find
     * the optimal value that maximises the likelihood of the alignment given the tree. The likelihood is calculated
     * according to equation 29 in <a href="https://doi.org/10.1371/journal.pcbi.1000172"> Rivas & Eddy, 2008</a>
     * @param min_val smallest value to search
     * @param max_val largest value to search
     * @param model_idx index of the substitution model to use
     * @param tree phylogenetic tree
     * @param geometricSeqLenParam geometric sequence length param for the alignment
     * @param aln the alignment
     * @return the optimal mu and lambda value
     * @throws IllegalArgumentException if the model is not supported
     */
    public static double optimiseMuLambda(double min_val, double max_val, int model_idx, Tree tree,
                                          double geometricSeqLenParam,
                                          EnumSeq.Alignment<Enumerable> aln) throws IllegalArgumentException {

        double[] F;
        double[][] IRM;
        Enumerable alpha;
        if (MODELS[model_idx].equals("JTT")) {
            F = JTT.F;
            IRM = JTT.Q;
            alpha = new Enumerable(JTT.S);
        } else if (MODELS[model_idx].equals("JC")) {
            F = JC.F(JC.S.length);
            IRM = JC.Q(1, JC.S.length);
            alpha = new Enumerable(JC.S);
        } else {
            throw new IllegalArgumentException(MODELS[model_idx] + "not implemented yet");
        }

        AlnLikelihood alnLikelihood = new AlnLikelihood(tree, aln, F, IRM, alpha, geometricSeqLenParam);

        return Minimise.brent(alnLikelihood, min_val, max_val);
    }

    /**
     * Function to calculate the likelihood of an alignment given a tree and gap augmented substitution model.
     */
    public static class AlnLikelihood implements Function {

        double[] F;
        double[][] IRM;
        Enumerable alpha;
        Tree tree;
        EnumSeq.Alignment<Enumerable> aln;
        double geometricSeqLenParam;

        public AlnLikelihood(Tree tree, EnumSeq.Alignment<Enumerable> aln, double[] F, double[][] IRM,
                             Enumerable alpha, double geometricSeqLenParam) {
            this.tree = tree;
            this.aln = aln;
            this.geometricSeqLenParam = geometricSeqLenParam;
            this.F = F;
            this.IRM = IRM;
            this.alpha = alpha;
        }

        /**
         * Get the log likelihood of a particular alignment given the tree and gap augmented substitution model.
         * This is used by a minimisation routine to find optimal params for the substitution model.
         * @param muLambda Assumes mu (deletion rate) and lambda (insertion rate) are equal
         * @return likelihood of the alignment given the tree + mu + lambda.
         */
        @Override
        public double f(double muLambda) {

            GapSubstModel newModel = new GapSubstModel(this.F, this.IRM, this.alpha, muLambda, muLambda);

            // trying to maximise the log likelihood
            return -1 * ThreadedPeeler.calcProbAlnGivenTree(tree, newModel, aln, geometricSeqLenParam,
                                                            alpha, NTHREADS);
        }
    }
}
