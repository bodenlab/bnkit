package asr;

import bn.ctmc.GapSubstModel;
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


        double geometric_seq_len_param = (double) 1 / aln.getAvgSeqLength();
        double[] rate_priors = {Math.log(0.25), Math.log(0.25), Math.log(0.25),Math.log(0.25)};
        double mu = 0.05;
        double lambda = 0.05;
        // This will be updated to a gap model later
        GapSubstModel MODEL;
        if (MODELS[MODEL_IDX].equals("JTT")) {
            MODEL = new JTTGap(mu, lambda);
        } else {
            throw new IllegalArgumentException(MODELS[MODEL_IDX] + "not implemented yet");
        }



        long startTime = System.currentTimeMillis();
        Double[][] columnPriors = computeColumnPriors(MODEL, tree, aln,
                MEAN_RATES.get(RATE_CATEGORY.HIGH), geometric_seq_len_param);

        long endTime = System.currentTimeMillis();
        long duration = endTime - startTime;
        System.out.println("Execution time: " + duration / 1000 + " s");

        double[][] prefix_sums = compute_prefix_sums(columnPriors);

        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;
        int[][] segments = assign_segments(columnPriors.length, rate_priors, prefix_sums, rho);
        System.out.println(Arrays.deepToString(segments));


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
        int nThreads = 4;
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
        ThreadedPeeler thread_pool = new ThreadedPeeler(peelers, nThreads);
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

}
