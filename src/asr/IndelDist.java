package asr;


import java.util.*;

public class IndelDist {

    private static final int SEG_START = 0;
    private static final int SEG_END = 1;
    private static final int SEG_RATE = 2;
    private static final int START = 0;
    private static final int RATE_ASSIGNED = 1;
    private static final double EXPECTED_INDEL_SEGMENT_LENGTH = 20.0;
    public static final double RHO = 1 / EXPECTED_INDEL_SEGMENT_LENGTH;
    public static final double[] RATE_PRIORS = {Math.log(0.25), Math.log(0.25), Math.log(0.25),Math.log(0.25)};


    public enum RATE_CATEGORY {
        LOW,
        MEDIUM,
        HIGH
    }

    public static final Map<RATE_CATEGORY, double[]> MEAN_RATES =
            Map.of(
                    RATE_CATEGORY.LOW, new double[]{0.05109101808171236, 0.13936648207205818, 0.2548798168737909, 0.5313293496383602},
                    RATE_CATEGORY.MEDIUM, new double[]{0.10822622872656855, 0.3236414936793391, 0.6232831770054337, 1.368582433909853},
                    RATE_CATEGORY.HIGH, new double[]{0.11819301147978692, 0.3823224943520077, 0.7695135589733401, 1.7645916247946944}
            );

    /**
     *
     * @param column_likelihoods array of rate x num_cols
     * @return an array of (num_cols + 1) x k with the likelihood
     * of observing up to that row with rate k
     */
    public static double[][] computePrefixSums(double[][] column_likelihoods) {

        int numberOfCols = column_likelihoods.length + 1;
        int numRates = column_likelihoods[0].length;

        double[][] partial_sums = new double[numberOfCols][numRates];
        for (int k = 0; k < numRates; k++) {
            double total = 0.0;
            for (int i = 1; i < numberOfCols; i++)  {
                total += column_likelihoods[k][i - 1];
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
     * Expand the segment order into an array of length equal to the number of columns
     * where each entry is the indel rate category assigned to that column.
     * @param segments list of segments
     * @return array of length equal to number of columns with indel rate category for each column
     */
    public static int[] expandSegmentOrder(int[][] segments) {

        int numCols = segments[segments.length - 1][SEG_END] + 1;
        int[] colRates = new int[numCols];
        for (int[] segment: segments) {
            int start = segment[SEG_START];
            int end = segment[SEG_END];
            int rate = segment[SEG_RATE];
            for (int col = start; col <= end; col++) {
                colRates[col] = rate;
            }
        }

        return colRates;
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

}
