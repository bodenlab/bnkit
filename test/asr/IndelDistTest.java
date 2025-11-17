package asr;

import bn.ctmc.GapSubstModel;
import bn.ctmc.matrix.JCGap;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;
import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class IndelDistTest {


    @Test
    public void testSingleColumnTwoRates() {
        int numCols = 1;
        double[] priors = {0.9, 0.1};  // two rates, rate 0 favored
        Double[][] ll = {{2.0, -1.0}};  // one column,
        double[][] prefix = IndelDist.computePrefixSums(ll);


        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;
        int[][] segs = IndelDist.assignSegments(numCols, priors, prefix, rho);

        assertEquals(1, segs.length);
        assertArrayEquals(new int[]{0,0,0}, segs[0]);  // start=0 end=0 rate=0
    }


    @Test
    public void testTwoColumnsSameRate() {
        double[] priors = {0.1, 0.9};
        Double[][] ll = {
                {-5.0, 5.0},   // col 1
                {-5.0, 5.0}    // col 2
        };
        double[][] prefix = IndelDist.computePrefixSums(ll);


        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;
        int[][] segs = IndelDist.assignSegments(2, priors, prefix, rho);

        // One long segment is optimal: length = 2
        assertEquals(1, segs.length);
        assertArrayEquals(new int[]{0,1,1}, segs[0]); // start=0 end=1 rate=0
    }

    @Test
    public void testTwoColumnsDifferentRates() {
        // rate0=good for col1, rate1=good for col2
        double[] priors = {0.8, 0.2};

        Double[][] ll = {
                {5.0, 0.0},  // col1 strongly prefers rate0
                {0.0, 4.0}   // col2 strongly prefers rate1
        };

        double[][] prefix = IndelDist.computePrefixSums(ll);


        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;
        int[][] segs = IndelDist.assignSegments(2, priors, prefix, rho);

        assertEquals(2, segs.length);

        // First segment: col0-col0 rate0
        assertArrayEquals(new int[]{0,0,0}, segs[0]);
        // Second segment: col1-col1 rate1
        assertArrayEquals(new int[]{1,1,1}, segs[1]);
    }

    @Test
    public void testMaxSegmentLength() {
        double[] priors = {1.0};
        Double[][] ll = {
                {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0},
                {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0},
                {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0},
                {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0},
                {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, {2.0}, // max segment length 50
                {2.0}, {2.0}
        };
        double[][] prefix = IndelDist.computePrefixSums(ll);

        double expected_segment_length = 50.0;
        double rho = 1 / expected_segment_length;

        int[][] segs = IndelDist.assignSegments(52, priors, prefix, rho);

        // Entire segment should be chosen
        assertEquals(2, segs.length);
        assertArrayEquals(new int[]{0,1,0}, segs[0]);
        assertArrayEquals(new int[]{2,51,0}, segs[1]);
    }

    @Test
    public void testThreeColumnsMixedRates() {
        double[] priors = {0.0, 0.0};

        // col1 prefers rate0, col2 prefers rate0, col3 prefers rate1
        Double[][] ll = {
                {3.0, 0.0},
                {3.0, 0.0},
                {0.0, 5.0}
        };
        double[][] prefix = IndelDist.computePrefixSums(ll);


        double expected_segment_length = 20.0;
        double rho = 1 / expected_segment_length;

        int[][] segs = IndelDist.assignSegments(3, priors, prefix, rho);

        assertEquals(2, segs.length);
        assertArrayEquals(new int[]{0, 1, 0}, segs[0]);  // cols 0–1 rate0
        assertArrayEquals(new int[]{2, 2, 1}, segs[1]);  // col 2 rate1
    }


    @Test
    public void testSingleColLikelihoodSingleDeletionCalcJC() {
        Tree tree;
        EnumSeq.Alignment<Enumerable> aln;
        try {
            tree = Utils.loadTree("test/resources/GapSubstModel_test_1.nwk");
            aln = Utils.loadAlignment("test/resources/GapSubstModel_test_1.aln", Enumerable.nacid);
        } catch (ASRException | IOException e) {
            throw new RuntimeException(e);
        }

        GapSubstModel MODEL = new JCGap(0.05, 0.05);
        double geometric_seq_len_param = (double) 1/ 400; // average seq length 400
        double prob = tree.logProbColGivenRate(aln, 1.0, MODEL, 0, geometric_seq_len_param);
        // expected value calculated by hand
        assertEquals(-13.465945445655642, prob, 1e-5);

    }

    @Test
    public void testSingleColLikelihoodCladeDeletionCalcJC() {
        Tree tree;
        EnumSeq.Alignment<Enumerable> aln;
        try {
            tree = Utils.loadTree("test/resources/GapSubstModel_test_1.nwk");
            aln = Utils.loadAlignment("test/resources/GapSubstModel_test_2.aln", Enumerable.nacid);
        } catch (ASRException | IOException e) {
            throw new RuntimeException(e);
        }

        GapSubstModel MODEL = new JCGap(0.05, 0.05);

        double geometric_seq_len_param = (double) 1/ 400; // average seq length 400
        double prob = tree.logProbColGivenRate(aln, 1.0, MODEL, 0, geometric_seq_len_param);
        // expected value calculated by hand
        assertEquals(-12.5800651892681382, prob, 1e-5);

    }
}
