package dat;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * Created by mikael on 13/10/16.
 */
public class IntervalSTTest  {
    private IntervalST<String> st0, st1, st1checkme, st1q, st2, st3, st4, st5, st6;
    private Set<Interval1D> ovlaps = new HashSet<>();
    private Set<Interval1D> queries = new HashSet<>();

    @Test
    public void get() throws Exception {
        Interval1D ival4 = new Interval1D(ivals2[4][0], ivals2[4][1]);
        Interval1D ival7 = new Interval1D(ivals2[7][0], ivals2[7][1]);
        Interval1D ival10 = new Interval1D(ivals2[10][0], ivals2[10][1]);
        Set<String> vals4 = st2.get(ival4);
        assertTrue(vals4.size() == 2);
        Set<String> vals7 = st2.get(ival7);
        assertTrue(vals7.size() == 2);
        Set<String> vals10 = st2.get(ival10);
        assertTrue(vals10.size() == 1);
    }

    @Test
    public void contains() throws Exception {
        for (int i = 0; i < ivals1.length; i ++)
            assertTrue(st1.contains(new Interval1D(ivals1[i][0], ivals1[i][1])));
        for (int i = 0; i < ivals2.length; i ++) {
            assertTrue(st2.contains(new Interval1D(ivals2[i][0], ivals2[i][1])));
            assertTrue(!st1.contains(new Interval1D(ivals2[i][0], ivals2[i][1])));
        }
        for (int i = 0; i < ivals3.length; i ++) {
            assertTrue(st3.contains(new Interval1D(ivals3[i][0], ivals3[i][1])));
            assertTrue(!st1.contains(new Interval1D(ivals3[i][0], ivals3[i][1])));
        }
    }

    @Test
    public void search() throws Exception {
        // search for an interval intersecting with a given interval
        for (int i = 0; i < ivals4.length; i ++) {
            Interval1D ival1 = new Interval1D(ivals4[i][0], ivals4[i][1]);
            //System.out.println("Searching for " + ival);
            Interval1D found = st4.search(ival1);
            //System.out.println("Finding " + found);
            assertTrue(ival1.equals(found));
        }
        assertTrue(st1.search(new Interval1D(1,4)) == null);
        assertTrue(st2.search(new Interval1D(29,33)) == null);
        assertTrue(st2.search(new Interval1D(800,80000000)) == null);
        assertTrue(st2.search(new Interval1D(231,238)) == null);
    }

    @Test
    public void search1() throws Exception {
        // search for an interval containing a given point
        for (int i = 0; i < ivals1.length; i ++) {
            for (int j = ivals1[i][0]; j < ivals1[i][1]; j ++)
                assertTrue(st1.search(j) != null);
        }
        assertTrue(st1.search(4) == null);
        assertTrue(st1.search(31) == null);
        assertTrue(st1.search(100) == null);
    }

    @Test
    public void searchAll() throws Exception {
        // find ALL intervals that overlap/intersect with that given
        Set<Interval1D> ref3 = new HashSet<>();
        for (int i = 0; i < ivals3.length; i++) {
            ref3.add(new Interval1D(ivals3[i][0], ivals3[i][1]));
        }
        for (Interval1D ival : ref3) {
            Set<Interval1D> all = st3.searchAll(ival);
            assertTrue(all != null);
            assertTrue(all.size() == 5);
        }
        for (int i = 0; i < ivals4.length; i++) {
            assertTrue(st4.searchAll(new Interval1D(ivals4[i][0], ivals4[i][1])).size() == 1);
            assertTrue(st5.searchAll(new Interval1D(ivals4[i][0], ivals4[i][1])).size() == 0);
        }
        System.out.println(ovlaps.size() + " contains-tests");
        for (Interval1D ov : ovlaps) {
            assertTrue(st6.search(ov) != null);
        }
    }

    @Test
    public void getClosest() throws Exception {
        Interval1D ival2030 = new Interval1D(20, 30);
        Interval1D ival3597 = new Interval1D(35, 97);
        Interval1D ival3234 = new Interval1D(32, 34);
        Interval1D ival_close1 = st1.getClosest(ival3234);
        assertTrue(ival_close1.equals(ival3597));
        Interval1D ival398_400 = new Interval1D(398, 400);
        Interval1D ival399_500 = new Interval1D(399, 500);
        Interval1D ival_close3 = st3.getClosest(ival398_400);
        assertTrue(ival_close3.equals(ival399_500));
        Interval1D ival498_500 = new Interval1D(498, 500);
        Interval1D ival498_499 = new Interval1D(498, 499);
        Interval1D ival500_600 = new Interval1D(500, 600);
        Interval1D ival602_605 = new Interval1D(602, 605);
        Interval1D ival_close4 = st3.getClosest(ival498_500);
        assertTrue(ival_close4.equals(ival498_499));
        Interval1D ival_close5 = st4.getClosest(ival498_500);
        assertTrue(ival_close5.equals(ival500_600));
        Interval1D ival_close6 = st5.getClosest(ival498_500);
        assertTrue(ival_close6.equals(ival602_605));
    }

    @Test
    public void getClosestOnSide() throws Exception {
        Interval1D ival2030 = new Interval1D(20, 30);
        Interval1D ival3597 = new Interval1D(35, 97);
        Interval1D ival3234 = new Interval1D(32, 34);
        Interval1D ival_close1a = st1.getClosest(ival3234, IntervalST.CLOSEST_BEFORE);
        Interval1D ans = getClosest(st1, ival3234, IntervalST.CLOSEST_BEFORE);
        assertTrue(ival_close1a.equals(ival2030));
        assertTrue(ival_close1a.equals(ans));
        Interval1D ival_close1b = st1.getClosest(ival3234, IntervalST.CLOSEST_AFTER);
        ans = getClosest(st1, ival3234, IntervalST.CLOSEST_AFTER);
        assertTrue(ival_close1b.equals(ival3597));
        assertTrue(ival_close1b.equals(ans));
        int cnt = 0;
        for (Interval1D i : st1q) {
            Interval1D before = st1checkme.getClosest(i, IntervalST.CLOSEST_BEFORE);
            ans = getClosest(st1checkme, i, IntervalST.CLOSEST_BEFORE);
            if (before == null) {
                assertTrue(before == null);
                cnt++;
            } else {
                assertTrue(before.signdist(i, true) >= 0);
                if (ans != null)
                    assertTrue(before.equals(ans));
            }
        }
        for (Interval1D i : st1q) {
            Interval1D after = st1checkme.getClosest(i, IntervalST.CLOSEST_AFTER);
            ans = getClosest(st1checkme, i, IntervalST.CLOSEST_AFTER);
            if (after == null) {
                assertTrue(after == null);
                cnt++;
            } else {
                assertTrue(after.signdist(i, true) <= 0);
                if (ans != null)
                    assertTrue(after.equals(ans));
            }
        }
    }

    public Interval1D getClosest(IntervalST tree, Interval1D query, int DIRECTION) {
        Set<Interval1D> all = tree.getAll();
        int best_dist = Integer.MAX_VALUE;
        Interval1D fav = null;
        boolean before = false;
        boolean after = false;
        for (Interval1D ival : all) {
            int dist = query.signdist(ival, true);
            if (dist < 0) {
                before = true; // ival is BEFORE query
                after = false; // ival is AFTER query
                dist = -dist;
            } else if (dist > 0) {
                before = false;
                after = true;
            } else { // dist == 0
                return null;
            }
            if (DIRECTION == IntervalST.CLOSEST_BEFORE_OR_AFTER || (DIRECTION == IntervalST.CLOSEST_BEFORE && before) || (DIRECTION == IntervalST.CLOSEST_AFTER && after)) {
                if (fav == null || dist < best_dist) {
                    best_dist = dist;
                    fav = ival;
                }
            }
        }
        return fav;
    }

    @Test
    public void getAll() throws Exception {
        Set<Interval1D> all1 = st1.getAll();
        assertTrue(all1.size() == ivals1.length);
        Set<Interval1D> all2 = st2.getAll();
        assertTrue(all2.size() == ivals2.length - 3); // 3 duplicates
    }

    @Test
    public void iterator() throws Exception {
        int prev = Integer.MIN_VALUE;
        for (Interval1D ival : st6) {
            assertTrue(ival.min >= prev);
            prev = ival.min;
        }
    }

    @Test
    public void flatten() throws Exception {
        Set<Interval1D> ivals = st1.flatten2Set();
        assertTrue(ivals.size() == 3);
        IntervalST<String> tst1 = st1.flatten2Tree();
        assertTrue(tst1.size() == 3);
        ivals = st0.flatten2Set(true);
        assertTrue(ivals.size() == 1);
        ivals = st0.flatten2Set(false);
        assertTrue(ivals.size() == 4);
        ivals = st3.flatten2Set();
        assertTrue(ivals.size() == 2);
        IntervalST<String> tst3 = st3.flatten2Tree();
        assertTrue(tst3.size() == 2);
/*        for (Interval1D ival : tst3) {
            System.out.println(ival);
            for (String v : tst3.get(ival))
                System.out.println("\t" + v);
        } */
    }

    // check integrity of count and interval fields
    @Test
    public void check() {
        assertTrue(st1.checkCount() && st1.checkMax());
        assertTrue(st2.checkCount() && st2.checkMax());
        assertTrue(st3.checkCount() && st3.checkMax());
        assertTrue(st4.checkCount() && st4.checkMax());
        assertTrue(st5.checkCount() && st5.checkMax());
    }

    @Test
    public void random() {
        for (Interval1D ival : ovlaps) {
            Interval1D match = st6.search(ival);
            assertTrue(match != null);
        }
        int nq = 0;
        for (Interval1D query : queries) {
            Interval1D tst_1 = st6.getClosest(query, IntervalST.CLOSEST_BEFORE_OR_AFTER);
            if (tst_1.dist(query, true) == 0)
                continue;
            nq ++;
            Interval1D ans_1 = getClosest(st6, query, IntervalST.CLOSEST_BEFORE_OR_AFTER);
            if (query.dist(tst_1, true) != query.dist(ans_1, true))
                System.out.println("Before or after: For query interval " + query + ": " + tst_1 + "(" + query.dist(tst_1, true) + ") is not equal to " + ans_1  + "(" + query.dist(ans_1, true) + ")");
            assertTrue(query.dist(tst_1, true) == query.dist(ans_1, true));
            Interval1D tst_2 = st6.getClosest(query, IntervalST.CLOSEST_BEFORE);
            Interval1D ans_2 = getClosest(st6, query, IntervalST.CLOSEST_BEFORE);
            if (tst_2.max != ans_2.max)
                System.out.println("Before: For query interval " + query + ": " + tst_2 + " is not equal to " + ans_2);
            assertTrue(tst_2.max == ans_2.max);
            Interval1D tst_3 = st6.getClosest(query, IntervalST.CLOSEST_AFTER);
            Interval1D ans_3 = getClosest(st6, query, IntervalST.CLOSEST_AFTER);
            if (tst_3.min != ans_3.min)
                System.out.println("After: For query interval " + query + ": " + tst_3 + " is not equal to " + ans_3);
            assertTrue(tst_3.min == ans_3.min);
        }
        System.out.println(nq + " closest-tests");
    }

    int[][] ivals0 = new int[][] {
            { 5, 21},
            {21, 50},
            {50, 99},
            {99, 100}};
    int[][] ivals1 = new int[][] {
            {10, 15},
            {11, 17},
            {20, 30},
            {21, 24},
            { 5, 21},
            {40, 50},
            {98, 99},
            {35, 97}};
    int[][] ivals1checkme = new int[][] {
            {4, 5},
            {6, 7},
            {8, 12},
            {10, 13},
            {15, 17},
            {21, 22},
            {35, 39},
            {41, 42},
            {109, 115}};
    int[][] ivals1q = new int[][] {
            {1, 3},
            {8, 9},
            {18, 19},
            {25, 27},
            {15, 16},
            {31, 32},
            {33, 34},
            {43, 44},
            {105, 106}};
    int[][] ivals2 = new int[][] {
            {  9,  15},
            {111, 117},
            {120, 130},
            {221, 224},
            {  6,  21},
            { 16,  17},
            {220, 230},
            {221, 224},
            {  6,  21},
            { 40,  51},
            { 98, 199},
            { 98, 199},
            { 35,  96}};
    String[] lab2 = new String[] {
            "noname",
            "noname",
            "noname",
            "rep1",
            "rep1",
            "noname",
            "noname",
            "rep2",
            "rep2",
            "noname",
            "rep1",
            "rep1",
            "noname"};
    int[][] ivals3 = new int[][] {
            {1000, 5000},
            {1100, 4700},
            {2000, 3000},
            {2100, 2400},
            { 503, 2100},
            { 400, 502},
            { 399, 501},
            { 399, 500},
            { 498, 499},
            { 498, 498}};
    int[][] ivals4 = new int[][]{ // all unique, non-intersecting
            {1, 2}, {30, 40}, {500, 600}, {7001, 8000}, {607, 6998}};
    int[][] ivals5 = new int[][]{ // all unique, non-intersecting, also not intersection with ivals4
            {0, 0}, {3, 29}, {6999,7000}, {602, 605}};

    @BeforeEach
    public void SetupTest() {
        st0 = new IntervalST<>();
        for (int i = 0; i < ivals0.length; i ++)
            st0.put(new Interval1D(ivals0[i][0], ivals0[i][1]), "I0_"+(i+1));
        st1 = new IntervalST<>();
        for (int i = 0; i < ivals1.length; i ++)
            st1.put(new Interval1D(ivals1[i][0], ivals1[i][1]), "I1_"+(i+1));
        st1checkme = new IntervalST<>(1);
        for (int i = 0; i < ivals1checkme.length; i ++)
            st1checkme.put(new Interval1D(ivals1checkme[i][0], ivals1checkme[i][1]), "I1c_"+(i+1));
        st1q = new IntervalST<>(1);
        for (int i = 0; i < ivals1q.length; i ++)
            st1q.put(new Interval1D(ivals1q[i][0], ivals1q[i][1]), "I1q_"+(i+1));
        st2 = new IntervalST<>(1);
        for (int i = 0; i < ivals2.length; i ++)
            st2.put(new Interval1D(ivals2[i][0], ivals2[i][1]), lab2[i]);;
        st3 = new IntervalST<>(1);
        for (int i = 0; i < ivals3.length; i ++)
            st3.put(new Interval1D(ivals3[i][0], ivals3[i][1]), "I3_"+(i+1));
        st4 = new IntervalST<>(1);
        for (int i = 0; i < ivals4.length; i ++)
            st4.put(new Interval1D(ivals4[i][0], ivals4[i][1]), "I4_"+(i+1));
        st5 = new IntervalST<>(1);
        for (int i = 0; i < ivals5.length; i ++)
            st5.put(new Interval1D(ivals5[i][0], ivals5[i][1]), "I5_"+(i+1));

        int N = 2000;
        Random rand = new Random(1);
        // generate N random intervals and insert into data structure
        st6 = new IntervalST<String>(1);
        for (int i = 0; i < N; i++) {
            int min = (int) (rand.nextDouble() * 1000000);
            int max = (int) (Math.abs(rand.nextGaussian()) * 500) + min;
            st6.put(new Interval1D(min, max), "I6_" + (i+1));
            int mid_tst = rand.nextInt(max - min + 1) + min;
            int rad_tst = (int)(rand.nextGaussian() * 1000);
            if (rad_tst < 0)
                ovlaps.add(new Interval1D(Math.max(0, mid_tst + rad_tst), mid_tst - rad_tst));
            else
                ovlaps.add(new Interval1D(Math.max(0, mid_tst - rad_tst), mid_tst + rad_tst));
        }
        for (int i = 0; i < N; i++) {
            int min = (int) (rand.nextDouble() * 1000000);
            int max = (int) (Math.abs(rand.nextGaussian()) * 100) + min;
            queries.add(new Interval1D(min, max));
        }
    }

}