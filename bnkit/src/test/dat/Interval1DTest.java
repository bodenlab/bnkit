package dat;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Created by mikael on 19/10/16.
 */
public class Interval1DTest {
    List<Interval1D> iset1 = new ArrayList<>();
    List<Interval1D> iset2 = new ArrayList<>();

    @BeforeEach
    public void setUp() throws Exception {
        for (int i = 0; i < ivals1.length; i ++)
            iset1.add(new Interval1D(ivals1[i][0], ivals1[i][1]));
        for (int i = 0; i < ivals2.length; i ++)
            iset2.add(new Interval1D(ivals2[i][0], ivals2[i][1]));
    }

    @Test
    public void clip() throws Exception {
        int[] correct_size1 = {0, 3, 5, 4, 1, 0, 0};
        int[] correct_size2 = {0, 2, 2, 2, 2, 0, 0};
        for (int n = 1; n < 7; n ++) {
//            System.out.println("N=" + n);
            Set<Interval1D> res1 = Interval1D.clip(iset1, n);
            assertTrue(res1.size() == correct_size1[n]);
//            for (Interval1D r : res1)
//                System.out.println("1: " + r);
            Set<Interval1D> res2 = Interval1D.clip(iset2, n);
            assertTrue(res2.size() == correct_size2[n]);
//            for (Interval1D r : res2)
//                System.out.println("2: " + r);
        }
    }

    @Test
    public void diff() throws Exception {
        Interval1D i2_5 = new Interval1D(2, 5);
        Interval1D i3_4 = new Interval1D(3, 4);
        Interval1D i4_6 = new Interval1D(4, 6);
        Interval1D i1_3 = new Interval1D(1, 3);
        Set<Interval1D> pair = Interval1D.diff(i2_5, i3_4);
        assertTrue(pair.size() == 2);
        Set<Interval1D> zero = Interval1D.diff(i3_4, i2_5);
        assertTrue(zero.size() == 0);
        Set<Interval1D> left = Interval1D.diff(i2_5, i4_6);
        for (Interval1D ival : left)
            assertTrue(ival.equals(new Interval1D(2, 4)));
        Set<Interval1D> right = Interval1D.diff(i4_6, i2_5);
        for (Interval1D ival : right)
            assertTrue(ival.equals(new Interval1D(5, 6)));
        left = Interval1D.diff(i1_3, i4_6);
        for (Interval1D ival : left)
            assertTrue(ival.equals(i1_3));
        right = Interval1D.diff(i4_6, i1_3);
        for (Interval1D ival : right)
            assertTrue(ival.equals(i4_6));
    }

    @Test
    public void diff1() throws Exception {
        Set<Interval1D> myset = new HashSet<>();
        myset.add(new Interval1D(1,4));
        myset.add(new Interval1D(6,10));
        myset.add(new Interval1D(12,15));
        myset.add(new Interval1D(16,18));
        Set<Interval1D> res1 = Interval1D.diff(new Interval1D(3, 17), myset);
        for (Interval1D ival : res1) {
            assertTrue(ival.equals(new Interval1D(15, 16)) || ival.equals(new Interval1D(4, 6)) || ival.equals(new Interval1D(10, 12)));
        }
        Set<Interval1D> res2 = Interval1D.diff(new Interval1D(7, 235), new HashSet<Interval1D>(iset1));
        for (Interval1D ival : res2) {
            assertTrue(ival.equals(new Interval1D(21, 35)) || ival.equals(new Interval1D(96, 98)) || ival.equals(new Interval1D(230, 235)));
        }
    }

    @Test
    public void union() throws Exception {
        Set<Interval1D> res = Interval1D.union(iset2);
        for (Interval1D ival : res) {
            assertTrue(ival.intersects(new Interval1D(500, 505)));
        }
    }

    int[][] ivals1 = new int[][] {
            {  9,  15},
            {211, 217},
            {120, 230},
            {221, 225},
            {  6,  21},
            { 16,  17},
            {220, 230},
            {221, 224},
            {  6,  21},
            { 40,  51},
            { 98, 199},
            { 98, 199},
            { 35,  96}};

    int[][] ivals2 = new int[][] {
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
/*
    N=1
            1: [6, 21]
            1: [98, 230]
            1: [35, 96]
            2: [503, 5000]
            2: [399, 502]
    N=2
            1: [6, 21]
            1: [98, 199]
            1: [211, 217]
            1: [220, 230]
            1: [40, 51]
            2: [1000, 4700]
            2: [399, 501]
    N=3
            1: [120, 199]
            1: [16, 17]
            1: [221, 225]
            1: [9, 15]
            2: [400, 500]
            2: [1100, 3000]
    N=4
            1: [221, 224]
            2: [2000, 2400]
            2: [498, 499] */
}