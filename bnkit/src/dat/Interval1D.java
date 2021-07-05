package dat;

/**
 * This code is based on the augmented interval binary search tree for interval with integer coordinates
 * by Robert Sedgewick and Kevin Wayne; see http://algs4.cs.princeton.edu/93intersection/
 */

import java.util.*;

/**
 * Define a one-dimensional interval.
 */
public class Interval1D implements Comparable<Interval1D> {
    public final int min;  // min endpoint
    public final int max;  // max endpoint

    /**
     * Create a new interval. Must have a min point smaller or equal to the max point.
     * @param min start point
     * @param max end point
     */
    public Interval1D(int min, int max) {
        if (min <= max) {
            this.min = min;
            this.max = max;
        } else
            throw new RuntimeException("Illegal interval");
    }

    /**
     * Create a new interval by performing the union between two intersecting intervals.
     * @param first
     * @param second
     * @return
     */
    public static Interval1D union(Interval1D first, Interval1D second) {
        if (first.intersects(second)) {
            int min = (first.min < second.min) ? first.min : second.min;
            int max = (first.max < second.max) ? second.max : first.max;
            return new Interval1D(min, max);
        } else
            throw new RuntimeException("Illegal union");
    }

    /**
     * Create a new interval by performing the intersection of two intervals, which must of course be intersecting.
     * @param first
     * @param second
     * @return
     */
    public static Interval1D intersection(Interval1D first, Interval1D second) {
        if (first.intersects(second)) {
            int min = (first.min > second.min) ? first.min : second.min;
            int max = (first.max > second.max) ? second.max : first.max;
            return new Interval1D(min, max);
        } else
            throw new RuntimeException("Illegal intersection");
    }

    /**
     * Calculate the Jaccard index between two intervals,
     * i.e. the size of their intersection divided by the size of their union
     * @param first interval
     * @param second interval
     * @return the Jaccard index
     */
    public static double jaccard(Interval1D first, Interval1D second) {
        if (first.intersects(second)) {
            int isect_min = (first.min > second.min) ? first.min : second.min;
            int isect_max = (first.max > second.max) ? second.max : first.max;
            int union_min = (first.min < second.min) ? first.min : second.min;
            int union_max = (first.max < second.max) ? second.max : first.max;
            int denom = union_max - union_min;
            if (denom > 0)
                return (double)(isect_max - isect_min) / (double)denom;
            return 0;
        } else
            return 0;
    }

    /**
     * Calculate the centre position (mid point) of this interval.
     * @return the mid point
     */
    public int getCentre() {
        return ((max - min) / 2) + min;
    }

    /**
     * Check if this interval intersects with another.
     * Note that this returns true if the intervals "touch" one another, i.e. they need not overlap
     * @see {@link Interval1D#overlaps(Interval1D)}
     * @param that other interval
     * @return true if this and that intervals intersect, false otherwise
     */
    public boolean intersects(Interval1D that) {
        if (that.max < this.min) return false;
        if (this.max < that.min) return false;
        return true;
    }

    /**
     * Check if this interval overlaps with another.
     * Note that this returns false if the intervals only "touch" one another,
     * instead they must share part of their intervals for this function to return true
     * @see {@link Interval1D#intersects(Interval1D)}
     * @param that other interval
     * @return true if this and that intervals overlap, false otherwise
     */
    public boolean overlaps(Interval1D that) {
        if (that.max <= this.min) return false;
        if (this.max <= that.min) return false;
        return true;
    }

    /**
     * Check if this interval contains a point
     * @param x point to investigate
     * @return true if the x is within the interval (including start and end points)
     */
    public boolean contains(int x) {
        return (min <= x) && (x <= max);
    }

    /**
     * Check if this interval contains another interval
     * @param other interval to investigate
     * @return true if the other interval is within this interval (including start and end points)
     */
    public boolean contains(Interval1D other) {
        return (min <= other.min) && (other.max <= max);
    }

    /**
     * Compare the intervals. If unequal, this operation imposes a strict ordering, meaning that
     * an interval that has an earlier start point will be BEFORE,
     * and only if the start points are the same in two intervals,
     * the interval with the earlier end point will be BEFORE.
     * @param that the other interval
     * @return -1 if this interval is before that; +1 if that interval is before this; 0 if equal.
     */
    public int compareTo(Interval1D that) {
        // check min value first: which interval comes first?
        return compareTo(this.min, this.max, that.min, that.max);
    }

    /**
     * Compare the intervals. If unequal, this operation imposes a strict ordering, meaning that
     * an interval that has an earlier start point will be BEFORE,
     * and only if the start points are the same in two intervals,
     * the interval with the earlier end point will be BEFORE.
     * @param this_min the start point of interval "this"
     * @param this_max the end point of interval "this"
     * @param that_min the start point of interval "that"
     * @param this_max the end point of interval "that"
     * @return -1 if this interval is before that; +1 if that interval is before this; 0 if equal.
     */
    public static int compareTo(int this_min, int this_max, int that_min, int that_max) {
        // check min value first: which interval comes first?
        if      (this_min < that_min) return -1;
        else if (this_min > that_min) return +1;
            // only check max value if min values are equal: which interval ends first?
        else if (this_max < that_max) return -1;
        else if (this_max > that_max) return +1;
        else                          return  0;
    }

    @Override
    public boolean equals(Object other) {
        if (this == other)
            return true;
        if (other == null)
            return false;
        if (getClass() != other.getClass())
            return false;
        return (this.min == ((Interval1D)other).min && this.max == ((Interval1D)other).max);
    }

    @Override
    public int hashCode() {
        return Objects.hash(this.min, this.max);
    }

    /**
     * Get the width of the interval
     * @return
     */
    public int getWidth() {
        return this.max - this.min;
    }

    /**
     * Calculate the distance between this and that other interval given
     * @param that other interval
     * @param minimum if true, then use minimum distance, if false, use centre to centre
     * @return the distance
     */
    public int dist(Interval1D that, boolean minimum) {
        if (minimum) {
            // that interval is BEFORE this
            if (this.min > that.max) return this.min - that.max;
            // that interval is AFTER this
            if (this.max < that.min) return that.min - this.max;
            return  0;
        } else { // centre-to-centre
            return Math.abs(that.getCentre() - this.getCentre());
        }
    }

    /**
     * Calculate the signed distance between this and that other interval given
     * @param that other interval
     * @param minimum if true, then use minimum distance, if false, use centre to centre
     * @return the signed distance (negative if that interval is before this, positive is that interval is after this)
     */
    public int signdist(Interval1D that, boolean minimum) {
        if (minimum) {
            // that interval is BEFORE this
            if (this.min > that.max) return  that.max - this.min;
            // that interval is AFTER this
            if (this.max < that.min) return that.min - this.max;
            return  0;
        } else { // centre-to-centre
            return that.getCentre() - this.getCentre();
        }
    }

    /**
     * Determine the interval or pair of intervals that define the (clipped) difference between the input intervals.
     * @param input1 interval 1
     * @param input2 interval 2
     * @return
     */
    public static Set<Interval1D> diff(Interval1D input1, Interval1D input2) {
        Set<Interval1D> ret = new HashSet<>();
        if (input1.min < input2.min) { // input1 begins before input2
            if (input1.max > input2.max) { // input2 ends before input1 so new interval is split
                // return pair
                Interval1D i1 = new Interval1D(input1.min, input2.min);
                Interval1D i2 = new Interval1D(input2.max, input1.max);
                ret.add(i1);
                ret.add(i2);
            } else { // input1 ends before or at the same time as input2
                Interval1D i1 = new Interval1D(input1.min, Math.min(input1.max, input2.min));
                ret.add(i1);
            }
        } else if (input1.min > input2.min) { // input2 begins before input1
            if (input1.max > input2.max) { // input2 ends before or at the same time as input1
                Interval1D i1 = new Interval1D(Math.max(input2.max, input1.min), input1.max);
                ret.add(i1);
            } else {
                // no interval
            }
        }
        return ret;
    }

    /**
     * Determine the set of intervals that define the (clipped) difference between the input1 interval and the input2 intervals.
     * @param foreground primary interval
     * @param background intervals
     * @return
     */
    public static Set<Interval1D> diff(Interval1D foreground, Set<Interval1D> background) {
        Set<Interval1D> ret = new HashSet<>();
        if (background == null) {
            ret.add(foreground);
            return ret;
        } else if (background.size() == 0) {
            ret.add(foreground);
            return ret;
        }
        int N = 1;
        int[] mins = new int[background.size()];
        int[] maxs = new int[background.size()];
        int idx = 0;
        for (Interval1D ival : background) {
            mins[idx] = ival.min;
            maxs[idx] = ival.max;
            idx += 1;
        }
        Arrays.sort(mins);
        Arrays.sort(maxs);
        int idx_min = 0;
        int idx_max = 0;
        int CNT = 0; // number of intervals currently entered but not yet exited
        boolean within = false; // inside an interval that is being created
        Interval1D current = null;
        if (foreground.min < mins[0])
            current = new Interval1D(foreground.min, foreground.max);
        while (idx_min < mins.length && idx_max < maxs.length) { // there are more elements in both the min and max lists
            if (mins[idx_min] < maxs[idx_max]) {   // background interval entered, which means that we may now be blocking the foreground
                CNT ++;
                if (current != null && CNT >= N) { // we are blocking the foreground, so complete the diff'ed foreground interval
                    // where does it start and end?
                    current = new Interval1D(current.min, Math.min(mins[idx_min], foreground.max));
                    ret.add(current);
                    current = null;
                }
                idx_min ++;
            } else if (mins[idx_min] > maxs[idx_max]) { // background interval exited, meaning that maybe the foreground is now accessible
                CNT --;
                if (current == null && CNT < N) { // create new tentative diff'ed foreground interval
                    current = new Interval1D(maxs[idx_max], foreground.max);
                }
                idx_max ++;
            } else {                                    // one entered, one exited
                idx_min ++;
                idx_max ++;
            }
        }
        if (idx_min < mins.length)
            throw new RuntimeException("Invalid intervals");
        if (idx_max < maxs.length) {
            if (current == null && foreground.max > maxs[CNT - N + idx_max]) {
                current = new Interval1D(maxs[CNT - N + idx_max], foreground.max);
                ret.add(current);
            }
        }
        return ret;
    }


    /**
     * Determine the set of intervals (zero, one or more) that define the union between the list of input intervals.
     * @param input list of intervals
     * @return
     */
    public static Set<Interval1D> union(List<Interval1D> input) {
        return clip(input, 1);
    }

    /**
     * Determine the set of intervals that define the N-fold intersection between the input intervals.
     * Note that zero-length intervals are NOT adding to the intersection count.
     * @param input
     * @param N
     * @return
     */
    public static Set<Interval1D> clip(List<Interval1D> input, int N) {
        if (N < 1)
            throw new RuntimeException("clip: N must be 1 or greater.");
        Set<Interval1D> ret = new HashSet<>();
        int[] mins = new int[input.size()];
        int[] maxs = new int[input.size()];
        int idx = 0;
        for (Interval1D ival : input) {
            mins[idx] = ival.min;
            maxs[idx] = ival.max;
            idx += 1;
        }
        Arrays.sort(mins);
        Arrays.sort(maxs);
        int idx_min = 0;
        int idx_max = 0;
        int CNT = 0; // number of intervals currently entered but not yet exited
        boolean within = false; // inside an interval that is being created
        Interval1D current = null;
        while (idx_min < mins.length && idx_max < maxs.length) {
            if (mins[idx_min] < maxs[idx_max]) {        // interval entered
                CNT ++;
                if (current == null && CNT >= N) {
                    current = new Interval1D(mins[idx_min], maxs[idx_max]);
                }
                idx_min ++;
            } else if (mins[idx_min] > maxs[idx_max]) { // interval exited
                CNT --;
                if (current != null && CNT < N) {
                    current = new Interval1D(current.min, maxs[idx_max]);
                    ret.add(current);
                    current = null;
                }
                idx_max ++;
            } else {                                    // one entered, one exited
                //CNT ++;
                if (current == null && CNT >= N) {
                    current = new Interval1D(mins[idx_min], maxs[idx_max]);
                }
                idx_min ++;
                //CNT --;
                if (current != null && CNT < N) {
                    current = new Interval1D(current.min, maxs[idx_max]);
                    ret.add(current);
                    current = null;
                }
                idx_max ++;
            }
        }
        if (idx_min < mins.length)
            throw new RuntimeException("Invalid intervals");
        if (idx_max < maxs.length && current != null) {
            current = new Interval1D(current.min, maxs[CNT - N + idx_max]);
            ret.add(current);
        }
        return ret;
    }

    public String toString() {
        return "[" + min + ", " + max + "]";
    }

    // test client
    public static void main(String[] args) {
        Interval1D a = new Interval1D(15, 20);
        Interval1D b = new Interval1D(25, 30);
        Interval1D c = new Interval1D(10, 28);
        Interval1D d = new Interval1D(40, 50);
        Interval1D bUc = Interval1D.union(b, c);
        Interval1D bIc = Interval1D.intersection(b, c);

        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("c = " + c);
        System.out.println("d = " + d);
        System.out.println("bUc = " + bUc);
        System.out.println("bIc = " + bIc);

        System.out.println("b intersects a = " + b.intersects(a));
        System.out.println("a intersects b = " + a.intersects(b));
        System.out.println("a intersects c = " + a.intersects(c));
        System.out.println("a intersects d = " + a.intersects(d));
        System.out.println("b intersects c = " + b.intersects(c));
        System.out.println("b intersects d = " + b.intersects(d));
        System.out.println("c intersects d = " + c.intersects(d));

    }

}