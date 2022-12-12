package dat;

/**
 * This code is based on the augmented interval binary search tree for interval with integer coordinates
 * by Robert Sedgewick and Kevin Wayne; see http://algs4.cs.princeton.edu/93intersection/
 * The tree inserts new nodes by traversing the tree, following the order that is defined on intervals,
 * and randomly inserts them as new nodes (with a probability proportional to the number of intervals below
 * the current node).
 * Duplicate intervals are (optionally) managed.
 *
 * FIXME:
 * centre-to-centre distances are not handled and a resolution is to have a separate binary search tree,
 * dividing intervals by the corresponding ordering...
 *
 */

import java.util.*;

/**
 * Binary search tree for storing int integer intervals, and for performing queries on them.
 * See https://en.wikipedia.org/wiki/Interval_tree, specifically the Augmented kind.
 * The present implementation balances the tree by using randomisation.
 * @param <Value> the data type which applies to values associated with intervals
 */
public class IntervalST<Value> implements Iterable<Interval1D> {

    public final int DISTANCE_MINIMUM = 1;
    public final int DISTANCE_CENTRE = 2; // centre-to-centre; is not implemented
    public final int SETUP_DISTANCE = DISTANCE_MINIMUM;

    public final int ISECT_STRATEGY_PICKONE = 1;
    public final int ISECT_STRATEGY_JACCARD = 2;
    public final int ISECT_STRATEGY_CENTRE = 3;
    public int SETUP_MULTIPLE_INTERSECTION = ISECT_STRATEGY_JACCARD; // handle overlap

    private Node root = null;   // root of the BST
    private Random rand = null;

    // BST helper node data type
    private class Node {
        Interval1D interval;      // key
        Set<Value> values = new HashSet<>();              // associated data
        Node left, right;         // left and right subtrees
        int N;                    // size of subtree rooted at this node
        int min, max;            // min and max endpoint in subtree rooted at this node

        Node(Interval1D interval, Value value) {
            this.interval = interval;
            this.values.add(value);
            this.N = 1;
            this.min = interval.min;
            this.max = interval.max;
        }
        public String toString() {
            return (left != null ? "o_":"x_") + min + interval.toString() + max + (right != null ? "_o":"_x");
        }
        boolean add(Value value) {
            return this.values.add(value);
        }
    }

    public IntervalST(int seed) {
        rand = new Random(seed);
    }

    public IntervalST() {
        this((int)System.currentTimeMillis());
    }

    @Override
    public Iterator<Interval1D> iterator() {
        return new Interval1DIterator();
    }

    private class Interval1DIterator implements Iterator<Interval1D> {
        final private Stack<Node> cursorstack = new Stack<>();
        private Node current;

        public Interval1DIterator() {
            current = root;
        }

        private Interval1DIterator(Node start) {
            current = start;
        }

        @Override
        public boolean hasNext() {
            return !cursorstack.isEmpty() || current != null;
        }

        @Override
        public Interval1D next() {
            while (current != null) {
                cursorstack.push(current);
                current = current.left;
            }
            current = cursorstack.pop();
            Node ret = current;
            current = current.right;
            return ret.interval;
        }
    }

    /***************************************************************************
     *  BST search
     ***************************************************************************/

    public boolean contains(Interval1D interval) {
        return (get(interval) != null);
    }

    // return value associated with the given key
    // if no such value, return null
    public Set<Value> get(Interval1D interval) {
        Node match = get(root, interval);
        if (match == null)
            return null;
        return match.values;
    }

    private Node get(Node x, Interval1D interval) {
        if (x == null) return null;
        int cmp = interval.compareTo(x.interval);
        if (cmp < 0) return get(x.left, interval);
        else if (cmp > 0) return get(x.right, interval);
        else return x;
    }

    private Set<Interval1D> getAll(Node x) {
        if (x == null)
            return Collections.EMPTY_SET;
        else {
            Set<Interval1D> ivals = new HashSet<>();
            if (x.left != null)
                ivals.addAll(getAll(x.left));
            if (x.right != null)
                ivals.addAll(getAll(x.right));
            ivals.add(x.interval);
            return ivals;
        }
    }

    public Set<Interval1D> getAll() {
        return getAll(root);
    }

    /***************************************************************************
     *  randomized insertion
     ***************************************************************************/
    public void put(Interval1D interval, Value value) {
        Node match = get(root, interval);
        if (match != null) { // duplicate interval
            match.add(value);
            // optionally:
            // remove(interval);
            // in which case we also need to insert new one as below...
        } else
            root = randomizedInsert(root, interval, value);
    }

    // make new node the root with uniform probability
    private Node randomizedInsert(Node x, Interval1D interval, Value value) {
        if (x == null) return new Node(interval, value);
        if (rand.nextDouble() * size(x) < 1.0) return rootInsert(x, interval, value);
        int cmp = interval.compareTo(x.interval);
        if (cmp < 0) x.left = randomizedInsert(x.left, interval, value);
        else x.right = randomizedInsert(x.right, interval, value);
        fix(x);
        return x;
    }

    private Node rootInsert(Node x, Interval1D interval, Value value) {
        if (x == null) return new Node(interval, value);
        int cmp = interval.compareTo(x.interval);
        if (cmp < 0) {
            x.left = rootInsert(x.left, interval, value);
            x = rotR(x);
        } else {
            x.right = rootInsert(x.right, interval, value);
            x = rotL(x);
        }
        return x;
    }


    /***************************************************************************
     *  deletion
     ***************************************************************************/
    private Node joinLR(Node a, Node b) {
        if (a == null) return b;
        if (b == null) return a;

        if (rand.nextDouble() * (size(a) + size(b)) < size(a)) {
            a.right = joinLR(a.right, b);
            fix(a);
            return a;
        } else {
            b.left = joinLR(a, b.left);
            fix(b);
            return b;
        }
    }

    public Value remove(Interval1D interval, Value value) {
        Node match = get(root, interval);
//        synchronized {
            if (match.values.contains(value)) {
                match.values.remove(value);
                if (match.values.isEmpty())
                    root = remove(root, interval);
                return value;
            }
        return null;
    }

    // remove and return all values associated with given interval;
    // if no such interval exists return null
    public Set<Value> remove(Interval1D interval) {
        Node match = get(root, interval);
        root = remove(root, interval);
        return match.values;
    }

    private Node remove(Node h, Interval1D interval) {
        if (h == null) return null;
        int cmp = interval.compareTo(h.interval);
        if (cmp < 0) h.left = remove(h.left, interval);
        else if (cmp > 0) h.right = remove(h.right, interval);
        else h = joinLR(h.left, h.right);
        fix(h);
        return h;
    }



    /***************************************************************************
     *  Interval searching
     ***************************************************************************/

    // return an interval in data structure that intersects the given interval;
    // return null if no such interval exists
    // running time is proportional to log N
    public Interval1D search(Interval1D interval) {
        return search(root, interval);
    }

    // return an interval in data structure that contains the given point;
    // return null if no such interval exists
    // running time is proportional to log N
    public Interval1D search(int point) {
        return search(root, point);
    }

    // look in subtree rooted at x
    public Interval1D search(Node x, Interval1D interval) {
        while (x != null) {
            if (interval.intersects(x.interval)) return x.interval;
            else if (x.left == null) x = x.right;
            else if (x.left.max < interval.min) x = x.right; // go right, because the interval is yet to be found, but needs to be greater than anything on the left
            else x = x.left;
        }
        return null;
    }

    // look in subtree rooted at x
    public Interval1D search(Node x, int point) {
        while (x != null) {
            if (x.interval.contains(point)) return x.interval;
            else if (x.left == null) x = x.right; // go right, because left does not exist
            else if (x.left.max < point)
                x = x.right;  // go right, because the interval is yet to be found, but needs to be greater than anything on the left
            else x = x.left;  // go left, because the interval is yet to be found, but needs to be smaller
        }
        return null;
    }

    /** Flag to set direction of interval search: either before or after the given interval */
    public static final int CLOSEST_BEFORE_OR_AFTER = 0;
    /** Flag to set direction of interval search: before the given interval */
    public static final int CLOSEST_BEFORE = 1;
    /** Flag to set direction of interval search: after the given interval */
    public static final int CLOSEST_AFTER = 2;

    /**
     * Retrieve the interval Y stored in the tree that is closest to the given interval X.
     * If the given interval overlaps with one or more stored intervals, one is returned:
     * if minimum distance--the interval Y with the greatest Jaccard index to X
     * if centre-to-centre distance--the interval Y with the smallest centre displacement from that of X
     * optionally, do not try to resolve overlaps, and return one of them arbitrarily.
     * @param query the interval for which the closest is sought
     * @return the interval closest to the given query interval
     */
    public Interval1D getClosest(Interval1D query) {
        return this.getClosest(query, CLOSEST_BEFORE_OR_AFTER);
    }

    /**
     * Retrieve the interval Y stored in the tree that is closest to the given interval X.
     * If the given interval overlaps with one or more stored intervals, one is returned:
     * if minimum distance--the interval Y with the greatest Jaccard index to X
     * if centre-to-centre distance--the interval Y with the smallest centre displacement from that of X
     * optionally, do not try to resolve overlaps, and return one of them arbitrarily.
     * @param query the interval for which the closest is sought
     * @param DIRECTION the direction of the search (0 is either, 1 is before, 2 is after)
     * @return the interval closest to the given query interval
     */
    public Interval1D getClosest(Interval1D query, int DIRECTION) {
        LinkedList<Interval1D> list = new LinkedList<Interval1D>();
        if (SETUP_MULTIPLE_INTERSECTION != ISECT_STRATEGY_PICKONE)
            searchAll(root, query, list); // check overlapping intervals first
        if (list.isEmpty()) { // no overlapping intervals are found so search for closest
            if (SETUP_DISTANCE == DISTANCE_MINIMUM) {
                switch (DIRECTION) {
                    case CLOSEST_BEFORE:
                        return getClosestOnSide(root, query, true);
                    case CLOSEST_AFTER:
                        return getClosestOnSide(root, query, false);
                    case CLOSEST_BEFORE_OR_AFTER:
                    default:
                        return getClosest(root, query);
                }
            } else
                throw new RuntimeException("Not implemented: DISTANCE=" + SETUP_DISTANCE);
                // return getClosest(root, query.getCentre());
        } else { // overlapping distances were found, so pick the one with greatest jaccard index
            Interval1D closestIval = null;
            if (SETUP_DISTANCE == DISTANCE_MINIMUM) {
                double best = 0;
                for (Interval1D x : list) {
                    double jacx = Interval1D.jaccard(query, x);
                    if (closestIval == null || jacx > best) {
                        best = jacx;
                        closestIval = x;
                    }
                }
            } else { // centre-to-centre
                throw new RuntimeException("Not implemented: DISTANCE=" + SETUP_DISTANCE);
                /*
                int closestDist = int.MAX_VALUE;
                for (Interval1D x : list) {
                    int distx = Math.abs(query.getCentre() - x.getCentre());
                    if (closestIval == null || distx < closestDist) {
                        closestDist = distx;
                        closestIval = x;
                    }
                }*/
            }
            return closestIval;
        }
    }

    /**
     * Recursively find the interval with the minimum distance to that given.
     * Does not guarantee that distances are sensible when overlapping intervals exist.
     * Essentially it assumes that overlaps have been eliminated prior.
     * @param cand node from which search starts
     * @param query interval
     * @return closest interval
     */
    private Interval1D getClosest(Node cand, Interval1D query) {
        Interval1D fav = null;
        int favDist = -1;
        while (cand != null) {
            if (query.compareTo(cand.interval) == 0) return cand.interval;
            int distx = query.dist(cand.interval, true);
            if (fav == null || distx <= favDist) {
                fav = cand.interval;
                favDist = distx;
            }
            if (cand.left == null) cand = cand.right;
            else if (cand.right == null) cand = cand.left;
            else if (cand.interval.min > query.max) cand = cand.left; // the smallest, indexed value (on left) is AFTER the query min
            else { // no way to choose without looking in the intervals below
                Interval1D favleft = null, favright = null;
                Interval1D baseleft = new Interval1D(cand.left.min, cand.left.max);
                Interval1D baseright = new Interval1D(cand.right.min, cand.right.max);
                int distleft = query.dist(baseleft, true);
                if (distleft < favDist) {
                    favleft = getClosest(cand.left, query);
                    distleft = (favleft != null ? query.dist(favleft, true) : Integer.MAX_VALUE);
                }
                int distright = query.dist(baseright, true);
                if (distright < favDist) {
                    favright = getClosest(cand.right, query);
                    distright = (favright != null ? query.dist(favright, true) : Integer.MAX_VALUE);
                }
                if (distleft < distright)
                    return (distleft < favDist ? favleft : fav);
                else
                    return (distright < favDist ? favright : fav);
            }
        }
        return fav;
    }

    /**
     * Recursively find the interval with the minimum distance to that given.
     * Assumes that the query interval does not overlap with any of the intervals.
     * @param cand node from which search starts
     * @param query interval
     * @param lowerSide if true, find closest on the lower side, else on the upper side
     * @return closest interval
     */
    private Interval1D getClosestOnSide(Node cand, Interval1D query, boolean lowerSide) {
        Interval1D favourite = null;
        int favDist = -1;
        int N = 0; // number of nodes looked at
        while (cand != null) { // cand is initially the root of the tree
            N ++;
            int distq = query.signdist(cand.interval, true);
            if (distq == 0) { // spot on
                favourite = cand.interval;
                break;
            }
            // if distq is negative, cand is BEFORE query
            // if distq is positive, cand is AFTER query
            if (lowerSide) { // we are looking for the closest NEGATIVE distance
                if (distq < 0) {
                    distq = -distq;
                    if (favourite == null || distq <= favDist) { // accept equal too, since the same interval beats intersecting
                        favourite = cand.interval;
                        favDist = distq;
                    }
                    if (cand.left == null) cand = cand.right;
                    else if (cand.right == null) cand = cand.left;
                    else if (cand.interval.min > query.min) cand = cand.left; // the smallest, indexed value (on left) is AFTER the query min
                    else if (cand.right.min > query.max) cand = cand.left; // the smallest possible value on right is AFTER the query max
                    else { // no way to choose without looking in the intervals below
                        Interval1D favleft = null, favright = null;
                        Interval1D baseleft = new Interval1D(cand.left.min, cand.left.max);
                        Interval1D baseright = new Interval1D(cand.right.min, cand.right.max);
                        int distleft = query.signdist(baseleft, true);
                        distleft = (distleft <= 0 ? -distleft : Integer.MAX_VALUE);
                        if (distleft < favDist) {
                            favleft = getClosestOnSide(cand.left, query, lowerSide);
                            distleft = (favleft != null ? query.dist(favleft, true) : Integer.MAX_VALUE);
                        }
                        int distright = query.signdist(baseright, true);
                        distright = (distright <= 0 ? -distright : Integer.MAX_VALUE);
                        if (distright < favDist) {
                            favright = getClosestOnSide(cand.right, query, lowerSide);
                            distright = (favright != null ? query.dist(favright, true) : Integer.MAX_VALUE);
                        }
                        if (distleft < distright)
                            return (distleft < favDist ? favleft : favourite);
                        else
                            return (distright < favDist ? favright : favourite);
                    }
                } else { // distq > 0
                    cand = cand.left;
                }
            } else { // lowerSide is false, looking for positive distance: cand AFTER query
                if (distq > 0) {
                    if (favourite == null || distq <= favDist) { // accept equal too, since the same interval beats intersecting
                        favourite = cand.interval;
                        favDist = distq;
                    }
                    if (cand.left == null) cand = cand.right;
                    else if (cand.right == null) cand = cand.left;
                    else if (query.max < cand.interval.min) cand = cand.left;
                    else if (query.min > cand.interval.min) cand = cand.right;
                    else cand = null;
                } else { // distq < 0
                    cand = cand.right;
                }
            }
        }
        //System.out.println("N=" + N);
        if (favourite == null)
            return null;
        return favourite;
    }

    /**
     * Recursively find the interval with the centre-to-centre distance to the point given.
     * Does not guarantee that distances are sensible when overlapping intervals exist.
     * @param x node from which search starts
     * @param ref point to search for
     * @return
     */
    private Interval1D getClosest(Node x, int ref) {
        Node closestNode = null;
        int closestDist = -1;
        while (x != null) {
            if (x.interval.getCentre() == ref) return x.interval;
            int distx = Math.abs(ref - x.interval.getCentre());
            if (closestNode == null) {
                closestNode = x;
                closestDist = distx;
            } else if (distx <= closestDist) { // accept equal too, since the same interval beats intersecting
                closestNode = x;
                closestDist = distx;
            }
            if (x.left == null) x = x.right;
            else if (distx > 0 && x.interval.min > ref) x = x.left;
            else if (x.left.max < ref) x = x.right;
            else x = x.left;
        }
        return closestNode.interval;
    }

    /**
     * Return *all* intervals in data structure that intersect the given interval.
     * Running time is proportional to R log N, where R is the number of intersections.
     * @param interval query interval
     * @return a set with all intervals that intersect with that given; empty if no intervals were found
     */
    public Set<Interval1D> searchAll(Interval1D interval) {
        LinkedList<Interval1D> list = new LinkedList<>();
        searchAll(root, interval, list);
        return new HashSet<Interval1D>(list);
    }

    /**
     * Recursively find all intervals in data structure under the given node that intersect the given interval.
     * Running time is proportional to R log N, where R is the number of intersections.
     * @param x node from which search is started
     * @param interval query interval
     * @param list the list which is populated while recursing the tree
     * @return a set with all intervals that intersect with that given
     */
    private boolean searchAll(Node x, Interval1D interval, LinkedList<Interval1D> list) {
        boolean found1 = false;
        boolean found2 = false;
        boolean found3 = false;
        if (x == null)
            return false;
        if (interval.intersects(x.interval)) {
            list.add(x.interval);
            found1 = true;
        }
        if (x.left != null && x.left.max >= interval.min)
            found2 = searchAll(x.left, interval, list);
        if (found2 || x.left == null || x.left.max < interval.min)
            found3 = searchAll(x.right, interval, list);
        return found1 || found2 || found3;
    }


    /***************************************************************************
     * Whole-tree operations
     ***************************************************************************/

    /**
     * In a union-like operation, extract the continuously intersecting intervals from the tree.
     * Flat intervals are created from intervals that simply "touch" one another;
     * they are not required to overlap.
     * @return the set of intervals that represent the (flat) union of the intervals in this tree
     */
    public Set<Interval1D> flatten2Set() {
        return flatten2Set(true);
    }
    /**
     * In a union-like operation, extract the continuously intersecting intervals from the tree.
     * @param FLATTEN_BY_TOUCH set to true so that flat intervals are created from intervals that simply "touch",
     *                         not extend over, one another; set to false, to require that they overlap
     * @return the set of intervals that represent the (flat) union of the intervals in this tree
     */
    public Set<Interval1D> flatten2Set(boolean FLATTEN_BY_TOUCH) {
        Set<Interval1D> ret = new HashSet<>();
        Interval1D current = null;
        for (Interval1D ival : this) {
            if (current == null) {
                current = new Interval1D(ival.min, ival.max);
            } else {
                if (ival.max > current.max) { // we will update
                    if (FLATTEN_BY_TOUCH ? (ival.min <= current.max) : (ival.min < current.max)) { // just extending the current
                        current = new Interval1D(current.min, ival.max);
                    } else { // ival.min > current.max, i.e. we need to create a new
                        ret.add(current);
                        current = new Interval1D(ival.min, ival.max);
                    }
                }
            }
        }
        if (current != null)
            ret.add(current);
        return ret;
    }

    /**
     * In a union-like operation, extract the intest continuous intervals from the tree.
     * @return the tree of intervals that represent the (flat) union of the intervals in this tree
     */
    public IntervalST<Value> flatten2Tree() {
        IntervalST<Value> ret = new IntervalST<Value>();
        List<Value> curr_vals = null;
        Interval1D current = null;
        for (Interval1D ival : this) {
            if (current == null) {
                curr_vals = new ArrayList<Value>();
                curr_vals.addAll(get(ival));
                current = new Interval1D(ival.min, ival.max);
            } else {
                if (ival.max > current.max) { // we will update
                    if (ival.min <= current.max) { // just extending the current
                        curr_vals.addAll(get(ival));
                        current = new Interval1D(current.min, ival.max);
                    } else { // ival.min > current.max, i.e. we need to stash old one, and create new
                        for (Value v : curr_vals)
                            ret.put(current, v);
                        current = new Interval1D(ival.min, ival.max);
                        curr_vals = new ArrayList<Value>();
                        curr_vals.addAll(get(ival));
                    }
                } else { // just stash values
                    curr_vals.addAll(get(ival));
                }
            }
        }
        if (current != null) {
            for (Value v : curr_vals)
                ret.put(current, v);
        }
        return ret;
    }

    /***************************************************************************
     *  useful binary tree functions
     ***************************************************************************/

    /**
     * Return number of intervals in tree.
     * Note that the number of values is only equal to this if no interval is associated with more than one value.
     * @return
     */
    public int size() {
        return size(root);
    }

    /**
     * Return number of intervals in subtree rooted at x
     * @return
     */
    private int size(Node x) {
        if (x == null) return 0;
        else return x.N;
    }

    /**
     * Return height of tree (empty tree height = 0)
     * @return
     */
    public int height() {
        return height(root);
    }

    /**
     * Return height of subtree under specified node (empty tree height = 0)
     * @param x the node from which height is calculated
     * @return
     */
    private int height(Node x) {
        if (x == null) return 0;
        return 1 + Math.max(height(x.left), height(x.right));
    }


    /***************************************************************************
     *  helper BST functions
     ***************************************************************************/

    // fix auxiliary information (subtree count, min and max fields)
    private void fix(Node x) {
        if (x == null) return;
        x.N = 1 + size(x.left) + size(x.right);
        x.min = min3(x.interval.min, min(x.left), min(x.right));
        x.max = max3(x.interval.max, max(x.left), max(x.right));
    }

    private int min(Node x) {
        if (x == null) return Integer.MAX_VALUE;
        return x.min;
    }

    private int max(Node x) {
        if (x == null) return Integer.MIN_VALUE;
        return x.max;
    }

    // precondition: a is not null
    private int min3(int a, int b, int c) {
        return Math.min(a, Math.min(b, c));
    }

    // precondition: a is not null
    private int max3(int a, int b, int c) {
        return Math.max(a, Math.max(b, c));
    }

    // right rotate
    private Node rotR(Node h) {
        Node x = h.left;
        h.left = x.right;
        x.right = h;
        fix(h);
        fix(x);
        return x;
    }

    // left rotate
    private Node rotL(Node h) {
        Node x = h.right;
        h.right = x.left;
        x.left = h;
        fix(h);
        fix(x);
        return x;
    }

    /***************************************************************************
     *  Debugging functions that test the integrity of the tree
     ***************************************************************************/

    public boolean checkCount() {
        return checkCount(root);
    }
    private boolean checkCount(Node x) {
        if (x == null) return true;
        return checkCount(x.left) && checkCount(x.right) && (x.N == 1 + size(x.left) + size(x.right));
    }
    public boolean checkMax() {
        return checkMax(root);
    }
    private boolean checkMax(IntervalST.Node x) {
        if (x == null) return true;
        return x.max == max3(x.interval.max, max(x.left), max(x.right));
    }

    public static void main(String[] args) {
        Interval1D a = new Interval1D(13, 20);
        Interval1D b = new Interval1D(25, 30);
        Interval1D c = new Interval1D(27, 34);
        Interval1D d = new Interval1D(40, 50);

        System.out.println("a = " + a);
        System.out.println("b = " + b);
        System.out.println("c = " + c);
        System.out.println("d = " + d);

        System.out.println("b intersects a = " + b.intersects(a));
        System.out.println("a intersects b = " + a.intersects(b));
        System.out.println("a intersects c = " + a.intersects(c));
        System.out.println("a intersects d = " + a.intersects(d));
        System.out.println("b intersects c = " + b.intersects(c));
        System.out.println("b intersects d = " + b.intersects(d));
        System.out.println("c intersects d = " + c.intersects(d));

        IntervalST<String> tree = new IntervalST<>();
        tree.put(a, "A");
        tree.put(b, "B");
        tree.put(c, "C");
        tree.put(d, "D");

        System.out.println(tree.search(16));
        System.out.println(tree.search(14));
        System.out.println(tree.getClosest(new Interval1D(21, 22)));
        System.out.println(tree.getClosest(new Interval1D(11, 12)));
        System.out.println(tree.getClosest(new Interval1D(55, 60)));
        System.out.println(tree.getClosest(new Interval1D(38, 39)));
        System.out.println(tree.getClosest(new Interval1D(29, 33)));
        System.out.println(tree.getClosest(new Interval1D(1, 30)));
        System.out.println(tree.getClosest(new Interval1D(1, 70)));
    }
}