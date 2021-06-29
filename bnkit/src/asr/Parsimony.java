package asr;

import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Class for inferring values on un-instantiated internal branch points in a TreeInstance.
 *
 * @author mikael
 */
public class Parsimony implements TreeDecor<List> {

    // TODO: re-write the constructor and decorate so that the tree instance is passed per TreeDecor,
    //  which in turn means that multi-threading can be implemented easily.
    // TODO: implement multi-threading in Prediction for INDEL inference

    private int nnodes;
    private final IdxTree tree;
    private Inference inf = null;
    private boolean recodeNull;

    /**
     * Set-up for performing parsimony
     */
    public Parsimony(IdxTree tree) {
        this(tree, false);
    }

    /**
     * Set-up for performing parsimony
     */
    public Parsimony(IdxTree tree, boolean recodeNull) {
        this.tree = tree;
        this.nnodes = tree.getSize();
        this.recodeNull = recodeNull;
    }

    /**
     * Exposes results of inference.
     * @param idx
     * @return
     */
    @Override
    public List getDecoration(int idx) {
        return getOptimal(idx);
    }

    /**
     * Completes inference.
     * @param ti
     */
    @Override
    public void decorate(TreeInstance ti) {
        if (ti.getTree() != this.tree)
            throw new ASRRuntimeException("Incompatible input to tree for parsimony");
        inf = infer(ti, recodeNull);
    }

    /**
     * Completes inference.
     * @param ti
     */
    public synchronized Inference infer(TreeInstance ti, boolean recodeNull) {
        Object[] key = ti.encode(recodeNull); // convert instances to standardised list
        Inference myinf = new Inference(ti);
        myinf.forward();
        myinf.backward();
        boolean nullAtAncestors = false;
        if (key[key.length - 1] == null) // an additional symbol was added to infer null at ancestors
            nullAtAncestors = true;
        for (int i = 0; i < nnodes; i ++) {
            if (myinf.optimal[i] != null) {
                boolean foundNull = false;
                List mylist = new ArrayList();
                for (Object sym : myinf.optimal[i]) {
                    try {
                        Integer mysym = (Integer) sym;
                        if (nullAtAncestors && mysym == key.length - 1)  // this has null at ancestors
                            foundNull = true;
                        else
                            mylist.add(key[mysym]);
                    } catch (ClassCastException e) {
                        throw new ASRRuntimeException("Invalid symbol inferred by parsimony: " + sym);
                    }
                }
                if (foundNull && mylist.size() == 0) // ancestor has only null as optimal value
                    myinf.optimal[i] = Collections.EMPTY_LIST;
                else if (mylist.size() > 0)
                    myinf.optimal[i] = mylist;
            }
        }
        return myinf;
    }

    /**
     * Performs the forward pass only in parsimony. Throws away the result that is required for backward pass.
     * Prefer {@link #decorate(TreeInstance)}.
     * @return
     */
    public double[] forward(TreeInstance ti) {
        Inference inf = new Inference(ti);
        return inf.forward();
    }

    public List getOptimal(int bpidx) {
        if (inf != null)
            return inf.optimal[bpidx];
        else
            return null;
    }

    /**
     * String representation of the instantiated tree on the Newick format.
     * @return string representation
     */
    public String toString() {
        return toString(0);
    }

    /**
     * String representation of the node and its children (recursively) on the Newick format.
     * @return string representation
     */
    public String toString(int bpidx) {
        StringBuilder sb = new StringBuilder();
        int[] children = tree.getChildren(bpidx);
        int cnt = 0;
        for (int child : children) {
            sb.append(toString(child));
            if (++ cnt < children.length)
                sb.append(",");
        }
        List assigns = getOptimal(bpidx);
        String istr = (TreeInstance.LABEL_INCLUDES_INDEX ? "#" + bpidx + "_" : "?");
        if (assigns != null) {
            StringBuilder sb2 = (TreeInstance.LABEL_INCLUDES_INDEX ? new StringBuilder("#" + bpidx + "_") : new StringBuilder());
            for (int i = 0; i < assigns.size(); i++)
                sb2.append(assigns.get(i) + (i < assigns.size() - 1 ? "|" : ""));
            istr = sb2.toString();
        }
        if (children.length < 1)
            return istr;
        else
            return "(" + sb.toString() + ")" + istr;
    }

    // inner class def ends
    private static Random random = new Random(System.currentTimeMillis());

    /**
     * Shuffles the elements in the array randomly.
     * Code based on methods in java.util.Collections.shuffle();
     */
    protected static void shuffle(int[] array) {
        int count = array.length;
        for (int i = count; i > 1; i--) {
            swap(array, i - 1, random.nextInt(i));
        }
    }

    /**
     * Helper function to shuffle.
     *
     * @param array
     * @param i
     * @param j
     */
    private static void swap(int[] array, int i, int j) {
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }

    /**
     * Generates an array with values from 0 to specified n.
     * The order is either ascending or shuffled, depending on parameter shuffled.
     *
     * @param n        number of elements, populating the array with 0 up to n - 1
     * @param shuffled if true, the array will be returned with elements in random order
     * @return the array
     */
    protected static int[] range(int n, boolean shuffled) {
        int[] ret = new int[n];
        for (int i = 0; i < n; i++)
            ret[i] = i;
        if (shuffled)
            shuffle(ret);
        return ret;
    }

    /**
     * If one or multiple solutions should be identified
     */
    public boolean SET_ONE_TARGET_PARSIMONY = false;
    /**
     * If solutions are listed in order or if they should be randomised
     */
    public boolean SET_RANDOM_PARSIMONY = false;

    /**
     * Inner class to hide some parameters used during inference; else, they could clog up memory if many Parsimony instances are active.
     */
    public class Inference {
        private final int[][][][] traceback;    // [node idx][parent value][child branch][best child value/s: 0, 1, 2, ...] = optimal child state,
        private final TreeInstance treeInstance;
        private final double[][] scores;        // [node idx][parent value] the optimal score for each parent value
        private final int nsym;
        private List[] optimal;                 // [node idx] list of (optimal) values at branch point


        public Inference(TreeInstance ti) {
            this.treeInstance = ti;
            this.nsym = ti.getNPossibleValues();
            this.scores = new double[nnodes][nsym];
            this.optimal = new List[nnodes];

            this.traceback = new int[nnodes][nsym][][];
            for (int bpidx = 0; bpidx < ti.getSize(); bpidx ++) {
                Object bpval = ti.getInstance(bpidx);
                if (bpval != null) {
                    this.optimal[bpidx] = Collections.singletonList(bpval);
                    for (int i = 0; i < nsym; i++)
                        this.scores[bpidx][i] = bpval.equals(ti.getValueByIndex(i)) ? 0 : (double) Integer.MAX_VALUE;
                } else {
                    this.optimal[bpidx] = null;
                }
            }
        }

        /**
         * Perform forward step of maximum parsimony.
         * Many values can be optimally parsimonious.
         *
         * @return the scores of the unique values at the root
         */
        public double[] forward() {
            return forward(0);
        }

        /**
         * Internal function that operates recursively to first initialise each node (forward),
         * stopping only once a value has been assigned to the node,
         * then to propagate scores from assigned nodes to root (backward).
         * Extended to deal with multiple values contributing to optimal scores
         *
         * @return the scores of the unique values at the root
         */
        protected double[] forward(int bpidx) {
            if (optimal[bpidx] != null) // has been instantiated so scores will have been set already
                return scores[bpidx];   // the above condition could also check instantiation.get(bpidx) but this would not apply to recursively determined values
            int[] children = tree.getChildren(bpidx);
            if (children.length == 0) { // un-instantiated leaf node
                for (int i = 0; i < nsym; i ++)
                    scores[bpidx][i] = 0; // all states are equally NOT penalised
                optimal[bpidx] = new ArrayList(); // this also flags we have/are about to resolve the scores of this node
                return scores[bpidx];
            }
            // This node is NOT instantiated, but has children, so prepare to find determine scores and optimal values
            // First, allocate space to hold optimal values; these are only added in "backward" BUT...
            optimal[bpidx] = new ArrayList(); // this also flags we have/are about to resolve the scores of this node
            // recurse into children nodes...
            // determine scores contributed by each child (cscores) BEFORE substitution penalties
            double[][] cscores = new double[children.length][];
            for (int c = 0; c < children.length; c++)
                cscores[c] = forward(children[c]); // one score from child c for each symbol
            // traceback array needs to hold all child symbol indices that contribute to (indexed) parent symbol score via (indexed) child
            for (int i = 0; i < nsym; i++)
                traceback[bpidx][i] = new int[children.length][];
            // need to work out best parent score for "scores"; they all start at 0 (init at allocation above)
            for (int c = 0; c < children.length; c++) {
                // loop through each possible parent assignment,
                // record what symbol in each child that best supports this (adding substitution penalties as we go)
                for (int i = 0; i < nsym; i++) {
                    double best_score = 9E9;
                    int best_cnt = 0; // this is how many symbols in child that need to be recorded (multiple are possible)
                    // next, loop through each possible value in this child to find best parent score (may be the score from multiple origins)
                    for (int j = 0; j < cscores[c].length; j++) {
                        double parent_score = cscores[c][j] + (i == j ? 0 : 1);
                        if (parent_score < best_score) { // new best, reset count to 1
                            best_score = parent_score;
                            best_cnt = 1;
                        } else if (parent_score == best_score) { // new equal best, add +1 to count
                            best_cnt ++;
                        } // else, let this parent score slip
                    }
                    // now we know what the best parent score is; work out all assignments in child c that give it (could be multiple)
                    traceback[bpidx][i][c] = new int[best_cnt]; // allocate space to hold their indices
                    int k = 0; // index for holding possible origins
                    // AGAIN, loop through each possible child symbol, again adding substitution penalties
                    for (int j = 0; j < cscores[c].length; j++) {
                        if (cscores[c][j] + (i == j ? 0 : 1) == best_score)
                            traceback[bpidx][i][c][k ++] = j; // save the kth optimal parent -> child (parent state i -> child state j)
                    }
                    scores[bpidx][i] += best_score; // the best we can do with parent symbol i, from child c
                } // finished the score for this parent each state i, for one child c
            } // finished all children here
            return scores[bpidx]; // done, return parent scores (recursing them as child scores up the tree)
        }

        /**
         * Performs the backward step of maximum parsimony:
         * tracing back all possible values that are optimally parsimonious.
         *
         * Two backward functions handle the "root" and internal nodes, respectively.
         * This handles the node as if it was root; goes by scores assigned to states/symbols.
         */
        public void backward() {
            int bpidx = 0; // root
            int best_index = 0; // find one index with the best score (could be many but one is enough)
            for (int i = 1; i < nsym; i++) {
                if (scores[bpidx][i] < scores[bpidx][best_index])
                    best_index = i;
            }
            // Go through each score and when it is "best", recurse into each child, propagating the state for the score
            // This iteration could be randomised (SET_RANDOM_PARSIMONY=true) so that the states are assigned in different order
            boolean butt_out = false;
            for (int parent_state : range(nsym, SET_RANDOM_PARSIMONY)) {
                if (scores[bpidx][parent_state] == scores[bpidx][best_index]) {
                    // now we know the index of (one of) the optimal parent state/s
                    optimal[bpidx].add(treeInstance.getValueByIndex(parent_state));
                    int[] children = tree.getChildren(bpidx);
                    for (int c = 0; c < children.length; c ++) {
                        int childidx = children[c];
                        // Go through each child state that (optimally) supports the current parent state
                        for (int best_in_child_index : range(traceback[bpidx][parent_state][c].length, SET_RANDOM_PARSIMONY)) {
                            int best_in_child_state = traceback[bpidx][parent_state][c][best_in_child_index]; // optimal child state
                            backward(childidx, best_in_child_state);
                            // if we are interested in only one optimal solution, we can butt out...
                            if (SET_ONE_TARGET_PARSIMONY) {
                                butt_out = true;
                                break;
                            }
                        }
                    }
                }
                if (butt_out)
                    break;
            }
        }

        /**
         * Two backwardParsimony functions handle the "root" and internal nodes, respectively.
         * This function handles internal nodes and should only be called by the former; goes by state assigned to parent.
         *
         * @param optimal_state the index of the state assigned to this node
         */
        private void backward(int bpidx, int optimal_state) {
            // now we know the index of the parent
            Object y = treeInstance.getValueByIndex(optimal_state);
            try {
                if (optimal[bpidx].contains(y)) // check so that the state is assigned only once...
                    return;                                 // ...avoid recurse since this value has been seen here before
            } catch (NullPointerException e) {
                e.printStackTrace();
            }
            optimal[bpidx].add(y);
            int[] children = tree.getChildren(bpidx);
            for (int c = 0; c < children.length; c ++) {
                int childidx = children[c];
                // For each child: choose the first optimal state, or go through each. Ordered randomly, or dictated by traceback/parent.
                for (int best_in_child_index : range(traceback[bpidx][optimal_state][c].length, SET_RANDOM_PARSIMONY)) {
                    int best_in_child_state = traceback[bpidx][optimal_state][c][best_in_child_index];
                    backward(childidx, best_in_child_state);
                    if (SET_ONE_TARGET_PARSIMONY)
                        break;
                }
            }
        }

        public double getScore(int parent_value_index) {
            Object val = treeInstance.getInstance(0);
            if (val != null ) // this node is set
                return (parent_value_index == treeInstance.getIndexByValue(val) ? 0 : Double.POSITIVE_INFINITY);
            double score = 0;
            int[] children = tree.getChildren(0);
            for (int chidx : children)
                score += getScore(chidx, parent_value_index);
            return score;
        }

        private double getScore(int bpidx, int parent_value_index) {
            Object val = treeInstance.getInstance(bpidx);
            if (val != null ) { // this node is set
                return (parent_value_index == treeInstance.getIndexByValue(val) ? 0 : 1);
            } else { // not instantiated
                double best = 9E9;
                for (int my_value_index = 0; my_value_index < nsym; my_value_index ++) {
                    double score = (parent_value_index == my_value_index ? 0 : 1);
                    int[] children = tree.getChildren(bpidx);
                    for (int chidx : children)
                        score += getScore(chidx, my_value_index);
                    if (score < best)
                        best = score;
                }
                return best;
            }
        }

        public List getOptimal(int idx) {
            return optimal[idx];
        }
    }

}
