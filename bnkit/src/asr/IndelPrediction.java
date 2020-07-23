package asr;

import dat.Interval1D;
import dat.phylo.IdxTree;
import dat.phylo.TreeInstance;
import dat.pog.Node;
import dat.pog.POGTree;
import dat.pog.POGraph;

import java.util.*;

/**
 * Container for indel data on ancestors, and factory methods to populate them.
 * Workflow:
 * Before calling a factory method, create/load data in POGs and tree; place these in POGTree.
 * Call factory method, in which POGTree is queried and an instance is extracted, suited to type of inference.
 */
public class IndelPrediction {
    private final POGTree pogTree; // input data contained in a POGTree
    private final Map<Object, POGraph> ancestors; // named ancestors

    /**
     * Basic constructor, not intended for use as it forgoes/assumes prior prediction.
     * @param pogTree reference input data
     * @param ancestors predicted ancestors as map, keyed by branch point label/ID, values are POGs
     */
    public IndelPrediction(POGTree pogTree, Map<Object, POGraph> ancestors) {
        this.pogTree = pogTree;
        this.ancestors = ancestors;
    }

    /**
     * Basic inference of gaps by position-specific parsimony.
     * @param pogTree POGs and phylogenetic tree
     * @return an instance of the class, incorporating ancestor POGs, named by their phylogenetic branch point label.
     */
    public static IndelPrediction PredictByParsimony(POGTree pogTree) {
        int nPos = pogTree.getPositions(); // find the number of indices that the POGs (input and ancestors) can use
        IdxTree tree = pogTree.getTree();  // indexed tree (quick access to branch points, no editing)
        Map<Object, POGraph> ancestors = new HashMap<>(); // prepare where predictions will go
        TreeInstance[] ti = pogTree.getNodeInstances(true);   // extract gap/no-gap (boolean) leaf instantiation for every position
        Parsimony[] pi = new Parsimony[nPos]; // prepare array where all inferred gap states will go
        // next stage should be multi-threaded... not so at the moment
        for (int i = 0; i < nPos; i ++) { // for every position...
            pi[i] = new Parsimony(ti[i]); // this is where the inferred states for an individual position goes
            pi[i].forward();  // inference part 1
            pi[i].backward(); // inference part 2
        }
        // unpack the results, branch point by branch point
        for (int j = 0; j < tree.getSize(); j ++) {         // for every (indexed) branch point (IdxTree defaults to depth-first order)
            if (tree.isLeaf(j))             // if leaf, ignore and
                continue;                   // jump to next branch point
            // else: ancestor branchpoint
            Object ancID = tree.getBranchPoint(j).getID();  // unique ancestor ID from BranchPoint class (modifiable only before assembled into IdxTree when creating Tree)
            POGraph pog = new POGraph(nPos);                // blank POG
            ancestors.put(ancID, pog);                      // put blank POG in place, to be modified below
            // each anchor-set (below) is a range of indices for the POG; the key is a position that MUST be entered, values contain admissible positions that follow
            Map<Integer, Set<Integer>> anchorsets = new HashMap<>();
            int current_anchor = -1; // the first jump always from the start position
            Set<Integer> anchorset = new HashSet<>();
            for (int i = 0; i < nPos; i ++) {   // now traverse all positions, consider if GAP, not-GAP or admissible
                if (pi[i].getOptimal(j).contains(Boolean.FALSE) && pi[i].getOptimal(j).contains(Boolean.TRUE)) { // admissible, but not required
                    anchorset.add(i);
                } else if (pi[i].getOptimal(j).contains(Boolean.FALSE)) { // ALWAYS character (i.e. not-GAP), so required position
                    anchorset.add(i);
                    anchorsets.put(current_anchor, anchorset);  // anchor set is linked to the position at which it started
                    anchorset = new HashSet<>();                // re-set anchor set
                    current_anchor = i;                         // next anchor set is headed by this, not-GAP required position
                }
                // else it is a GAP only (so NOT added to admissible set, NOT an anchor)
            }
            anchorset.add(nPos);    // finish-up last anchor set
            anchorsets.put(current_anchor, anchorset);
            // next, peruse anchor sets, adding Nodes for all anchored or admissible positions
            for (Map.Entry<Integer, Set<Integer>> entry : anchorsets.entrySet()) {
                int anchor = entry.getKey();
                if (anchor >= 0)
                    pog.addNode(anchor, new Node());
                for (int to : entry.getValue()) { // link anchored Node to each of the admissible nodes in the anchor set
                    if (to < nPos)
                        pog.addNode(to, new Node());
                    pog.addEdge(anchor, to);
                }
                for (int from : entry.getValue()) { // all possible pairs of admissible Nodes are linked
                    for (int to : entry.getValue()) {
                        if (from < to)
                            pog.addEdge(from, to);
                    }
                }
            }
        }
        return new IndelPrediction(pogTree, ancestors);
    }


    /**
     * Simple Indel Code based on Simmons and Ochoterena "Gaps as characters..." 2000
     *
     * Seq.A    GG---1---CCTT------3-----GG
     * Seq.B    GG---1---CCTT------3-----GG
     * Seq.C    GGAAA---2--TT-4-AC-5--AAAGG
     * Seq.D    GGAAA---2--TT-4-AC-5--AAAGG
     * Seq.E    GGAAACCCCCCTTCAAACCCCAAAAGG
     *
     *          1 2 3 4 5
     *          ---------
     * Seq.A    1 0 1 - -
     * Seq.B    1 0 1 - -
     * Seq.C    0 1 0 1 1
     * Seq.D    0 1 0 1 1
     * Seq.E    0 0 0 0 0
     *
     * @param pogTree
     * @return instance of IndelPrediction
     */
    public static IndelPrediction PredictByIndelParsimony(POGTree pogTree) {
        int nPos = pogTree.getPositions(); //
        IdxTree tree = pogTree.getTree();
        Map<Object, POGraph> ancestors = new HashMap<>();
        // Retrieve an instance for each indel across the whole alignment (ordered by interval tree)
        // (state for each etant-indel: present, absent or permissible/neutral, as per Simmons and Ochoterena, 2000)
        TreeInstance[] ti = pogTree.getIndelInstances();
        Parsimony[] pi = new Parsimony[ti.length];
        // next stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) { // for each "indel"
            pi[i] = new Parsimony(ti[i]);
            pi[i].forward();
            pi[i].backward();

        }
        // extra table
        int i = 0; // interval index
        System.out.println("Indels---------");
        for (Interval1D ival : pogTree.getIntervalTree())
            System.out.println(i++ + "\t" + ival);
            // now decorate the ancestors
        System.out.println("Sequences---------");
        i = 0; // interval index
        for (Interval1D ival : pogTree.getIntervalTree())
            System.out.print("\t" + i++);
        System.out.println();
        for (int j = 0; j < tree.getSize(); j++) {
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length == 0) { // not an ancestor
                POGraph pog = pogTree.getExtant(ancID);
                System.out.print(ancID + "\t");
                if (pog != null) {
                    i = 0;
                    for (Interval1D ival : pogTree.getIntervalTree()) {
                        StringBuilder sb = new StringBuilder();
                        List<Boolean> calls = pi[i].getOptimal(j);
                        for (Boolean b : calls)
                            sb.append(b.toString().substring(0,1));
                        System.out.print(sb.toString() + "\t");
                        i ++;
                    }
                }
                System.out.println();
                continue; // skip the code below
            }
            // else: ancestor branch point
            // (1) Find ancestor STATE for each INDEL
            // all unambiguously TRUE indels, as they will preclude "contained" indels
            Set<Interval1D> definitively_true = new HashSet<>();
            // all INDELs which are OPTIONAL
            Set<Interval1D> optionally_true = new HashSet<>();
            System.out.print(ancID + "\t");
            i = 0; // interval index
            for (Interval1D ival : pogTree.getIntervalTree()) {
                StringBuilder sb = new StringBuilder();
                List<Boolean> calls = pi[i].getOptimal(j);
                for (Boolean b : calls)
                    sb.append(b.toString().substring(0,1));
                System.out.print(sb.toString() + "\t");
                if (calls.contains(Boolean.TRUE)) { // INDEL can be TRUE
                    optionally_true.add(ival);
                    if (calls.size() == 1) // the ONLY value is TRUE so DEFINITIVELY so
                        definitively_true.add(ival);
                }
                i++;
            }
            System.out.println();
            int[][] edges_optional = new int[optionally_true.size()][2];
            i = 0;
            for (Interval1D edge : optionally_true)
                edges_optional[i++] = new int[]{edge.min, edge.max};
            // Second, use only indels that are not contained within an unambiguously TRUE indel
            Set<Interval1D> discarded = new HashSet<>();
            i = 0; // interval index
            for (Interval1D ival : pogTree.getIntervalTree()) {
                int count_contains = 0;
                for (Interval1D definitive : definitively_true)
                    count_contains += definitive.contains(ival) ? 1 : 0;
                if (count_contains > 1)
                    discarded.add(ival);
                i += 1;
            }
            definitively_true.removeAll(discarded);
            int[][] edges_definitive = new int[definitively_true.size()][];
            i = 0;
            for (Interval1D edge : definitively_true)
                edges_definitive[i++] = new int[]{edge.min, edge.max};
            POGraph pog = POGraph.createFromEdgeIndicesWithoutDeadends(nPos, edges_optional, edges_definitive);
            ancestors.put(ancID, pog);
        }
        return new IndelPrediction(pogTree, ancestors);
    }

    public static IndelPrediction PredictBySimpleIndelParsimony(POGTree pogTree) {
            int nPos = pogTree.getPositions(); //
            IdxTree tree = pogTree.getTree();
            Map<Object, POGraph> ancestors = new HashMap<>();
            TreeInstance[] ti = pogTree.getIndelInstances();
            Parsimony[] pi = new Parsimony[ti.length];
            // next stage should be multi-threaded... not so at the moment
            for (int i = 0; i < ti.length; i ++) { // for each "indel"
                pi[i] = new Parsimony(ti[i]);
                pi[i].forward();
                pi[i].backward();
            }
        // establish containment relationships between indels
        Map<Integer, Set<Integer>> contain_map = new HashMap<>();
        int ibig = 0;
        for (Interval1D ival_big : pogTree.getIntervalTree()) {
            Set<Integer> small_ivals = new HashSet<>();
            contain_map.put(ibig, small_ivals);
            int ismall = 0;
            for (Interval1D ival_small : pogTree.getIntervalTree()) {
                if (!ival_big.equals(ival_small)) { // could instead check ibig != ismall
                    if (ival_big.contains(ival_small))
                        small_ivals.add(ismall);
                }
                ismall += 1;
            }
            ibig += 1;
        }
        // now decorate the ancestors
        for (int j = 0; j < tree.getSize(); j ++) {
            if (tree.getChildren(j).length == 0) // not an ancestor
                continue; // skip the code below
            // else: ancestor branchpoint
            Object ancID = tree.getBranchPoint(j).getID();
            POGraph pog = new POGraph(nPos);
            ancestors.put(ancID, pog);

            // First identify all unambiguously TRUE indels, as they will preclude "contained" indels
            Set<Interval1D> definitively_true = new HashSet<>();
            int i = 0; // interval index
            for (Interval1D ival : pogTree.getIntervalTree()) {
                List<Boolean> calls = pi[i].getOptimal(j);
                if (calls.size() == 1 && calls.contains(Boolean.TRUE)) { // has ONLY gap
                    definitively_true.add(ival);
                }
                i ++;
            }
            // Second, use only indels that are not contained within an unambiguously TRUE indel
            Set<Interval1D> discarded = new HashSet<>();
            i = 0; // interval index
            for (Interval1D ival : pogTree.getIntervalTree()) {
                int count_contains = 0;
                for (Interval1D definitive : definitively_true)
                    count_contains += definitive.contains(ival) ? 1 : 0;
                if (count_contains > 1)
                    discarded.add(ival);
                i += 1;
            }
            definitively_true.removeAll(discarded);
            // Third, distinguish predictions that have not been precluded, and build POG
            i = 0; // interval index
            if (ancID.toString().equals("10"))
                i = 0;
            for (Interval1D ival : pogTree.getIntervalTree()) {
                boolean precluded = false; // assume this indel is not precluded (by containment)
                if (!definitively_true.contains(ival)) { // the indel is NOT one that is unambiguously optimal
                    for (Interval1D definitive : definitively_true) { // check if it is contained and therefore precluded by one
                        if (definitive.contains(ival)) {
                            precluded = true;
                            break;
                        }
                    }
                }
                if (!precluded) {
                    List<Boolean> calls = pi[i].getOptimal(j);
                    if (calls.contains(Boolean.TRUE)) { // has gap
                        if (pog.isIndex(ival.min))
                            pog.addNode(ival.min, new Node());
                        if (pog.isIndex(ival.max))
                            pog.addNode(ival.max, new Node());
                        pog.addEdge(ival.min, ival.max);
                    }
                    if (calls.contains(Boolean.FALSE)) { // "contained" indels that are TRUE should be incorporated
                        for (int ismall : contain_map.get(i)) {
                            List<Boolean> calls_small = pi[ismall].getOptimal(j);
                            if (calls_small.contains(Boolean.TRUE)) { // has gap
                                if (pog.isIndex(ival.min))
                                    pog.addNode(ival.min, new Node());
                                if (pog.isIndex(ival.max))
                                    pog.addNode(ival.max, new Node());
                                pog.addEdge(ival.min, ival.max);
                            }
                        }
                    }
                }
                i ++;
            }
        }
        return new IndelPrediction(pogTree, ancestors);
    }

    /**
     * Find the TRUE indels that will complete this region
     * @param start
     * @param end
     * @return
     */
    private Set<Interval1D> patchThisRegion(Set<Interval1D> optionally_true, int start, int end) {
        // FIXME: not implemented
        return null;
    }

    public POGraph getAncestor(Object ancID) {
        POGraph ancestor = ancestors.get(ancID);
        return ancestor;
    }
}
