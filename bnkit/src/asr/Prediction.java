package asr;

import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.GLOOME1;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.Interval1D;
import dat.phylo.IdxTree;
import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;
import dat.pog.*;

import java.util.*;

/**
 * Container class for
 * - data on phylogenetic tree and ancestors, linked by name and position, to indexed, position specific instance trees
 * - factory methods to predict indels (required to instantiate this class)
 * - constructor explicitly linking ancestor POGs to trees
 * - methods to infer character content
 *
 * Workflow:
 * 1.  Create/load data in (extant) POGs/alignments and tree; place these in POGTree.
 * 2.  Call factory methods to infer insertions and deletions (selection of algorithms is available), in which POGTree is queried and template POGs are created
 * 3a. Set evolutionary model
 * 3. Call methods to decorate ancestors with character states (joint reconstruction) and/or probability distributions (marginal reconstruction)
 *
 */
public class Prediction {
    public static boolean DEBUG = GRASP.VERBOSE;    // print out various information
    private final POGTree pogTree;                  // input data contained in a POGTree
    private final IdxTree phylotree;                // the input phylogenetic tree (is also accessible via POGTree)
    private int[] ancidxs = null;                   // store the sub-set of indices that are used for internal nodes (ancestors)
    //private final Map<Object, POGraph> ancestors; // named ancestors
    private final POGraph[] ancarr;                 // ancestors by branchpoint index
    private final IdxTree[] positrees;              // position-specific tree, or an edited form of the original tree for the purpose of inferring content
    // Indexing of ancestors by ID to branch point index in phylogenetic tree and position-specific trees
    private final int[][] positidxs;                // position-specific tree indices [aligned pos]["global" branchpoint idx]
    private final EnumDistrib[][] distribs;         // Probability distributions of ancestor states by branchpoint and position index
    private Object[][] states = null;               // Actual values at inferred branch points, indexed by branchpoint index and position

    /**
     * Basic constructor, not intended for use as it forgoes/assumes prior prediction.
     * @param pogTree reference input data
     * @param ancestors predicted ancestors as map, keyed by branch point label/ID, values are POGs
     */
    public Prediction(POGTree pogTree, Map<Object, POGraph> ancestors) {
        this.pogTree = pogTree;
        this.phylotree = pogTree.getTree();
        //this.ancestors = ancestors;
        this.ancarr = new POGraph[pogTree.getTree().getSize()];
        for (Map.Entry<Object, POGraph> entry : ancestors.entrySet()) { // make sure names are set
            entry.getValue().setName(entry.getKey().toString());
            this.ancarr[getBranchpointIndex(entry.getKey())] = entry.getValue();
        }
        this.positrees = new IdxTree[pogTree.getPositions()]; // by default there's one tree for each index in the alignment/POG
        this.positidxs = new int[pogTree.getPositions()][];   // by default there's one tree for each index in the alignment/POG
        this.distribs = new EnumDistrib[phylotree.getSize()][];
    }

    /**
     * Retrieve number of positions in reconstruction, which defines the bounds of indices that can be accessed in any of the composite POGs.
     * @return number of positions in reconstruction
     */
    public int getPositions() {
        return pogTree.getPositions();
    }

    private int[] getAncestorIndices() {
        if (ancidxs == null)
            ancidxs = phylotree.getAncestors();
        return ancidxs;
    }

    /**
     * Retrieve the original phylogenetic tree, unspecific to position.
     * @return index tree
     */
    public IdxTree getTree() {
        return phylotree;
    }

    /**
     * Retrieve position-specific index tree that is based on the original phylogenetic tree but has only indices for ancestor nodes,
     * which are not marked as absent indels.
     * The index tree is cached inside this class instance, so that it can be quickly retrieved when required again.
     * @param position index in alignment/POG
     * @return index tree for specified position in alignment/POG
     */
    public IdxTree getTree(int position) {
        if (positrees[position] == null) {
            // we need to construct this from ancestor POGs...
            IdxTree phylo = pogTree.getTree();
            Set<Integer> pruneMe = new HashSet<>();
            for (int idx : phylo) {
                Object ancid = phylo.getLabel(idx);
                if (!phylo.isLeaf(idx)) { // ancestor
                    POGraph pog = getAncestor(ancid);
                    if (pog == null)
                        throw new ASRRuntimeException("Invalid ancestor ID " + ancid);
                    if (!pog.isNode(position))
                        pruneMe.add(idx);
                }
            }
            int[] indices = phylo.getPrunedIndex(pruneMe);
            // save indices for quick re-retrieval later
            positidxs[position] = indices;
            // save tree for quick re-retrieval later
            positrees[position] = IdxTree.createPrunedTree(phylo, indices);
        } // else the index tree was already cached...
        return positrees[position];
    }

    /**
     * Access what ancestors that are available
     * @return
    public Map<Object, POGraph> getAncestors() {
        return ancestors;
    }
     */

    /**
     * Map identifier to a branch point index for direct access to content arrays
     * @param id
     * @return branch point index, or -1 if not found
     */
    public int getBranchpointIndex(Object id) {
        return pogTree.getTree().getIndex(id);
    }

    /**
     * Get the ancestor POG for a particular ancestor ID
     * @param ancID
     * @return
     */
    public POGraph getAncestor(Object ancID) {
        int bpidx = getBranchpointIndex(ancID);
        return bpidx != -1 ? ancarr[bpidx] : null;
    }

    /**
     * Get the extant POG for a particular extant sequence ID
     * @param extID sequence name
     * @return the POG, null if not available
     */
    public POGraph getExtant(Object extID) {
        return pogTree.getExtant(extID);
    }

    /**
     * Retrieve all predictions by instantiating POGs by the requested type (e.g. states for joint reconstructions, distributions for marginal).
     * Ancestors that have not been inferred are NOT included in the returned map.
     * @param mode inference mode, currently JOINT and MARGINAL are supported (see enum defined in GRASP class)
     * @return map from identifier to POGraph
     */
    public Map<Object, POGraph> getAncestors(GRASP.Inference mode) {
        Map<Object, POGraph> ancestors = new HashMap<>();
        // iterate through all ancestors, and extracting states from those that have been inferred
        for (int idx : getAncestorIndices()) {
            POGraph pog = ancarr[idx];
            if (mode == GRASP.Inference.JOINT) {
                if (states != null)
                    pog.decorateNodes(states[idx]);
            } else if (mode == GRASP.Inference.MARGINAL) {
                if (distribs[idx] != null)
                    pog.decorateNodes(distribs[idx]);
            }
            ancestors.put(phylotree.getLabel(idx), pog);
        }
        return ancestors;
    }

    /**
     * Retrieve one prediction by instantiating the POG with the requested inference type (e.g. states for joint reconstructions, distributions for marginal)
     * @param ancID identifier of ancestor
     * @param mode inference mode, currently JOINT and MARGINAL are supported (see enum defined in GRASP class)
     * @return requested POG, or null if not available
     */
    public POGraph getAncestor(Object ancID, GRASP.Inference mode) {
        int bpidx = getBranchpointIndex(ancID);
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID (not found in tree) " + ancID);
        POGraph pog0 = ancarr[bpidx];
        if (pog0 == null)
            throw new ASRRuntimeException("Invalid ancestor ID (not inferred) " + ancID);
        if (mode == GRASP.Inference.JOINT) {
            pog0.decorateNodes(states[bpidx]);
        } else if (mode == GRASP.Inference.MARGINAL) {
            pog0.decorateNodes(distribs[bpidx]);
        }
        return pog0;
    }

    /**
     * Get the marginal distributions for a specified ancestor and substitution model.
     * Note that distributions are null for positions that are not part of the ancestor (i.e. the result of a deletion in an earlier ancestor,
     * or parts that precede an insertion in a later ancestor)
     * @param ancestorID the ancestor ID
     * @param MODEL the substitution model
     * @return
     */
    public EnumDistrib[] getMarginal(Object ancestorID, SubstModel MODEL) {
        int bpidx = getBranchpointIndex(ancestorID);                            // the index of the ancestor as it appears in the phylogenetic tree
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID (not found in tree) " + ancestorID);
        if (distribs[bpidx] == null) {                                          // the ancestor has not yet been inferred, so DO it...
            IdxTree[] trees = new IdxTree[getPositions()];                      // this is how many position-specific trees we are dealing with
            MaxLhood.Marginal[] inf = new MaxLhood.Marginal[getPositions()];    // which is also how many inferences we will carry out
            for (int pos = 0; pos < getPositions(); pos ++) {                   // for each position...
                trees[pos] = getTree(pos);                                      //   this is the tree with indels imputed
                int ancidx = positidxs[pos][bpidx];                             //   index for sought ancestor in the position-specific tree
                if (ancidx >= 0)                                                //   which may not exist, i.e. part of an indel, but if it is real...
                    inf[pos] = new MaxLhood.Marginal(ancidx, trees[pos], MODEL);//     set-up the inference
            }
            distribs[bpidx] = new EnumDistrib[pogTree.getPositions()];
            for (int pos = 0; pos < getPositions(); pos ++) {                   // for each position...
                int specidx = positidxs[pos][bpidx];                            //   index for sought ancestor in the position-specific tree
                if (specidx >= 0) {                                             //   which may not exist, i.e. part of an indel, but if it is real...
                    TreeInstance ti = pogTree.getNodeInstance(pos, trees[pos]); //     get the instances at the leaves at that position, and...
                    inf[pos].decorate(ti);                                      //     perform inference
                    distribs[bpidx][pos] = inf[pos].getDecoration(specidx);     //     extract distribution of marginal prob
                }
            }
        }
        // at this point we know the ancestor has been inferred
        return distribs[bpidx];
    }

    /**
     * Perform joint reconstruction across all ancestors, and all positions
     * @param MODEL evolutionary model
     * @return the states at all ancestors that assign the greatest likelihood to the observed states at extant sequences
     */
    public Object[][] getJoint(SubstModel MODEL) {
        this.states = new Object[getTree().getSize()][getPositions()];
        IdxTree[] trees = new IdxTree[getPositions()];              // this is how many position-specific trees we are dealing with
        MaxLhood.Joint[] inf = new MaxLhood.Joint[getPositions()];  // which is also how many inferences we will carry out
        for (int pos = 0; pos < inf.length; pos++) {                // for each position...
            trees[pos] = getTree(pos);                              //   this is the tree with indels imputed
            inf[pos] = new MaxLhood.Joint(trees[pos], MODEL);       //     set-up the inference
        }
        for (int idx : getAncestorIndices()) {
//      for  (Map.Entry<Object, POGraph> entry : getAncestors().entrySet()) {
            // Object ancID = entry.getKey();
            for (int pos = 0; pos < getPositions(); pos ++) {       // for each position...
                //int idx = pogTree.getTree().getIndex(ancID);        //   index for sought ancestor in the phylogenetic tree
                int ancidx = positidxs[pos][idx];                   //   index for sought ancestor in the position-specific tree
                if (ancidx >= 0) {                                  //   which may not exist, i.e. part of an indel, but if it is real...
                    TreeInstance ti = pogTree.getNodeInstance(pos, trees[pos]); //     get the instances at the leaves at that position, and...
                    inf[pos].decorate(ti);                                      //     perform inference
                    states[idx][pos] = inf[pos].getDecoration(ancidx);          //     extract state
                }
            }
        }
        return states;
    }

    /**
     * Retrieve the joint reconstruction for an ancestor; this is one ancestor out of many which need to be computed jointly.
     * If none has been determined before, all ancestors will be computed.
     * @param ancestorID the ancestor
     * @param MODEL evolutionary model
     * @return the states at the ancestor (together with all the others) that assign the greatest likelihood to the observed states at extant sequences
     */
    public Object[] getJoint(Object ancestorID, SubstModel MODEL) {
        if (states == null)   // the ancestors has not yet been inferred
            getJoint(MODEL);
        int bpidx = getBranchpointIndex(ancestorID);                            // the index of the ancestor as it appears in the phylogenetic tree
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID (not found in tree) " + ancestorID);
        return states[bpidx];
    }

    /**
     * Retrieve the inferred ancestor but as a sequence; requires that it has been inferred already
     * @param ancID ancestor ID
     * @param mode inference mode, currently JOINT and MARGINAL are supported (see enum defined in GRASP class)
     * @param gappy whether gaps should be included (as they appear from the POG index)
     * @return the sequence
     */
    public EnumSeq getSequence(Object ancID, GRASP.Inference mode, boolean gappy) {
        int bpidx = getBranchpointIndex(ancID);                            // the index of the ancestor as it appears in the phylogenetic tree
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID: " + ancID);
        int[] idxs = getConsensus(bpidx);
        int N = getPositions();
        Object[] elems = new Object[gappy ? N : idxs.length];
        EnumSeq seq = gappy ? new EnumSeq.Gappy(pogTree.getDomain()) : new EnumSeq(pogTree.getDomain());
        if (mode == GRASP.Inference.JOINT) {
            if (states == null) // not inferred yet
                throw new ASRRuntimeException("Joint inference has not been performed: " + ancID);
            for (int i = 0; i < idxs.length; i ++)
                elems[gappy ? idxs[i] : i] = states[bpidx][idxs[i]];
        } else if (mode == GRASP.Inference.MARGINAL) {
            if (distribs[bpidx] == null) // not inferred yet
                throw new ASRRuntimeException("Marginal inference has not been performed: " + ancID);
            for (int i = 0; i < idxs.length; i ++)
                elems[gappy ? idxs[i] : i] = distribs[bpidx][idxs[i]].getMax();
        }
        seq.set(elems);
        seq.setName(ancID.toString());
        return seq;
    }

    public int[] getConsensus(Object ancID) {
        int bpidx = getBranchpointIndex(ancID);                            // the index of the ancestor as it appears in the phylogenetic tree
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID: " + ancID);
        return getConsensus(bpidx);
    }

    public int[] getConsensus(int bpidx) {
        POGraph pog = this.ancarr[bpidx];
        if (pog == null)
            throw new ASRRuntimeException("Ancestor has not been inferred: index is " + bpidx);
        // collect info to make decisions...
        int[] leaves = phylotree.getLeaves(bpidx); // determine all branch points of leaves (i.e. extants) under this ancestor
        // go through POG nodes, to determine the transition "weights"
        PriorityQueue<Integer> queue = new PriorityQueue<>();
        for (int idx : pog.getForward()) // to start us off: add all indices emanating from start
            queue.add(idx);
        // iterate through a queue, to which nodes are added if linked from "current" node
        while (queue.size() > 0) { // until empty...
            int curr = queue.poll();
            int[] nexts = pog.getForward(curr);
            if (pog.isEndNode(curr)) { // check if next node can be terminal
                // if so, add to nexts
                int[] nnexts = nexts;
                nexts = new int[nnexts.length + 1];
                for (int i = 0; i < nnexts.length; i ++)
                    nexts[i] = nnexts[i];
                nexts[nnexts.length] = pog.size(); // add the end terminus
            }
            double[] rates = pogTree.getEdgeRates(curr, nexts, leaves);
            // set the weights
            for (int i = 0; i < nexts.length; i ++) {
                POGraph.StatusEdge edge = pog.getEdge(curr, nexts[i]);
                edge.setWeight(-Math.log(rates[i])); // neg log of prob; so P=1 means zero weight, low P means high weight
            }
            // add indices of all nodes that can be visited next
            for (int idx : nexts)
                if (idx != pog.maxsize())
                    queue.add(idx);
        }
        // with weights set, we find the optimal path (pick one if several; else need to query individual edges and assemble)
        int[] consensus = pog.getMostSupported();
        return consensus;
    }

    /**
     * Basic inference of gaps by position-specific parsimony.
     * @param pogTree POGs and phylogenetic tree
     * @return an instance of the class, incorporating ancestor POGs, named by their phylogenetic branch point label.
     */
    public static Prediction PredictByParsimony(POGTree pogTree) {
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
        return new Prediction(pogTree, ancestors);
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
    public static Prediction PredictByIndelParsimony(POGTree pogTree) {
        int nPos = pogTree.getPositions(); //
        IdxTree tree = pogTree.getTree();
        Map<Object, POGraph> ancestors = new HashMap<>();
        // Retrieve an instance for each indel across the whole alignment (ordered by interval tree)
        // (state for each extant-indel: present, absent or permissible/neutral, as per Simmons and Ochoterena, 2000)
        // initially "permissible/neutral" is encoded as null (the variable is uninstantiated); see POGraph.getSimpleGapCode)
        // inference will infer true, false, or accept that both true and false can be correct
        TreeInstance[] ti = pogTree.getIndelInstances();
        Parsimony[] pi = new Parsimony[ti.length];
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) { // for each "indel"
            pi[i] = new Parsimony(ti[i]);
            pi[i].forward();
            pi[i].backward();
        }
        if (DEBUG) {
            // print out tables...
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
        }
        // the code below
        // (1) regardless, if an ancestor or extant, we can pull out what the INDEL states are: absent (false), present (true) or permissible (true/false)
        // (2) if an ancestor, an ancestor POG is created, using the info from (1)
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branch point, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length == 0) { // not an ancestor
                if (DEBUG) {
                    POGraph pog = pogTree.getExtant(ancID);
                    System.out.print(ancID + "\t");
                    if (pog != null) {
                        int i = 0;
                        for (Interval1D ival : pogTree.getIntervalTree()) {
                            StringBuilder sb = new StringBuilder();
                            List calls = pi[i].getOptimal(j);
                            for (Object b : calls) // each "b" is a Boolean
                                sb.append(b.toString().substring(0, 1)); // this converts each value to "t" or "f"
                            System.out.print(sb.toString() + "\t");
                            i++;
                        }
                    }
                    System.out.println();
                }
                continue; // skip the code below, only predictions for ancestors are used to compose POGs
            }
            // else: ancestor branch point
            // (1) Find ancestor STATE for each INDEL, but to resolve (2) what the POG looks like we determine...
            // all unambiguously TRUE indels, as they will preclude "contained" indels
            Set<Interval1D> definitively_true = new HashSet<>();
            // all INDELs which are OPTIONAL
            Set<Interval1D> optionally_true = new HashSet<>();
            if (DEBUG) System.out.print(ancID + "\t");
            int i = 0; // interval index
            for (Interval1D ival : pogTree.getIntervalTree()) {
                List<Boolean> calls = pi[i].getOptimal(j);
                if (DEBUG) {
                    StringBuilder sb = new StringBuilder();
                    for (Boolean b : calls)
                        sb.append(b.toString().substring(0, 1));
                    System.out.print(sb.toString() + "\t");
                }
                if (calls.contains(Boolean.TRUE)) { // INDEL can be TRUE
                    optionally_true.add(ival);
                    if (calls.size() == 1) // the ONLY value is TRUE so DEFINITIVELY so
                        definitively_true.add(ival);
                }
                i ++;
            }
            if (DEBUG) System.out.println();
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
            // finally put the info into a POG
            POGraph pog = POGraph.createFromEdgeIndicesWithoutDeadends(nPos, edges_optional, edges_definitive);
            ancestors.put(ancID, pog);
        }
        return new Prediction(pogTree, ancestors);
    }

    /**
     * Simple Indel Code based on Simmons and Ochoterena "Gaps as characters..." (2000).
     * Inference with ML, using a gain/loss model similar to FastML (based on Cohen and Pupko 2010 M1).
     * Note that FastML uses a mixture of two models (termed M1 and M2) but both have loss rates which dominate gain.
     * @param pogTree
     * @return instance of IndelPrediction
     */
    public static Prediction PredictByIndelMaxLhood(POGTree pogTree) {
        return Prediction.PredictByIndelMaxLhood(pogTree, new GLOOME1());
    }

    /**
     * Simple Indel Code based on Simmons and Ochoterena "Gaps as characters..." (2000).
     * Inference with ML.
     * @param pogTree
     * @param gain_loss_model model for gain and loss events
     * @return instance of IndelPrediction
     */
    public static Prediction PredictByIndelMaxLhood(POGTree pogTree, SubstModel gain_loss_model) {
        int nPos = pogTree.getPositions(); // number of positions in alignment/POG
        IdxTree tree = pogTree.getTree();  // the tree that represents phylogenetic relationships between extants
        Map<Object, POGraph> ancestors = new HashMap<>();
        // Retrieve an instance for each indel across the whole alignment (ordered by interval tree).
        // This is a "deletion" as it skips positions. At the moment, zero-distance skips are also included, since when viewed
        // across a full sequence they may correspond to "gaps". (There is a potential to resolve some of them without inference though.)
        // State for each extant-indel: present, absent or permissible/neutral, as per Simmons and Ochoterena (2000).
        // Initially "permissible/neutral" is encoded as null (the variable is uninstantiated); see POGraph.getSimpleGapCode.
        // ML inference will infer the most likely combination of true and false across the full tree;
        // that is, maximise the probability of the observed extant states GIVEN the combination of the indel states at the ancestors.
        // We use a "joint" inference approach, which is NOT what FastML does; Ashkenazy et al (2012) reports that FastML outputs the
        // "posterior probability for each indel site at each ancestral node of the phylogeny.
        // Most likely character states in the ancestral nodes are reported only in positions
        // that are inferred to be non-gapped with a probability 􏰀0.5."
        // To do this would require a switch to marginal inference, then thresholding for 0.5.
        TreeInstance[] ti = pogTree.getIndelInstances(); // instantiate a tree for each "indel", assigning leaf states as per extants
        MaxLhood.Joint[] ji = new MaxLhood.Joint[ti.length];
        for (int i = 0; i < ji.length; i++) // for each "indel" we need to infer either gain or loss, so set-up inference
            ji[i] = new MaxLhood.Joint(tree, gain_loss_model);
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) { // for each "indel"
            ji[i].decorate(ti[i]);
        }
        if (DEBUG) {
            // print out tables...
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
        }
        // the code below
        // (1) regardless, if an ancestor or extant, we can pull out what the INDEL states are:
        // gain (no "skip", encoded false) or loss ("skip", encoded true)
        // (2) if an ancestor, an ancestor POG is created, using the info from (1)
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branchpoint, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length == 0) { // not an ancestor, so we just print out debug info below before continuing with next
                if (DEBUG) {
                    POGraph pog = pogTree.getExtant(ancID);
                    System.out.print(ancID + "\t");
                    if (pog != null) {
                        int i = 0;
                        for (Interval1D ival : pogTree.getIntervalTree()) {
                            StringBuilder sb = new StringBuilder();
                            Boolean call = (Boolean)ji[i].getDecoration(j);
                            sb.append(call ? "L" : "G"); // Loss and Gain, respectively
                            System.out.print(sb.toString() + "\t");
                            i++;
                        }
                    }
                    System.out.println();
                }
                continue; // skip the code below, only predictions for ancestors are used to compose POGs
            }
            // else: ancestor branch point
            // (1) Find ancestor STATE for each INDEL, but to resolve (2) what the POG looks like we determine...
            // all TRUE indels ("skips"), as they will preclude "contained" indels;
            // there is provision BUT this is not needed in this implementation for OPTIONAL skips
            // TODO: The code below is based on the code for parsimony inference may deal with situations that cannot occur
            Set<Interval1D> definitively_true = new HashSet<>();
            // all INDELs which are OPTIONAL
            Set<Interval1D> optionally_true = new HashSet<>();
            if (DEBUG) System.out.print(ancID + "\t");
            int i = 0; // interval index
            for (Interval1D ival : pogTree.getIntervalTree()) {
                Boolean call = (Boolean)ji[i].getDecoration(j);
                if (DEBUG) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(call ? "L" : "G"); // Loss and Gain, respectively
                    System.out.print(sb.toString() + "\t");
                }
                if (call) {
                    optionally_true.add(ival);
                    definitively_true.add(ival);
                }
                i ++;
            }
            if (DEBUG) System.out.println();
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
            // finally put the info into a POG
            POGraph pog = POGraph.createFromEdgeIndicesWithoutDeadends(nPos, edges_optional, edges_definitive);
            ancestors.put(ancID, pog);
        }
        return new Prediction(pogTree, ancestors);
    }

    /**
     * Bi-directional edge parsimony for inference of indel states in ancestor POGs.
     * @param pogTree
     * @return instance of IndelPrediction
     */
    public static Prediction PredictByBidirEdgeParsimony(POGTree pogTree) {
        int nPos = pogTree.getPositions(); //
        IdxTree tree = pogTree.getTree();
        Map<Object, POGraph> ancestors = new HashMap<>();
        // Retrieve an instance for each indel across the whole alignment (ordered by interval tree)
        // (state for each extant-indel: present, absent or permissible/neutral, as per Simmons and Ochoterena, 2000)
        // initially "permissible/neutral" is encoded as null (the variable is uninstantiated); see POGraph.getSimpleGapCode)
        // inference will infer true, false, or accept that both true and false can be correct
        TreeInstance[] tif = new TreeInstance[nPos + 2]; // forward
        TreeInstance[] tib = new TreeInstance[nPos + 2]; // backward
        for (int i = -1; i <= nPos; i ++) {
            tif[i+1] = pogTree.getEdgeInstance(i, POGTree.EDGE_FORWARD);
            tib[i+1] = pogTree.getEdgeInstance(i, POGTree.EDGE_BACKWARD);
        }
        Parsimony[] pif = new Parsimony[tif.length];
        Parsimony[] pib = new Parsimony[tib.length];
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = -1; i <= nPos; i++) { // for each position, FIXME: some unnecessary inferences at termini below
            pif[i+1] = new Parsimony(tif[i+1]);
            pib[i+1] = new Parsimony(tib[i+1]);
            pif[i+1].forward();
            pif[i+1].backward();
            pib[i+1].forward();
            pib[i+1].backward();
        }
        // inference done, now assemble...
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branchpoint, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length > 0) { // an ancestor
                EdgeMap emap = new EdgeMap();
                for (int i = -1; i <= nPos; i ++) {
                    if (i != nPos) {
                        List solutsf = pif[i + 1].getOptimal(j);
                        for (Object s : solutsf) {
                            int next = ((Integer) s).intValue();
                            emap.add(i, next);
                        }
                    }
                    if (i != -1) {
                        List solutsb = pib[i + 1].getOptimal(j);
                        for (Object s : solutsb) {
                            int prev = ((Integer) s).intValue();
                            emap.add(prev, i);
                        }
                    }
                }
                POGraph pog = POGraph.createFromEdgeMap(nPos, emap);
                ancestors.put(ancID, pog);
            }
        }
        return new Prediction(pogTree, ancestors);
    }

    /**
     * Bi-directional edge max likelihood inference of indel states in ancestor POGs.
     * @param pogTree
     * @return instance of Prediction
     */
    public static Prediction PredictByBidirEdgeMaxLHood(POGTree pogTree) {
        int nPos = pogTree.getPositions(); //
        IdxTree tree = pogTree.getTree();
        Map<Object, POGraph> ancestors = new HashMap<>();
        TreeInstance[] tif = new TreeInstance[nPos + 2]; // forward
        TreeInstance[] tib = new TreeInstance[nPos + 2]; // backward
        for (int i = -1; i <= nPos; i ++) {
            tif[i+1] = pogTree.getEdgeInstance(i, POGTree.EDGE_FORWARD);
            tib[i+1] = pogTree.getEdgeInstance(i, POGTree.EDGE_BACKWARD);
        }
        // ML inference will maximise the JOINT probability of the observed extant edge states GIVEN the combination of the edge states at the ancestors.
        TreeDecor[] jif = new TreeDecor[tif.length];
        TreeDecor[] jib = new TreeDecor[tib.length];
        for (int i = 0; i < jif.length; i++) {
            Object[] possible = tif[i].getPossible();
            if (possible.length < 1) { // nothing to infer, not used at all
                jif[i] = null;
            } else if (possible.length < 2) { // nothing to infer, can only take one value, so create blanket output
                jif[i] = new TreeInstance.BlanketTreeDecor(tif[i].getSize(), possible[0]);
            } else {
                SubstModel substmodel = new JC(1, possible); // need to know the alphabet...
                jif[i] = new MaxLhood.Joint(tree, substmodel);
            }
        }
        for (int i = 0; i < jib.length; i++) {
            Object[] possible = tib[i].getPossible();
            if (possible.length < 1) { // nothing to infer, not used at all
                jib[i] = null;
            } else if (possible.length < 2) { // nothing to infer, can only take one value
                jib[i] = new TreeInstance.BlanketTreeDecor(tib[i].getSize(), possible[0]);
            } else {
                SubstModel substmodel = new JC(1, possible); // need to know the alphabet...
                jib[i] = new MaxLhood.Joint(tree, substmodel);
            }
        }
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < tif.length; i++) {
            if (jif[i] != null)
                jif[i].decorate(tif[i]);
        }
        for (int i = 0; i < tib.length; i++) {
            if (jib[i] != null)
                jib[i].decorate(tib[i]);
        }
        // inference done, now assemble...
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branchpoint, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length > 0) { // an ancestor
                EdgeMap emap = new EdgeMap();
                for (int i = -1; i <= nPos; i ++) {
                    if (i != nPos) {
                        Object solutsf = jif[i + 1].getDecoration(j);
                        int next = ((Integer) solutsf).intValue();
                        emap.add(i, next);
                    }
                    if (i != -1) {
                        Object solutsb = jib[i + 1].getDecoration(j);
                        int prev = ((Integer) solutsb).intValue();
                        emap.add(prev, i);
                    }
                }
                POGraph pog = POGraph.createFromEdgeMap(nPos, emap);
                ancestors.put(ancID, pog);
            }
        }
        return new Prediction(pogTree, ancestors);
    }

}
