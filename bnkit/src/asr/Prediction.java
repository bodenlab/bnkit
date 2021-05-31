package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.GLOOME1;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Interval1D;
import dat.IntervalST;
import dat.file.Newick;
import dat.phylo.IdxTree;
import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;
import dat.pog.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

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
    public int NTHREADS = 4;                        // how many threads to utilise
    private final POGTree pogTree;                  // input data contained in a POGTree
    private final IdxTree phylotree;                // the input phylogenetic tree (is also accessible via POGTree)
    private int[] ancidxs = null;                   // store the sub-set of indices that are used for internal nodes (ancestors)
    //private final Map<Object, POGraph> ancestors; // named ancestors
    private final POGraph[] ancarr;                 // ancestors by branchpoint index
    private final IdxTree[] positrees;              // position-specific tree, or an edited form of the original tree for the purpose of inferring content
    private TreeInstance[]  treeinstances;          // position-specific tree instances, which contain instantiated (by extants) and inferred content (duplicating the content in POGs)
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

    /**
     * Retrieve all indices in the original phylogenetic tree that represent ancestors
     * @return ancestor branch point indices
     */
    private int[] getAncestorIndices() {
        if (ancidxs == null)
            ancidxs = phylotree.getAncestors();
        return ancidxs;
    }

    /**
     * Map a local, position-specific branch point index to the global index, applicable in the original phylogenetic tree.
     * Computationally costly...
     * @param pos position in POG/alignment
     * @param local_idx the index at the given position
     * @return global index, or -1 if the local index is not used
     */
    private int local2global(int pos, int local_idx) {
        for (int global_idx = 0; global_idx < positidxs[pos].length; global_idx ++) {
            if (positidxs[pos][global_idx] == local_idx)
                return global_idx;
        }
        return -1;
    }

    /**
     * Map a global branch point index to a local, position-specific index
     * @param pos position in alignment or POG
     * @param global_idx branch point index in phylogenetic tree
     * @return the index in the local, position-specific tree
     */
    private int global2local(int pos, int global_idx) {
        return positidxs[pos][global_idx];
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
            // pruneMe contains indices that SHOULD BE REMOVED
            int[] indices = phylo.getPrunedIndex(pruneMe);
            // save indices for quick re-retrieval later
            positidxs[position] = indices;
            // save tree for quick re-retrieval later
            positrees[position] = IdxTree.createPrunedTree(phylo, indices);
        } // else the index tree was already cached...
        return positrees[position];
    }

    /**
     * Retrieve all necessary values to instantiate the tree at a specified position. These values are based on:
     * Each leaf (branch point with no children) is assigned a state as input to inference.
     * After inference, each ancestor branch point is assigned a state GIVEN leaf states.
     *
     * @param position the position in the input alignment/or POG
     * @param mode inference mode, currently only GRASP.Inference.JOINT is supported
     * @return a TreeInstance with states from input data and from inference
     */
    public TreeInstance getTreeInstance(int position, GRASP.Inference mode) {
        if (mode == GRASP.Inference.JOINT) {
            if (states != null) { // states[bpidx][pos]
                IdxTree tree = getTree(position);
                Object[] assigned = new Object[tree.getSize()];
                for (int i = 0; i < assigned.length; i ++) {
                    int global_idx = local2global(position, i);
                    if (global_idx != -1) { // exists in "local", position-specific tree
                        POGraph pog = pogTree.getExtant(global_idx);
                        if (pog != null) {
                            SymNode node = (SymNode) pog.getNode(position);
                            if (node != null)
                                assigned[i] = node.get();
                        } else
                            assigned[i] = states[global_idx][position];
                    }
                }
                TreeInstance ti = new TreeInstance(getTree(position), assigned);
                return ti;
            }
            throw new ASRRuntimeException("Ancestor states have not been inferred");
        } else if (mode == GRASP.Inference.MARGINAL) {
            throw new ASRRuntimeException("Not implemented: requires all ancestors present in nominated position to have been inferred");
        }
        throw new ASRRuntimeException("Unknown inference mode: " + mode);
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
            if (pog != null) {
                if (mode == GRASP.Inference.JOINT) {
                    if (states != null) {
                        pog.decorateNodes(states[idx]);
                        ancestors.put(phylotree.getLabel(idx), pog);
                    }
                } else if (mode == GRASP.Inference.MARGINAL) {
                    if (distribs[idx] != null) {
                        pog.decorateNodes(distribs[idx]);
                        ancestors.put(phylotree.getLabel(idx), pog);
                    }
                }
            }
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
            // FIXME: create an index map for "inf" to enable generics <EnumDistrib>?
            TreeDecor[] inf = new TreeDecor[getPositions()];    // which is also how many inferences we will carry out
            for (int pos = 0; pos < getPositions(); pos ++) {                   // for each position...
                trees[pos] = getTree(pos);                                      //   this is the tree with indels imputed
                int ancidx = positidxs[pos][bpidx];                             //   index for sought ancestor in the position-specific tree
                if (ancidx >= 0)                                                //   which may not exist, i.e. part of an indel, but if it is real...
                    inf[pos] = new MaxLhoodMarginal(ancidx, trees[pos], MODEL);//     set-up the inference
            }

            distribs[bpidx] = new EnumDistrib[pogTree.getPositions()];
            treeinstances = new TreeInstance[pogTree.getPositions()];
            for (int pos = 0; pos < getPositions(); pos ++) {                   // for each position...
                int specidx = positidxs[pos][bpidx];                            //   index for sought ancestor in the position-specific tree
                if (specidx >= 0) {                                             //   which may not exist, i.e. part of an indel, but if it is real...
                    treeinstances[pos] = pogTree.getNodeInstance(pos, trees[pos], positidxs[pos]); //     get the instances at the leaves at that position, and...
                }
            }
            ThreadedDecorators threadpool = new ThreadedDecorators(inf, treeinstances, GRASP.NTHREADS);
            try {
                Map<Integer, TreeDecor> ret = threadpool.runBatch();
                for (int pos = 0; pos < getPositions(); pos ++) {                   // for each position...
                    int specidx = positidxs[pos][bpidx];                            //   index for sought ancestor in the position-specific tree
                    if (specidx >= 0) {                                             //   which may not exist, i.e. part of an indel, but if it is real...
                        distribs[bpidx][pos] = (EnumDistrib)ret.get(pos).getDecoration(specidx);     //     extract distribution of marginal prob
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
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
        // FIXME: create an index map for "inf" to enable generics <EnumDistrib>?
        TreeDecor[] inf = new TreeDecor[getPositions()];  // which is also how many inferences we will carry out
        for (int pos = 0; pos < inf.length; pos++) {                // for each position...
            trees[pos] = getTree(pos);                              //   this is the tree with indels imputed
            inf[pos] = new MaxLhoodJoint(trees[pos], MODEL);       //     set-up the inference
        }
        treeinstances = new TreeInstance[getPositions()];
        for (int pos = 0; pos < getPositions(); pos ++) {        // for each position...
            treeinstances[pos] = pogTree.getNodeInstance(pos, trees[pos], positidxs[pos]); // get the instances at the leaves at that position, and...
        }
        ThreadedDecorators threadpool = new ThreadedDecorators(inf, treeinstances, GRASP.NTHREADS);
        try {
            Map<Integer, TreeDecor> ret = threadpool.runBatch();for (int pos = 0; pos < getPositions(); pos ++) {       // for each position...
                for (int idx : getAncestorIndices()) {
                    int ancidx = positidxs[pos][idx];                   //   index for sought ancestor in the position-specific tree
                    if (ancidx >= 0)                                   //   which may not exist, i.e. part of an indel, but if it is real...
                        states[idx][pos] = ret.get(pos).getDecoration(ancidx);          //     extract state
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
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
     * @param ancID ancestor ID (as named internally in GRASP, e.g. "3" for "N3" which then links to a specific branch point index, say 5
     * @param mode inference mode, currently JOINT and MARGINAL are supported (see enum defined in GRASP class)
     * @param gappy whether gaps should be included (as they appear from the POG index)
     * @return the sequence
     */
    public EnumSeq getSequence(Object ancID, GRASP.Inference mode, boolean gappy) {
        int bpidx = getBranchpointIndex(ancID);                            // the index of the ancestor as it appears in the phylogenetic tree
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID: " + ancID);
        int[] idxs = getConsensus(bpidx);
        if (idxs == null)
            throw new ASRRuntimeException("Failed to find optimal path for ancestor ID: " + ancID);
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
        String name = ancID.toString().startsWith("N") ? ancID.toString() : "N" + ancID.toString();
        seq.setName(name); // FIXME: internal label is here re-named to have an "N" in-front: make naming strategy more principled?
        return seq;
    }

    /**
     * Retrieve the most-supported sequence for a given ancestor.
     * Currently, it is the most probable path as estimated by looking at all sequences descendant to the ancestor,
     * using a dynamic programming algorithm
     * @param ancID the identifier/label for the ancestor
     * @return the positions of the corresponding ancestor POG that make up the "consensus" path
     */
    public int[] getConsensus(Object ancID) {
        int bpidx = getBranchpointIndex(ancID);                            // the index of the ancestor as it appears in the phylogenetic tree
        if (bpidx == -1)
            throw new ASRRuntimeException("Invalid ancestor ID: " + ancID);
        return getConsensus(bpidx);
    }

    /**
     * Retrieve the most-supported sequence for a given ancestor.
     * Currently, it is the most probable path as estimated by looking at all sequences descendant to the ancestor,
     * using a dynamic programming algorithm
     * @param bpidx the branch point index of the ancestor, in the phylogenetic tree
     * @return the positions of the corresponding ancestor POG that make up the "consensus" path
     */
    public int[] getConsensus(int bpidx) {
        POGraph pog = this.ancarr[bpidx];
        if (pog == null)
            throw new ASRRuntimeException("Ancestor has not been inferred: index is " + bpidx);
        // collect info to make decisions...
        int[] leaves = phylotree.getLeaves(bpidx); // determine all branch points of leaves (i.e. extants) under this ancestor
        // go through POG nodes, to determine the transition "weights"
        PriorityQueue<Integer> queue = new PriorityQueue<>();
        Set<Integer> visiting = new HashSet<>();
        for (int idx : pog.getForward()) { // to start us off: add all indices emanating from start
            queue.add(idx);
            visiting.add(idx);
        }
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
                nexts[nnexts.length] = pog.maxsize(); // add the end terminus
            }
            double[] rates = pogTree.getEdgeRates(curr, nexts, leaves);
            // set the weights
            for (int i = 0; i < nexts.length; i ++) {
                POGraph.StatusEdge edge = pog.getEdge(curr, nexts[i]);
                if (edge != null) {
                    double w = -Math.log(rates[i]);
                    edge.setWeight(w > 10000 ? 10000 : w); // neg log of prob; so P=1 means zero weight, low P means high weight
                } else {
                    throw new ASRRuntimeException("Invalid POG with missing edge: \"" + pog.getName() + "\" edge=" + curr + "-" +nexts[i]);
                }
            }
            // add indices of all nodes that can be visited next
            for (int idx : nexts) {
                if (idx != pog.maxsize() && !visiting.contains(idx)) {
                    queue.add(idx);
                    visiting.add(idx);
                }
            }
        }
        // with weights set, we find the optimal path (pick one if several; else need to query individual edges and assemble)
        int[] consensus = pog.getMostSupported();
        return consensus;
    }

    /**
     * Save all position-specific trees, as instantiated by either the input data (extants) or inference
     * @param directory the name of the directory in which all individual files will be saved
     * @throws IOException
     * @throws ASRException
     */
    public void saveTreeInstances(String directory) throws IOException, ASRException {
        File file = new File(directory);
        StringBuilder sb = new StringBuilder();
        if (file.mkdirs()) { // true if the directory was created, false otherwise
        } else {
            System.err.println("Directory " + directory + " already exists");
            //throw new ASRException("Directory " + directory + " already exists");
        }
        for (int pos = 0; pos < getPositions(); pos ++) {
            String name = directory + "/T" + Integer.toString(pos + 1) + ".nwk";
            Newick.save(getTreeInstance(pos, GRASP.Inference.JOINT), name);
        }
    }


    // --------------------------------------------------------------------------------------------------------------- //
    // static "factory" methods for constructing Prediction instances
    // including indel inference by parsimony and maximum likelihood
    // --------------------------------------------------------------------------------------------------------------- //

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
                    if (pog.getNode(anchor) == null)
                        pog.addNode(anchor, new Node());
                for (int to : entry.getValue()) { // link anchored Node to each of the admissible nodes in the anchor set
                    if (to < nPos)
                        if (pog.getNode(to) == null)
                            pog.addNode(to, new Node());
                    pog.addEdge(anchor, to, new POGraph.StatusEdge(true));
                }
                for (int from : entry.getValue()) { // all possible pairs of admissible Nodes are linked
                    for (int to : entry.getValue()) {
                        if (from < to)
                            pog.addEdge(from, to, new POGraph.StatusEdge(true));
                    }
                }
            }
        }
        return new Prediction(pogTree, ancestors);
    }

    public static Prediction PredictbyMaxLhood(POGTree pogTree){

        int nPos = pogTree.getPositions(); // find the number of indices that the POGs (input and ancestors) can use
        IdxTree tree = pogTree.getTree();  // indexed tree (quick access to branch points, no editing)
        Map<Object, POGraph> ancestors = new HashMap<>(); // prepare where predictions will go
        TreeInstance[] ti = pogTree.getNodeInstances(true);   // extract gap/no-gap (boolean) leaf instantiation for every position

        MaxLhoodJoint[] ji = new MaxLhoodJoint[ti.length];
        for (int i = 0; i < ji.length; i++) {
            if (i == 1) {
            }
            Object[] possible = {true, false};
            SubstModel substmodel = new JC(1, possible); // need to know the alphabet...
            ji[i] = new MaxLhoodJoint(tree, substmodel);
        }
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) {
            ji[i].decorate(ti[i]);
        }

        for (int j = 0; j < tree.getSize(); j ++) {         // for every (indexed) branch point (IdxTree defaults to depth-first order)
            if (tree.isLeaf(j))             // if leaf, ignore and
                continue;

            Object ancID = tree.getBranchPoint(j).getID();  // unique ancestor ID from BranchPoint class (modifiable only before assembled into IdxTree when creating Tree)
        POGraph pog = new POGraph(nPos);                // blank POG
        ancestors.put(ancID, pog);                      // put blank POG in place, to be modified below
        // each anchor-set (below) is a range of indices for the POG; the key is a position that MUST be entered, values contain admissible positions that follow
        Map<Integer, Set<Integer>> anchorsets = new HashMap<>();
        int current_anchor = -1; // the first jump always from the start position
        Set<Integer> anchorset = new HashSet<>();

        for (int i = 0; i < nPos; i ++) {   // now traverse all positions, consider if GAP, not-GAP or admissible
             if (ji[i].getDecoration(j) == Boolean.FALSE) { // ALWAYS character (i.e. not-GAP), so required position
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
                if (pog.getNode(anchor) == null)
                    pog.addNode(anchor, new Node());
            for (int to : entry.getValue()) { // link anchored Node to each of the admissible nodes in the anchor set
                if (to < nPos)
                    if (pog.getNode(to) == null)
                        pog.addNode(to, new Node());
                pog.addEdge(anchor, to, new POGraph.StatusEdge(true));
            }
            for (int from : entry.getValue()) { // all possible pairs of admissible Nodes are linked
                for (int to : entry.getValue()) {
                    if (from < to)
                        pog.addEdge(from, to, new POGraph.StatusEdge(true));
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
    public static Prediction PredictByIndelParsimony(POGTree pogTree, Boolean forceLinear) {
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

            // Add in all of the linear neighbours that are partially ordered

            if (forceLinear) {
                for (Interval1D poval : pogTree.getPOVals()) {
                    optionally_true.add(poval);
                    int count_contains = 0;
                    for (Interval1D definitive : definitively_true)
                        count_contains += definitive.contains(poval) ? 1 : 0;
                    if (count_contains > 1)
                        discarded.add(poval);
                }
            }

            int[][] edges_optional = new int[optionally_true.size()][2];
            i = 0;
            for (Interval1D edge : optionally_true)
                edges_optional[i++] = new int[]{edge.min, edge.max};

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
     * Simple Indel Code based on Simmons and Ochoterena "Gaps as characters..." 2000
     * This version does incorporate multiple, optional INDELs
     * It DOES generate POGs, but they are single-path POGs.
     * @param pogTree
     * @return instance of Prediction
     */
    public static Prediction PredictBySICP(POGTree pogTree) {
        int nPos = pogTree.getPositions(); //
        Random rand = new Random(nPos); // random seed set here
        IdxTree tree = pogTree.getTree();
        Map<Object, POGraph> ancestors = new HashMap<>();
        // Retrieve an instance for each indel across the whole alignment (ordered by interval tree)
        // (state for each extant-indel: present, absent or permissible/neutral, as per Simmons and Ochoterena, 2000)
        // initially "permissible/neutral" is encoded as null (the variable is uninstantiated); see POGraph.getSimpleGapCode)
        TreeInstance[] ti = pogTree.getIndelInstances(); // a pogTree has a list of "indels"; here leaves are instantiated with applicable indels
        // inference will infer true, false, or accept that both true and false can be correct
        Parsimony[] pi = new Parsimony[ti.length];
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) { // for each "indel"
            pi[i] = new Parsimony(ti[i]);
            pi[i].forward();
            pi[i].backward();
        }
        // the code below
        // (1) regardless, if an ancestor or extant, we can pull out what the INDEL states are: absent (false), present (true) or permissible (true/false)
        // (2) if an ancestor, an ancestor POG is created, using the info from (1)
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branch point, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length == 0) // not an ancestor
                continue; // skip the code below, only predictions for ancestors are used to compose POGs
            // else: ancestor branch point
            // (1) Find ancestor STATE for each INDEL, and
            // (2) resolve what the POG looks like...

            // First, construct a list to include all unambiguously true indels, some of which are
            // rendered inapplicable (due to being precluded by others)
            List<Interval1D> unambiguous = new ArrayList<>();
            List<Interval1D> ambiguous = new ArrayList<>();
            int i = 0; // interval index; this order is decided above when indels are instantiated and inferred
            for (Interval1D ival : pogTree.getIntervalTree()) { // order specific to pogTree, and linked with ti and pi
                List<Boolean> calls = pi[i].getOptimal(j); // for ancestor index j
                if (calls.contains(Boolean.TRUE)) { // INDEL can be TRUE
                    if (calls.size() == 1) // the ONLY value is TRUE so DEFINITIVELY include
                        unambiguous.add(ival);
                    else
                        ambiguous.add(ival);
                }
                i ++;
            }
            // the order in which the intervals are considered is important: sorted by first start-index, within-which end-index
            Collections.sort(unambiguous);
            // Second, construct an interval tree definitive, with INDELs that are not contained within a TRUE INDEL
            IntervalST<Boolean> definitive = new IntervalST<>();    // to hold all unambiguously TRUE and not-precluded INDELs
            Set<Integer> valididx = new HashSet<>();                // the set of indices that are used to hold all unambiguous calls
            Interval1D precluder = null;                            // the interval that is the last to have been added, when considered "in order"
            for (int cnt = 0; cnt < unambiguous.size(); cnt ++) {
                Interval1D current = unambiguous.get(cnt);          // "current" interval under consideration...
                if (cnt < unambiguous.size() - 1) {                 // there is at least one more after this...
                    Interval1D next = unambiguous.get(cnt + 1);     // so look-ahead to the "next" interval
                    if (!next.contains(current)) {                  // next is not precluding current...
                        if (precluder != null) {                    // consider if the last-addition does
                            if (!precluder.contains(current)) {     // last-addition does NOT preclude the current one either, so...
                                definitive.put(current, true);// add current interval to interval tree, true indicates that it is unambiguous
                                valididx.add(current.min);
                                valididx.add(current.max);
                                precluder = current;                // update last-addition
                            }
                        } else {                                    // there isn't a "last-addition", so...
                            definitive.put(current, true);    // add current
                            valididx.add(current.min);
                            valididx.add(current.max);
                            precluder = current;                    // update last-addition to current
                        }
                    }                                               // else: next interval precludes current, so can ignore current
                } else { // none after so include...
                    if (precluder != null) {                        // consider if the last-addition does
                        if (!precluder.contains(current)) {         // last-addition does NOT preclude the current one either, so...
                            definitive.put(current, true);    // add current interval to interval tree, the cnt is not relevant at this stage
                            valididx.add(current.min);
                            valididx.add(current.max);
                        }
                    } else {                                        // there isn't a "last-addition", so...
                        definitive.put(current, true);
                        valididx.add(current.min);
                        valididx.add(current.max);
                    }
                }
            }
            // optionally, add ambiguous calls, i.e. indels that are optimally both true and false.
            // With SICP, an unambiguous call for an indel A precludes other calls, say B, if B is contained in A,
            // regardless of B being unambiguous or ambiguous.
            // However, an ambiguous call for an indel C does NOT preclude calls for other ambiguous calls.
            // To incorporate ambiguous calls, we thus (a) refrain from adding those which are contained by unambiguous
            // indels (which are not themselves contained), and (b) add ambiguous calls that have start and end points that
            // are supported by unambiguous calls.
            // the code below is an altered (mostly extended) version of the above strategy.
            // Sorting is important: ambiguous indels added to the definitive interval tree in-order,
            // will never contain those that follow.
            // Last statement is un-true!!!! This is what needs to be fixed!!!
            Collections.sort(ambiguous);
            EdgeMap emap = new EdgeMap();
            List<Interval1D> ambigedges = new ArrayList<>();
            for (int cnt = 0; cnt < ambiguous.size(); cnt ++) {
                Interval1D current = ambiguous.get(cnt);        // "current" interval under consideration...
                if (valididx.contains(current.min) && valididx.contains(current.max)) { // require both indices to be included from unambiguous calls
                    boolean not_contained = true;
                    for (Interval1D overlap : definitive.searchAll(current))
                        if (overlap.contains(current)) {
                            not_contained = false;
                            break;
                        }
                    if (not_contained)
                        ambigedges.add(current);
                        //emap.add(current.min, current.max);
                        //definitive.put(current, false);   // add current interval to interval tree, false means it is ambiguous
                }
            }
            for (Interval1D ival : ambigedges)
                definitive.put(ival, false);
            // finally, we are now in a position to create edges for a POG, including edges that are just linkers,
            // representing discontinuous sequence without decision
            int prev = -1;
            for (Interval1D edge : definitive) {
                if (edge.min > prev) { // linker required, so time to patch the way we anticipate sequence based implementations do...
                    for (int ptr = prev; ptr < edge.min; ptr++)
                        emap.add(ptr, ptr + 1);
                }
                emap.add(edge.min, edge.max);
                if (definitive.get(edge).contains(true))    // possibly test if it is unambiguous or ambiguous, before deciding to...
                    emap.add(edge.min, edge.max);           // label the edge as "reciprocated"
                prev = edge.max;
            }
            // finally put the info into a POG
            POGraph pog = POGraph.createFromEdgeMap(nPos, emap);
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
    public static Prediction PredictBySICML(POGTree pogTree) {
        Object[] possible = {true, false};
        SubstModel substmodel = new JC(1, possible);
        return Prediction.PredictBySICML(pogTree, substmodel);
    }

    /**
     * Simple Indel Code based on Simmons and Ochoterena "Gaps as characters..." (2000).
     * Inference with ML.
     * @param pogTree
     * @param gain_loss_model model for gain and loss events
     * @return instance of IndelPrediction
     */
    public static Prediction PredictBySICML(POGTree pogTree, SubstModel gain_loss_model) {
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
        MaxLhoodJoint[] ji = new MaxLhoodJoint[ti.length];
        for (int i = 0; i < ji.length; i++) { // for each "indel" we need to infer either gain or loss, so set-up inference
            ji[i] = new MaxLhoodJoint(tree, gain_loss_model);
        }
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) { // for each "indel" we need to infer either gain or loss, so set-up inference
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
        // (1) regardless, if an ancestor or extant, we can pull out what the INDEL states are: absent (false), present (true)
        // (2) if an ancestor, an ancestor POG is created, using the info from (1)
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branch point, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length == 0) { // not an ancestor, so we just print out debug info below before continuing with next
                if (DEBUG) {
                    POGraph pog = pogTree.getExtant(ancID);
                    System.out.print(ancID + "\t");
                    if (pog != null) {
                        int i = 0;
                        for (Interval1D ival : pogTree.getIntervalTree()) {
                            System.out.print(((Boolean)ji[i].getDecoration(j) ? "L" : "G") + "\t");
                            i++;
                        }
                    }
                    System.out.println();
                }
                continue; // skip the code below, only predictions for ancestors are used to compose POGs
            }
            // else: ancestor branch point
            // (1) Find ancestor STATE for each INDEL, and
            // (2) resolve what the POG looks like...
            // First, construct a list to include all unambiguously true indels, some of which are
            // rendered inapplicable (due to being precluded by others)
            List<Interval1D> unambiguous = new ArrayList<>();
            if (DEBUG) System.out.print(ancID + "\t");
            int i = 0; // interval index; this order is decided above when indels are instantiated and inferred
            for (Interval1D ival : pogTree.getIntervalTree()) { // order specific to pogTree, and linked with ti and pi
                Boolean call = (Boolean)ji[i].getDecoration(j);
                if (DEBUG)
                    System.out.print((call ? "L" : "G") + "\t");
                if (call)
                    unambiguous.add(ival);
                i ++;
            }
            if (DEBUG) System.out.println();
            // the order in which the intervals are considered is important: sorted by first start-index, within-which end-index
            Collections.sort(unambiguous);
            // Second, construct an interval tree definitive, with INDELs that are not contained within a TRUE INDEL
            IntervalST<Boolean> definitive = new IntervalST<>();    // to hold all unambiguously TRUE and not-precluded INDELs
            Set<Integer> valididx = new HashSet<>();                // the set of indices that are used to hold all unambiguous calls
            Interval1D precluder = null;                            // the interval that is the last to have been added, when considered "in order"
            for (int cnt = 0; cnt < unambiguous.size(); cnt ++) {
                Interval1D current = unambiguous.get(cnt);          // "current" interval under consideration...
                if (cnt < unambiguous.size() - 1) {                 // there is at least one more after this...
                    Interval1D next = unambiguous.get(cnt + 1);     // so look-ahead to the "next" interval
                    if (!next.contains(current)) {                  // next is not precluding current...
                        if (precluder != null) {                    // consider if the last-addition does
                            if (!precluder.contains(current)) {     // last-addition does NOT preclude the current one either, so...
                                definitive.put(current, true);// add current interval to interval tree, true indicates that it is unambiguous
                                valididx.add(current.min);
                                valididx.add(current.max);
                                precluder = current;                // update last-addition
                            }
                        } else {                                    // there isn't a "last-addition", so...
                            definitive.put(current, true);    // add current
                            valididx.add(current.min);
                            valididx.add(current.max);
                            precluder = current;                    // update last-addition to current
                        }
                    }                                               // else: next interval precludes current, so can ignore current
                } else { // none after so include...
                    if (precluder != null) {                        // consider if the last-addition does
                        if (!precluder.contains(current)) {         // last-addition does NOT preclude the current one either, so...
                            definitive.put(current, true);    // add current interval to interval tree, the cnt is not relevant at this stage
                            valididx.add(current.min);
                            valididx.add(current.max);
                        }
                    } else {                                        // there isn't a "last-addition", so...
                        definitive.put(current, true);
                        valididx.add(current.min);
                        valididx.add(current.max);
                    }
                }
            }
            // finally, we are now in a position to create edges for a POG, including edges that are just linkers,
            // representing discontinuous sequence without decision
            EdgeMap emap = new EdgeMap();
            int prev = -1;
            for (Interval1D edge : definitive) {
                if (edge.min > prev) { // linker required, so time to patch the way we anticipate sequence based implementations do...
                    for (int ptr = prev; ptr < edge.min; ptr++)
                        emap.add(ptr, ptr + 1);
                }
                emap.add(edge.min, edge.max);
                if (definitive.get(edge).contains(true))    // possibly test if it is unambiguous or ambiguous, before deciding to...
                    emap.add(edge.min, edge.max);           // label the edge as "reciprocated"
                prev = edge.max;
            }
            // finally put the info into a POG
            POGraph pog = POGraph.createFromEdgeMap(nPos, emap);
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
    public static Prediction PredictByIndelMaxLhood(POGTree pogTree, Boolean forceLinear) {
        Object[] possible = {true, false};
        SubstModel substmodel = new JC(1, possible);
        return Prediction.PredictByIndelMaxLhood(pogTree, forceLinear, substmodel);
    }

    /**
     * Simple Indel Code based on Simmons and Ochoterena "Gaps as characters..." (2000).
     * Inference with ML.
     * @param pogTree
     * @param gain_loss_model model for gain and loss events
     * @return instance of IndelPrediction
     */
    public static Prediction PredictByIndelMaxLhood(POGTree pogTree, Boolean forceLinear, SubstModel gain_loss_model) {
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
        MaxLhoodJoint[] ji = new MaxLhoodJoint[ti.length];
        for (int i = 0; i < ji.length; i++) { // for each "indel" we need to infer either gain or loss, so set-up inference
            ji[i] = new MaxLhoodJoint(tree, gain_loss_model);
        }
        // Below is where the main inference occurs
        // this stage should be multi-threaded... not so at the moment
        for (int i = 0; i < ti.length; i++) { // for each "indel" we need to infer either gain or loss, so set-up inference
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

            // Add in all of the linear neighbours that are partially ordered

            if (forceLinear) {

                for (Interval1D poval : pogTree.getPOVals()) {
                    optionally_true.add(poval);
                    int count_contains = 0;
                    for (Interval1D definitive : definitively_true)
                        count_contains += definitive.contains(poval) ? 1 : 0;
                    if (count_contains > 1)
                        discarded.add(poval);
                }
            }

            definitively_true.removeAll(discarded);

            int[][] edges_optional = new int[optionally_true.size()][2];
            i = 0;
            for (Interval1D edge : optionally_true)
                edges_optional[i++] = new int[]{edge.min, edge.max};

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
                jif[i] = new MaxLhoodJoint(tree, substmodel);
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
                jib[i] = new MaxLhoodJoint(tree, substmodel);
            }
        }
        ThreadedDecorators fpool = new ThreadedDecorators(jif, tif, (GRASP.NTHREADS + 1)/ 2);
        ThreadedDecorators bpool = new ThreadedDecorators(jib, tib, (GRASP.NTHREADS + 1)/ 2);
        // Below is where the main inference occurs
        try {
            Map<Integer, TreeDecor> fret = fpool.runBatch();
            Map<Integer, TreeDecor> bret = bpool.runBatch();
        } catch (Exception e) {
            e.printStackTrace();
        }
        // inference done, now assemble...
        for (int j = 0; j < tree.getSize(); j++) { // we look at each branchpoint, corresponding to either an extant or ancestor sequence
            Object ancID = tree.getBranchPoint(j).getID();
            if (tree.getChildren(j).length > 0) { // an ancestor
                EdgeMap emap = new EdgeMap();
                for (int i = -1; i <= nPos; i ++) {
                    if (i != nPos) {


                        if (DEBUG) {
                            System.out.println("i " + i);
                            System.out.println("j " + j);
                            System.out.println(jif[i + 1].getDecoration(j));
                        }
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
