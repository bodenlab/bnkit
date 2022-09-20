package dat.pog;

import asr.ASRRuntimeException;
import dat.EnumSeq;
import dat.Enumerable;
import dat.Interval1D;
import dat.IntervalST;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;
import json.JSONArray;
import json.JSONObject;

import java.io.IOException;
import java.util.*;

/**
 * A collection of POGs, linked up with a phylogenetic tree.
 * This class is intended to draw information from and integrate the POGs and the tree,
 * to present information necessary for analysis and inference.
 * The class TreeInstance is the vehicle for a unit of information to be transferred
 * to inference algorithms.
 * @author mikael
 */
public class POGTree {

    private IdxTree phylotree;              // the phylogenetic tree for the sequences
    private int[] leafidxs = null;          // store the sub-set of indices that are used for leaf nodes (extants)
    // private Map<String, POGraph> extants;   // POGs for extant sequences
    // FIXED: refactoring extants to be in an array, see "extarr"
    private POGraph[] extarr;               // POGs for extant sequences, indexed by branchpoint index
    private Map<Object, Integer> id2bpidx;  // the indices of extant sequences
    private IntervalST<Integer> ivals;      // Aggregation of all indels indicated by extant sequences (refactored to map-to bpidx Dec 2020)
    private final int nNodes;

    /**
     * Construct a POGTree from an alignment of sequences;
     * each extant POG is named to match 1-1 against nodes in a phylogenetic tree.
     * We also collect all "indels" as intervals on the sequence indices.
     *
     * @param aln a multiple sequence alignment
     * @param tree a phylogenetic tree
     */
    public POGTree(dat.EnumSeq.Alignment<Enumerable> aln, IdxTree tree) {
        this.phylotree = tree;
        this.id2bpidx = new HashMap<>();
        this.extarr = new POGraph[tree.getSize()]; // store extant POGs in array, indexed by branchpoint index
//        this.domain = aln.getDomain();
        this.nNodes = aln.getWidth();
        this.ivals = new IntervalST<>();
        for (int j = 0; j < aln.getHeight(); j++) {
            EnumSeq.Gappy<Enumerable> gseq = aln.getEnumSeq(j);
            POGraph pog = new POGraph(aln.getWidth());
            pog.setName(gseq.getName());
            int bpidx = tree.getIndex(gseq.getName());
            extarr[bpidx] = pog;
            id2bpidx.put(gseq.getName(), bpidx);
            int start = 0, from = -1, to = -1;
            Object[] syms = gseq.get();
            for (int i = start; i < aln.getWidth(); i++) {
                if (syms[i] == null)
                    continue;
                else if (to == -1) {
                    to = i;
                    Node n = new SymNode(syms[i]);
                    n.setXLabel(i + 1);
                    pog.addNode(to, n);
                    pog.addEdge(from, to);
                    ivals.put(new Interval1D(from, to), bpidx);
                } else {
                    from = to;
                    to = i;
                    Node n = new SymNode(syms[i]);
                    n.setXLabel(i + 1);
                    pog.addNode(to, n);
                    pog.addEdge(from, to);
                    ivals.put(new Interval1D(from, to), bpidx);
                }
            }
            from = to;
            to = pog.maxsize();
            pog.addEdge(from, to);
            ivals.put(new Interval1D(from, to), bpidx);
        }
    }

    /**
     * Construct a POGTree from a collection of POGs, represented by a map keyed by sequence name,
     * to match 1-1 against nodes in a phylogenetic tree.
     *
     * @param extants
     * @param tree
     */
    public POGTree(Map<String, POGraph> extants, IdxTree tree) {
        this.phylotree = tree;
        this.id2bpidx = new HashMap<>();
        if (extants.size() != tree.getNLeaves())
            throw new ASRRuntimeException("Invalid combination of extant POGs and tree");
        this.extarr = new POGraph[tree.getSize()]; // store extant POGs in array, indexed by branchpoint index
        Enumerable mydomain = null;
        int mynNodes = -1;
        this.ivals = new IntervalST<>();
        int j = 0;
        for (Map.Entry<String, POGraph> entry : extants.entrySet()) {
            String name = entry.getKey();
            POGraph pog = entry.getValue();
            if (mynNodes == -1 || mynNodes == pog.nNodes)
                mynNodes = pog.nNodes;
            else
                throw new ASRRuntimeException("Invalid POG size: " + pog.nNodes + " should be " + mynNodes);
            int bpidx = tree.getIndex(name);
            extarr[bpidx] = pog;
            id2bpidx.put(name, bpidx);
            for (int i = 0; i < pog.nNodes; i ++) {
                if (pog.nodes[i] == null)
                    continue;
                else {
                    Node n = pog.getNode(i);
                    //new SymNode();
//                    n.setXLabel(i + 1);
//                    pog.addNode(to, n);
//                    pog.addEdge(from, to);
//                    ivals.put(new Interval1D(from, to), bpidx);
                }
            }
            j += 1;
        }
        this.nNodes = mynNodes;
    }


//    @Override
//    public boolean equals(Object o) {
//        if (this == o) return true;
//        if (!(o instanceof POGTree)) return false;
//        POGTree pogTree = (POGTree) o;
//        return nNodes == pogTree.nNodes && Objects.equals(phylotree, pogTree.phylotree) && Arrays.equals(leafidxs, pogTree.leafidxs) && Arrays.equals(extarr, pogTree.extarr) && Objects.equals(domain, pogTree.domain);
//    }

    @Override
    public int hashCode() {
        int result = Objects.hash(phylotree, nNodes);
        result = 31 * result + Arrays.hashCode(leafidxs);
        result = 31 * result + Arrays.hashCode(extarr);
        return result;
    }

    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        json.put("Hashcode", hashCode());
        json.put("Tree", phylotree.toJSON());
        List<POGraph> pogs = new ArrayList<>();
        for (int i : phylotree.getLeaves())
            pogs.add(extarr[i]);
        json.put("Extants", POGraph.toJSONArray(pogs));
        return json;
    }

    public static POGTree fromJSON(JSONObject json) {
        JSONArray jexts = json.optJSONArray("Extants");
        if (jexts == null)
            throw new ASRRuntimeException("Failed to load extants: no field \"Extants\" in JSON.");
        Map<String, POGraph> extmap = new HashMap<>();
        for (int i = 0; i < jexts.length(); i ++) {
            JSONObject obj = jexts.getJSONObject(i);
            POGraph pog = POGraph.fromJSON(obj);
            extmap.put(pog.getName(), pog);
        }
        JSONObject jtree = json.optJSONObject("Tree");
        if (jtree == null)
            throw new ASRRuntimeException("Failed to load tree: no field \"Extants\" in JSON.");
        IdxTree tree = IdxTree.fromJSON(jtree);
        return new POGTree(extmap, tree);
    }


    /**
     * Retrieve the data type (domain) of the characters on the alignment,
     * nodes in the POG and target states of the phylogenetic inference.
     * @return the domain
     */
//    public Enumerable getDomain() {
//        return domain;
//    }

    /**
     * Retrieve the phylogenetic tree, as an index tree
     * @return the tree
     */
    public IdxTree getTree() {
        return phylotree;
    }

    /**
     * Return the indices that are used for extants
     * @return indices as an array
     */
    private int[] getLeafIndices() {
        if (leafidxs == null)
            leafidxs = phylotree.getLeaves();
        return leafidxs;
    }

    /**
     * Get the name of the (extant) sequence that is positioned at the given branchpoint index.
     * Only use for extant sequences, and sparingly, as it is not efficient.
     * @param index branchpoint in the tree
     * @return the name
     */
    private Object getID(int index) {
        return getTree().getLabel(index);
    }

    /**
     * Get the POG for a given extant sequence by its label
     * @param id sequence name
     * @return the POG
     */
    public POGraph getExtant(Object id) {
        int bpidx = id2bpidx.get(id);
        return bpidx != -1 ? extarr[bpidx] : null;
    }

    /**
     * Get the POG for a given extant sequence by its branchpoint index
     * @param bpidx branchpoint index
     * @return the POG
     */
    public POGraph getExtant(int bpidx) {
        return extarr[bpidx];
    }

    /**
     * Get the index of the branchpoint that this extant sequence is placed
     * @param name name of extant sequence
     * @return the index in the tree, which is the leaf node with the named sequence
     */
    public int getIndex(String name) {
        Integer bpidx = id2bpidx.get(name);
        return bpidx != null ? bpidx.intValue() : -1;
    }

    /**
     * Get the data structure holding all index as collected from extant sequences.
     * @return the interval tree that contains all indels (as intervals)
     */
    public IntervalST<Integer> getIntervalTree() {
        return ivals;
    }

// MB removed: see above
//    /**
//     * Get the data structure holding all the linear neighbours that are partially ordered
//     * @return the interval tree that contains all partially ordered neighbours (as intervals)
//     */
//    public IntervalST<Integer> getPOVals() {
//        return povals;
//    }

    /**
     * Get the number of positions that comprise the full collection of sequences, and indicate their homology
     * @return number of positions, indices that are used to retrieve content in extant and ancestral sequences
     */
    public int getPositions() {
        return nNodes;
    }

    /**
     * Determine a conditional probability for each nominated target transition from a given source node in a POG,
     * based on a nominated set of extant sequences.
     * @param from source node index in POG
     * @param targets target node indices in POG
     * @param extants the sequences on which the calculation is based, indexed by branch point
     * @return
     */
    public double[] getEdgeRates(int from, int[] targets, int[] extants) {
        double[] rates = new double[targets.length];
        Arrays.fill(rates, 1); // pseudo count
        double denom = targets.length; // consider pseudo count
        for (int subidx : extants) {
            if (phylotree.isLeaf(subidx)) {
                int len = extarr[subidx].size();
                denom += len; // normalisation factor
                for (int i = 0; i <targets.length; i ++) {
                    int to = targets[i];
                    if (extarr[subidx].isEdge(from, to)) {
                        rates[i] += len; // count each extant in proportion to its length
                    }
                }
            }
        }
        // Above: count each sequence making a given jump, based on its length: longer sequences have greater weight
        // Below: normalise the count, based on the overall sequence lengths
        for (int i = 0; i <targets.length; i ++)
            rates[i] /= denom; // normalise counts
        return rates;
    }

    /**
     * Retrieve the indices of all extants that have a given edge
     * @param from  source node for edge
     * @param to target node for edge
     * @return indices of extants that match
     */
    public int[] getExtantsWithEdge(int from, int to) {
        int[] subtree = phylotree.getLeaves();
        Set<Integer> matched = new HashSet<>();
        for (int subidx : subtree) {
            if (extarr[subidx].isEdge(from, to))
                matched.add(subidx);
        }
        int[] arr = new int[matched.size()];
        int i = 0;
        for (int y : matched) arr[i ++] = y;
        return arr;
    }

    /**
     * Retrieve the indices of all extants under the given branch point index that have a given edge
     * @param from  source node for edge
     * @param to target node for edge
     * @param bpidx the branch point of the subtree under which extants are checked
     * @return indices of extants that match
     */
    public int[] getExtantsWithEdge(int from, int to, int bpidx) {
        Set<Integer> subtree = phylotree.getSubtreeIndices(bpidx);
        Set<Integer> matched = new HashSet<>();
        for (int subidx : subtree) {
            if (phylotree.isLeaf(subidx)) {
                if (extarr[subidx].isEdge(from, to))
                    matched.add(subidx);
            }
        }
        int[] arr = new int[matched.size()];
        int i = 0;
        for (int y : matched) arr[i ++] = y;
        return arr;
    }

    /**
     * Determine TreeInstances for all node indices in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * Importantly, this is by reference to the phylogenetic tree provided at construction, NOT one which has been processed
     * for indels, say.
     * @return an array with TreeInstances indexed by position
     */
    public TreeInstance[] getNodeInstances() {
        return getNodeInstances(false);
    }

    /**
     * Determine TreeInstances for all node indices in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * Importantly, this is by reference to the phylogenetic tree provided at construction, NOT one which has been processed
     * for indels, say.
     * @param nullNotnull convert absence of assignment of leaf node to TRUE, presence of assignment to FALSE
     * @return an array with TreeInstances indexed by position
     */
    public TreeInstance[] getNodeInstances(boolean nullNotnull) {
        int nPos = this.getPositions();
        TreeInstance[] ti = new TreeInstance[nPos];
        for (int i = 0; i < nPos; i ++)
            ti[i] = getNodeInstance(i, nullNotnull);
        return ti;
    }

    /**
     * Determine the TreeInstance for a given node index in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * Importantly, this is by reference to the phylogenetic tree provided at construction, NOT one which has been processed
     * for indels, say.
     *
     * @param index index of column of alignment, or position in POG
     * @return a tree instance with states set according to extant sequences' character content
     */
    public TreeInstance getNodeInstance(int index) {
        return getNodeInstance(index, false);
    }

    /**
     * Determine the TreeInstance for a given node index in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * Importantly, this is by reference to the phylogenetic tree provided at construction, NOT one which has been processed
     * for indels, say.
     *
     * @param index index of column of alignment, or position in POG
     * @param nullNotnull convert absence of assignment of leaf node to TRUE, presence of assignment to FALSE
     * @return a tree instance with states set according to extant sequences
     */
    public TreeInstance getNodeInstance(int index, boolean nullNotnull) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i : getLeafIndices()) {
            POGraph pog = extarr[i];
            if (pog != null) {
                SymNode node = (SymNode) pog.getNode(index);
                if (node != null && !nullNotnull)
                    instarr[i] = node.get();
                else if (node == null && nullNotnull)
                    instarr[i] = Boolean.TRUE;  // gap
                else if (node != null)
                    instarr[i] = Boolean.FALSE; // not gap
            }
        }
        return new TreeInstance(phylotree, instarr);
    }

    /**
     * Determine a TreeInstance for a given node index in the POGs, for a specified tree, based on symbols at leaf nodes;
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     *
     * @param index index of column of alignment, which is the same as the position in the POG
     * @param tree tree that specifies the relationship between leaf nodes; note that it the tree instance uses
     *             the index of the provided (position-specific) tree, NOT the phylogenetic tree of the POGTree
     * @param idxmap  map global branch point index to position-specific branch point index; this is used to map
     *                character states of the extant sequences (which use the phylogenetic tree of the POGTree)
     *                to the position-specific tree
     * @return a tree instance with states set according to extant sequences' character states
     */
    public TreeInstance getNodeInstance(int index, IdxTree tree, int[] idxmap) {
        Object[] instarr = new Object[tree.getSize()];
        int[] leaves = getLeafIndices();
        for (int i : leaves) {
            POGraph pog = extarr[i];
            if (pog != null) {
                try {
                    SymNode node = (SymNode) pog.getNode(index);
                    if (node != null)
                        // FIXME: below is a hack, ensuring that a string is converted into a character (one of the Enumerable values)
                        // FIXME continued: it could be a string if it has been loaded from JSON...
                        instarr[idxmap[i]] = node.get().toString().charAt(0);
                } catch (ClassCastException e) {
                }
            }
        }
        return new TreeInstance(tree, instarr);
    }

    public static boolean EDGE_FORWARD = true;
    public static boolean EDGE_BACKWARD = false;

    /**
     * Determine the TreeInstance for edges to or from a specified index in the POGs, based on edges extracted from input alignment;
     * the return value will represent the tree with nodes representing the edges that each extant takes, and
     * un-instantiated slots intended for the edges that the ancestors use.
     *
     * @param index index of column of alignment, or position in POG
     * @param EDGE_STATUS_FORWARD set to true, if seeking the set of forward looking edges, set to false if backward-looking
     * @return a tree instance with states set according to extant sequences' edges (as represented by their resp POG)
     */
    public TreeInstance getEdgeInstance(int index, boolean EDGE_STATUS_FORWARD) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i : getLeafIndices()) {
            POGraph pog = extarr[i];
            if (pog != null) {
                if (index == -1 || index == pog.maxsize()) {
                    if (index == -1 && EDGE_STATUS_FORWARD) {
                        int[] edges = pog.getForward(); // get start node
                        if (edges != null)
                            if (edges.length > 0)
                                instarr[i] = edges[0];
                    } else if (index == pog.maxsize() && !EDGE_STATUS_FORWARD) {
                        int[] edges = pog.getBackward(); // get terminal node
                        if (edges != null)
                            if (edges.length > 0)
                                instarr[i] = edges[0];
                    }
                } else {
                    if (pog.isStartNode(index) && !EDGE_STATUS_FORWARD) { // looking backward from the 1st node in the POG
                        instarr[i] = -1;
                    } else if (pog.isEndNode(index) && EDGE_STATUS_FORWARD) { // looking forward from last node in the POG
                        instarr[i] = pog.maxsize();
                    } else if (pog.isNode(index)) {
                        int[] edges = (EDGE_STATUS_FORWARD == EDGE_FORWARD ? pog.getForward(index) : pog.getBackward(index));
                        if (edges != null)
                            if (edges.length > 0)
                                instarr[i] = edges[0];
                    }
                }
            }
        }
        return new TreeInstance(phylotree, instarr);
    }

    /**
     * Determine TreeInstances for all edges in the POGs (collectively) based on Simple Gap Coding (Simmons and Ochoterena, 2000)
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * @return an array with TreeInstances indexed by edge in order of interval tree
     */
    public TreeInstance[] getIndelInstances() {
        List<TreeInstance> instances = new ArrayList<>();
        for (Interval1D ival : ivals) {
            if (ival.getWidth() > 1 || ival.min == -1 || ival.max == this.nNodes) // exclude non-gaps
                instances.add(getIndelInstance(ival));
        }
        TreeInstance[] insts = new TreeInstance[instances.size()];
        instances.toArray(insts);
        return insts;
    }

    /**
     * Determine TreeInstance for one edge across all the POGs based on Simple Gap Coding (Simmons and Ochoterena, 2000)
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * @return one TreeInstance
     */
    public TreeInstance getIndelInstance(Interval1D indel) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i : getLeafIndices()) {
            POGraph pog = extarr[i];
            if (pog != null)
                instarr[i] = pog.getSimpleGapCode(indel.min, indel.max);
        }
        return new TreeInstance(phylotree, instarr);
    }

    /**
     * Retrieve the set of identifiers for all extants
     * @return set of identifiers
     * @see POGTree#getExtant(Object)
     *
     */
    public Set<Object> getExtantIDs() {
        return id2bpidx.keySet();
    }

    public static void main(String[] args) {
        try {
            EnumSeq.Alignment aln = new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/dp16_poag5.aln", Enumerable.aacid));
            Tree tree = Tree.load("/Users/mikael/simhome/ASR/dp16_poag5.nwk", "newick");
            POGTree pogt = new POGTree(aln, tree);
            for (Object name : pogt.getExtantIDs()) {
                POGraph pog = (POGraph) pogt.getExtant(name);
                System.out.println(name + "\t" + pog);
                //pog.saveToDOT("/Users/mikael/simhome/ASR/" + name + ".dot");
            }
        } catch (IOException e) {
            System.err.println(e);
        }
    }


}
