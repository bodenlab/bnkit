package dat.pog;

import dat.EnumSeq;
import dat.Enumerable;
import dat.Interval1D;
import dat.IntervalST;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;

import java.io.IOException;
import java.util.*;

/**
 * A collection of POGs, linked up with a phylogenetic tree.
 * This class is intended to draw information from and integrate the POGs and the tree,
 * to present information necessary for analysis and inference.
 * The class TreeInstance is the vehicle for a unit of information to be transferred
 * to inference algorithms.
 */
public class POGTree {

    private IdxTree phylotree;              // the phylogenetic tree for the sequences
    // TODO: consider refactor extants to be in an array?
    private Map<String, POGraph> extants;   // POGs for extant sequences
    private IntervalST<String> ivals;       // Aggregation of all indels indicated by extant sequences
    private Enumerable domain;              // the alphabet
    private final int nNodes;

    /**
     * Construct a POGTree from an alignment of sequences;
     * each extant POG is named to match 1-1 against nodes in a phylogenetic tree.
     * We also collect all "indels" as intervals on the sequence indices.
     *
     * @param aln the muultiple sequence alignment
     * @param tree the phylogenetic tree
     */
    public POGTree(dat.EnumSeq.Alignment<Enumerable> aln, IdxTree tree) {
        this.phylotree = tree;
        this.extants = new HashMap<>();
        this.domain = aln.getDomain();
        this.nNodes = aln.getWidth();
        this.ivals = new IntervalST<>();
        for (int j = 0; j < aln.getHeight(); j++) {
            EnumSeq.Gappy<Enumerable> gseq = aln.getEnumSeq(j);
            POGraph pog = new POGraph(aln.getWidth());
            extants.put(gseq.getName(), pog);
            int start = 0, from = -1, to = -1;
            Object[] syms = gseq.get();
            for (int i = start; i < aln.getWidth(); i++) {
                if (syms[i] == null)
                    continue;
                else if (to == -1) {
                    to = i;
                    pog.addNode(to, new SymNode(syms[i]));
                    pog.addEdge(from, to);
                    ivals.put(new Interval1D(from, to), gseq.getName());
                } else {
                    from = to;
                    to = i;
                    pog.addNode(to, new SymNode(syms[i]));
                    pog.addEdge(from, to);
                    ivals.put(new Interval1D(from, to), gseq.getName());
                }
            }
            from = to;
            to = pog.maxsize();
            pog.addEdge(from, to);
            ivals.put(new Interval1D(from, to), gseq.getName());
        }
    }

    /**
     * Retrieve the data type (domain) of the characters on the alignment,
     * nodes in the POG and target states of the phylogenetic inference.
     * @return the domain
     */
    public Enumerable getDomain() {
        return domain;
    }

    /**
     * Retrieve the phylogenetic tree, as an index tree
     * @return the tree
     */
    public IdxTree getTree() {
        return phylotree;
    }

    /**
     * Get the data structure holding all index as collected from extant sequences.
     * @return the interval tree that contains all indels (as intervals)
     */
    public IntervalST<String> getIntervalTree() {
        return ivals;
    }

    /**
     * Get the number of positions that comprise the full collection of sequences, and indicate their homology
     * @return number of positions, indices that are used to retrieve content in extant and ancestral sequences
     */
    public int getPositions() {
        return nNodes;
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

        throw new RuntimeException("Not implemented");
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
     * @param index index of POG, or column of alignment
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
     * @param index       index of POG, or column of alignment
     * @param nullNotnull convert absence of assignment of leaf node to TRUE, presence of assignment to FALSE
     * @return a tree instance with states set according to extant sequences
     */
    public TreeInstance getNodeInstance(int index, boolean nullNotnull) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i = 0; i < phylotree.getSize(); i ++) {
            if (phylotree.getChildren(i).length == 0) { // leaf node, see if we can instantiate
                Object label = phylotree.getLabel(i);
                if (label != null) {
                    POGraph pog = extants.get(label);
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
            } // else parent... ignore
        }
        return new TreeInstance(phylotree, instarr);
    }

    /**
     * Determine a TreeInstance for a given node index in the POGs, for a specified tree, based on symbols at leaf nodes;
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     *
     * @param index       index of POG, or column of alignment
     * @param tree        tree that specifies the relationship between leaf nodes; note that the original tree in the POGTree is overridden
     * @return a tree instance with states set according to extant sequences' character states
     */
    public TreeInstance getNodeInstance(int index, IdxTree tree) {
        Object[] instarr = new Object[tree.getSize()];
        for (int i = 0; i < tree.getSize(); i ++) {
            if (tree.getChildren(i).length == 0) { // leaf node, see if we can instantiate
                Object label = tree.getLabel(i);
                if (label != null) {
                    POGraph pog = extants.get(label);
                    if (pog != null) {
                        SymNode node = (SymNode) pog.getNode(index);
                        if (node != null)
                            instarr[i] = node.get();
                    }
                }
            } // else parent... ignore
        }
        return new TreeInstance(tree, instarr);
    }

    public static boolean EDGE_FORWARD = true;
    public static boolean EDGE_BACKWARD = false;

    /**
     * Determine the TreeInstance for edges to or from a specified index in the POGs, based on edges extracted from input alignment;
     * the return value will represent the tree with nodes representing the edges that the extant takes, and
     * uninstantiated slots intended for the edges that the ancestors use.
     *
     * @param index index of POG, or column of alignment
     * @param EDGE_STATUS_FORWARD set to true, if seeking the set of forward looking edges, set to false if backward-looking
     * @return a tree instance with states set according to extant sequences' edges (as represented by their resp POG)
     */
    public TreeInstance getEdgeInstance(int index, boolean EDGE_STATUS_FORWARD) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i = 0; i < phylotree.getSize(); i++) {
            if (phylotree.getChildren(i).length == 0) { // leaf node, see if we can instantiate
                Object label = phylotree.getLabel(i);
                if (label != null) {
                    POGraph pog = extants.get(label);
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
            } // else parent... ignore
        }
        return new TreeInstance(phylotree, instarr);
    }

    /**
     * Determine TreeInstances for all edges in the POGs (collectively) based on Simple Gap Coding (Simmons and Ochoterena, 2000)
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * @return an array with TreeInstances indexed by edge in order of interval tree
     */
    public TreeInstance[] getIndelInstances() {
        TreeInstance[] insts = new TreeInstance[ivals.size()];
        int i = 0;
        for (Interval1D ival : ivals)
            insts[i ++] = getIndelInstance(ival);
        return insts;
    }

    /**
     * Determine TreeInstance for one edge across all the POGs based on Simple Gap Coding (Simmons and Ochoterena, 2000)
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * @return one TreeInstance
     */
    public TreeInstance getIndelInstance(Interval1D indel) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i = 0; i < phylotree.getSize(); i ++) {
            if (phylotree.isLeaf(i)) { // leaf node, see if we can instantiate
                Object label = phylotree.getLabel(i);
                if (label != null) {
                    POGraph pog = extants.get(label);
                    if (pog != null)
                        instarr[i] = pog.getSimpleGapCode(indel.min, indel.max);
                }
            } // else parent... ignore
        }
        return new TreeInstance(phylotree, instarr);
    }

    /**
     * Get the POG for a given extant sequence by its label
     * @param label sequence name
     * @return the POG
     */
    public POGraph getExtant(Object label) {
        return extants.get(label);
    }

    public static void main(String[] args) {
        try {
            EnumSeq.Alignment aln = new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/dp16_poag5.aln", Enumerable.aacid));
            Tree tree = Tree.load("/Users/mikael/simhome/ASR/dp16_poag5.nwk", "newick");
            POGTree pogt = new POGTree(aln, tree);
            for (Object name : pogt.extants.keySet()) {
                POGraph pog = (POGraph) pogt.extants.get(name);
                System.out.println(name + "\t" + pog);
                pog.saveToDOT("/Users/mikael/simhome/ASR/" + name + ".dot");
            }
        } catch (IOException e) {
            System.err.println(e);
        }
    }


}
