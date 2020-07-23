package dat.pog;

import dat.EnumSeq;
import dat.Enumerable;
import dat.Interval1D;
import dat.IntervalST;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.phylo.TreeInstance;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * A collection of POGs, linked up with a phylogenetic tree.
 * This class is intended to draw information from and integrate the POGs and the tree,
 * to present information necessary for analysis and inference.
 * The class TreeInstance is the vehicle for a unit of information to be transferred
 * to inference algorithms.
 */
public class POGTree {

    private IdxTree phylotree;              // the phylogenetic tree for the sequences
    private IdxTree[] positrees;            // position-specific tree, or an edited form of the original tree for the purpose of inferring content
    private Map<String, POGraph> extants;   // POGs for extant sequences
    private POGraph poag;                   // POG of alignment
    private IntervalST<String> ivals;       // Aggregation of all indels indicated by extant sequences
    private Enumerable domain;              // the alphabet
    private final int nNodes;

    /**
     * Construct a POGTree from an alignment of sequences, each named
     * to match 1-1 against nodes in a phylogenetic tree.
     *
     * @param aln
     * @param tree
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
                    pog.addNode(to, new SymNode(getDomain(), syms[i]));
                    pog.addEdge(from, to);
                    ivals.put(new Interval1D(from, to), gseq.getName());
                } else {
                    from = to;
                    to = i;
                    pog.addNode(to, new SymNode(getDomain(), syms[i]));
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

    public Enumerable getDomain() {
        return domain;
    }

    public IdxTree getTree() {
        return phylotree;
    }

    public IntervalST<String> getIntervalTree() {
        return ivals;
    }

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
    public POGTree(Map<String, POGraph> extants, Tree tree) {
        this.phylotree = tree;

        throw new RuntimeException("Not implemented");
    }

    /**
     * Determine TreeInstances for all node indices in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     * @return an array with TreeInstances indexed by position
     */
    public TreeInstance[] getNodeInstances() {
        int nPos = this.getPositions();
        TreeInstance[] ti = new TreeInstance[nPos];
        for (int i = 0; i < nPos; i ++)
            ti[i] = getNodeInstance(i, false);
        return ti;
    }

    /**
     * Determine TreeInstances for all node indices in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
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
     *
     * @param index index of POG, or column of alignment
     * @return
     */
    public TreeInstance getNodeInstance(int index) {
        return getNodeInstance(index, false);
    }

    /**
     * Determine the TreeInstance for a given node index in the POGs, based on symbols assigned to leaf nodes (class SymNode);
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     *
     * @param index       index of POG, or column of alignment
     * @param nullNotnull convert absence of assignment of leaf node to TRUE, presence of assignment to FALSE
     * @return
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

    public static int EDGE_TO_NODE = 1;
    public static int EDGE_FROM_NODE = 2;

    /**
     * Determine the TreeInstance for edges to or from a specified index in the POGs, based on edges extracted from input alignment;
     * the return value will represent the tree with nodes representing ancestors to be populated with content.
     *
     * @param index index of POG, or column of alignment
     * @return
     */
    public TreeInstance getEdgeInstance(int index, int EDGE_STATUS) {
        Object[] instarr = new Object[phylotree.getSize()];
        for (int i = 0; i < phylotree.getSize(); i++) {
            if (phylotree.getChildren(i).length == 0) { // leaf node, see if we can instantiate
                Object label = phylotree.getLabel(i);
                if (label != null) {
                    // ...

                }
            } // else parent... ignore
        }
        throw new RuntimeException("Not implemented");
//        return new TreeInstance(tree, instarr);
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
