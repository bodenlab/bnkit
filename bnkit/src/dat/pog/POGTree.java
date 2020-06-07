package dat.pog;

import dat.EnumSeq;
import dat.Enumerable;
import dat.PhyloTree;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * A collection of POGs, linked up with a phylogenetic tree.
 */
public class POGTree {

    private PhyloTree tree;
    private Map<String, POGraph<Enumerable>> extants;
    private Enumerable domain;

    /**
     * Construct a POGTree from an alignment of sequences, each named
     * to match 1-1 against nodes in a phylogenetic tree.
     * @param aln
     * @param tree
     */
    public POGTree(dat.EnumSeq.Alignment<Enumerable> aln, dat.PhyloTree tree) {
        this.tree = tree;
        this.extants = new HashMap<>();
        this.domain = aln.getDomain();
        for (int j = 0; j < aln.getHeight(); j ++) {
            EnumSeq.Gappy<Enumerable> gseq = aln.getEnumSeq(j);
            POGraph<Enumerable> pog = new POGraph<>(aln.getWidth());
            extants.put(gseq.getName(), pog);
            int start = 0, from = -1, to = -1;
            for (int i = start; i < aln.getWidth(); i ++) {
                Object sym = gseq.get(i);
                if (sym == null)
                    continue;
                else if (to == -1) {
                    to = i;
                    pog.addNode(to, new SymNode(getDomain(), sym));
                    pog.addEdge(from, to);
                } else {
                    from = to;
                    to = i;
                    pog.addNode(to, new SymNode(getDomain(), sym));
                    pog.addEdge(from, to);
                }
            }
            from = to;
            to = pog.maxsize();
            pog.addEdge(from, to);
        }
    }

    public Enumerable getDomain() {
        return domain;
    }

    /**
     * Construct a POGTree from a collection of POGs, represented by a map keyed by sequence name,
     * to match 1-1 against nodes in a phylogenetic tree.
     * @param extants
     * @param tree
     */
    public POGTree(Map<String, POGraph<Enumerable>> extants, dat.PhyloTree tree) {
        this.tree = tree;

    }

    /**
     * Determine the tree for a given index of the sequences at leaf nodes;
     * the tree will contain nodes representing ancestors with content.
     * @param index
     * @return
     */
    PhyloTree getTree(int index) {
        throw new RuntimeException("Not implemented");
    }

    PhyloTree getTree() {
        return tree;
    }


    public static void main(String[] args) {
        try {
            EnumSeq.Alignment aln = new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/dp16_poag5.aln", Enumerable.aacid));
            PhyloTree tree = PhyloTree.load("/Users/mikael/simhome/ASR/dp16_poag5.nwk", "newick");
            POGTree pogt = new POGTree(aln, tree);
            for (Object name : pogt.extants.keySet()) {
                POGraph pog = (POGraph)pogt.extants.get(name);
                System.out.println(name + "\t" + pog);
                pog.saveToDOT("/Users/mikael/simhome/ASR/" + name + ".dot");
            }
        } catch (IOException e) {
            System.err.println(e);
        }
    }


}
