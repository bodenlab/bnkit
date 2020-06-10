package dat.pog;

import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.Tree;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;

class POGTreeTest {

    static EnumSeq.Alignment aln = null;
    static POGTree pogt1 = null;
    static Tree tree = null;

    @BeforeAll
    static void setPogt1() {
        try {
            aln = new EnumSeq.Alignment(EnumSeq.Gappy.loadClustal("src/test/resources/default.aln", Enumerable.aacid));
            tree = Tree.load("bnkit/src/test/resources/default.nwk", "newick");
            pogt1 = new POGTree(aln, tree);
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    @Test
    void getTree() {
        System.out.println(pogt1.getTree().getRoot());
    }


}