package dat.phylo;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Collections;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class IdxTreeTest {

    @Test
    void createPrunedTree() {
        Set<Integer> s = IdxTree.getIndices(false, new Object[] {true, false, true, true, false, true});
        IdxTree pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 2, pruned.getSize());
        assertEquals(5, pruned.getRoots().length); // confirmed by manual inspection (MB)
        s = IdxTree.getIndices(false, new Object[] {false, true}); // remove root only, two subtrees should form
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 1, pruned.getSize());
        assertEquals(2, pruned.getRoots().length);
        s = IdxTree.getIndices(false, new Object[] {true, true, true, false}); // remove a leaf node only, one tree still
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 1, pruned.getSize());
        assertEquals(1, pruned.getRoots().length);
        s = IdxTree.getIndices(false, new Object[] {true, true, true, true, true, true, true, false, false, false}); // remove two leafs and their ancestor, one tree still
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 3, pruned.getSize());
        assertEquals(1, pruned.getRoots().length);
        s = IdxTree.getIndices(false, new Object[] {true, true, true, true, true, true, true, true, false, false}); // remove two leafs, one tree still
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 2, pruned.getSize());
        assertEquals(1, pruned.getRoots().length);
        s = IdxTree.getIndices(false, new Object[] {true, true, true, true, true, true, true, false, true, true}); // remove the ancestor of two leafs, three trees...
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 1, pruned.getSize());
        assertEquals(3, pruned.getRoots().length);
        IdxTree pruned_identical = IdxTree.createPrunedTree(defaultTree, defaultTree.getPrunedIndex(s, true));
        assertEquals(pruned.getSize(), pruned_identical.getSize());
        s = IdxTree.getIndices(false, new Object[] {true, true, false, true, true, false, true, false, false, true, true, true, true, false});
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 5, pruned.getSize());
        assertEquals(6, pruned.getRoots().length); // did not check manually!
        pruned_identical = IdxTree.createPrunedTree(defaultTree, defaultTree.getPrunedIndex(s, true));
        assertEquals(pruned.getSize(), pruned_identical.getSize());
        s = IdxTree.getIndices(false, new Object[] {false, true, true, true, false, true, true, false, true, false, false, true, true, true, true, false});
        pruned = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(defaultTree.getSize() - 6, pruned.getSize());
        assertEquals(7, pruned.getRoots().length); // did not check manually!
        pruned_identical = IdxTree.createPrunedTree(defaultTree, defaultTree.getPrunedIndex(s, true));
        assertEquals(pruned.getSize(), pruned_identical.getSize());
    }

    @Test
    void createSubtree() {
        IdxTree subtree = IdxTree.createSubtree(defaultTree, 0);
        assertEquals(subtree.getSize(), defaultTree.getSize());
        for (int bp = 0; bp < defaultTree.getSize(); bp ++) {
            IdxTree pruned = IdxTree.createPrunedTree(defaultTree, Collections.singleton(bp));
            System.out.println(pruned);
            int[] roots = pruned.getRoots();
            //System.out.print(roots.length + " roots: ");
            int total = 0;
            for (int rootidx : roots) {
                subtree = IdxTree.createSubtree(pruned, rootidx);
                total += subtree.getSize();
                //System.out.print(subtree.getSize() + " ");
            }
            //System.out.println("= " + pruned.getSize());
            assertEquals(pruned.getSize(), total);
        }
    }
    static IdxTree defaultTree = null;

    @BeforeAll
    static void setThingsUp() {
        Tree tree;
        try {
            tree = Tree.load("src/test/resources/default.nwk", "newick");
            defaultTree = (IdxTree)tree;
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    int nSEEDS = 10;
    int[] nextants      = new int[] {5, 25, 100};               // number of extants
    Tree[] trees = new Tree[nextants.length * nSEEDS];

    void setUp() {
        double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
        double GAMMA_SCALE = 0.2;
        for (int i = 0; i < nextants.length; i ++) {
            int N = nextants[i]; // number of leaves
            for (int SEED = 0; SEED < nSEEDS; SEED ++)
                trees[nSEEDS * i + SEED] = Tree.Random(N, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, 2, 2);
        }
        int z = 2;
    }

    @Test
    void toJSON() {
        setUp();
        assertTrue(defaultTree.toJSON().toString().equals(IdxTree.fromJSON(defaultTree.toJSON()).toJSON().toString()));
        assertEquals(defaultTree.hashCode(), IdxTree.fromJSON(defaultTree.toJSON()).hashCode());
        for (IdxTree t : trees)
            assertTrue(t.toJSON().toString().equals(IdxTree.fromJSON(t.toJSON()).toJSON().toString()));
        for (IdxTree t : trees)
            assertEquals(t.hashCode(), IdxTree.fromJSON(t.toJSON()).hashCode());
    }

    @Test
    void getIndicesOfOrphanedTrees1() {
        Set<Integer> s = IdxTree.getIndices(false, new Object[] {true, true /*N1*/, false /*N2*/, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false /*N9*/, true, true, false /*N11*/, true, true});
        for (int i : s) {
            System.out.println("Delete branch point: " + i  + "[" + defaultTree.getLabel(i) + "]");
        }
        IdxTree t = IdxTree.createPrunedTree(defaultTree, s);
        assertEquals(7, t.getRoots().length); // checked manually (MB)
        for (int i : t.getRoots()) {
            System.out.println("Root: " + i + "[" + t.getLabel(i) + "]");
            for (int j : t.getSubtreeIndices(i))
                System.out.println("\t" + j + "[" + t.getLabel(j) + "]");
        }
        Set<Integer> pruned = t.getIndicesOfOrphanedTrees();
        for (int i : pruned) {
            System.out.println("Remove: " + i + "[" + t.getLabel(i) + "]");
        }
        assertEquals(0, pruned.size()); // checked manually (MB)
    }

    @Test
    void getIndicesOfOrphanedTrees2() {
        Set<Integer> s = IdxTree.getIndices(false, new Object[] {true, true /*N1*/, false /*N2*/, true, true /*N3*/, true, true, true, true, true, true, true, true, true, true, false /*N8*/, true, true, true, false /*N10*/, true,
                true /*N11*/, true, true, true, true, true, true, true, true /*N15*/, true, true, true, true, true, true, false /*N18*/});
        for (int i : s) {
            System.out.println("Delete branch point: " + i  + "[" + defaultTree.getLabel(i) + "]");
        }
        IdxTree t = IdxTree.createPrunedTree(defaultTree, s);
        for (int i : t.getRoots()) {
            System.out.println("Root: " + i + "[" + t.getLabel(i) + "]");
            for (int j : t.getSubtreeIndices(i))
                System.out.println("\t" + j + "[" + t.getLabel(j) + "]");
        }
        Set<Integer> pruned = t.getIndicesOfOrphanedTrees();
        for (int i : pruned) {
            System.out.println("Remove: " + i + "[" + t.getLabel(i) + "]");
        }
        assertEquals(9, t.getRoots().length); // checked manually (MB)
        assertEquals(3, pruned.size()); // checked manually (MB)
        IdxTree clean = IdxTree.createPrunedTree(t, pruned);
        assertEquals(8, clean.getRoots().length); // one root fewer now
    }

    @Test
    void getIndicesOfOrphanedTrees3() {
        Set<Integer> s = IdxTree.getIndices(false, new Object[] {true, true /*N1*/, false /*N2*/, true, true /*N3*/, true, false /*N5*/, true, true, true, true, true, false /*N7*/, true, true, false /*N8*/, true, true, false /*N9*/, true /*N10*/, false /*s18*/,
                false /*N11*/, true, true, true, true, true, true, true, true /*N15*/, true, true, true, true, true, true, false /*N18*/});
        for (int i : s) {
            System.out.println("Delete branch point: " + i  + "[" + defaultTree.getLabel(i) + "]");
        }
        IdxTree t = IdxTree.createPrunedTree(defaultTree, s);
        for (int i : t.getRoots()) {
            System.out.println("Root: " + i + "[" + t.getLabel(i) + "]");
            for (int j : t.getSubtreeIndices(i))
                System.out.println("\t" + j + "[" + t.getLabel(j) + "]");
        }
        Set<Integer> pruned = t.getIndicesOfOrphanedTrees();
        for (int i : pruned) {
            System.out.println("Remove: " + i + "[" + t.getLabel(i) + "]");
        }
        assertEquals(14, t.getRoots().length); // checked manually (MB)
        assertEquals(3, pruned.size()); // checked manually (MB)
        IdxTree clean = IdxTree.createPrunedTree(t, pruned);
        assertEquals(12, clean.getRoots().length); // two roots fewer
    }

    @Test
    void getPrunedTree() {
        Set<Integer> s = IdxTree.getIndices(false, new Object[] {true, true /*N1*/, false /*N2*/, true, true /*N3*/, true, false /*N5*/, true, true, true, true, true, false /*N7*/, true, true, false /*N8*/, true, true, false /*N9*/, true /*N10*/, false /*s18*/,
                false /*N11*/, true, true, true, true, true, true, true, true /*N15*/, true, true, true, true, true, true, false /*N18*/});
        for (int i : s) {
            System.out.println("Delete branch point: " + i  + "[" + defaultTree.getLabel(i) + "]");
        }
        int[] indices = defaultTree.getPrunedIndex(s, true);
        IdxTree clean = IdxTree.createPrunedTree(defaultTree, indices);
        for (int i : clean.getRoots()) {
            System.out.println("Root: " + i + "[" + clean.getLabel(i) + "]");
            for (int j : clean.getSubtreeIndices(i))
                System.out.println("\t" + j + "[" + clean.getLabel(j) + "]");
        }
        assertEquals(12, clean.getRoots().length); // two roots fewer

    }
}