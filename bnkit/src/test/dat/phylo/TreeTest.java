package dat.phylo;

import dat.file.Newick;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class TreeTest {

    String s1 = "((A:0.3,B:0.4):0.1,C:0.4):0.2,(D:0.5,E:0.1):0.2"; // first comma is that preceding D
    String s2 = "(E:0.1,D:0.5):0.2,(C:0.4,(B:0.4,A:0.3):0.1):0.2"; // first comma is that preceding C
    String ss1 = "(" + s1 + ");";
    String ss2 = "(" + s2 + ");";
    String ss3 = "((E:0.1,'D@ --':0.5)[&label=\"N_42\"]:0.2,(C:0.4,(B:0.4,A:0.3)'_mysubtree is cool':0.1):0.2)"; // try formatting, little NEXUS incorporated
    String ss4 = "((E:0.1,'D':0.5)[&label=\"N_42\",!hilight={38,2.203965965632127,#ffffcc}]:0.2,(C:0.4,(B:0.4,A:0.3)'mysubtree':0.1):0.2)"; // try formatting, more NEXUS

    @Test
    void getDescendants() {
        Tree t1 = Newick.parse(ss1);
        Tree t2 = Newick.parse(ss2);
        Collection<BranchPoint> bps1 = t1.getRoot().getDescendants();
        Collection<BranchPoint> bps2 = t2.getRoot().getDescendants();
        assertEquals(bps1.size(), bps2.size());
        for (BranchPoint bp1 : bps1) {
            boolean found = false;
            for (BranchPoint bp2 : bps2) {
                if (bp1.getLabel().equals(bp2.getLabel())) {
                    if (bp1.isLeaf())
                        assertEquals(bp1.getDistance(), bp2.getDistance());
                    found = true;
                    break;
                } else {
                }
            }
            assertEquals(true, found);
        }
    }

    @Test
    void newick1() {
        int pos = Newick.getComma(s1);
        assertEquals(29, pos);
        pos = Newick.getComma(s2);
        assertEquals(17, pos);
        int errors = 0;
        for (int i = 0; i < ss1.length(); i ++) {
            try {
                StringBuilder ns = new StringBuilder(ss1);
                ns.setCharAt(i, 'x');
                Tree t = Newick.parse(ns.toString());
            } catch (RuntimeException e) {
                ++errors;
//                System.err.println(errors + ": " + e.getMessage());
            }
        }
        assertEquals(36, errors);
        errors = 0;
        char[] testme = new char[] {' ', '\n', '\t'};
        for (char ch : testme) {
            for (int i = 0; i < ss2.length(); i++) {
                try {
                    StringBuilder ns = new StringBuilder(ss2);
                    ns.setCharAt(i, ch);
                    Tree t = Newick.parse(ns.toString());
//                    System.out.println("Fine: " + ns.toString() + "\t" + t.toString());
                } catch (RuntimeException e) {
                    ++errors;
//                    System.err.println(errors + ": " + e.getMessage());
                }
            }
        }
        assertEquals(19 * testme.length, errors);
        Tree t3 = null;
        try {
            t3 = Newick.parse(ss3);
        } catch (RuntimeException e) {
        }
        assertNotNull(t3);
//        System.out.println(t3);
    }

    @Test
    void load() {
        try {
            tree = Tree.load("bnkit/src/test/resources/small.nwk", "newick");
        } catch (IOException e) {
            assertTrue(false);
        }
        assertNotNull(tree);
    }

    @Test
    void save() {
        Tree t3 = null;
        Tree t3l = null;
        try {
            t3 = Newick.parse(ss3);
            int n = t3.getRoot().getDescendants().size() + 1;
            t3.save("bnkit/src/test/resources/saveme.nwk", "newick");
            t3.save("bnkit/src/test/resources/saveme.anwk", "ancestor");
            Newick.save(t3, "bnkit/src/test/resources/saveme.snwk", Newick.MODE_STRIPPED);
            t3l = Tree.load("bnkit/src/test/resources/saveme.nwk", "newick");
            assertNotNull(t3l);
            assertNotNull(t3l.find("_mysubtree is cool"));
            BranchPoint root1 = t3l.getRoot();
            assertEquals(n, root1.getSubtree().size());
            t3l = Tree.load("bnkit/src/test/resources/saveme.anwk", "newick");
            assertNotNull(t3l);
            assertNull(t3l.find("_mysubtree is cool"));
            BranchPoint root2 = t3l.getRoot();
            assertEquals(root1.getSubtree().size(), root2.getSubtree().size());
            t3l = Tree.load("bnkit/src/test/resources/saveme.snwk", "newick");
            assertNotNull(t3l);
            assertNull(t3l.find("_mysubtree is cool"));
            Tree dtree = Newick.load("bnkit/src/test/resources/default.nwk");
            dtree.save("bnkit/src/test/resources/default.anwk", "ancestor");
        } catch (IOException e) {
            assertTrue(false);
        }

        assertNotNull(t3);
    }
    @Test
    void find() {
    }

    @Test
    void iterator() {
        Tree t1 = Newick.parse(ss1);
        int cnt = 0;
        Set<Object> all = new HashSet<>();
        for (int i : t1) {
            BranchPoint bp = t1.getBranchPoint(i);
            all.add(bp.getID());
            cnt += 1;
        }
        assertEquals(t1.getSize(), cnt);
        assertEquals(t1.getSize(), all.size());
    }

    @Test
    void iterator2() {
        try {
            tree = Tree.load("bnkit/src/test/resources/small.nwk", "newick");
            Iterator<Integer> iter = tree.getBreadthFirstIterator();
            Set<Integer> all = new HashSet<>();
            int prevdepth = 0;
            while (iter.hasNext()) {
                Integer idx = iter.next();
                all.add(idx);
                int currdepth = tree.getDepth(idx);
                assertTrue(prevdepth <= currdepth);
                prevdepth = currdepth;
            }
            int cnt = 0;
            iter = tree.getDepthFirstIterator();
            while (iter.hasNext()) {
                Integer idx = iter.next();
                assertTrue(all.contains(idx));
                cnt += 1;
            }
            assertEquals(all.size(), cnt);
        } catch (IOException e) {
            assertTrue(false);
        }
    }


    static Tree tree = null;
    static Tree random_tree = null;

    @BeforeAll
    static void setThingsUp() {
        try {
            tree = Tree.load("bnkit/src/test/resources/default.nwk", "newick");
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    @Test
    void random() {
        for (int n = 6; n < 20; n ++) {
            String[] leaves = new String[n];
            for (int i = 0; i < n; i ++)
                leaves[i] = "A" + (i + 1);
            Tree random_tree1 = Tree.Random(leaves, n, 1, 2, 2, 2);
            Tree random_tree2 = Tree.Random(leaves, n + 1, 1, 2, 3, 3);
            for (String leaf : leaves) {
                BranchPoint bp1 = random_tree1.find(leaf);
                BranchPoint bp2 = random_tree2.find(leaf);
                assertTrue(bp1.getLabel().equals(bp2.getLabel()));
            }
            int sum1 = 0, sum2 = 0;
            for (int idx : random_tree1.getLeaves())
                sum1 += random_tree1.getDepth(idx);
            for (int idx : random_tree2.getLeaves())
                sum2 += random_tree2.getDepth(idx);
            if (!(sum1 > sum2)) {
                System.out.println(random_tree1.toString() + "\t" + sum1);
                System.out.println(random_tree2.toString() + "\t" + sum2);
            }
            assertTrue(sum1 > sum2);
        }
    }
}