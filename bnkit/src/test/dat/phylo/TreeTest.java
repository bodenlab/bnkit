package dat.phylo;

import dat.file.Newick;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Collection;

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
        Collection<Tree.BranchPoint> bps1 = t1.getRoot().getDescendants();
        Collection<Tree.BranchPoint> bps2 = t2.getRoot().getDescendants();
        assertEquals(bps1.size(), bps2.size());
        for (Tree.BranchPoint bp1 : bps1) {
            boolean found = false;
            for (Tree.BranchPoint bp2 : bps2) {
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
            Tree.BranchPoint root1 = t3l.getRoot();
            assertEquals(n, root1.getSubtree().size());
            t3l = Tree.load("bnkit/src/test/resources/saveme.anwk", "newick");
            assertNotNull(t3l);
            assertNull(t3l.find("_mysubtree is cool"));
            Tree.BranchPoint root2 = t3l.getRoot();
            assertEquals(root1.getSubtree().size(), root2.getSubtree().size());
            t3l = Tree.load("bnkit/src/test/resources/saveme.snwk", "newick");
            assertNotNull(t3l);
            assertNull(t3l.find("_mysubtree is cool"));
        } catch (IOException e) {
            assertTrue(false);
        }

        assertNotNull(t3);
    }
    @Test
    void find() {
    }

    static Tree tree = null;

    @BeforeAll
    static void setThingsUp() {
        try {
            tree = Tree.load("bnkit/src/test/resources/default.nwk", "newick");
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

}