package dat.phylo;

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

    @Test
    void getDescendants() {
        Tree t1 = Tree.parseNewick(ss1);
        Tree t2 = Tree.parseNewick(ss2);
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
        int pos = Tree.getComma(s1);
        assertEquals(29, pos);
        pos = Tree.getComma(s2);
        assertEquals(17, pos);
        int errors = 0;
        for (int i = 0; i < ss1.length(); i ++) {
            try {
                StringBuilder ns = new StringBuilder(ss1);
                ns.setCharAt(i, 'x');
                Tree t = Tree.parseNewick(ns.toString());
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
                    Tree t = Tree.parseNewick(ns.toString());
                    System.out.println("Fine: " + ns.toString() + "\t" + t.toString());
                } catch (RuntimeException e) {
                    ++errors;
//                    System.err.println(errors + ": " + e.getMessage());
                }
            }
        }
        assertEquals(19 * testme.length, errors);
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