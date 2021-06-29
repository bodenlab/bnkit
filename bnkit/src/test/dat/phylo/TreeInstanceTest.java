package dat.phylo;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

class TreeInstanceTest {

    IdxTree defaultTree;
    TreeInstance ti;
    Map<Integer, Object> assign = new HashMap<>();

    @BeforeEach
    void setUp() {
        Tree tree;
        try {
            tree = Tree.load("bnkit/src/test/resources/default.nwk", "newick");
            defaultTree = (IdxTree)tree;
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
        assign.put(defaultTree.getIndex("s22"), 11);
        assign.put(defaultTree.getIndex("s20"), 11);
        assign.put(defaultTree.getIndex("s19"), 22);
        assign.put(defaultTree.getIndex("s1"), 11);
        assign.put(defaultTree.getIndex("s2"), 33);
        assign.put(defaultTree.getIndex("s8"), 33);
        assign.put(defaultTree.getIndex("s7"), 33);
        ti = new TreeInstance(defaultTree, assign);
    }

    @Test
    @DisplayName("encode(true)")
    void encode1() {
        int cnt_notnull = 0;
        for (int i : defaultTree) {
            Object val = assign.get(i);
            if (val != null) {
                assertEquals(val, ti.getInstance(i));
                cnt_notnull += 1;
            }
        }
        System.out.println("TreeInstance with " + defaultTree.getSize() + " nodes, " + cnt_notnull + " leaves assigned non-NULL");
        Object[] key= ti.encode();
        int unique = new HashSet(assign.values()).size();
        int cnt_check = 0;
        for (int i : defaultTree) {
            Object val = assign.get(i);
            if (val != null) {
                System.out.println(i + " used to be " + val + "; now " + ti.getInstance(i));
                assertNotEquals(val, ti.getInstance(i));
                cnt_check += 1;
            } else if (ti.getInstance(i) != null) { // newly assigned node
                System.out.println(i + " used to be NULL; now " + ti.getInstance(i));
                assertEquals(true, defaultTree.isLeaf(i));
                assertEquals(unique, ti.getInstance(i));
            }
        }
        assertEquals(cnt_notnull, cnt_check);
    }
    @Test
    @DisplayName("encode(false)")
    void encode2() {
        int cnt_notnull = 0;
        for (int i : defaultTree) {
            Object val = assign.get(i);
            if (val != null) {
                assertEquals(val, ti.getInstance(i));
                cnt_notnull += 1;
            }
        }
        System.out.println("TreeInstance with " + defaultTree.getSize() + " nodes, " + cnt_notnull + " leaves assigned non-NULL");
        Object[] key= ti.encode(false);
        int unique = new HashSet(assign.values()).size();
        int cnt_check = 0;
        for (int i : defaultTree) {
            Object val = assign.get(i);
            if (val != null) {
                System.out.println(i + " used to be " + val + "; now " + ti.getInstance(i));
                assertNotEquals(val, ti.getInstance(i));
                cnt_check += 1;
            } else if (ti.getInstance(i) != null) { // newly assigned node
                System.out.println(i + " used to be NULL; now " + ti.getInstance(i));
                assertEquals(true, defaultTree.isLeaf(i));
                assertEquals(unique, ti.getInstance(i));
            }
        }
        assertEquals(cnt_notnull, cnt_check);
    }

    @Test
    void decode1() {
        Object[] key= ti.encode(true);
        int cnt_null = 0;
        for (int i : defaultTree) {
            if (defaultTree.isLeaf(i) && ti.getInstance(i) == null)
                cnt_null += 1;
        }
        System.out.println("TreeInstance with " + defaultTree.getSize() + " nodes, " + defaultTree.getNLeaves() + " leaves, of which " + cnt_null + " are assigned NULL");
        ti.decode(key);
        int cnt_check = 0;
        for (int i : defaultTree) {
            Object val = assign.get(i);
            if (val != null) {
                System.out.println(i + " used to be " + val + "; now back to " + ti.getInstance(i));
                assertEquals(val, ti.getInstance(i));
            }
            if (defaultTree.isLeaf(i) && ti.getInstance(i) == null)
                cnt_check += 1;
        }
        assertNotEquals(cnt_null, cnt_check);
    }
    @Test
    void decode2() {
        Object[] key= ti.encode(false);
        int cnt_null = 0;
        for (int i : defaultTree) {
            if (defaultTree.isLeaf(i) && ti.getInstance(i) == null)
                cnt_null += 1;
        }
        System.out.println("TreeInstance with " + defaultTree.getSize() + " nodes, " + defaultTree.getNLeaves() + " leaves, of which " + cnt_null + " are assigned NULL");
        ti.decode(key);
        int cnt_check = 0;
        for (int i : defaultTree) {
            Object val = assign.get(i);
            if (val != null) {
                System.out.println(i + " used to be " + val + "; now back to " + ti.getInstance(i));
                assertEquals(val, ti.getInstance(i));
            }
            if (defaultTree.isLeaf(i) && ti.getInstance(i) == null)
                cnt_check += 1;
        }
        assertEquals(cnt_null, cnt_check);
    }
}