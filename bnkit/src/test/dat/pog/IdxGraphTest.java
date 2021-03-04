package dat.pog;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;

import static org.junit.jupiter.api.Assertions.*;

class IdxGraphTest {

    IdxEdgeGraph.DefaultGraph g1 = null;
    IdxEdgeGraph.DefaultGraph g2 = null;
    IdxEdgeGraph.DefaultGraph g3 = null;

    @Test
    void addNode() {
        for (int i = -5; i < g1.maxsize() + 5; i ++) {
            assertFalse(g1.isNode(i));
            try {
                g1.addNode(i);
            } catch (RuntimeException ignored) { }
            assertEquals(g1.isIndex(i), g1.isNode(i));
        }
    }

    @Test
    void removeNode() {
        for (int i = 0; i < g1.maxsize() - 5; i ++)
            g1.addNode(i);
        for (int i = g1.maxsize() - 1; i >= 0; i --) {
            assertEquals(i < g1.maxsize() - 5, g1.isNode(i));
            try {
                g1.removeNode(i);
            } catch (RuntimeException e) { }
            assertFalse(g1.isNode(i));
        }
    }

    private class MyNode extends Node {
        final String label;
        public MyNode(String label) {
            this.label = label;
        }
        public String toString() { return label; }
    }

    void createGraphs() {
        for (int i = 0; i < g1.maxsize(); i++) {
            g1.addNode(i);
            g2.addNode(i);
            g3.addNode(i, new MyNode("Node[" + i + "]"));
        }
        for (int i = 0; i < g1.maxsize() - 1; i += 2)
            g3.addEdge(-1, i);
        for (int from = 0; from < g1.maxsize(); from += 2) {
            for (int to = 1; to < g1.maxsize(); to += 2) {
                g1.addEdge(from, to); // un-directed
                g2.addEdge(from, to); // directed
                g3.addEdge(from, to);
            }
        }
        for (int i = 1; i < g1.maxsize(); i += 3)
            g3.addEdge(g1.maxsize() - i, g3.maxsize());
    }

    @Test
    void addEdge() {
        createGraphs();
        for (int from = 0; from < g1.maxsize(); from += 1) {
            for (int to = 0; to < g1.maxsize(); to += 1) {
                if (to % 2 == 1 && from % 2 == 0) {
                    assertTrue(g1.isEdge(from, to));
                    assertTrue(g1.isEdge(to, from));
                    assertFalse(g2.isEdge(to, from));
                    assertTrue(g2.isEdge(from, to));
                } else {
                    assertEquals((from % 2 == 1 && to % 2 == 0), g2.isEdge(to, from));
                    assertFalse(g2.isEdge(from, to));
                }
            }
        }
    }

    void crippleGraphs() {
        int count = 0;
        for (int from = 0; from < g1.maxsize(); from += 1) {
            for (int to = 0; to < g1.maxsize(); to += 1) {
                count += 1;
                if (count % 3 == 0) {
                    try {
                        g1.removeEdge(from, to);
                    } catch (RuntimeException e) {
                    }
                    try {
                        g2.removeEdge(from, to);
                    } catch (RuntimeException e) {
                    }
                }
            }
        }
    }

    int removeEdges(int every, IdxGraph g) {
        int cnt = 0, rcnt = 0;
        for (int from = 0; from < g.maxsize(); from ++) {
            if (g.isNode(from)) {
                for (int to = 0; to < g.maxsize(); to++) {
                    if (g.isNode(to)) {
                        if (g.isEdge(from, to)) {
                            cnt += 1;
                            if (cnt % every == 0) {
                                g.removeEdge(from, to);
                                rcnt += 1;
                            }
                        }
                    }
                }
            }
        }
        return rcnt;
    }

    @Test
    void getEdgeCount() {
        createGraphs();
        int expect_g1_g2 = (g1.maxsize() / 2) * (g1.maxsize() / 2);
        int expect_g3 = g1.maxsize() / 2 + g1.maxsize() / 3 + expect_g1_g2;
        assertEquals(expect_g1_g2, g1.getEdgeCount());
        assertEquals(expect_g1_g2, g2.getEdgeCount());
        assertEquals( expect_g3 , g3.getEdgeCount());
        int rcnt = removeEdges(5, g1);
        assertEquals(expect_g1_g2 - rcnt, g1.getEdgeCount());
        rcnt = removeEdges(6, g2);
        assertEquals(expect_g1_g2 - rcnt, g2.getEdgeCount());
        rcnt = removeEdges(7, g3);
        assertEquals( expect_g3 - rcnt, g3.getEdgeCount());
    }

    @Test
    void removeEdge() {
        createGraphs();
        crippleGraphs();
        int count = 0;
        for (int from = 0; from < g1.maxsize(); from += 1) {
            for (int to = 0; to < g1.maxsize(); to += 1) {
                count += 1;
                if (count % 3 == 0) {
                    assertEquals(false, g1.isEdge(to, from));
                    assertEquals(false, g2.isEdge(from, to));
                } else if (to % 2 == 1 && from % 2 == 0) {
                    assertEquals(true, g1.isEdge(from, to));
                    assertEquals(true, g1.isEdge(to, from));
                    assertEquals(false, g2.isEdge(to, from));
                    assertEquals(true, g2.isEdge(from, to));
                } else {
                    assertEquals((from % 2 == 1 && to % 2 == 0), g2.isEdge(to, from));
                    assertEquals(false, g2.isEdge(from, to));
                }
            }
        }
    }

    @Test
    void getNodeIndices() {
        createGraphs();
        int count = 0;
        HashMap<Integer, int[]> m1 = new HashMap<>();
        HashMap<Integer, int[]> m2 = new HashMap<>();
        for (int from = 0; from < g1.maxsize(); from += 1) {
            count += 1;
            int[] idx1 = g1.getNodeIndices(from);
            m1.put(count, idx1);
            int[] idx2 = g2.getNodeIndices(from);
            m2.put(count, idx2);
            if (count % 3 == 0) {
                if (idx1.length > 0)
                    g1.removeEdge(from, idx1[0]);
                if (idx2.length > 0)
                    g2.removeEdge(from, idx2[0]);
            }
        }
        count = 0;
        for (int from = 0; from < g1.maxsize(); from += 1) {
            count += 1;
            int[] idx1 = g1.getNodeIndices(from);
            int[] idx1x = m1.get(count);
            int[] idx2 = g2.getNodeIndices(from);
            int[] idx2x = m2.get(count);
            if (count % 3 == 0) {
                for (int i = 0; i < idx1.length; i ++)
                    assertEquals(true, idx1[i] == idx1x[i + 1]);
                for (int i = 0; i < idx2.length; i ++)
                    assertEquals(true, idx2[i] == idx2x[i + 1]);
            } else {
                for (int i = 0; i < idx2.length; i ++)
                    assertEquals(true, idx2[i] == idx2x[i]);
            }
        }
    }

    @BeforeEach
    void setUp() {
        g1 = new IdxEdgeGraph.DefaultGraph(10, true, false);
        g2 = new IdxEdgeGraph.DefaultGraph(100, false, false);
        g3 = new IdxEdgeGraph.DefaultGraph(200, false, true);
    }
}