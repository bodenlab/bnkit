package dat.pog;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

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

    @Test
    void getStarts() {
        IdxGraph g0 = new IdxGraph(10, false, true);
        for (int i = 0; i < 10; i ++)
            g0.addNode(i, new Node());
        g0.addEdge(-1, 0);
        g0.addEdge(0, 1);
        g0.addEdge(0, 2);
        g0.addEdge(-1, 3);
        g0.addEdge(3, 4);
        g0.addEdge(2, 5);
        g0.addEdge(4, 5);
        g0.addEdge(5, 6);
        g0.addEdge(6, 7);
        g0.addEdge(7, 8);
        g0.addEdge(8, 9);
        g0.addEdge(6, 10);
        g0.addTerminalEdge(9);
        int[] s = g0.getStarts();
        int[] e = g0.getEnds();
        assertTrue(s.length == 2);
        assertTrue(e.length == 2);
    }


    @Test
    void getOrdered() {
        IdxGraph g0 = new IdxGraph(10, false, true);
        for (int i = 0; i < 10; i ++)
            g0.addNode(i, new Node());
        g0.addEdge(-1, 0);
        g0.addEdge(0, 1);
        g0.addEdge(0, 2);
        g0.addEdge(1, 3);
        g0.addEdge(3, 4);
        g0.addEdge(2, 5);
        g0.addEdge(4, 5);
        g0.addEdge(5, 6);
        g0.addEdge(6, 7);
        g0.addEdge(7, 8);
        g0.addEdge(8, 9);
        g0.addTerminalEdge(9);
        int[] ord = g0.getOrdered(1, true);
        for (int i = 0; i < ord.length; i ++) {
            System.out.print(ord[i] + "\t");
            if (i > 0)
                assertTrue(ord[i-1] < ord[i]);
        }
        System.out.println();
        assertEquals(3, ord.length);
        ord = g0.getOrdered(2, true);
        for (int i = 0; i < ord.length; i ++) {
            System.out.print(ord[i] + "\t");
        }
        System.out.println();
        assertEquals(1, ord.length);
        ord = g0.getOrdered(9, false);
        for (int i = 0; i < ord.length; i ++) {
            System.out.print(ord[i] + "\t");
            if (i > 0)
                assertTrue(ord[i-1] > ord[i]);
        }
        System.out.println();
        assertEquals(4, ord.length);
        ord = g0.getOrdered(1, false);
        for (int i = 0; i < ord.length; i ++)
            System.out.print(ord[i] + "\t");
        System.out.println();
        assertEquals(1, ord.length);
        ord = g0.getOrdered(0, false);
        for (int i = 0; i < ord.length; i ++)
            System.out.print(ord[i] + "\t");
        System.out.println();
        assertEquals(0, ord.length);
    }

    @Test
    void getOrderedWrap() {
        IdxGraph g0 = new IdxGraph(11, false, true);
        for (int i = 0; i < 11; i ++)
            g0.addNode(i, new Node());
        g0.addEdge(-1, 0);
        g0.addEdge(0, 1);
        g0.addEdge(0, 2);
        g0.addEdge(-1, 6);
        g0.addEdge(1, 3);
        g0.addEdge(3, 4);
        g0.addEdge(2, 5);
        g0.addEdge(4, 5);
        g0.addEdge(5, 6);
        g0.addEdge(6, 7);
        g0.addEdge(7, 8);
        g0.addEdge(8, 9);
        g0.addEdge(6, 10);
        g0.addTerminalEdge(9);
        g0.addTerminalEdge(10);
        int idx = 0;
        int[] wrap = g0.getOrderedWrap(idx, true); //
        for (int i = 0; i < wrap.length; i ++)
            System.out.println(idx + "\t+" + (i + 1) + "\t" + wrap[i]);
        assertTrue(wrap.length == 4);
        idx = 6;
        wrap = g0.getOrderedWrap(idx, true); //
        for (int i = 0; i < wrap.length; i ++)
            System.out.println(idx + "\t+" + (i + 1) + "\t" + wrap[i]);
        assertTrue(wrap.length == 4);

    }

    @Test
    void isIndexTopologicalOrder() {
        IdxGraph g0 = new IdxGraph(11, false, true);
        for (int i = 0; i < 11; i ++)
            g0.addNode(i, new Node());
        g0.addEdge(-1, 0);
        g0.addEdge(0, 1);
        g0.addEdge(0, 2);
        g0.addEdge(-1, 6);
        g0.addEdge(1, 3);
        g0.addEdge(3, 4);
        g0.addEdge(2, 5);
        g0.addEdge(4, 5);
        g0.addEdge(5, 6);
        g0.addEdge(6, 7);
        g0.addEdge(7, 8);
        g0.addEdge(8, 9);
        g0.addEdge(6, 10);
        g0.addTerminalEdge(9);
        g0.addTerminalEdge(10);
        assertTrue(g0.isIndexTopologicalOrder());
        g0 = new IdxGraph(11, false, true);
        for (int i = 10; i >= 0; i --)
            g0.addNode(i, new Node());
        g0.addEdge(-1, 0);
        g0.addEdge(-1, 6);
        g0.addEdge(4, 5);
        g0.addEdge(5, 7);
        g0.addEdge(7, 8);
        g0.addEdge(8, 9);
        g0.addEdge(6, 10);
        g0.addTerminalEdge(9);
        assertTrue(g0.isIndexTopologicalOrder());
        g0.addTerminalEdge(10);
        g0.addEdge(0, 1);
        g0.addEdge(0, 2);
        g0.addEdge(1, 3);
        g0.addEdge(3, 4);
        g0.addEdge(2, 5);
        g0.addEdge(7, 6);
        assertFalse(g0.isIndexTopologicalOrder());
    }

    @Test
    void toJSON() {
        createGraphs();
        System.out.println(g1.toJSON().toString());
        System.out.println(g2.toJSON().toString());
        System.out.println(g3.toJSON().toString());
        HashSet<IdxGraph> gs = new HashSet<>();
        gs.add(g1);
        gs.add(g2);
        gs.add(g3);
        System.out.println(IdxGraph.toJSONArray(gs).toString());
    }

    @Test
    void fromJSON() {
        createGraphs();
        assertTrue(g1.equals(IdxGraph.fromJSON(g1.toJSON())));
        assertTrue(g2.equals(IdxGraph.fromJSON(g2.toJSON())));
        assertTrue(g3.equals(IdxGraph.fromJSON(g3.toJSON())));
    }

    @Test
    void getEdgeIndex() {
        createGraphs();
        IdxGraph[] gs = {g1, g2, g3};
        for (IdxGraph g : gs) {
            System.out.println("Graph maxsize: " + g.maxsize() + "\tDirected: " + g.isDirected() + "\tTerminated: " + g.isTerminated());
            Set<Integer> used = new HashSet<>();
            for (int from = -1; from <= g.maxsize(); from += 1) {
                for (int to = -1; to <= g.maxsize(); to += 1) {
                    if (g.isEdge(from, to)) {
                        int idx = g.getEdgeIndex(from, to);
                        if (g.isDirected())
                            assertFalse(used.contains(idx)); // indices MUST be unique if edges are directed (e.g. 0 -> 1 is different from 1 -> 0)
                        used.add(idx);
                        System.out.println("Edge index: " + idx + "\tFrom: " + from + "\tTo: " + to);
                        if (g.isDirected()) {
                            assertEquals(from, g.getFrom(idx));
                            assertEquals(to, g.getTo(idx));
                        } else {
                            assertEquals(Math.min(from, to), g.getFrom(idx));
                            assertEquals(Math.max(from, to), g.getTo(idx));
                        }
                    }
                }
            }
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