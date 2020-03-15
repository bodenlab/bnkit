package dat.pog;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;

import static org.junit.jupiter.api.Assertions.*;

class IdxGraphTest {

    IdxGraph g1 = null;
    IdxGraph g2 = null;

    @Test
    void addNode() {
        for (int i = -5; i < g1.size() + 5; i ++) {
            assertEquals(false, g1.isNode(i));
            try {
                g1.addNode(i);
            } catch (RuntimeException e) { }
            assertEquals(g1.isIndex(i), g1.isNode(i));
        }
    }

    @Test
    void removeNode() {
        for (int i = 0; i < g1.size() - 5; i ++)
            g1.addNode(i);
        for (int i = g1.size() - 1; i >= 0; i --) {
            assertEquals(i < g1.size() - 5, g1.isNode(i));
            try {
                g1.removeNode(i);
            } catch (RuntimeException e) { }
            assertEquals(false, g1.isNode(i));
        }
    }

    void createGraphs() {
        for (int i = 0; i < g1.size(); i++) {
            g1.addNode(i);
            g2.addNode(i);
        }
        for (int from = 0; from < g1.size(); from += 2) {
            for (int to = 1; to < g1.size(); to += 2) {
                g1.addEdge(from, to); // un-directed
                g2.addEdge(from, to); // directed
            }
        }
    }

    @Test
    void addEdge() {
        createGraphs();
        for (int from = 0; from < g1.size(); from += 1) {
            for (int to = 0; to < g1.size(); to += 1) {
                if (to % 2 == 1 && from % 2 == 0) {
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

    void crippleGraphs() {
        int count = 0;
        for (int from = 0; from < g1.size(); from += 1) {
            for (int to = 0; to < g1.size(); to += 1) {
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

    @Test
    void removeEdge() {
        createGraphs();
        crippleGraphs();
        int count = 0;
        for (int from = 0; from < g1.size(); from += 1) {
            for (int to = 0; to < g1.size(); to += 1) {
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
        for (int from = 0; from < g1.size(); from += 1) {
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
        for (int from = 0; from < g1.size(); from += 1) {
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
        g1 = new IdxGraph(10, true, false);
        g2 = new IdxGraph(100, false, false);
    }
}