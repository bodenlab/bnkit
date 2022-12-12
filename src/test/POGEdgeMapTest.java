import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import reconstruction.POGEdgeMap;

import static org.junit.jupiter.api.Assertions.*;

class POGEdgeMapTest {

    static int[][] pairs = {{1, 2}, {2, 3}, {1, 2}, {3, 5}, {2, 4}, {4, 6}, {6, 7}, {5, 7}};
    static int[][] pairs2 = {{1, 2}, {2, 1}, {1, 2}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {3, 4}, {5, 4}};

    @BeforeAll
    public static void setUp() throws Exception {
    }

    @Test
    void add() {
        POGEdgeMap m = new POGEdgeMap();
        for (int i = 0; i < pairs.length; i ++)
            m.add(pairs[i][0], pairs[i][1]);
        assertEquals(pairs.length - 1, m.size());
        assertEquals(1, m.size(true));
        assertEquals(3, m.size(2));
        POGEdgeMap m2 = new POGEdgeMap();
        for (int i = 0; i < pairs2.length; i ++)
            m2.add(pairs2[i][0], pairs2[i][1]);
        assertEquals(2, m2.size(true));
    }

    @Test
    void getLinkedIndices() {
        POGEdgeMap m = new POGEdgeMap();
        for (int i = 0; i < pairs.length; i ++)
            m.add(pairs[i][0], pairs[i][1]);
        assertEquals(true, m.getLinkedIndices(2).contains(1));
        assertEquals(true, m.getLinkedIndices(2).contains(4));
        assertEquals(false, m.getLinkedIndices(2).contains(2));
    }

    @Test
    void isEdge() {
        POGEdgeMap m = new POGEdgeMap();
        for (int i = 0; i < pairs.length; i ++)
            m.add(pairs[i][0], pairs[i][1]);
        assertEquals(true, m.isEdge(2, 3));
        assertEquals(false, m.isEdge(2, 3, true));
        assertEquals(false, m.isEdge(4, 7));
        assertEquals(true, m.isEdge(4, 6));
    }
}