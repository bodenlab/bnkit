package dat.pog;

import dat.Enumerable;
import dat.Interval1D;
import dat.IntervalST;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class POGraphTest {

    static int N = 30;
    POGraph pog = null;
    Set<Integer> allNodes = new HashSet<>();

    @BeforeEach
    void setupPOG() {
        pog = new POGraph(N);
        for (int i = 0; i < N; i ++) {
            pog.addNode(i, new EnumNode(Enumerable.aacid));
        }
        int last = -1;
        for (int from = -1; from < N-10; from += 2) {
            int to = Math.min(from + 2, N);
            pog.addEdge(from, to);
            allNodes.add(to);
            last = to;
        }
        if (last < N)
            pog.addEdge(last, N);
        for (int from = -1; from < N-10; from += 7) {
            int to = Math.min(from + 7, N);
            pog.addEdge(from, to);
            allNodes.add(to);
            last = to;
        }
        if (last < N)
            pog.addEdge(last, N);
        for (int from = -1; from < N-10; from += 9) {
            int to = Math.min(from + 9, N);
            pog.addEdge(from, to);
            allNodes.add(to);
            last = to;
        }
        if (last < N)
            pog.addEdge(last, N);
        try {
            pog.saveToDOT("bnkit/src/test/resources/pogtest.dot");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void getForward() {
        assertEquals(3, recurseForward(-1, 0));
    }

    int recurseForward(int idx, int nrounds) {
        int[] next = pog.getForward(idx);
        int shortest = N;
        assertTrue(nrounds < N);
        if (next.length == 0) {
            return nrounds;
        }
        for (int i = 0; i < next.length; i ++) {
            int n = recurseForward(next[i], nrounds + 1);
            if (n < shortest)
                shortest = n;
        }
        return shortest;
    }


    @Test
    void visitForward() {
        Set<Integer> found = recurseForward(-1, new HashSet<>());
        assertTrue(found.containsAll(allNodes));
    }

    Set<Integer> recurseForward(int idx, Set<Integer> mySet) {
        int[] next = pog.getForward(idx);
        if (next.length == 0)
            return mySet;
        for (int i = 0; i < next.length; i ++) {
            mySet.add(next[i]);
            recurseForward(next[i], mySet);
        }
        return mySet;
    }



    @Test
    void getBackward() {
        assertEquals(3, recurseBackward(N, 0));
    }

    int recurseBackward(int idx, int nrounds) {
        int[] before = pog.getBackward(idx);
        int shortest = N;
        assertTrue(nrounds < N);
        if (before.length == 0) {
            return nrounds;
        }
        for (int i = 0; i < before.length; i ++) {
            int n = recurseBackward(before[i], nrounds + 1);
            if (n < shortest)
                shortest = n;
        }
        return shortest;
    }


    @Test
    void visitBackward() {
        Set<Integer> found = recurseBackward(N, new HashSet<>());
        assertTrue(found.containsAll(allNodes));
    }

    Set<Integer> recurseBackward(int idx, Set<Integer> mySet) {
        int[] before = pog.getBackward(idx);
        if (before.length == 0)
            return mySet;
        for (int i = 0; i < before.length; i ++) {
            mySet.add(before[i]);
            recurseBackward(before[i], mySet);
        }
        return mySet;
    }


    @Test
    void getIndels() {
        Set<int[]> ivset = pog.getIndels();
        for (int[] ival : ivset) {
            System.out.println("<" + ival[0] + ", " + ival[1] + ">");
        }
    }

    @Test
    void getSimpleGapCode() {

        POGraph pog1 = new POGraph(20);
        int[] nodes = new int[] {0, 1, 4, 5, 12, 13, 18, 19};
        for (int n : nodes) pog1.addNode(n, new Node());
        pog1.addEdge(-1, nodes[0]);
        for (int i = 0; i < nodes.length - 1; i ++)
            pog1.addEdge(nodes[i], nodes[i + 1]);
        pog1.addEdge(nodes[nodes.length - 1], nodes.length);

        assertEquals(POGraph.GAP_STATUS_PRESENT, pog1.getSimpleGapCode( 1,  4));
        assertEquals(POGraph.GAP_STATUS_PRESENT, pog1.getSimpleGapCode( 5, 12));
        assertEquals(POGraph.GAP_STATUS_PRESENT, pog1.getSimpleGapCode(13, 18));
        assertEquals(POGraph.GAP_STATUS_UNKNOWN, pog1.getSimpleGapCode( 6, 12));
        assertEquals(POGraph.GAP_STATUS_UNKNOWN, pog1.getSimpleGapCode( 1,  3));
        assertEquals(POGraph.GAP_STATUS_UNKNOWN, pog1.getSimpleGapCode( 7, 10));
        assertEquals(POGraph.GAP_STATUS_ABSENT,  pog1.getSimpleGapCode( 5, 13));
        assertEquals(POGraph.GAP_STATUS_ABSENT,  pog1.getSimpleGapCode( 4, 12));
        assertEquals(POGraph.GAP_STATUS_ABSENT,  pog1.getSimpleGapCode( 1, 18));

        assertEquals(POGraph.GAP_STATUS_PRESENT, pog.getSimpleGapCode(13, 20));
        assertEquals(POGraph.GAP_STATUS_PRESENT, pog.getSimpleGapCode( 6, 13));
        assertEquals(POGraph.GAP_STATUS_PRESENT, pog.getSimpleGapCode(21, pog.maxsize()));
        assertEquals(POGraph.GAP_STATUS_PRESENT, pog.getSimpleGapCode(-1,  1));
        assertEquals(POGraph.GAP_STATUS_ABSENT,  pog.getSimpleGapCode( 9, 19));
        assertEquals(POGraph.GAP_STATUS_UNKNOWN,  pog.getSimpleGapCode( 8, 14));

    }

}