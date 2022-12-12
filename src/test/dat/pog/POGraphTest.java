package dat.pog;

import asr.ASRException;
import dat.Enumerable;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

class POGraphTest {

    static int N = 30;
    POGraph pog = null;
    POGraph[] pogs = new POGraph[N];
    POGraph dijk = null;
    POGraph ex1 = null;

    Set<Integer> allNodes = new HashSet<>();

    @BeforeEach
    void setupPOG() {
        int[] start = new int[] {0,2};
        int[] end = new int[] {7,9};
        int[][] adj = new int[][] {
            {0,1,0,0,0,0,1,0,0,0},
            {1,0,1,1,0,1,0,0,0,0},
            {0,1,0,0,1,0,0,0,0,0},
            {0,0,0,0,0,1,0,1,0,0},
            {0,0,1,0,0,0,1,0,0,0},
            {0,1,0,0,0,0,0,0,1,0},
            {1,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,1,0,0,1},
            {0,0,0,0,0,1,0,0,0,1},
            {0,0,0,0,0,0,0,1,0,0},
        };
        ex1 = POGraph.createFromAdjacency(start, end, adj);
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
        /*
        try {
            //ex1.saveToDOT("bnkit/src/test/resources/ex1.dot");
            pog.saveToDOT("bnkit/src/test/resources/pogtest.dot");
        } catch (IOException e) {
            e.printStackTrace();
        }
         */
    }

    void setupPOGS() {
        Random rand = new Random(N);
        // we will create N POGs
        for (int i = 0; i < N; i ++) {
            int myN = (i + 2) * 3; // max size of this POG is 3 * 2, 3, ..., N + 1
            pogs[i] = new POGraph(myN);
            pogs[i].setName("P" + i);
            int addN = rand.nextInt(myN) + 1;       // how many nodes we'll add
            Set<Integer> nodes = new HashSet<>();   // nodes will be kept in a set
            for (int j = 0; j < addN; j ++) {
                int n = rand.nextInt(addN);         // the index we'll give it; note NOT ordered
                pogs[i].addNode(n, new Node());     // add it
                nodes.add(n);                       // add the index to a set
            }
            // we are going to assume that the array below is the topological sort of indices
            Integer[] nodearr = new Integer[nodes.size()];
            nodes.toArray(nodearr);
            for (int i1 = 0; i1 < nodearr.length - 1; i1 ++) {  // look at each node in topological order
                int nforw = nodearr.length - i1 - 1;            // how many nodes FORWARD that can be connected
                int npick = Math.min(rand.nextInt(nforw), rand.nextInt(nforw)) + 1;         // pick a number of them, at least 1, at most all
                for (int j = 0; j < npick; j ++) {
                    int i2 = rand.nextInt(nforw) + i1 + 1;      // pick a node index that is AFTER i1
                    boolean recip = rand.nextBoolean();
                    // add edge: even-numbered graphs will have "real" weights (0-1),
                    // odd-numbered, will have either 1, 2 or 3 as weights (ints),
                    // which increases the prob of multiple paths to be optimal
//                    pogs[i].addEdge(nodearr[i1], nodearr[i2], new POGraph.StatusEdge(recip, rand.nextDouble()));
                    pogs[i].addEdge(nodearr[i1], nodearr[i2], new POGraph.StatusEdge(recip, i % 2 == 0 ? rand.nextDouble() + 1: rand.nextInt(3) + 1));
                }
            }
            int onestart = -1;
            for (int j = 0; j < rand.nextInt(nodearr.length) + 1; j ++) { // this is how many start nodes we'll pick (maybe), at least 1, at most all
                int i1 = Math.min(rand.nextInt(nodearr.length), rand.nextInt(nodearr.length));      // this is a random start node
                onestart = i1;
                boolean recip = rand.nextBoolean();
//                pogs[i].addEdge(-1, nodearr[i1], new POGraph.StatusEdge(recip, rand.nextDouble()));
                pogs[i].addEdge(-1, nodearr[i1], new POGraph.StatusEdge(recip, i % 2 == 0 ? rand.nextDouble() + 1 : rand.nextInt(3) + 1));
                if (rand.nextBoolean()) // a chance that we finish prematurely, to bias POGs to have a small number of terminal nodes
                    break;
            }
            for (int j = 0; j < rand.nextInt(nodearr.length - onestart) + 1; j ++) { // this is how many end nodes we'll pick (maybe), at least 1, at most all after last startnode
                int i2 = Math.max(rand.nextInt(nodearr.length - onestart), rand.nextInt(nodearr.length - onestart)) + onestart;         // this is a random end node
                boolean recip = rand.nextBoolean();
//                pogs[i].addEdge(nodearr[i2], myN, new POGraph.StatusEdge(recip, rand.nextDouble()));
                pogs[i].addEdge(nodearr[i2], myN, new POGraph.StatusEdge(recip, i % 2 == 0 ? rand.nextDouble() + 1 : rand.nextInt(3) + 1));
                if (rand.nextBoolean()) // a chance that we finish prematurely, to bias POGs to have a small number of terminal nodes
                    break;
            }
            // Finally ensure that is at least one path from start to end (here: to end from start)
            int current = myN; // the node index used for edge to link a source to destination node
            while (current != -1) {
                int prev = current; // remember the destination
                int[] before = pogs[i].getBackward(current);
                if (before.length > 0) { // there was one or more edges already, pick one
                    current = before[rand.nextInt(before.length)];
                } else {
                    // new edge is needed
                    // find the current node in the topologically sorted list...
                    int k = 0;
                    for (; k < nodearr.length; k ++) {
                        if (nodearr[k] == current) // found it,
                            break;
                    }
                    current = nodearr[rand.nextInt(k + 1)]; // pick one BEFORE k, OR start node
                    if (current == prev)
                        current = -1;
                    boolean recip = rand.nextBoolean();
                    pogs[i].addEdge(current, prev, new POGraph.StatusEdge(recip, i % 2 == 0 ? rand.nextDouble() + 1 : rand.nextInt(3) + 1));
                }
            }
        }
        /*
        try {
            POGraph.saveToDOT("bnkit/src/test/resources/pogstest", pogs);
        } catch (ASRException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
         */
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
    void getDepth() {
        int[] depths = ex1.getDepths(true);
        assertTrue(depths[9] == 1);
        assertTrue(depths[8] == 2);
        assertTrue(depths[5] == 3);
        assertTrue(depths[0] == 3);
        assertEquals(4, depths[2]);
        assertEquals(3, depths[1]);
        assertEquals(2, depths[3]);
        depths = ex1.getDepths(false);
        assertEquals(2, depths[6]);
        assertEquals(2, depths[1]);
        assertEquals(4, depths[8]);
        assertEquals(4, depths[9]);
    }

    @Test
    void getDepth2() {
        int[] depths = pog.getDepths(true);
        // for each node, at least ONE of the next nodes must be at exactly ONE less level of depth, all others one less OR GREATER
        for (int i = 0; i < depths.length; i ++) {
            boolean oneless = false;
            for (int j : pog.getForward(i)) {
                if (depths[j] == depths[i] - 1)
                    oneless = true;
                assertTrue(depths[j] >= depths[i] - 1);
            }
            assertTrue(oneless || pog.getForward(i).length == 0);
        }
        depths = pog.getDepths(false);
        // for each node, at least ONE of the next nodes must be at exactly ONE less level of depth, all others one less OR GREATER
        for (int i = 0; i < depths.length; i ++) {
            boolean oneless = false;
            for (int j : pog.getBackward(i)) {
                if (depths[j] == depths[i] - 1)
                    oneless = true;
                assertTrue(depths[j] >= depths[i] - 1);
            }
            assertTrue(oneless || pog.getBackward(i).length == 0);
        }
    }

    @Test
    void isPath() {
        assertTrue(pog.isPath(6,19));
        assertTrue(pog.isPath(7,26));
        assertFalse(pog.isPath(8,20));
        assertFalse(pog.isPath(19,26));
        assertFalse(pog.isPath(26,7));
        assertTrue(pog.isPath(7,7));
        assertTrue(pog.isPath(-1,pog.size()));
        assertFalse(pog.isPath(0,pog.size()));
    }

    @Test
    void getTopologicalOrder() {
        setupPOGS();
        for (POGraph p : pogs) {
            int[] order = p.getTopologicalOrder();
            assertTrue(p.size() + 1 == order.length);
            for (int i = 0; i < order.length; i ++) {
                for (int j = i + 1; j < order.length; j ++)
                    assertFalse(p.isPath(order[j], order[i]));
            }
        }
    }

    /*
    @Test
    void getMostSupported() {
        setupPOGS();
        Random rand = new Random(N);
        for (int i = 0; i < pogs.length; i ++) {
            POGraph p = pogs[i];
            // 1. calc best path S -> E
            // 2. remove one edge in best path
            // 3. recalc best path, repeat from 1 until no path from S to E can be found, ensure that best path score decreases monotonically each time
            //System.out.println(p);
            double score = -1;
            int cnt = 0;
            while (true) {
                POGSearch astar = p.getMostSupported();
                int[] path = astar.getOnePath();
                if (path == null)
                    break;
                double prev = score;
                score = astar.getCost();
                if (path.length < 2 || score < 0)
                    break;
                int start = rand.nextInt(path.length - 1); // decide on what to remove from the best path
                if (score < prev || p.getName().equals("P12") || p.getName().equals("P14")|| p.getName().equals("P16")|| p.getName().equals("P19")) {
                    System.out.print(cnt + ": Failed with " + p.getName() + " at  " + p.getEdgeCount() + " when " + p + " at cost " + prev + " -> " + score + ": ");
                    for (int j = 0; j < path.length; j++)
                        System.out.print(path[j] + " ");
                    System.out.print("removing " + path[start] + "->" + path[start + 1] + " " + p.getEdge(path[start], path[start + 1]) + "[" + astar.minW + "]");
                    for (int j = 0; j < astar.toposort.length; j++)
                        System.out.print(astar.toposort[j] + " ");
                    System.out.println();
                    try { p.saveToDOT("bnkit/src/test/resources/pogstest/" + p.getName() + "_" + cnt + ".dot"); } catch (IOException e) {}
                    cnt ++;
                }
                p.disableEdge(path[start], path[start + 1]); // cripple the graph
                //assertTrue(score >= prev); // cost should go up (or stay the same) by crippling the graph
            }
            assertTrue(score > -1); // at least the original POG was solved
        }
    }

    @Test
    void getMostSupported2() {
        setupPOGS();
        Random rand = new Random(N);
        for (int i = 0; i < pogs.length; i ++) {
            POGraph p = pogs[i];
            // 1. calc best path S -> E
            // 2. add a better edge to best path
            // 3. recalc best path, repeat from 1 until no path from S to E can be found, ensure that best path score increases monotonically each time
            double score = Double.POSITIVE_INFINITY;
            while (true) {
                POGSearch astar = p.getMostSupported();
                int[] path = astar.getOnePath();
                if (path == null)
                    break;
                double prev = score;
                score = astar.getCost();
                if (path.length < 2 || score < 0)
                    break;
                assertTrue(score < prev); // cost should go down
                int startidx = rand.nextInt(path.length - 1);
                int startnode = path[startidx];
                int middlenode = path[startidx + 1];
                int endnode = startidx >= path.length - 2 ? p.maxsize() : path[startidx + 2];
                POGraph.StatusEdge mimicme = rand.nextBoolean() ? p.getEdge(startnode, middlenode) : p.getEdge(middlenode, endnode);
                p.addEdge(startnode, endnode, new POGraph.StatusEdge(mimicme.getReciprocated(), mimicme.getWeight() * 0.99));
            }
            assertTrue(score < Double.POSITIVE_INFINITY); // at least the original POG was solved
        }
    }
*/

    @Test
    void dijkstra() {

    }

    @Test
    void getOptimal() {
        setupPOGS();

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

    @Test
    void toJSON() {
        setupPOGS();
        try {
            POGraph.saveToJSON("/tmp", Arrays.asList(pogs));
            for (int i = 0; i < pogs.length; i ++) {
                //System.out.println(pogs[i].toJSON());
                //System.out.println(POGraph.fromJSON(pogs[i].toJSON()));
                assertTrue(pogs[i].equals(POGraph.fromJSON(pogs[i].toJSON())));
            }

            List<POGraph> pogs_loaded = POGraph.loadFromJSON("/tmp/pogs.json");
            assertEquals(pogs.length, pogs_loaded.size());
            for (int i = 0; i < pogs.length; i ++) {
                //System.out.println(pogs[i].toJSON());
                //System.out.println(pogs_loaded.get(i).toJSON());
                assertTrue(pogs[i].equals(pogs_loaded.get(i)));
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}