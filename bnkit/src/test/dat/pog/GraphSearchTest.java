package dat.pog;

import asr.ASRException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertTrue;

class GraphSearchTest {

    IdxEdgeGraph<POGraph.StatusEdge> dijk = null;
    IdxEdgeGraph<POGraph.StatusEdge> pog1 = null;

    /*
       /-B-6--E-\
     10  8    3  11
    A-3--C-9--F--8-H
      8  4    1  5
       \-D-7--G-/
     */
    @BeforeEach
    void setupDijkstraExample() {
        dijk = new IdxEdgeGraph<POGraph.StatusEdge>(8, true, false);
        String[] nodenames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H"};
        for (int i = 0; i < nodenames.length; i ++) {
            Node n = new POGraph.StatusNode();
            n.setLabel(nodenames[i]);
            dijk.addNode(i, n);
        }
        dijk.addEdge(-1, 0 /*A*/, new POGraph.StatusEdge(true, 0));
        dijk.addEdge(0 /*A*/, 1 /*B*/, new POGraph.StatusEdge(true, 10));
        dijk.addEdge(0 /*A*/, 2 /*C*/, new POGraph.StatusEdge(true, 3));
        dijk.addEdge(0 /*A*/, 3 /*D*/, new POGraph.StatusEdge(true, 8));
        dijk.addEdge(2 /*C*/, 1 /*B*/, new POGraph.StatusEdge(true, 8));
        dijk.addEdge(1 /*B*/, 4 /*E*/, new POGraph.StatusEdge(true, 6));
        dijk.addEdge(2 /*C*/, 5 /*F*/, new POGraph.StatusEdge(true, 9));
        dijk.addEdge(2 /*C*/, 3 /*D*/, new POGraph.StatusEdge(true, 4));
        dijk.addEdge(3 /*D*/, 6 /*G*/, new POGraph.StatusEdge(true, 7));
        dijk.addEdge(6 /*G*/, 5 /*F*/, new POGraph.StatusEdge(true, 1));
        dijk.addEdge(5 /*F*/, 7 /*H*/, new POGraph.StatusEdge(true, 8));
        dijk.addEdge(5 /*F*/, 4 /*E*/, new POGraph.StatusEdge(true, 3));
        dijk.addEdge(6 /*G*/, 7 /*H*/, new POGraph.StatusEdge(true, 5));
        dijk.addEdge(7 /*H*/, 4 /*E*/, new POGraph.StatusEdge(true, 11));
        dijk.addEdge(4 /*E*/, 8 /*end*/, new POGraph.StatusEdge(true, 0));
    }

    @BeforeEach
    void setupPOGExample1() {
        pog1 = new IdxEdgeGraph<POGraph.StatusEdge>(8, false, true);
        String[] nodenames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H"};
        for (int i = 0; i < nodenames.length; i ++) {
            Node n = new POGraph.StatusNode();
            n.setLabel(nodenames[i]);
            pog1.addNode(i, n);
        }
        pog1.addEdge(-1, 0 /*A*/, new POGraph.StatusEdge(true));
        pog1.addEdge(-1, 2 /*C*/, new POGraph.StatusEdge(false));
        pog1.addEdge(0 /*A*/, 1 /*B*/, new POGraph.StatusEdge(true));
        pog1.addEdge(0 /*A*/, 2 /*C*/, new POGraph.StatusEdge(false));
        pog1.addEdge(0 /*A*/, 3 /*D*/, new POGraph.StatusEdge(true));
        pog1.addEdge(1 /*C*/, 2 /*B*/, new POGraph.StatusEdge(true));
        pog1.addEdge(1 /*B*/, 4 /*E*/, new POGraph.StatusEdge(false));
        pog1.addEdge(2 /*C*/, 4 /*F*/, new POGraph.StatusEdge(false));
        pog1.addEdge(2 /*C*/, 3 /*D*/, new POGraph.StatusEdge(false));
        pog1.addEdge(3 /*D*/, 6 /*G*/, new POGraph.StatusEdge(true));
        pog1.addEdge(5 /*G*/, 6 /*F*/, new POGraph.StatusEdge(false));
        pog1.addEdge(5 /*F*/, 7 /*H*/, new POGraph.StatusEdge(true));
        pog1.addEdge(4 /*F*/, 5 /*E*/, new POGraph.StatusEdge(true));
        pog1.addEdge(6 /*G*/, 7 /*H*/, new POGraph.StatusEdge(false));
        pog1.addEdge(4 /*H*/, 7 /*E*/, new POGraph.StatusEdge(false));
        pog1.addEdge(4 /*E*/, 8 /*end*/, new POGraph.StatusEdge(false));
        pog1.addEdge(7 /*E*/, 8 /*end*/, new POGraph.StatusEdge(true));
    }

    @Test
    void dijkstra1() {
        DijkstraSearch<POGraph.StatusEdge> search = new DijkstraSearch<>(dijk, 0, 4);
        System.out.print("Path: ");
        for (int n : search.getOnePath()) {
            System.out.print(n + " " + search.getCost(n) + "\t");
        }
        System.out.println("\tTotal cost: " + search.getCost());
    }
    @Test
    void dijkstra2() {
        for (POGraph.StatusEdge e : dijk.getEdges().values())
            e.setWeight(e.getWeight() <=7 ? 1 : 2);
    /*
       /-B-1--E-\
      4  2    1  2
    A-1--C-2--F--2-H
      2  1    1  1
       \-D-1--G-/
     */
        DijkstraSearch<POGraph.StatusEdge> search = new DijkstraSearch<>(dijk, 0, 4);
        System.out.println("Paths are made up of: ");
        for (POGEdge e : search.getOptimal()) {
            System.out.println("\t" + e);
        }
        System.out.println();
        dijk.getEdge(0, 1).setWeight(4); // A -- B = 4
        dijk.getEdge(3, 6).setWeight(0); // D -- G = 0
        dijk.getEdge(2, 5).setWeight(3); // C -- F = 3
    /*
       /-B-1--E-\
      4  2    1  2
    A-1--C-3--F--2-H
      2  1    1  1
       \-D-0--G-/
     */
        search = new DijkstraSearch<>(dijk, 0, 4);
        System.out.println("Paths are made up of: ");
        for (POGEdge e : search.getOptimal()) {
            System.out.println("\t" + e);
        }
        System.out.println();
    }

    @Test
    void pog1() {
        DijkstraSearch<POGraph.StatusEdge> search = new DijkstraSearch<>(pog1, -1, 8, DijkstraSearch.PRIORITY_RECIPROCATED);
        System.out.print("POG1 Path: ");
        for (int n : search.getOnePath()) {
            System.out.print(n + " " + search.getCost(n) + "\t");
        }
        System.out.println("\tTotal cost: " + search.getCost());
        System.out.println("POG1 Paths are made up of: ");
        for (POGEdge e : search.getOptimal()) {
            System.out.println("\t" + e);
        }
        System.out.println();

    }

    int N = 30;
    POGraph[] pogs = new POGraph[N];

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
    }

    @Test
    void getMostSupported() {
        setupPOGS();
        Random rand = new Random(N);
        int cnt = 0;
        for (int i = 0; i < pogs.length; i ++) {
            POGraph p = pogs[i];
            // 1. calc best path S -> E
            // 2. remove one edge in best path
            // 3. recalc best path, repeat from 1 until no path from S to E can be found, ensure that best path score decreases monotonically each time
            //System.out.println(p);
            double score = -1;
            while (true) {
                DijkstraSearch dsearch = new DijkstraSearch(p);
                int[] path = dsearch.getOnePath();
                if (path == null)
                    break;
                double prev = score;
                score = dsearch.getCost();
                if (path.length < 2 || score < 0)
                    break;
                int start = rand.nextInt(path.length - 1); // decide on what to remove from the best path
                p.disableEdge(path[start], path[start + 1]); // cripple the graph
                assertTrue(score >= prev); // cost should go up (or stay the same) by crippling the graph
                cnt += 1;
            }
            assertTrue(score > -1); // at least the original POG was solved
        }
        System.out.println(cnt + " tests across " + pogs.length + " POGs");
    }

    @Test
    void getMostSupported2() {
        setupPOGS();
        Random rand = new Random(N);
        int cnt = 0;
        for (int i = 0; i < pogs.length; i ++) {
            POGraph p = pogs[i];
            // 1. calc best path S -> E
            // 2. add a _better_ edge to best path
            // 3. recalc best path, repeat from 1 until no path from S to E can be found, ensure that best path score increases monotonically each time
            double score = Double.POSITIVE_INFINITY;
            while (true) {
                DijkstraSearch dsearch = new DijkstraSearch(p);
                int[] path = dsearch.getOnePath();
                if (path == null)
                    break;
                double prev = score;
                score = dsearch.getCost();
                if (path.length < 3 || score < 0)
                    break;
                assertTrue(score < prev); // cost should go down
                cnt += 1;
                int startidx = rand.nextInt(path.length - 2);
                int startnode = path[startidx];
                int middlenode = path[startidx + 1];
                int endnode = startidx >= path.length - 2 ? p.maxsize() : path[startidx + 2];
                POGraph.StatusEdge mimicme = rand.nextBoolean() ? p.getEdge(startnode, middlenode) : p.getEdge(middlenode, endnode);
                p.addEdge(startnode, endnode, new POGraph.StatusEdge(mimicme.getReciprocated(), mimicme.getWeight() * 0.99));
            }
            assertTrue(score < Double.POSITIVE_INFINITY); // at least the original POG was solved
        }
        System.out.println(cnt + " tests across " + pogs.length + " POGs");
    }

}