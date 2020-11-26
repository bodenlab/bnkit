package dat.pog;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class GraphSearchTest {

    IdxEdgeGraph<POGraph.StatusEdge> dijk = null;

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

    @Test
    void dijkstra1() {
        DijkstraSearch<POGraph.StatusEdge> search = new DijkstraSearch<>(dijk, 0, 4);
        System.out.print("Path: ");
        for (int n : search.getOnePath()) {
            System.out.print(n + " " + search.getCost(n) + "\t");
        }
        System.out.println();
    }
    @Test
    void dijkstra2() {
        DijkstraSearch<POGraph.StatusEdge> search = new DijkstraSearch<>(dijk, 0, 4);
        System.out.println("Paths are made up of: ");
        for (POGEdge e : search.getOptimal()) {
            System.out.println("\t" + e);
        }
        System.out.println();
    }
}