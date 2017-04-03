import api.PartialOrderGraph;
import org.junit.jupiter.api.*;

import java.util.HashMap;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Functions to test the partial order graph implementation (api and code).
 *
 * Created by marnie on 29/3/17.
 */
public class POGraphTests {
    @Test
    @DisplayName("Small PO Graph test")
    public void smallPOGraph(){
        String filepath = "src/test/resources/testPOGraphSmall.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        System.out.println("Consensus sequence: " + apiGraph.getConsensusSequence());
        assertEquals("PNAR", apiGraph.getConsensusSequence());
    }

    @Test
    @DisplayName("Small PO Graph test edge weights")
    public void smallPOGraphEdgeWeights(){
        String filepath = "src/test/resources/testPOGraphSmall.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        for (Integer nodeId : apiGraph.getNodeIDs()) {
            System.out.println(Integer.toString(nodeId) + ":" + apiGraph.getLabel(nodeId));
            HashMap<Integer, Double> edgeWeights = new HashMap<>(apiGraph.getOutEdgeWeights(nodeId));
            for (Integer nextId : edgeWeights.keySet())
                System.out.println(Integer.toString(nodeId) + " -> " + Integer.toString(nextId) + ":" + Double.toString
                        (edgeWeights.get(nextId)));
        }
    }

    @Test
    @DisplayName("Sort small PO Graph")
    public void smallPOGraphSort(){
        String filepath = "src/test/resources/testPOGraphSmall.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        Integer[] ids = apiGraph.sort();
        for (Integer id : ids)
            System.out.println(id);
        assertArrayEquals(new Integer[]{0,1,2,3,4,5}, ids);
    }

}