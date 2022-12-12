import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;

import api.PartialOrderGraph;
import dat.POGraph;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.Consensus;
import vis.POAGJson;

/**
 * Functions to test the partial order graph implementation (api and code).
 *
 * Created by marnie on 29/3/17.
 */
public class POGraphTests {
    @Test
    @DisplayName("Small Partial Order Graph test")
    public void smallPartialOrderGraph() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphSmall.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        System.out.println("Consensus sequence: " + apiGraph.getConsensusSequence());
        assertEquals("PNAR", apiGraph.getConsensusSequence());
        POAGJson json = new POAGJson(apiGraph);

        /* We want to check that we get the same result reading back in the JSON object */
        JSONObject jsonObject = json.toJSON();
        Consensus c = new Consensus(jsonObject);
        String regeneratedConsensus = c.getSupportedSequence(false);
        assertEquals("PNAR", regeneratedConsensus);

    }

    @Test
    @DisplayName("Small PO Graph test")
    public void smallPOGraph() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphSmall.dot";
        POGraph graph = new POGraph(filepath);
        System.out.println("Consensus sequence: " + graph.getSupportedSequence(false));
        assertEquals("PNAR", graph.getSupportedSequence(false));
        POAGJson jsongraph = new POAGJson(new PartialOrderGraph(graph));
        System.out.println(jsongraph.toJSON().toString());
        Map<Integer, List<Integer>> lists = graph.getSequenceNodeMapping();
        for (Integer nodeId : graph.getNodeIDs()) {
            graph.setCurrent(nodeId);
            System.out.println(nodeId + ":" + graph.getCurrentBase());
        }
        for (Integer seqId : lists.keySet()){
            System.out.println("Seq: " + seqId);
            for (Integer id : lists.get(seqId))
                System.out.println(id);
        }
    }

    @Test
    @DisplayName("MSA PO Graph test")
    public void MSAPOGraph() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphMSAMed.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        for (Integer id : apiGraph.getNodeIDs())
            System.out.println(apiGraph.getLabel(id));
    }

    @Test
    @DisplayName("Load Graph test")
    public void LoadPOGraph() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphMSAFourLevels.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        System.out.println(apiGraph.toString());
    }

    @Test
    @DisplayName("Small PO Graph test edge weights")
    public void smallPOGraphEdgeWeights() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphSmall.dot";
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
    public void smallPOGraphSort() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphSmall.dot";
        PartialOrderGraph apiGraph = new PartialOrderGraph(filepath);
        Integer[] ids = apiGraph.sort();
        for (Integer id : ids)
            System.out.println(id);
        assertArrayEquals(new Integer[]{0,1,2,3,4,5}, ids);
    }

    @Test
    @DisplayName("Convert PO Graph to API PartialOrderGraph")
    public void graphToAPITest() throws IOException {
        String filepath = "bnkit/src/test/resources/small.aln";
        POGraph graph = new POGraph(filepath);
        graph.saveSequences("bnkit/src/text/resources/small_graph.aln", "clustal");
        PartialOrderGraph apiGraph = new PartialOrderGraph(graph);
        System.out.println(graph.toString());
        System.out.println(apiGraph.toString());
        assertEquals(graph.toString(), apiGraph.toString());
    }

    @Test
    @DisplayName("Convert PO Graph to JSON")
    public void graphToJSONTest() throws IOException {
        String filepath = "bnkit/src/test/resources/testPOGraphMSAEightLevels.dot";
        PartialOrderGraph graph = new PartialOrderGraph(filepath);
        POAGJson json = new POAGJson(graph);
        JSONObject obj = json.toJSON();
        System.out.println(obj.toString());
    }


}
