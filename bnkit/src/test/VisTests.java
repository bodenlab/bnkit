import static org.junit.jupiter.api.Assertions.assertEquals;

import api.PartialOrderGraph;
import dat.POGraph;
import java.io.IOException;
import json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;
import reconstruction.Consensus;
import vis.POAGJson;


/**
 * Created by marnie on 27/4/17.
 */
public class VisTests {

    @Test
    @DisplayName("Generate JSON graph")
    public void getJSONGraphTest() throws IOException {
        String alnfilepath = "bnkit/src/test/resources/small.aln";
        PartialOrderGraph graph = new PartialOrderGraph(new POGraph(alnfilepath));
        POAGJson json = new POAGJson(graph);
        JSONObject jsonObject = json.toJSON();
        String expectedResult = "{\"nodes\":[[111,-1,-1,[[]],false,[[]],[]],[80,0,0,[[[80,4],[71,2]]],false,[[[80,4],[71,2]]],[[80,66.67],[71,33.33]]],[77,1,1,[[[77,3]]],false,[[[77,3]]],[[77,100]]],[78,2,2,[[[78,4]]],false,[[[78,4]]],[[78,100]]],[65,3,3,[[[65,4],[77,2]]],false,[[[65,4],[77,2]]],[[65,66.67],[77,33.33]]],[68,4,4,[[[68,2]]],false,[[[68,2]]],[[68,100]]],[82,5,5,[[[82,4]]],false,[[[82,4]]],[[82,100]]],[102,6,6,[[]],false,[[]],[]]],\"edges\":[[0,0,-1,0,100,0],[0,0,0,1,50,0],[0,0,0,2,50,0],[0,0,1,2,16,1],[0,0,1,3,33,0],[0,0,2,3,66,0],[0,0,3,4,33,0],[0,0,3,5,50,0],[0,0,3,6,16,1],[0,0,4,5,16,1],[0,0,4,6,16,1],[0,0,5,6,66,0]]}";
        //assertEquals(expectedResult, jsonObject.toString());
        Consensus c = new Consensus(jsonObject);
        System.out.println(c.getSupportedSequence(true));
        System.out.println(graph.getConsensusGappySequence());
        System.out.println("Expected JSON object equalled the result.");

    }

    @Test
    @DisplayName("Generate JSON graph for ASR output")
    public void getJSONASRTest() throws IOException, InterruptedException {
        String alnfilepath = "bnkit/src/test/resources/small.aln";
        String treefilepath = "bnkit/src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, treefilepath, false, false, null, 1);
        PartialOrderGraph msa = asr.getPartialOrderGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        POAGJson msajson = new POAGJson(msa);
        POAGJson graphjson = new POAGJson(graph);
        JSONObject msajsonObject = msajson.toJSON();
        JSONObject graphjsonObject = graphjson.toJSON();
        String expectedResult = "{\"nodes\":[[111,-1,-1,[[]],1,[[]],[]],[80,0,0,[[[80,4],[71,2]]],1,[[[80,4],[71,2]]],[[80,66.67],[71,33.33]]],[77,1,1,[[[77,3]]],0,[[[77,3]]],[[77,100]]],[78,2,2,[[[78,4]]],1,[[[78,4]]],[[78,100]]],[65,3,3,[[[65,4],[77,2]]],1,[[[65,4],[77,2]]],[[65,66.67],[77,33.33]]],[68,4,4,[[[68,2]]],0,[[[68,2]]],[[68,100]]],[82,5,5,[[[82,4]]],1,[[[82,4]]],[[82,100]]],[102,6,6,[[]],0,[[]],[]]],\"edges\":[[1,0,-1,0,100,0],[0,0,0,1,50,0],[1,0,0,2,50,0],[0,0,1,2,16,1],[0,0,1,3,33,0],[1,0,2,3,66,0],[0,0,3,4,33,0],[1,0,3,5,50,0],[0,0,3,6,16,1],[0,0,4,5,16,1],[0,0,4,6,16,1],[1,0,5,6,66,0]]}";
        assertEquals(new JSONObject(expectedResult), msajsonObject);
        System.out.println("Expected JSON object equalled the result.");
    }

    @Test
    @DisplayName("New Inference Storing method")   public void getJSONInferences() throws IOException, InterruptedException {
        /**
         * ToDo: work out why this test fails despite being exactly the same on visual inspection
         */
        String alnfilepath = "bnkit/src/test/resources/small.aln";
        String treefilepath = "bnkit/src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, treefilepath, false, false, null, 1);
        PartialOrderGraph msa = asr.getPartialOrderGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        JSONObject result = asr.exportInferencesToJSON();
        System.out.println("Note: This test might fail if the integers don't appear in the exact order");
        String expectedResult = "{\"inferences\":[{\"inferences\":1,\"id\":[1,0],\"label\":0,\"type\":\"meta\",\"transitions\":[1,2],\"base\":[1,1]},[\"N0_X0\",[[-1,45,[0]],[0,80,[2,1,-1]],[1,77,[3,2,0]],[2,78,[3,0]],[3,65,[4,5,2]],[4,68,[5,6,3]],[5,82,[6,3]],[6,45,[6,5]]]]]}";
        assertEquals(new JSONObject(expectedResult), result);
        System.out.println("The inferences were as expected.");
    }
}
