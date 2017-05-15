import api.PartialOrderGraph;
import json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;
import vis.POAGJson;

import java.io.IOException;

/**
 * Created by marnie on 27/4/17.
 */
public class ASRTests {

    @Test
    @DisplayName("Get ancestral partial order graph")
    public void getAncestralGraphTest() throws IOException {
        String alnfilepath = "src/test/resources/large.aln";
        String nwkfilepath = "src/test/resources/large.nwk";
        ASRPOG asr = new ASRPOG(null, nwkfilepath, alnfilepath, "N3_72.0", true);
        PartialOrderGraph msa = asr.getMSAGraph();
        PartialOrderGraph graph = asr.getGraph("N3_72.0");
        System.out.println(msa.toString());
        System.out.println(graph.toString());

        POAGJson msajson = new POAGJson(msa);
        JSONObject msaObj = msajson.toJSON();
        POAGJson json = new POAGJson(graph);
        JSONObject jsonObj = json.toJSON();
        System.out.println(jsonObj.toString());
    }
}
