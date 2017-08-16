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
        String alnfilepath = "src/test/resources/2U1_aligned_trimmed.aln";
        String nwkfilepath = "src/test/resources/2U1_final.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true);
        asr.saveALN("src/test/cyp2u1", "fasta");
        PartialOrderGraph msa = asr.getMSAGraph();
        asr.saveMSAGraph("src/test/msa");
        asr.saveGraph("src/test/cyp2u1", "root");
        System.out.println(msa.toString());

        POAGJson msajson = new POAGJson(msa);
        JSONObject msaObj = msajson.toJSON();
        System.out.println(msaObj.toString());
    }

    @Test
    @DisplayName("Small ASR")
    public void performSmallASR() throws IOException {
        String alnfilepath = "src/test/resources/edge1.aln";
        String nwkfilepath = "src/test/resources/edge1.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true);
        PartialOrderGraph msa = asr.getMSAGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        System.out.println(msa.toString());
        System.out.println(graph.toString());
    }
}
