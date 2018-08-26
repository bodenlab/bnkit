package test;


import api.PartialOrderGraph;
import dat.POGraph;
import json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;
import vis.POAGJson;

import java.io.IOException;


/**
 * Created by marnie on 27/4/17.
 */
public class VisTests {

    @Test
    @DisplayName("Generate JSON graph")
    public void getJSONGraphTest() throws IOException {
        String alnfilepath = "src/test/resources/small.aln";
        PartialOrderGraph graph = new PartialOrderGraph(new POGraph(alnfilepath));
        POAGJson json = new POAGJson(graph);
        JSONObject jsonObject = json.toJSON();
        System.out.println(jsonObject.toString());
    }

    @Test
    @DisplayName("Generate JSON graph for ASR output")
    public void getJSONASRTest() throws IOException, InterruptedException {
        String alnfilepath = "src/test/resources/small.aln";
        String treefilepath = "src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, treefilepath, false, "None", false, null, 1);
        PartialOrderGraph msa = asr.getPartialOrderGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        POAGJson msajson = new POAGJson(msa);
        POAGJson graphjson = new POAGJson(graph);
        JSONObject msajsonObject = msajson.toJSON();
        JSONObject graphjsonObject = graphjson.toJSON();
        System.out.println(msajsonObject.toString());
    }
}
