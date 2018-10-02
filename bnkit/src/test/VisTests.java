

import static org.junit.jupiter.api.Assertions.assertEquals;

import api.PartialOrderGraph;
import dat.POGraph;
import java.io.IOException;
import json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;
import vis.POAGJson;


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
        String expectedResult = "{\"nodes\":[{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"P\",\"value\":4},{\"label\":\"G\",\"value\":2}]},\"x\":0,\"y\":0,\"id\":0,\"label\":\"PG\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"P\",\"value\":66.66666666666666},{\"x_label\":\"G\",\"value\":33.33333333333333}]},\"seq\":{\"chars\":[{\"label\":\"P\",\"value\":4},{\"label\":\"G\",\"value\":2}]}},{\"consensus\":false,\"mutants\":{\"chars\":[]},\"x\":-1,\"y\":0,\"id\":-1,\"label\":\"initial\",\"class\":\"\",\"lane\":0,\"seq\":{\"chars\":[]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"M\",\"value\":3}]},\"x\":1,\"y\":0,\"id\":1,\"label\":\"M\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"M\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"M\",\"value\":3}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"N\",\"value\":4}]},\"x\":2,\"y\":0,\"id\":2,\"label\":\"N\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"N\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"N\",\"value\":4}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"A\",\"value\":4},{\"label\":\"M\",\"value\":2}]},\"x\":3,\"y\":0,\"id\":3,\"label\":\"AM\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"A\",\"value\":66.66666666666666},{\"x_label\":\"M\",\"value\":33.33333333333333}]},\"seq\":{\"chars\":[{\"label\":\"A\",\"value\":4},{\"label\":\"M\",\"value\":2}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"D\",\"value\":2}]},\"x\":4,\"y\":0,\"id\":4,\"label\":\"D\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"D\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"D\",\"value\":2}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"R\",\"value\":4}]},\"x\":5,\"y\":0,\"id\":5,\"label\":\"R\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"R\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"R\",\"value\":4}]}},{\"consensus\":false,\"mutants\":{\"chars\":[]},\"x\":6,\"y\":0,\"id\":6,\"label\":\"final\",\"class\":\"\",\"lane\":0,\"seq\":{\"chars\":[]}}],\"max_depth\":0,\"edges\":{\"edges_1:3\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":1,\"y2\":0,\"weight\":33,\"from\":1,\"x2\":3,\"to\":3},\"edges_1:2\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":1,\"y2\":0,\"weight\":16,\"from\":1,\"x2\":2,\"to\":2},\"edges_-1:0\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":-1,\"y2\":0,\"weight\":100,\"from\":-1,\"x2\":0,\"to\":0},\"edges_2:3\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":2,\"y2\":0,\"weight\":66,\"from\":2,\"x2\":3,\"to\":3},\"edges_3:5\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":3,\"y2\":0,\"weight\":50,\"from\":3,\"x2\":5,\"to\":5},\"edges_3:4\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":3,\"y2\":0,\"weight\":33,\"from\":3,\"x2\":4,\"to\":4},\"edges_4:6\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":4,\"y2\":0,\"weight\":16,\"from\":4,\"x2\":6,\"to\":6},\"edges_3:6\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":3,\"y2\":0,\"weight\":16,\"from\":3,\"x2\":6,\"to\":6},\"edges_4:5\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":4,\"y2\":0,\"weight\":16,\"from\":4,\"x2\":5,\"to\":5},\"edges_5:6\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":5,\"y2\":0,\"weight\":66,\"from\":5,\"x2\":6,\"to\":6},\"edges_0:2\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":0,\"y2\":0,\"weight\":50,\"from\":0,\"x2\":2,\"to\":2},\"edges_0:1\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":0,\"y2\":0,\"weight\":50,\"from\":0,\"x2\":1,\"to\":1}}}";
        assertEquals(expectedResult, jsonObject.toString());
        System.out.println("Expected JSON object equalled the result.");

    }

    @Test
    @DisplayName("Generate JSON graph for ASR output")
    public void getJSONASRTest() throws IOException, InterruptedException {
        String alnfilepath = "src/test/resources/small.aln";
        String treefilepath = "src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, treefilepath, false, false, null, 1);
        PartialOrderGraph msa = asr.getPartialOrderGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        POAGJson msajson = new POAGJson(msa);
        POAGJson graphjson = new POAGJson(graph);
        JSONObject msajsonObject = msajson.toJSON();
        JSONObject graphjsonObject = graphjson.toJSON();
        String expectedResult = "{\"nodes\":[{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"P\",\"value\":4},{\"label\":\"G\",\"value\":2}]},\"x\":0,\"y\":0,\"id\":0,\"label\":\"PG\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"P\",\"value\":66.66666666666666},{\"x_label\":\"G\",\"value\":33.33333333333333}]},\"seq\":{\"chars\":[{\"label\":\"P\",\"value\":4},{\"label\":\"G\",\"value\":2}]}},{\"consensus\":false,\"mutants\":{\"chars\":[]},\"x\":-1,\"y\":0,\"id\":-1,\"label\":\"initial\",\"class\":\"\",\"lane\":0,\"seq\":{\"chars\":[]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"M\",\"value\":3}]},\"x\":1,\"y\":0,\"id\":1,\"label\":\"M\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"M\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"M\",\"value\":3}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"N\",\"value\":4}]},\"x\":2,\"y\":0,\"id\":2,\"label\":\"N\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"N\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"N\",\"value\":4}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"A\",\"value\":4},{\"label\":\"M\",\"value\":2}]},\"x\":3,\"y\":0,\"id\":3,\"label\":\"AM\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"A\",\"value\":66.66666666666666},{\"x_label\":\"M\",\"value\":33.33333333333333}]},\"seq\":{\"chars\":[{\"label\":\"A\",\"value\":4},{\"label\":\"M\",\"value\":2}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"D\",\"value\":2}]},\"x\":4,\"y\":0,\"id\":4,\"label\":\"D\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"D\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"D\",\"value\":2}]}},{\"consensus\":false,\"mutants\":{\"chars\":[{\"label\":\"R\",\"value\":4}]},\"x\":5,\"y\":0,\"id\":5,\"label\":\"R\",\"class\":\"\",\"lane\":0,\"graph\":{\"bars\":[{\"x_label\":\"R\",\"value\":100}]},\"seq\":{\"chars\":[{\"label\":\"R\",\"value\":4}]}},{\"consensus\":false,\"mutants\":{\"chars\":[]},\"x\":6,\"y\":0,\"id\":6,\"label\":\"final\",\"class\":\"\",\"lane\":0,\"seq\":{\"chars\":[]}}],\"max_depth\":0,\"edges\":{\"edges_1:3\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":1,\"y2\":0,\"weight\":33,\"from\":1,\"x2\":3,\"to\":3},\"edges_1:2\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":1,\"y2\":0,\"weight\":16,\"from\":1,\"x2\":2,\"to\":2},\"edges_-1:0\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":-1,\"y2\":0,\"weight\":100,\"from\":-1,\"x2\":0,\"to\":0},\"edges_2:3\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":2,\"y2\":0,\"weight\":66,\"from\":2,\"x2\":3,\"to\":3},\"edges_3:5\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":3,\"y2\":0,\"weight\":50,\"from\":3,\"x2\":5,\"to\":5},\"edges_3:4\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":3,\"y2\":0,\"weight\":33,\"from\":3,\"x2\":4,\"to\":4},\"edges_4:6\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":4,\"y2\":0,\"weight\":16,\"from\":4,\"x2\":6,\"to\":6},\"edges_3:6\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":3,\"y2\":0,\"weight\":16,\"from\":3,\"x2\":6,\"to\":6},\"edges_4:5\":{\"consensus\":false,\"single\":true,\"reciprocated\":false,\"y1\":0,\"x1\":4,\"y2\":0,\"weight\":16,\"from\":4,\"x2\":5,\"to\":5},\"edges_5:6\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":5,\"y2\":0,\"weight\":66,\"from\":5,\"x2\":6,\"to\":6},\"edges_0:2\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":0,\"y2\":0,\"weight\":50,\"from\":0,\"x2\":2,\"to\":2},\"edges_0:1\":{\"consensus\":false,\"single\":false,\"reciprocated\":false,\"y1\":0,\"x1\":0,\"y2\":0,\"weight\":50,\"from\":0,\"x2\":1,\"to\":1}}}";
        assertEquals(msajsonObject.toString(), expectedResult);
        System.out.println("Expected JSON object equalled the result.");
    }

    @Test
    @DisplayName("New Inference Storing method")
    public void getJSONInferences() throws IOException, InterruptedException {
        String alnfilepath = "src/test/resources/small.aln";
        String treefilepath = "src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, treefilepath, false, false, null, 1);
        PartialOrderGraph msa = asr.getPartialOrderGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        JSONObject result = asr.exportInferencesToJSON();
        JSONObject expectedResult = new JSONObject("{\"inferences\":[{\"inferences\":1,\"id\":[1,0],\"label\":0,\"type\":\"meta\",\"transitions\":[1,2],\"base\":[1,1]},[\"N0_X0\",[[-1,45,[0]],[0,80,[1,2,-1]],[1,77,[3,2,0]],[2,78,[3,0]],[3,65,[5,4,2]],[4,68,[5,6,3]],[5,82,[6,3]],[6,45,[6,5]]]]]}");
        assertEquals(result, expectedResult);
        System.out.println("The inferences were as expected.");
    }
}
