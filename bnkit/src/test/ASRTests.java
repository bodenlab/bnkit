import api.PartialOrderGraph;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;

import java.io.IOException;

/**
 * Created by marnie on 27/4/17.
 */
public class ASRTests {

    @Test
    @DisplayName("Get ancestral partial order graph")
    public void getAncestralGraphTest() throws IOException {
        String alnfilepath = "src/test/resources/small.aln";
        String nwkfilepath = "src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true);
        PartialOrderGraph graph = asr.getGraph("root");
      //  asr.saveMSAGraph();
        System.out.println(graph.toString());
    }
}
