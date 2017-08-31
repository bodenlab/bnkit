import api.PartialOrderGraph;
import dat.POGraph;
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
        String alnfilepath = "src/test/resources/tawfik.aln";
        String nwkfilepath = "src/test/resources/tawfik.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true, false, null, 1);
        //asr.saveALN("src/test/cyp2u1", "fasta");
        POGraph msa = asr.getMSAGraph();
        PartialOrderGraph root = asr.getGraph("root");
        //PartialOrderGraph msa = asr.getPartialOrderGraph();
        POAGJson jsongraph = new POAGJson(root);
        System.out.println(root.getConsensusSequence());
       // asr.saveMSAGraph("src/test/msa");
       // asr.saveGraph("src/test/cyp2u1", "root");
        asr.saveSupportedAncestors("src/test/root");
        //System.out.println(msa.toString());
    }


    @Test
    @DisplayName("Small ASR")
    public void performSmallASR() throws IOException {
        //String alnfilepath = "src/test/resources/270717_2U1_var_region1.aln";
        //String nwkfilepath = "src/test/resources/270717_2U1.nwk";
        String alnfilepath = "src/test/resources/small.aln";
        String nwkfilepath = "src/test/resources/small.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true, false, null, 1);
        asr.saveGraph("src/test/cyp2u1", "root");
        PartialOrderGraph msa = asr.getPartialOrderGraph();
        PartialOrderGraph graph = asr.getGraph("root");
        System.out.println(msa.toString());
        System.out.println(graph.toString());
        System.out.println(graph.getConsensusSequence());
        System.out.println(graph.getConsensusGappySequence());
        for (String node : asr.getAncestralDict().keySet()) {
            PartialOrderGraph anc = asr.getGraph(node);
            System.out.println(node + ": " + anc.getConsensusGappySequence());
        }
    }
}
