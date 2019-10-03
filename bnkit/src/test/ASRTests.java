import static org.junit.jupiter.api.Assertions.assertEquals;

import api.PartialOrderGraph;
import dat.EnumSeq;
import dat.Enumerable;
import dat.POGraph;
import java.io.IOException;
import java.util.List;
import json.JSONObject;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import reconstruction.ASRPOG;
import reconstruction.Consensus;
import reconstruction.Inference;
import vis.POAGJson;

/**
 * Created by marnie on 27/4/17.
 */
public class ASRTests {

    @Test
    @DisplayName("Get ancestral partial order graph")
    public void getAncestralGraphTest() throws IOException, InterruptedException {
        System.out.println("======================= Get ancestral partial " +
                "order graph ====================\n");
        String alnfilepath = "bnkit/src/test/resources/tawfik.aln";
        String nwkfilepath = "bnkit/src/test/resources/tawfik.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true, false, null, 1);
        PartialOrderGraph root = asr.getGraph("root");

        POGraph pog = asr.getAncestor("N22_68");
        String expectedGappyConsensus = "GYVMTVLG-PLPINEMGVTLMHEHILLDASGKWVPPCCCSDRHLAEMPVKMENLGELSLNPLMSRDNCQLFDVDVAIDELTKYRALGGETVVDPTNIGIGRDPKALARIARLT-GLNIIMGTGLYLEPSHPEWVKIS-----------SVEQLTERLIYDLGGAEEKPEVLAGLIGEIGISS-RFTPDEEKSLRAAGRASAATGVPIEVHLPGWERLGHRVLDILE-QEGADLRHTVLCHMNPSFADKRYQRELAQRGAFLEYDMIGMSYYYADES------AQSPSDEENARAIRELIDDGYIQQILLSQDVFLKTM----------LTRYGGHGYGYILKHFVPRLRRHGVSGEQLETLMIGNPQRVFGG-----";

        String expectedSupportedSeq =
                "SQVQTVTG-PIDVEQLGKTLVHEHVFVLGE-----------EFRQNYQAEWD----------------EEERIADAVEKLTELKSLGIDTIVDPTVIGLGRYIPRIQRIAEQV-DLNIVVATGIYTYNEVPFQFHYSGPGL----LFDGPEPMVEMFVKDIEDGIAGTGVRAGFL-KCAIEEQGLTPGVERVMRAVAQAHVRTGAPITVHTHAHSESGLEAQRVLA-EEGADLTKVVIGHSG-DSTDLDYLCELADAGSYLGMDRF-----GLDV---------LLPFEERVDTVAELCRRGYADRMVLAHDASCFID---WFPPEARAAAVPNWNYRHISEDVLPALRERGVTEEQIQTMLVDNPRRYFGS-----";
        String gappySeq = root.getConsensusGappySequence();
        String supportedSeq = pog.getSupportedSequence(true);
        assertEquals(gappySeq, expectedGappyConsensus);
        System.out.println("     Testing the expected gappy consensus " +
                "sequence \n");
        System.out.println("Expected: " + expectedGappyConsensus);
        System.out.println("Received: " + gappySeq);
        System.out.println("\n****    PASSED\n");
        System.out.println("     Testing the expected supported consensus" +
                " " +
                "sequence \n");
        assertEquals(supportedSeq, expectedSupportedSeq);
        System.out.println("Expected: " + supportedSeq);
        System.out.println("Received: " + expectedSupportedSeq);
        System.out.println("\n****    PASSED\n");

        System.out.println("     Testing we are able to re-generate the consensus from the inference ROOT: \n");
        POAGJson json = new POAGJson(root);

        /* We want to check that we get the same result reading back in the JSON object */
        JSONObject jsonObject = json.toJSON();
        Consensus c = new Consensus(jsonObject);
        String regeneratedConsensus = c.getSupportedSequence(true);

        assertEquals(expectedGappyConsensus, regeneratedConsensus);


        assertEquals(expectedGappyConsensus, regeneratedConsensus);
        System.out.println("Expected: " + expectedGappyConsensus);
        System.out.println("Received: " + regeneratedConsensus);
        System.out.println("\n****    PASSED\n");

        PartialOrderGraph n22 = asr.getGraph("N22_68");
        json = new POAGJson(n22);

        System.out.println("     Testing we are able to re-generate the consensus from the inference N22_68: \n");

        /* We want to check that we get the same result reading back in the JSON object */
        jsonObject = json.toJSON();
        c = new Consensus(jsonObject);
        regeneratedConsensus = c.getSupportedSequence(true);
        assertEquals(supportedSeq, regeneratedConsensus);
        System.out.println("Expected: " + expectedSupportedSeq);
        System.out.println("Received: " + regeneratedConsensus);
        System.out.println("\n****    PASSED\n");
        assertEquals(expectedSupportedSeq, regeneratedConsensus);

    }


    @Test
    @DisplayName("huge ASR")
    public void performSmallASR() throws IOException, InterruptedException {
        System.out.println("======================= Perform Huge ASR " +
                "=========================\n");
        String alnfilepath = "bnkit/src/test/resources/tawfik.aln";
        String nwkfilepath = "bnkit/src/test/resources/tawfik.nwk";
        ASRPOG asr = new ASRPOG(alnfilepath, nwkfilepath, true, false, null, 1);
        String gappySeq = asr.getGraph("root").getConsensusGappySequence();
        String expectedGappyConsensus = "GYVMTVLG-PLPINEMGVTLMHEHILLDASGKWVPPCCCSDRHLAEMPVKMENLGELSLNPLMSRDNCQLFDVDVAIDELTKYRALGGETVVDPTNIGIGRDPKALARIARLT-GLNIIMGTGLYLEPSHPEWVKIS-----------SVEQLTERLIYDLGGAEEKPEVLAGLIGEIGISS-RFTPDEEKSLRAAGRASAATGVPIEVHLPGWERLGHRVLDILE-QEGADLRHTVLCHMNPSFADKRYQRELAQRGAFLEYDMIGMSYYYADES------AQSPSDEENARAIRELIDDGYIQQILLSQDVFLKTM----------LTRYGGHGYGYILKHFVPRLRRHGVSGEQLETLMIGNPQRVFGG-----";
        assertEquals(gappySeq, expectedGappyConsensus);
        System.out.println("     Testing the expected gappy consensus " +
                "sequence \n");
        System.out.println("Expected: " + expectedGappyConsensus);
        System.out.println("Received: " + gappySeq);
        System.out.println("\n****    PASSED\n");
    }

    @Test
    @DisplayName("Run reconstruction from input data")
    public void runFromData() throws IOException, InterruptedException {
        System.out.println("======================= Run reconstruction from " +
                "input data =========================\n");
        String alnfilepath = "bnkit/src/test/resources/tawfik.aln";
        String nwkfilepath = "bnkit/src/test/resources/tawfik.nwk";
        String treeNwk = "(RTXKPRP:0.00000014,RTXKlebvar:0.00623058,(RTX_Pseudo:0.10125108,(RTX_3K2g:0.52474419,((Symbiobact:0.37069573,(PHP_Escher:0.14236022,(PHP_Yersin:0.27406260,(PHP_Photor:0.13809403,PHP_Xenorh:0.42798709)59:0.07439548)79:0.11321042)100:0.66251453)98:0.28009990,((PLLDeinoco:0.42937975,(PLLGeoKaus:0.07205125,PLLGeobThe:0.04452138)99:0.28466264)100:0.89834731,(((1HZY_pte:0.00000001,PTEFlavob:0.00302678)97:0.07465645,(2R1N_opd:0.00323286,PTEAgrobac:0.00332231)81:0.02820201)100:1.19982396,((PLLSulAcid:0.13201170,(SisPox_a:0.04040092,ssopoxmo:0.05938749)100:0.15659953)100:0.61438202,((PLLRhodoco:0.39323398,((PLLAhIA:0.00324601,PLLQsdA:0.00000023)100:0.19348514,(PLLBreviba:0.17059149,PLLDermaco:0.24217329)68:0.09748923)100:0.15423775)88:0.12323455,(PLLStrepto:0.57408811,(PLLMycsubs:0.03654787,(PLLPPH:0.00000001,(PLLMycbovi:0.00321720,PLLMycobCD:0.00324499)22:0.00000022)100:0.05401624)99:0.14298798)94:0.09766462)99:0.50935379)82:0.20681095)94:0.37463577)95:0.33701264)100:0.83757149)92:0.27920519)100:0.21425280);";
        List<EnumSeq.Gappy<Enumerable>> seqs = EnumSeq.Gappy.loadClustal(alnfilepath, Enumerable.aacid_ext);

        ASRPOG asr = null;
        System.out.println("running joint... ");
        asr =  new ASRPOG(alnfilepath, nwkfilepath, true, false, null, 1);

        JSONObject infs = asr.exportInferencesToJSON();
        ASRPOG test  = new ASRPOG(null, 1, infs, seqs, treeNwk);


        System.out.println("\n\nTesting differences in the saved " +
                "inferences...\n");
        for (String a : test.getAncestralInferences().keySet()) {
            System.out.println(a + "--------");
            for (int i = 0; i < test.getAncestralInferences().get(a).size(); i++) {
                Inference it = test.getAncestralInferences().get(a).get(i);
                Inference im = asr.getAncestralInferences().get(a).get(i);
                if (!it.toString().equalsIgnoreCase(im.toString())) {
                    System.out.println(it.toString());
                    System.out.println(im.toString());
                }
            }
        }

        System.out.println("\n\nTesting differences in the graph " +
                "structures...\n");

        if (!test.getMSAGraph().toString().equalsIgnoreCase(asr.getMSAGraph().toString())) {
            System.out.println("MSA Graphs (template graphs) are not equal... ");
        }

        POGraph testPO = test.getAncestor("root");
        POGraph asrPO = asr.getAncestor("root");

        if (!testPO.toString().equalsIgnoreCase(asrPO.toString())) {

            System.out.println("POGraphs are not equal...");

            for (Integer nodeId : testPO.getNodeIDs()) {
                System.out.println("Node ID: " + nodeId + " -----------------");
                testPO.setCurrent(nodeId);
                asrPO.setCurrent(nodeId);
                for (Integer outE : testPO.getNextIDs())
                    if (!asrPO.getNextIDs().contains(outE))
                        System.out.println("    edge " + outE + " does not match");
                if (testPO.getCurrentConsensusFlag() != asrPO.getCurrentConsensusFlag())
                    System.out.println("    consensus does not match for node");
                if (!testPO.getCurrentLabel().equalsIgnoreCase(asrPO.getCurrentLabel()))
                    System.out.println("    label " + testPO.getCurrentLabel() + " | " + asrPO.getCurrentLabel() + " does not match for node");
            }


        }

        PartialOrderGraph testG = test.getGraph("root");
        PartialOrderGraph oldG = asr.getGraph("root");

        if (!testG.toString().equalsIgnoreCase(oldG.toString())) {

            System.out.println("PartialOrderGraphs are not equal...");

            for (Integer nodeId : testG.getNodeIDs()) {
                System.out.println("Node ID: " + nodeId + " -----------------");
                for (Integer outE : testG.getOutEdgeWeights(nodeId).keySet())
                    if (!oldG.getOutEdgeWeights(nodeId).containsKey(outE))
                        System.out.println("    edge " + outE + " does not match");
                if (testG.getConsensusMembership(nodeId) != oldG.getConsensusMembership(nodeId))
                    System.out.println("    consensus does not match for node");
                if (!testG.getLabel(nodeId).equalsIgnoreCase(oldG.getLabel(nodeId)))
                    System.out.println("    label " + testG.getLabel(nodeId) + " | " + oldG.getLabel(nodeId) + " does not match for node");
            }


        }


        POAGJson json = new POAGJson(testG);
        POAGJson jsonOld = new POAGJson(oldG);

        // make sure node IDs line up with the correct positioning in the MSA graph
        JSONObject graph =  json.toJSON();
        JSONObject graphOld = jsonOld.toJSON();

    }
}
