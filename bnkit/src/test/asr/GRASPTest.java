package asr;

import bn.ctmc.SubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.IdxTree;
import dat.phylo.Tree;
import dat.pog.POGTree;
import dat.pog.POGraph;
import org.junit.jupiter.api.*;

import java.io.IOException;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertTrue;

class GRASPTest {

    static POGTree pogTree;
    static SubstModel jtt;
    static GRASP.Inference joint;



    @BeforeAll
    static void setup() {
        try {

            EnumSeq.Alignment aln = Utils.loadAlignment("src/test/resources/indelible_generated/700_3_2232.fasta", Enumerable.aacid);
            Tree tree = Utils.loadTree("src/test/resources/indelible_generated/700_3_2232.nwk");


            pogTree = new POGTree(aln, tree);

            jtt = SubstModel.createModel("JTT");
            joint = GRASP.Inference.JOINT;


        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    @Test
    void test_Indelible_Data_Returns_Prediction() {

        Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree);
        assertTrue(indelpred instanceof Prediction);
        indelpred.getJoint(jtt);
        Map<Object, POGraph> pogs = indelpred.getAncestors(joint);


        POGraph[] ancestors = new POGraph[pogs.size()];

        int i = 0;
        for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
            ancestors[i ++] = entry.getValue();
        EnumSeq[] ancseqs = new EnumSeq[pogs.size()];

        int ii = 0;

        for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
            ancseqs[ii++] = indelpred.getSequence(entry.getKey(), joint, true);
        }

        for (int s = 0; s < ancseqs.length; s++) {
            Object[] str = ancseqs[s].get();

            String anc_str = "";

            for (i = 0; i < str.length; i++) {
                try {
                    anc_str += str[i].toString();
                } catch (NullPointerException npe) {
                    anc_str += "-";
                }
            }
            System.out.println(anc_str);
        }
    }



    // Basic tests to check that each indel method is working

    @Test
    @Disabled
    @DisplayName("Check that BEP works on a basic example")
    void testBEP()  {

        Prediction prediction = Prediction.PredictByBidirEdgeParsimony(pogTree);
        assertTrue(prediction instanceof Prediction);
        prediction.getJoint(jtt);
        Map<Object, POGraph> pogs = prediction.getAncestors(joint);
        POGraph[] ancestors = new POGraph[pogs.size()];

    }

    @Test
    @Disabled
    @DisplayName("Check that BEML works on a basic example")
    void testBEML()  {

        Prediction prediction = Prediction.PredictByBidirEdgeMaxLHood(pogTree);
        assertTrue(prediction instanceof Prediction);
        prediction.getJoint(jtt);

    }

    @Test
    @DisplayName("Check that SICP works on a basic example")
    void testSICP()  {

        Prediction prediction = Prediction.PredictByIndelParsimony(pogTree);
        assertTrue(prediction instanceof Prediction);
        prediction.getJoint(jtt);

    }

    @Test
    @DisplayName("Check that SICML works on a basic example")
    void testSICML()  {

        Prediction prediction = Prediction.PredictByIndelMaxLhood(pogTree);
        assertTrue(prediction instanceof Prediction);
        prediction.getJoint(jtt);

    }

    @Test
    @DisplayName("Check that PSP works on a basic example")
    void testPSP()  {


        Prediction prediction = Prediction.PredictByParsimony(pogTree);
        assertTrue(prediction instanceof Prediction);
        prediction.getJoint(jtt);
        Map<Object, POGraph> pogs = prediction.getAncestors(joint);
        POGraph[] ancestors = new POGraph[pogs.size()];
        int ii = 0;
        for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
            ancestors[ii ++] = entry.getValue();
        ii = 0;

        EnumSeq[] ancseqs = new EnumSeq[pogs.size()];
        for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
            ancseqs[ii ++] = prediction.getSequence(entry.getKey(), joint, true);




    }


}
