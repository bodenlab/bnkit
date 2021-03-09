package asr;

import bn.ctmc.SubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.Utils;
import dat.phylo.Tree;
import dat.pog.POGTree;
import dat.pog.POGraph;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;


class FailingReconstructions {

    static SubstModel jtt;
    static GRASP.Inference joint;

    @BeforeAll
    @Disabled
    static void setup() {
        jtt = SubstModel.createModel("JTT");
        joint = GRASP.Inference.JOINT;
    }


    // This test is currently failing
    @Test
    @DisplayName("Path inferred in root ancestor that leads to nodes / edges that are not in root ancestor")
    void inferred_Path_Lead_To_Non_Inferred_Path() throws IOException, ASRException {

        try {

            EnumSeq.Alignment aln = Utils.loadAlignment("src/test/resources/currently_failing/SICML_FailedPath_1_6.aln", Enumerable.aacid);
            Tree tree = Utils.loadTree("src/test/resources/currently_failing/SICML_FailedPath_1_6.nwk");


            POGTree pogTree = new POGTree(aln, tree);


            Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree);
            assertTrue(indelpred instanceof Prediction);
            indelpred.getJoint(jtt);

            Map<Object, POGraph> pogs = indelpred.getAncestors(joint);


            POGraph[] ancestors = new POGraph[pogs.size()];

            int i = 0;
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
                ancestors[i++] = entry.getValue();
            EnumSeq[] ancseqs = new EnumSeq[pogs.size()];

            int ii = 0;

            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                ancseqs[ii++] = indelpred.getSequence(entry.getKey(), joint, true);
            }

        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    // This test is currently failing
    @Test
    @DisplayName("A certain node in an ancestor has only two possible edges leaving from it, neither of which get inferred")
    void neither_Of_Two_Possible_Edges_Inferred() throws IOException, ASRException {

        try {

            EnumSeq.Alignment aln = Utils.loadAlignment("src/test/resources/currently_failing/700_3_2232_cutdown.fasta", Enumerable.aacid);
            Tree tree = Utils.loadTree("src/test/resources/currently_failing/700_3_2232.nwk");

            POGTree pogTree = new POGTree(aln, tree);


            Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree);
            assertTrue(indelpred instanceof Prediction);
            indelpred.getJoint(jtt);

            Map<Object, POGraph> pogs = indelpred.getAncestors(joint);


            POGraph[] ancestors = new POGraph[pogs.size()];

            int i = 0;
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet())
                ancestors[i++] = entry.getValue();
            EnumSeq[] ancseqs = new EnumSeq[pogs.size()];

            int ii = 0;

            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                ancseqs[ii++] = indelpred.getSequence(entry.getKey(), joint, true);
            }

        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }


    }
}


