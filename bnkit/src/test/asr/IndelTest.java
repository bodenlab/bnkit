package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.GLOOME1;
import dat.EnumSeq;
import dat.Enumerable;
import dat.Interval1D;
import dat.IntervalST;
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


class IndelTest {

    static SubstModel jtt;
    static GRASP.Inference joint;

    @BeforeAll
    static void setup() {
        jtt = SubstModel.createModel("JTT");
        joint = GRASP.Inference.JOINT;
    }

    /**
     * Helper method to print ancestors from a prediction
     * @param indelpred Prediction to print ancestors from
     * @param pogTree   POGTree representing the tree
     * @param printPOVals Boolean for whether to print out if an added partially ordered edge is present at an ancestor
     */
    static void printAncestors(Prediction indelpred, POGTree pogTree, Boolean printPOVals) {

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

        IntervalST<Integer> povals = pogTree.getPOVals();

        if (printPOVals) {

            for (POGraph ancgraph : ancestors) {


                for (Interval1D poval : povals) {
                    System.out.println(poval);
                    System.out.println("This edge is present in this ancestor - " + ancgraph.isPath(poval.min, poval.max)+ "\n");
                }
            }
        }

        for (EnumSeq ancseq : ancseqs) {
            Object[] str = ancseq.get();

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


    // This test recreates Figure 6.7 in Gabe's thesis
    // It doesn't fail but Position Specific methods infer a gap in the ultimate ancestor where there shouldn't be one
    @Test
    @DisplayName("Position specific methods infer a gap in ultimate ancestor")
    void PSP_FAIL() throws IOException, ASRException {

        try {

            Tree tree = Utils.loadTree("src/test/resources/indels/PSP_FAIL_6.nwk");
            EnumSeq.Alignment aln = Utils.loadAlignment("src/test/resources/indels/PSP_FAIL_6.aln", Enumerable.aacid);
            POGTree pogTree = new POGTree(aln, tree);

            boolean forceLinear = false;

//            Prediction indelpred = Prediction.PredictByParsimony(         pogTree);
            Prediction indelpred = Prediction.PredictbyMaxLhood(pogTree);

//            Prediction indelpred = Prediction.PredictByIndelParsimony(pogTree, forceLinear);
//            Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree, forceLinear);

//            Prediction indelpred = Prediction.PredictByBidirEdgeParsimony(pogTree);
//            Prediction indelpred = Prediction.PredictByBidirEdgeMaxLHood(         pogTree);

            printAncestors(indelpred, pogTree, forceLinear);

        }

            catch (IOException | ASRException e) {
                System.err.println(e.getMessage());
                System.exit(1);
            }
    }

    // This test recreates Figure 6.8 in Gabe's thesis
    // It currently fails for both SIC methods unless columns are forced to be linear (using the forceLinear paramater)
    @Test
    @DisplayName("Simple Indel Coding fails due to lack of path in ancestral graph")
    void SIC_FAIL() throws IOException, ASRException {

        try {

            Tree tree = Utils.loadTree("src/test/resources/indels/SIC_FAIL_6.nwk");
            EnumSeq.Alignment aln = Utils.loadAlignment("src/test/resources/indels/SIC_FAIL_6.aln", Enumerable.aacid);
            POGTree pogTree = new POGTree(aln, tree);

            Boolean forceLinear = false;

//            Prediction indelpred = Prediction.PredictByParsimony(         pogTree);
//            Prediction indelpred = Prediction.PredictbyMaxLhood(         pogTree);

            Prediction indelpred = Prediction.PredictByIndelParsimony(pogTree, forceLinear);
//            Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree, forceLinear);

//            Prediction indelpred = Prediction.PredictByBidirEdgeParsimony(         pogTree);
//            Prediction indelpred = Prediction.PredictByBidirEdgeMaxLHood(         pogTree);


            printAncestors(indelpred, pogTree, forceLinear);

        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }

    }
}