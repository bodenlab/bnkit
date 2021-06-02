package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import dat.EnumSeq;
import dat.Enumerable;
import dat.Interval1D;
import dat.IntervalST;
import dat.file.FastaWriter;
import dat.file.Utils;
import dat.phylo.Tree;
import dat.pog.IdxGraph;
import dat.pog.POGTree;
import dat.pog.POGraph;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;


class IndelTest {

    static SubstModel jtt;
    static GRASP.Inference joint;

    @BeforeAll
    static void setup() {
        GRASP.VERBOSE = true;
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
            Prediction indelpred = Prediction.PredictBySICP(pogTree);
//            Prediction indelpred = Prediction.PredictByIndelParsimony(pogTree, forceLinear);
//            Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree, forceLinear);

//            Prediction indelpred = Prediction.PredictByBidirEdgeParsimony(         pogTree);
//            Prediction indelpred = Prediction.PredictByBidirEdgeMaxLHood(         pogTree);


            printAncestors(indelpred, pogTree, forceLinear);

        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }

    }

    // This test performs a reconstruction of INDELs with SICP, into POGs.
    // It then compares the predictions with the gappy reconstruction of FastML using parsimony for INDELs
    @Test
    @DisplayName("Compare SICP POGs against FastML sequences")
    void SICP_v_FastML() throws IOException, ASRException {

        String OUTPUT = "/Users/mikael/simhome/ASR/SICPtests/100";
        try {

            Tree tree = Utils.loadTree("src/test/resources/indels/FastML_test_100.nwk");
            EnumSeq.Alignment input = Utils.loadAlignment("src/test/resources/indels/FastML_test_100.aln", Enumerable.aacid);
            EnumSeq.Alignment output = Utils.loadAlignment("src/test/resources/indels/FastMLp_ancestors_extants_100.aln", Enumerable.aacid);
            Map<String, Integer> name2idx = new HashMap<>();
            String[] names = output.getNames();
            for (int i = 0; i < names.length; i ++)
                name2idx.put(names[i], i);
            POGTree pogTree = new POGTree(input, tree);
            Prediction indelpred = Prediction.PredictBySICP(pogTree);
//            Prediction indelpred = Prediction.PredictByIndelParsimony(pogTree, false);
            Map<String, List<String>> logerr = new HashMap<>();
            int maxincorr = 0;
            int maxidx = -1;
            for (Integer idx : tree) {
                if (!tree.isLeaf(idx)) { // is ancestor
                    POGraph pog = indelpred.getAncestor(tree.getLabel(idx));
                    if (pog != null) {
                        int fastml_idx = name2idx.get("N" + tree.getLabel(idx));
                        EnumSeq.Gappy fastml_seq = output.getEnumSeq(fastml_idx);
                        int correct = 0;
                        int incorrect = 0;
                        int prev = -1;
                        Object prevx = null;
                        List<String> errs = new ArrayList<>();
                        for (int i = 0; i < fastml_seq.length(); i ++) {
                            Object x = fastml_seq.get(i);
                            if (x != null) {
                                if (pog.isEdge(prev, i))
                                    correct += 1;
                                else {
                                    incorrect += 1;
                                    if (prevx != null)
                                        errs.add(String.format("%s%d-%s%d", prevx, prev, x, i));
                                    else
                                        errs.add(String.format("%d-%c%s", prev, x, i));
                                }
                                prev = i;
                                prevx = x;
                            }
                        }
                        System.out.println(idx + "\tN" + tree.getLabel(idx) + "\t" + correct + "\t" + (correct + incorrect) + "\t" + incorrect + "\t" + tree.getDepth(idx));
                        if (incorrect > 0) {
                            logerr.put("N" + tree.getLabel(idx), errs);
                            if (incorrect > maxincorr) {
                                maxidx = idx;
                                maxincorr = incorrect;
                            }
                        }
                    }
                }
            }
            if (maxidx > -1)
                System.out.println("Branchpoint " + maxidx + " failed on " + maxincorr);
            indelpred.getJoint(jtt);
            Map<Object, POGraph> pogs = indelpred.getAncestors(joint);
            POGraph[] ancestors = new POGraph[pogs.size()];
            EnumSeq[] ancseqs = new EnumSeq[pogs.size()];
            int ii = 0;
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                ancestors[ii ] = entry.getValue();
                ancseqs[ii ++] = indelpred.getSequence(entry.getKey(), joint, true);
            }
            FastaWriter fw = new FastaWriter(new File(OUTPUT, "GRASP_ancestors.fasta"));
            fw.save(ancseqs);
            fw.close();
            IdxGraph.saveToDOT(OUTPUT, ancestors);
            assertEquals(0, logerr.size());
        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }

    }
    // This test performs a reconstruction of INDELs with SICML, into POGs.
    // It then compares the predictions with the gappy reconstruction of FastML using ML for INDELs
    @Test
    @DisplayName("Compare SICML POGs against FastML sequences")
    void SICML_v_FastML() throws IOException, ASRException {
        String OUTPUT = "/Users/mikael/simhome/ASR/SICMLtests/old50";
        try {
            Tree tree = Utils.loadTree("src/test/resources/indels/FastML_50.nwk");
            EnumSeq.Alignment input = Utils.loadAlignment("src/test/resources/indels/FastML_test_50.aln", Enumerable.aacid);
            EnumSeq.Alignment output = Utils.loadAlignment("src/test/resources/indels/FastMLml_ancestors_extants_50.aln", Enumerable.aacid);
            Map<String, Integer> name2idx = new HashMap<>();
            String[] names = output.getNames();
            for (int i = 0; i < names.length; i ++)
                name2idx.put(names[i], i);
            POGTree pogTree = new POGTree(input, tree);
            Object[] possible = {true, false};
            SubstModel substmodel = new JC(1, possible);
            Prediction indelpred = Prediction.PredictBySICML(pogTree, substmodel);
            //Prediction indelpred = Prediction.PredictByIndelMaxLhood(pogTree, false);
            Map<String, List<String>> logerr = new HashMap<>();
            for (Integer idx : tree) {
                if (!tree.isLeaf(idx)) { // is ancestor
                    POGraph pog = indelpred.getAncestor(tree.getLabel(idx));
                    if (pog != null) {
                        int fastml_idx = name2idx.get("N" + tree.getLabel(idx));
                        EnumSeq.Gappy fastml_seq = output.getEnumSeq(fastml_idx);
                        int correct = 0;
                        int incorrect = 0;
                        int prev = -1;
                        Object prevx = null;
                        List<String> errs = new ArrayList<>();
                        for (int i = 0; i < fastml_seq.length(); i ++) {
                            Object x = fastml_seq.get(i);
                            if (x != null) {
                                if (pog.isEdge(prev, i))
                                    correct += 1;
                                else {
                                    incorrect += 1;
                                    if (prevx != null)
                                        errs.add(String.format("%s%d-%s%d", prevx, prev, x, i));
                                    else
                                        errs.add(String.format("%d-%c%s", prev, x, i));
                                }
                                prev = i;
                                prevx = x;
                            }
                        }
                        System.out.println(idx + "\tN" + tree.getLabel(idx) + "\t" + correct + "\t" + (correct + incorrect));
                        if (incorrect > 0)
                            logerr.put("N" + tree.getLabel(idx), errs);
                    }
                }
            }
            indelpred.getJoint(jtt);
            Map<Object, POGraph> pogs = indelpred.getAncestors(joint);
            POGraph[] ancestors = new POGraph[pogs.size()];
            EnumSeq[] ancseqs = new EnumSeq[pogs.size()];
            int ii = 0;
            for (Map.Entry<Object, POGraph> entry : pogs.entrySet()) {
                ancestors[ii ] = entry.getValue();
                ancseqs[ii ++] = indelpred.getSequence(entry.getKey(), joint, true);
            }
            FastaWriter fw = new FastaWriter(new File(OUTPUT, "GRASP_ancestors.fasta"));
            fw.save(ancseqs);
            fw.close();
            IdxGraph.saveToDOT(OUTPUT, ancestors);
        } catch (IOException | ASRException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }

    }
}