package asr;

import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JTT;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.FastaWriter;
import dat.phylo.Tree;
import dat.pog.POAGraph;
import dat.pog.POGTree;
import dat.pog.POGraph;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

class TrAVISTest {
    boolean DEBUG = false;
    String OUTPUT = "/Users/mikael/simhome/ASR/travis/";
    String[] INDELS = new String[] {"BEP", "BEML", "SICP", "SICML", "PSP", "PSML"}; // same as in asr.GRASP
    //int[] methods = new int[] {0, 1, 2, 3, 4, 5};                   // run these INDEL methods
    int[] methods = new int[] {0, 1, 2, 4};                   // run these INDEL methods
    //int[] nextants      = new int[] {5, 10, 50, 100};               // number of extants
    int[] nextants      = new int[] {5};               // number of extants
    double[] scaledist  = new double[] {0.1, 0.2, 0.3, 0.5};   // scale evolutionary distances so that tree stretches this

    Map<Integer, TrAVIS.TrackTree> trackMap = new HashMap<>();  // map with trackers for each configuration
    int nSEEDS = 10;
    double GAMMA_SHAPE = 1.1; // setting to 1.0 will introduce values very close to zero
    double GAMMA_SCALE = 0.2;
    //Object[] seq = new Object[] {'A','A','A','A','A', 'Q','Q','Q','Q','Q', 'P','P','P','P','P', 'W','W','W','W','W'};
    Object[] seq = new Object[] {'M', 'A', 'Q', 'P', 'W'};
    EnumSeq ancseq = new EnumSeq(Enumerable.aacid);
    SubstModel model = new JTT();

    int getKey(int i, int j, int seed) {
        return i + j * nextants.length + seed * (nextants.length * scaledist.length);
    }
    int getI(int key) {
        int mkey = key % (nextants.length * scaledist.length);
        return mkey % nextants.length; // the remainder will give us i, e.g. for i=2, j=1, seed=3, the key is 2 + 1x4 + 3x(4x4) = 2 + 4 + 48 = 54, and 54 % 16 = 6 then 6 % 4 = 2
    }
    int getJ(int key) {
        int mkey = key % (nextants.length * scaledist.length);
        return mkey / nextants.length; // the remainder will give us i, e.g. for i=3, j=1, seed=2, the key is 3 + 1x4 + 2x(4x4) = 3 + 4 + 32 = 39, and 39 % 16 = 7 then 7 / 4 = 1
    }

    @BeforeEach
    void setUp() {
        ancseq.set(seq);
        for (int i = 0; i < nextants.length; i ++) {
            int N = nextants[i];
            for (int j = 0; j < scaledist.length; j ++) {
                double D = scaledist[j];
                for (int SEED = 0; SEED < nSEEDS; SEED ++) {
                    //             TrackTree tracker = new TrackTree(tree, ancseq, MODEL, SEED, RATESGAMMA==null?-1:RATESGAMMA);
                    Tree tree = Tree.Random(N, SEED, GAMMA_SHAPE, 1.0 / GAMMA_SCALE, 2, 2);
                    tree.adjustDistances(D); // scale distances
                    TrAVIS.TrackTree t = new TrAVIS.TrackTree(tree, ancseq, model, SEED, 0.7); // gamma a=0.7 typical protein rate variation
                    int key = getKey(i, j, SEED);
                    trackMap.put(key, t);
                }
            }
        }
    }

    @Test
    void setTrackMap() {
        System.out.println("Have access to " + trackMap.size() + " trackers");
    }


    @Test
    void compareRecon() {
        int maxkey = nSEEDS * nextants.length * scaledist.length; // number of tests
        int[][] n0missed = new int[methods.length][maxkey];
        int[][] n0indels = new int[methods.length][maxkey];
        // for each RECONSTRUCTION do ...
        for (Map.Entry<Integer, TrAVIS.TrackTree> entry : trackMap.entrySet()) {
            TrAVIS.TrackTree t = entry.getValue();
            int key = entry.getKey();
            EnumSeq[] a = t.getAlignment();
            List<EnumSeq> aext = new ArrayList<>();       // extants as generated
            Map<Integer, EnumSeq> aanc = new HashMap<>(); // ancestors as generated
            Tree tree = t.tree;
            int n0flips_total = 0;  // number of missed indel flips at N0 in this reconstruction
            int n0insertions = 0;   // number of insertions that have been introduced under N0 (regardless of size)
            int n0deletions = 0;    // number of deletions that have been introduced under N0 (regardless of size)
            for (int i = 0; i < a.length; i++) { // inspect the alignment, to create an index-to-sequence mapping (simulated extants and ancestors)
                int idx = tree.getIndex(a[i].getName());
                if (idx == -1 && a[i].getName().startsWith("N")) { // possibly an ancestor, so internal name does NOT begin with "N"
                    try {
                        Integer id = Integer.parseInt(a[i].getName().substring(1)); // look for numeric ID after "N"
                        idx = tree.getIndex(id);
                    } catch (NumberFormatException e) {
                        System.err.println("Could not find index for label " + a[i].getName());
                    }
                }
                if (idx != -1) {
                    if (tree.isLeaf(idx))
                        aext.add(a[i]);
                    aanc.put(idx, a[i]);
                } else
                    System.err.println("Could not find index for label " + a[i].getName());
            }
            EnumSeq.Alignment aln = new EnumSeq.Alignment(aext);
            if (DEBUG)
                System.out.println("Tracked reconstruction key=" + key + "\t" + a.length + " sequences, of which " + aext.size() + " are extants\n");
            POGTree pogtree = new POGTree(aln, tree);
            // for each METHOD do ...
            for (int midx = 0; midx < methods.length; midx++) {
                int METHOD = methods[midx];
                Prediction indelpred = null;
                switch (METHOD) {
                    case 0: // "BEP", "BEML", "SICP", "SICML", "PSP", "PSML"
                        indelpred = Prediction.PredictByBidirEdgeParsimony(pogtree);
                        break;
                    case 1:
                        indelpred = Prediction.PredictByBidirEdgeMaxLhood(pogtree);
                        break;
                    case 2:
                        indelpred = Prediction.PredictBySICP(pogtree);
                        break;
                    case 3:
                        indelpred = Prediction.PredictBySICML(pogtree);
                        break;
                    case 4:
                        indelpred = Prediction.PredictByParsimony(pogtree);
                        break;
                    case 5:
                        indelpred = Prediction.PredictByMaxLhood(pogtree);
                        break;
                }
                indelpred.getJoint(model);
                Map<Object, POGraph> pogs = indelpred.getAncestors(GRASP.Inference.JOINT);
                EnumSeq[] ancseqs = new EnumSeq[tree.getSize()];
                for (Map.Entry<Object, POGraph> e : pogs.entrySet()) {
                    Object name = e.getKey();
                    int idx = tree.getIndex(name);
                    if (idx == -1 && name.toString().startsWith("N")) { // possibly an ancestor, and internal name does NOT begin with "N"
                        try {
                            Integer id = Integer.parseInt(name.toString().substring(1)); // look for numeric ID after "N"
                            idx = tree.getIndex(id);
                        } catch (NumberFormatException e2) {
                            System.err.println("Could not find index for label " + name.toString());
                        }
                    }
                    Object[] oseq = indelpred.getSequence(name, GRASP.Inference.JOINT, true);
                    ancseqs[idx] = new EnumSeq(aln.getDomain());
                    ancseqs[idx].setName(name.toString());
                    ancseqs[idx].set(oseq);
                }
                for (int idx = 0; idx < ancseqs.length; idx++) {
                    EnumSeq predict = ancseqs[idx];
                    if (predict == null)
                        continue;
                    Object[] predict_nogap = predict.getStripped();
                    if (!tree.isLeaf(idx)) {
                        EnumSeq actual = aanc.get(idx);
                        Object[] actual_nogap = actual.getStripped();
                        assertEquals(actual.getName(), predict.getName());
                        if (DEBUG) {
                            StringBuilder sb = new StringBuilder();
                            for (int chidx : tree.getChildren(idx))
                                sb.append(tree.getBranchPoint(chidx).getLabel() + ";");
                            System.out.println(actual.getName() + "\t <-- " + sb);
                            System.out.println("Actual: \t" + actual);
                            System.out.println("Predict:\t" + predict);
                        }
                        int[] indel = new int[actual.length()]; // state of the sequence traversal::
                        //  0 means A and P are Match, or both are Gap
                        //  1 means P is Gap, A is character
                        // -1 means A is Gap, P is character
                        int flips = 0; // times flipped state
                        int prev = 0;
                        for (int i = 0; i < actual.length(); i++) {
                            if (actual.get(i) != predict.get(i)) // different
                                indel[i] = (actual.get(i) == null) ? -1 : (predict.get(i) == null ? +1 : 0);
                            if (prev != indel[i] && indel[i] != 0)
                                flips += 1;
                            prev = indel[i];
                        }
                        for (int chidx : tree.getChildren(idx)) {
                            int[] deletions = (int[]) t.getTreeWithDeletions().getInstance(chidx);
                            int[] insertions = (int[]) t.getTreeWithInsertions().getInstance(chidx);
                            int ndeletions = 0;
                            int ninsertions = 0;
                            int nmissed = 0;
                            for (int i = 0; i < deletions.length; i++) {
                                if (deletions[i] > 0) {
                                    n0deletions += 1;
                                    ndeletions += deletions[i];
                                    if (DEBUG)
                                        System.out.println("\tDELETE " + "@" + i + ":" + deletions[i]);
                                }
                            }
                            for (int i = 0; i < insertions.length; i++) {
                                if (insertions[i] > 0) {
                                    n0insertions += 1;
                                    ninsertions += insertions[i];
                                    if (DEBUG)
                                        System.out.println("\tINSERT " + "@" + i + ":" + insertions[i]);
                                }
                            }
                            EnumSeq chseq = null;
                            if (tree.isLeaf(chidx)) {
                                chseq = aanc.get(chidx);
                                if (DEBUG)
                                    System.out.println("Extant: \t" + chseq + "\t(" + tree.getBranchPoint(chidx).getLabel() + ")");
                            } else {
                                chseq = ancseqs[chidx];
                                if (DEBUG)
                                    System.out.println("Reconst:\t" + ancseqs[chidx] + "\t(" + tree.getBranchPoint(chidx).getLabel() + ")");
                            }
                        }
                        int diffseq = (predict_nogap.length - actual_nogap.length);
                        if (flips > 0 && DEBUG)
                            System.out.println("\tIndels missed:\t" + flips + "\tDifference:\t" + (diffseq));
                        if (idx == 0)
                            n0flips_total += flips;
                    }
                }
                if (DEBUG)
                    System.out.println("Missed indels at N0:\t" + n0flips_total + "\tSimulation: " + n0insertions + " insertions and " + n0deletions + " deletions were injected\n" + "============================================");
                n0missed[midx][key] = n0flips_total;
                n0indels[midx][key] = n0insertions + n0deletions;
            }
            // done with all methods

        }
        // done with all reconstructions
        // assemble a summary table...

        double[][] means = new double[methods.length][nextants.length * scaledist.length];
        System.out.print("#Tips\tScale\t");
        for (int midx = 0; midx < methods.length; midx ++)
            System.out.print(INDELS[methods[midx]] + " \t");
        System.out.println();
        for (int nidx = 0; nidx < nextants.length; nidx ++) {
            for (int sidx = 0; sidx < scaledist.length; sidx ++) {
                System.out.print(String.format("%4d\t%-4.2f\t", nextants[nidx], scaledist[sidx]));
                for (int midx = 0; midx < methods.length; midx ++) {
                    for (int SEED = 0; SEED < nSEEDS; SEED ++)
                        means[midx][nidx + sidx * nextants.length] += (n0missed[midx][getKey(nidx, sidx, SEED)] / (double)n0indels[midx][getKey(nidx, sidx, SEED)]) / nSEEDS;
                    System.out.print(String.format("%-5.3f\t", means[midx][nidx + sidx * nextants.length]));
                }
                System.out.println();
            }
        }
    }

    @Test
    void compareRecon2() {
        int maxkey = nSEEDS * nextants.length * scaledist.length; // number of tests
        int[][] missed = new int[methods.length][maxkey];
        int[][] notmissed = new int[methods.length][maxkey];
        double[][] nfactors = new double[methods.length][maxkey]; // normalisation factors
        // for each RECONSTRUCTION do ...
        for (Map.Entry<Integer, TrAVIS.TrackTree> entry : trackMap.entrySet()) {
            TrAVIS.TrackTree t = entry.getValue();
            int key = entry.getKey();
            EnumSeq[] a = t.getAlignment();
            List<EnumSeq> aext = new ArrayList<>();       // extants as generated
            Map<Integer, EnumSeq> aanc = new HashMap<>(); // ancestors as generated
            Tree tree = t.tree;
            for (int i = 0; i < a.length; i++) {
                int idx = tree.getIndex(a[i].getName());
                if (idx == -1 && a[i].getName().startsWith("N")) { // possibly an ancestor, so internal name does NOT begin with "N"
                    try {
                        Integer id = Integer.parseInt(a[i].getName().substring(1)); // look for numeric ID after "N"
                        idx = tree.getIndex(id);
                    } catch (NumberFormatException e) {
                        System.err.println("Could not find index for label " + a[i].getName());
                    }
                }
                if (idx != -1) {
                    if (tree.isLeaf(idx))
                        aext.add(a[i]);
                    aanc.put(idx, a[i]);
                } else
                    System.err.println("Could not find index for label " + a[i].getName());
            }
            EnumSeq.Alignment aln = new EnumSeq.Alignment(aext);
            if (DEBUG)
                System.out.println("Tracked reconstruction key=" + key + "\t" + a.length + " sequences, of which " + aext.size() + " are extants\n");
            POGTree pogtree = new POGTree(aln, tree);
            // for each METHOD do ...
            for (int midx = 0; midx < methods.length; midx++) {
                int METHOD = methods[midx];
                Prediction indelpred = null;
                switch (METHOD) {
                    case 0: // "BEP", "BEML", "SICP", "SICML", "PSP", "PSML"
                        indelpred = Prediction.PredictByBidirEdgeParsimony(pogtree);
                        break;
                    case 1:
                        indelpred = Prediction.PredictByBidirEdgeMaxLhood(pogtree);
                        break;
                    case 2:
                        indelpred = Prediction.PredictBySICP(pogtree);
                        break;
                    case 3:
                        indelpred = Prediction.PredictBySICML(pogtree);
                        break;
                    case 4:
                        indelpred = Prediction.PredictByParsimony(pogtree);
                        break;
                    case 5:
                        indelpred = Prediction.PredictByMaxLhood(pogtree);
                        break;
                }
                indelpred.getJoint(model);
                Map<Object, POGraph> pogs = indelpred.getAncestors(GRASP.Inference.JOINT);
                EnumSeq[] ancseqs = new EnumSeq[tree.getSize()];
                int correct_total = 0;
                int incorrect_total = 0;
                double nfactor = 0;
                for (Map.Entry<Object, POGraph> e : pogs.entrySet()) {
                    Object name = e.getKey();
                    int idx = tree.getIndex(name);
                    if (idx == -1 && name.toString().startsWith("N")) { // probably an ancestor, internal name does NOT begin with "N"
                        try {
                            Integer id = Integer.parseInt(name.toString().substring(1)); // look for numeric ID after "N"
                            idx = tree.getIndex(id);
                        } catch (NumberFormatException e2) {
                            System.err.println("Could not find index for label " + name.toString());
                        }
                    }
                    // do the thing...
                    POGraph pog = e.getValue();
                    EnumSeq seq = aanc.get(idx);
                    int correct = 0;
                    int incorrect = 0;
                    int prev = -1;
                    Object prevx = null;
                    List<String> errs = new ArrayList<>();
                    for (int i = 0; i < seq.length(); i++) {
                        Object x = seq.get(i);
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
                    correct_total += correct;
                    incorrect_total += incorrect;
                    nfactor += pog.getEdgeCount() - (double) pog.size() + 1; // number of "redundant" edges
                    //System.out.println(idx + "\tN" + tree.getLabel(idx) + "\t" + correct + "\t" + (correct + incorrect) + "\t" + incorrect + "\t" + tree.getDepth(idx));
                }
                if (DEBUG)
                    System.out.println("Missing edges\t" + incorrect_total);
                notmissed[midx][key] = correct_total;
                missed[midx][key] = incorrect_total;
                nfactors[midx][key] += nfactor / pogs.size();
            }

            // done with all methods

        }
        // done with all reconstructions
        // assemble a summary table...

        double[][] means = new double[methods.length][nextants.length * scaledist.length];
        System.out.print("#Tips\tScale\t");
        for (int midx = 0; midx < methods.length; midx ++)
            System.out.print(INDELS[methods[midx]] + " \t");
        System.out.println();
        for (int nidx = 0; nidx < nextants.length; nidx ++) {
            for (int sidx = 0; sidx < scaledist.length; sidx ++) {
                System.out.print(String.format("%4d\t%-4.2f\t", nextants[nidx], scaledist[sidx]));
                for (int midx = 0; midx < methods.length; midx ++) {
                    for (int SEED = 0; SEED < nSEEDS; SEED ++)
                        means[midx][nidx + sidx * nextants.length] +=
                                nfactors[midx][getKey(nidx, sidx, SEED)] *
                                        (missed[midx][getKey(nidx, sidx, SEED)] / (double)(missed[midx][getKey(nidx, sidx, SEED)] + notmissed[midx][getKey(nidx, sidx, SEED)])) / nSEEDS;
                    System.out.print(String.format("%-5.3f\t", means[midx][nidx + sidx * nextants.length]));
                }
                System.out.println();
            }
        }
    }

    @Test
    void compareRecon3() {
        int maxkey = nSEEDS * nextants.length * scaledist.length; // number of tests
        int[][] missed = new int[methods.length][maxkey];
        int[][] notmissed = new int[methods.length][maxkey];
        double[][] nfactors = new double[methods.length][maxkey]; // normalisation factors
        // for each RECONSTRUCTION do ...
        for (Map.Entry<Integer, TrAVIS.TrackTree> entry : trackMap.entrySet()) {
            TrAVIS.TrackTree t = entry.getValue();
            int key = entry.getKey();
            EnumSeq[] a = t.getAlignment();
            List<EnumSeq> aext = new ArrayList<>();       // extants as generated
            Map<Integer, EnumSeq> aanc = new HashMap<>(); // ancestors as generated
            Tree tree = t.tree;
            for (int i = 0; i < a.length; i++) {
                int idx = tree.getIndex(a[i].getName());
                if (idx == -1 && a[i].getName().startsWith("N")) { // possibly an ancestor, so internal name does NOT begin with "N"
                    try {
                        Integer id = Integer.parseInt(a[i].getName().substring(1)); // look for numeric ID after "N"
                        idx = tree.getIndex(id);
                    } catch (NumberFormatException e) {
                        System.err.println("Could not find index for label " + a[i].getName());
                    }
                }
                if (idx != -1) {
                    if (tree.isLeaf(idx))
                        aext.add(a[i]);
                    aanc.put(idx, a[i]);
                } else
                    System.err.println("Could not find index for label " + a[i].getName());
            }
            EnumSeq.Alignment aln = new EnumSeq.Alignment(aext);
            if (DEBUG)
                System.out.println("Tracked reconstruction key=" + key + "\t" + a.length + " sequences, of which " + aext.size() + " are extants\n");
            POGTree pogtree = new POGTree(aln, tree);
            // for each METHOD do ...
            Prediction[] indelpreds = new Prediction[methods.length];
            for (int midx = 0; midx < methods.length; midx++) {
                int METHOD = methods[midx];
                ;
                switch (METHOD) {
                    case 0: // "BEP", "BEML", "SICP", "SICML", "PSP", "PSML"
                        indelpreds[midx] = Prediction.PredictByBidirEdgeParsimony(pogtree);
                        break;
                    case 1:
                        indelpreds[midx] = Prediction.PredictByBidirEdgeMaxLhood(pogtree);
                        break;
                    case 2:
                        indelpreds[midx] = Prediction.PredictBySICP(pogtree);
                        break;
                    case 3:
                        indelpreds[midx] = Prediction.PredictBySICML(pogtree);
                        break;
                    case 4:
                        indelpreds[midx] = Prediction.PredictByParsimony(pogtree);
                        break;
                    case 5:
                        indelpreds[midx] = Prediction.PredictByMaxLhood(pogtree);
                        break;
                }
                Prediction indelpred = indelpreds[midx];
                indelpred.getJoint(model);
                Map<Object, POGraph> pogs = indelpred.getAncestors(GRASP.Inference.JOINT);
                EnumSeq[] ancseqs = new EnumSeq[tree.getSize()];
                int correct_total = 0;
                int incorrect_total = 0;
                double nfactor = 0;
                for (Map.Entry<Object, POGraph> e : pogs.entrySet()) {
                    Object name = e.getKey();
                    int idx = tree.getIndex(name);
                    if (idx == -1 && name.toString().startsWith("N")) { // probably an ancestor, internal name does NOT begin with "N"
                        try {
                            Integer id = Integer.parseInt(name.toString().substring(1)); // look for numeric ID after "N"
                            idx = tree.getIndex(id);
                        } catch (NumberFormatException e2) {
                            System.err.println("Could not find index for label " + name.toString());
                        }
                    }
                    // do the thing...
                    POGraph pog = e.getValue();
                    EnumSeq seq = aanc.get(idx);
                    int correct = 0;
                    int incorrect = 0;
                    int prev = -1;
                    Object prevx = null;
                    List<String> errs = new ArrayList<>();
                    for (int i = 0; i < seq.length(); i++) {
                        Object x = seq.get(i);
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
                    correct_total += correct;
                    incorrect_total += incorrect;
                    nfactor += pog.getEdgeCount() - (double) pog.size() + 1; // number of "redundant" edges
                    //System.out.println(idx + "\tN" + tree.getLabel(idx) + "\t" + correct + "\t" + (correct + incorrect) + "\t" + incorrect + "\t" + tree.getDepth(idx));
                }
                if (DEBUG)
                    System.out.println("Missing edges\t" + incorrect_total);
                notmissed[midx][key] = correct_total;
                missed[midx][key] = incorrect_total;
                nfactors[midx][key] += nfactor / pogs.size();
            }
            // done with all methods
            // now check if they resulted in different predictions
            int nCorrect = 0;
            for (int midx = 0; midx < methods.length; midx++) {
                nCorrect += missed[midx][key] == 0 ? 1 : 0;
            }
            if (nCorrect > 0 && nCorrect < methods.length) {
                /********************/
                try {
                    System.out.println("Extants=" + nextants[getI(key)] + "\tScale=" + scaledist[getJ(key)] + "\tKey=" + key + " (saved)");
                    FastaWriter fw = new FastaWriter(OUTPUT + key + ".fa");
                    EnumSeq[] aextarr = new EnumSeq[aext.size()];
                    aext.toArray(aextarr);
                    fw.save(aextarr);
                    fw.close();
                    fw = new FastaWriter(OUTPUT + key + "_sim.fa");
                    fw.save(t.getAlignment());
                    fw.close();
                    t.tree.save(OUTPUT + key + ".nwk", "nwk");
                } catch (IOException e) {
                    System.out.println("Extants=" + nextants[getI(key)] + "\tScale=" + scaledist[getJ(key)] + "\tKey=" + key + " (failed saving)");
                }
                System.out.println("\tMethod\tMissed\tNot missed");
                for (int midx = 0; midx < methods.length; midx++) {
                    System.out.println("\t" + String.format("%5s", INDELS[methods[midx]]) + "\t" + String.format("%5d", missed[midx][key]) + "\t" + String.format("%5d", notmissed[midx][key]));
                }
            }

        }
        // done with all reconstructions
        // assemble a summary table...

        double[][] means = new double[methods.length][nextants.length * scaledist.length];
        System.out.print("#Tips\tScale\t");
        for (int midx = 0; midx < methods.length; midx ++)
            System.out.print(INDELS[methods[midx]] + " \t");
        System.out.println();
        for (int nidx = 0; nidx < nextants.length; nidx ++) {
            for (int sidx = 0; sidx < scaledist.length; sidx ++) {
                System.out.print(String.format("%4d\t%-4.2f\t", nextants[nidx], scaledist[sidx]));
                for (int midx = 0; midx < methods.length; midx ++) {
                    for (int SEED = 0; SEED < nSEEDS; SEED ++)
                        means[midx][nidx + sidx * nextants.length] +=
                                nfactors[midx][getKey(nidx, sidx, SEED)] *
                                        (missed[midx][getKey(nidx, sidx, SEED)] / (double)(missed[midx][getKey(nidx, sidx, SEED)] + notmissed[midx][getKey(nidx, sidx, SEED)])) / nSEEDS;
                    System.out.print(String.format("%-5.3f\t", means[midx][nidx + sidx * nextants.length]));
                }
                System.out.println();
            }
        }
    }

    @Test
    void checkGeneratedAlignments() {
        double[][] lengths = new double[nextants.length][scaledist.length];
        double[][] nindels = new double[nextants.length][scaledist.length];
        for (Map.Entry<Integer, TrAVIS.TrackTree> entry : trackMap.entrySet()) {
            TrAVIS.TrackTree t = entry.getValue();
            int key = entry.getKey();
            EnumSeq[] a = t.getAlignment();
            lengths[getI(key)][getJ(key)] += (a[0].length()/(double)nSEEDS);
            POAGraph poag = t.getPOAG();
            nindels[getI(key)][getJ(key)] += ((poag.getEdgeCount() - a[0].length())/(double)nSEEDS);
        }
        System.out.println("Alignment width");
        for (int j = 0; j < scaledist.length; j ++)
            System.out.print(String.format("\t%-3.1f  ", scaledist[j]));
        System.out.println();
        for (int i = 0; i < nextants.length; i ++) {
            int N = nextants[i];
            System.out.print(String.format("%-3d\t",N));
            for (int j = 0; j < scaledist.length; j++)
                System.out.print(String.format("%-3.1f \t", lengths[i][j]));
            System.out.println();
        }
        System.out.println("\nNumber of indels");
        for (int j = 0; j < scaledist.length; j ++)
            System.out.print(String.format("\t%-3.1f  ", scaledist[j]));
        System.out.println();
        for (int i = 0; i < nextants.length; i ++) {
            int N = nextants[i];
            System.out.print(String.format("%-3d\t",N));
            for (int j = 0; j < scaledist.length; j++)
                System.out.print(String.format("%-3.1f \t", nindels[i][j]));
            System.out.println();
        }
    }

    @Test
    void compareInference() {

    }
}