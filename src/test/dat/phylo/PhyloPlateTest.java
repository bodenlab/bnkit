package dat.phylo;

import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.node.CPT;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;
import dat.file.FastaReader;
import json.JSONObject;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

class PhyloPlateTest {

    static IdxTree mblTree = null;
    static EnumSeq.Alignment mblAln = null;

    @BeforeAll
    static void setThingsUp() {
        Tree tree;
        try {
            tree = Tree.load("data/master_relabel.nwk", "newick");
            mblTree = (IdxTree)tree;
            //Newick.save(mblTree, "data/master_relabel.nwk", Newick.MODE_STRIPPED);
            FastaReader r = new FastaReader("data/master_t10.aln", Enumerable.aacid, '-');
            EnumSeq.Gappy[] eseqs = r.loadGappy();
            mblAln = new EnumSeq.Alignment(eseqs);
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    @Test
    void templateToJSON() {
        System.out.println(Enumerable.bool.toJSON());
        System.out.println(Enumerable.nacid.toJSON());
        PhyloPlate.Modes template = new PhyloPlate.Modes(new Enumerable[] {Enumerable.bool, Enumerable.nacid});
        PhyloPlate.Plate plate = new PhyloPlate.Plate("Test", template);
        plate.addNode(new CPT(new EnumVariable(Enumerable.nacidwn, "DNA1wN"), plate.getParents(new int[] {0})));
        plate.addNode(new CPT(new EnumVariable(Enumerable.nacid, "DNA2"), plate.getParents(new int[] {0,1})));
        JSONObject j1 = plate.toJSON();
        System.out.println(j1);
        PhyloPlate.Plate p1 = PhyloPlate.Plate.fromJSON(j1);
        System.out.println(p1.toJSON());
        assertTrue(plate.toJSON().toString().equals(PhyloPlate.Plate.fromJSON(plate.toJSON()).toJSON().toString()));
    }

    @Test
    void mblmotifdata() {

        /* Metal binding at (positions in alignment, starting at column 1)
        Alpha site: 184, 186, 306
        Beta site: 188, 189, 567
        Bridging residue: 332
         */
        int[] sites = new int[] {183, 185, 305,  187, 188, 566};
        // int NSTATES = 3;
        long SEED = 3;
        Enumerable state1 = new Enumerable(new Object[] {'A','B','C'});
        Enumerable state2 = new Enumerable(new Object[] {'a','b','c'});
        PhyloPlate.Modes template = new PhyloPlate.Modes(new Enumerable[] {state1, state2});
        PhyloPlate phybn = new PhyloPlate(mblTree, template);
        boolean isMaster = true;
        CPT[] masters = new CPT[sites.length];
        for (int idx : mblTree) {
            if (mblTree.isLeaf(idx)) {
                PhyloPlate.Plate plate1 = phybn.getPlate(idx);
                for (int s = 0; s < sites.length; s ++) {
                    CPT pos = new CPT(new EnumVariable(Enumerable.aacid, "Pos" + (sites[s] + 1)), plate1.getParents(new int[]{s / (sites.length / 2)}));
                    pos.randomize(SEED + s);
                    plate1.addNode(pos);
                    if (isMaster)
                        masters[s] = pos;
                    else
                        pos.tieTo(masters[s]);
                }
                if (isMaster) {
                    phybn.setMaster(plate1);
                    isMaster = false;
                }
                phybn.getBN().add(plate1.bnodes);
            }
        }

        String[] seqnames = mblAln.getNames();
        Object[][] data = new Object[seqnames.length][sites.length];
        Map<String, Object[]> seqmap = new HashMap<>();
        for (int i = 0; i < seqnames.length; i ++) {
            for (int s = 0; s < sites.length; s ++)
                data[i][s] = mblAln.getEnumSeq(i).get(sites[s]);
            seqmap.put(seqnames[i], data[i]);
        }

        phybn.trainEM(seqnames, data, SEED);
        System.out.println(phybn.getMasterJSON());

        for (Object key : state1.getValues()) {
            System.out.println(">" + key);
            for (int i = 0; i < 100; i ++) {
                for (int b = 0; b < phybn.getMaster().bnodes.length / 2; b ++) {
                    BNode bnode = phybn.getMaster().bnodes[b];
                    Object ch = bnode.getDistrib(new Object[]{key}).sample();
                    System.out.print(ch.toString());
                }
                System.out.println();
            }
        }
        for (Object key : state2.getValues()) {
            System.out.println(">" + key);
            for (int i = 0; i < 100; i ++) {
                for (int b = phybn.getMaster().bnodes.length / 2; b < phybn.getMaster().bnodes.length; b ++) {
                    BNode bnode = phybn.getMaster().bnodes[b];
                    Object ch = bnode.getDistrib(new Object[]{key}).sample();
                    System.out.print(ch.toString());
                }
                System.out.println();
            }
        }

        for (int idx : mblTree) {
            if (mblTree.isLeaf(idx)) {
                BNode[] bnodes = phybn.getBNodes(idx);
                PhyloPlate.Plate plate = phybn.getPlate(idx);
                String seqname = plate.name;
                Object[] motif = seqmap.get(seqname);
                if (motif == null)
                    System.out.println("Did not find motif for " + seqname);
                else {
                    if (plate != null) { // not hidden, so can be instantiated and inferred
                        int i = 0;
                        for (BNode bnode : plate.bnodes) {
                            bnode.setInstance(motif[i]);
                            i++;
                        }
                    }
                }
            }
        }
        VarElim ve = new VarElim();
        ve.instantiate(phybn.getBN());
        String[] motifs = new String[mblTree.getSize()];
        Object[][] infval = new Object[mblTree.getSize()][2];
        for (int idx : mblTree) {
//            if (mblTree.isParent(idx)) {
            BNode[] querynodes = phybn.getBNodes(idx);
            for (int j = 0; j < querynodes.length; j ++) {
                Query q = ve.makeQuery(querynodes[j].getVariable());
                CGTable r1 = (CGTable) ve.infer(q);
                Object value = r1.query(querynodes[j].getVariable());
                infval[idx][j] = value.toString().replace(' ', '_');
            }
        }

        String[] infdists = new String[mblTree.getSize()];
        for (int idx : mblTree) {
            if (mblTree.isLeaf(idx)) {
                // System.out.println(mblTree.getLabel(idx) + "\t" + value);
                Object[] motif = seqmap.get(mblTree.getLabel(idx));
                StringBuilder sb = new StringBuilder();
                for (Object ch : motif)
                    sb.append(ch == null ? "-" : ch.toString());
                motifs[idx] = sb.toString();
                infdists[idx] = motifs[idx];
            } else
                infdists[idx] = infval[idx][0].toString() + infval[idx][1].toString();
        }
        TreeInstance t_marg = new TreeInstance(mblTree, infdists);

        Object[][] infval2 = new Object[mblTree.getSize()][2];
        Query q = ve.makeMPE();
        CGTable r1 = (CGTable) ve.infer(q);
        Variable.Assignment[] assign = r1.getMPE();
        for (Variable.Assignment assign1 : assign) {
            String name = assign1.var.toString();
            int whichstate = 0;
            int shift = name.indexOf("_1.");
            if (shift >= 0)
                name = name.substring(0, shift);
            else {
                shift = name.indexOf("_2.");
                if (shift >= 0) {
                    name = name.substring(0, shift);
                    whichstate = 1;
                }
            }
            int idx = mblTree.getIndex(name);
            if (idx >= 0) {
                infval2[idx][whichstate] = assign1.val;
            } else
                System.out.println("No value for " +name + " (" + assign1.var.toString() + ")");
        }
        String[] infstates = new String[mblTree.getSize()];
        for (int i = 0; i < infval2.length; i ++)
            if (infval2[i][0] != null && infval2[i][1] != null)
                infstates[i] = infval2[i][0].toString() + infval2[i][1].toString();
        TreeInstance t_joint = new TreeInstance(mblTree, infstates);
        try {
            t_marg.save("data/master_infm6x2_s3.nwk");
            t_joint.save("data/master_infj6x2_s3.nwk");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}