package bn.ctmc;

import bn.BNode;
import bn.Distrib;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.matrix.JC;
import bn.ctmc.matrix.JTT;
import bn.factor.CachedFactor;
import bn.factor.FactorCache;
import bn.node.CPT;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;
import dat.colourschemes.Clustal;
import dat.file.AlnWriter;
import dat.file.FastaWriter;
import dat.file.Newick;
import dat.phylo.*;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.TimeUnit;

import static dat.Enumerable.aacid;
import static org.junit.jupiter.api.Assertions.*;
import util.MilliTimer;

/**
 * Many of the tests below are written to debug the implementations, evaluate efficiency and accuracy.
 */
class SubstNodeTest {

    SubstModel model = SubstModel.createModel("JTT");
    FactorCache cache = new FactorCache();
    boolean CACHE_FACTORS = true;

    /**
     * @param N number of nodes
     * @param NDIST number of time distances
     * @param TIME time from root to leaf
     */
    PhyloBN getBarebone(int N, int NDIST, double TIME) {
        SubstNode[] snarr = new SubstNode[N];
        double[] tarr = new double[N];
        EnumVariable[] vararr = new EnumVariable[N];
        long START_TIME = System.currentTimeMillis();
        Random rand = new Random(N);
        double delta_t = TIME / (double)N;
        double sum_t = 0;
        Map<Double, Integer> tmap = new HashMap<>();
        for (int i = 0; i < snarr.length; i ++) {
            int NDIV = NDIST / 2 + 1;
            double adjust = ((rand.nextInt(NDIV - 1) + 1)/(double)NDIV) * delta_t * (rand.nextBoolean() ? -1 : +1);
            double t = delta_t + adjust;
            tarr[i] = t;
            if (!tmap.containsKey(t))
                tmap.put(t, 0);
            tmap.put(t, tmap.get(t) + 1);
            vararr[i] = new EnumVariable(aacid, "N" + i);
            if (i == 0)
                snarr[i] = new SubstNode(vararr[i], model);
            else
                snarr[i] = new SubstNode(vararr[i], vararr[i - 1], model, t);
            sum_t += t;
        }
        PhyloBN pbn = PhyloBN.createBarebone(snarr);
        if (CACHE_FACTORS)
            System.out.println("Creating BN stretching " + sum_t + " time units, over " + N + " nodes, of which " + pbn.setCache(cache) + " are cache-enabled");
        System.out.println("Number of different time distances: " + tmap.size() + " \t");
        for (Map.Entry<Double, Integer> entry : tmap.entrySet()) {
            System.out.print(entry.getKey() + ": " + entry.getValue() + ", ");
        }
        System.out.println();
        long ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
        System.out.println(String.format("Creating %d nodes in %d min, %d s or %d ms", N, TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME)), TimeUnit.MILLISECONDS.toMillis(ELAPSED_TIME)));
        assertTrue(pbn.isValid());
        BNode root = pbn.getBNode(0);
        Map<Object, Double> selfmap = new HashMap<>();
        for (Object sym : model.getDomain().getValues())
            selfmap.put(sym, 1.0);
        // loop through all nodes, starting with root (0), moving towards the leaf
        for (int i = 0; i < pbn.getBN().getNodes().size(); i ++) {
            BNode node = pbn.getBNode(i);
            for (Object sym : model.getDomain().getValues()) {
                double p;
                if (node.isRoot())
                    p = node.get(sym); // P(sym)
                else
                    p = node.get(sym, sym); // P(sym | sym)
                if (i >= 2) {
                    if (tarr[i] <= tarr[i - 1]) // decreasing distance (or staying the same)
                        assertTrue(selfmap.get(sym) <= p); // make sure the prob of substitution to itself is decreasing with greater distance
                    else
                        assertTrue(selfmap.get(sym) > p); // make sure the prob of substitution to itself is increasing with greater distance
                }
                double sum = 0;
                for (Object other : model.getDomain().getValues()) {
                    if (other != sym) {
                        double pp;
                        if (node.isRoot())
                            pp = node.get(other); // P(other)
                        else
                            pp = node.get(other, sym); // P(other | sym)
                        assertTrue(pp <= 1.0 && pp >= 0.0);
                        sum += pp;
                    }
                }
                assertTrue((sum + p) > 0.99 && (sum + p) < 1.01 ); // the probs should sum roughly to 1
                selfmap.put(sym, p);
            }
        }
        return pbn;
    }

    Tree getTree(int NLEAVES, long seed) {
        String[] labels = new String[NLEAVES];
        for (int i = 0; i < NLEAVES; i ++)
            labels[i] = "S" + String.format("%03d", i);
        return Tree.Random(labels, seed, 1.1, 5.0, 2, 2);
    }

    @Test
    void getProb() {
        int N = 200, NDIST = 10;
        double TIME = 5.0;
        double[] tarr = new double[N];
        PhyloBN pbn = getBarebone(N, NDIST, TIME);
        BNode leaf = pbn.getBNode(N - 1);
        leaf.print();
    }

    @Test
    void getProb_Marginal1() {
        int N = 200, NDIST = 10;
        double TIME = 5.0;
        MilliTimer timer = new MilliTimer();
        PhyloBN pbn = getBarebone(N, NDIST, TIME);
        BNode root = pbn.getBNode(0);
        BNode leaf = pbn.getBNode(N - 1);
        //for (Object sym : new Object[] {'A'}) {
        int cnt = 0;
        for (Object sym : /*model.getDomain().getValues()*/ new Object[] {'A','C','D','R','P'}) {
            timer.start("Reset");
            pbn.getBN().resetNodes();
            timer.stopStart("Reset", "setInst");
            leaf.setInstance(sym);
            timer.stopStart("setInst","VE inst");
            VarElim ve = new VarElim();
            ve.instantiate(pbn.getBN());
            timer.stopStart("VE inst","makeQuery");
            Query q = ve.makeQuery(root.getVariable());
            timer.stopStart("makeQuery","Infer");
            CGTable qr = (CGTable)ve.infer(q);
            timer.stopStart("Infer", "Distrib");
            Distrib d = qr.query(root.getVariable());
            timer.stop("Distrib");
            System.out.println(root + " is " + d);
            System.out.println("Cache size: " + cache.size());
            cnt += 1;
        }
        timer.report();
        cache.reportCache();
    }
    @Test
    void getProb_Joint1() {
        int N = 200, NDIST = 10;
        double TIME = 5.0;
        MilliTimer timer = new MilliTimer();
        PhyloBN pbn = getBarebone(N, NDIST, TIME);
        BNode root = pbn.getBNode(0);
        BNode leaf = pbn.getBNode(N - 1);
        //for (Object sym : new Object[] {'A'}) {
        int cnt = 0;
        for (Object sym : /*model.getDomain().getValues()*/ new Object[] {'W','C','D','R','P'}) {
            timer.start("Reset");
            pbn.getBN().resetNodes();
            timer.stopStart("Reset", "setInst");
            leaf.setInstance(sym);
            timer.stopStart("setInst","VE inst");
            VarElim ve = new VarElim();
            ve.instantiate(pbn.getBN());
            timer.stopStart("VE inst","makeMPE");
            Query q = ve.makeMPE();
            timer.stopStart("makeMPE","Infer");
            CGTable qr = (CGTable)ve.infer(q);
            ve.timer.report();
            timer.stopStart("Infer", "Assign");
            Variable.Assignment[] as = qr.getMPE();
            timer.stop("Assign");
            for (Variable.Assignment a : as) {
                if (a.var.getName().equals("N0"))
                    System.out.println("\t" + a);
            }
            System.out.println("Cache size: " + cache.size());
            cnt += 1;
        }
        timer.report();
        cache.reportCache();
    }

    @Test
    void getProb_Marginal2() {
        int NLEAVES = 20;
        long SEED = 1L;
        MilliTimer timer = new MilliTimer();
        Random rand = new Random(SEED);
        Object[] vals = model.getDomain().getValues();
        IdxTree t = getTree(NLEAVES, SEED);
        try {
            Newick.save(t, "/Users/mikael/simhome/ASR/infer2022/marg20m.nwk", Newick.MODE_ANCESTOR);
        } catch (IOException e) {
            System.err.println(e);
            System.exit(1);
        }
        PhyloBN pbn2 = PhyloBN.create(t, model);
        System.out.println("Cached nodes: " + pbn2.setCache(cache) + " of " + pbn2.getBN().getNodes().size());
        for (int pos = 0; pos < 100; pos ++) {
            timer.start("Reset");
            pbn2.getBN().resetNodes();
            timer.stopStart("Reset", "setInst");
            int[] chidx = t.getLeaves();
            for (int i : chidx) {
                String label = t.getLabel(i).toString();
                BNode leaf = pbn2.getBN().getNode(label);
                Object sym = vals[rand.nextInt(vals.length)];
                leaf.setInstance(sym);
            }
            timer.stopStart("setInst","VE inst");
            VarElim ve = new VarElim();
            ve.instantiate(pbn2.getBN());
            timer.stopStart("VE inst","makeQuery");
            BNode root = pbn2.getBNode(0);
            Query q = ve.makeQuery(root.getVariable()); // root
            timer.stopStart("makeQuery","Infer");
            CGTable qr = (CGTable)ve.infer(q);
            timer.stopStart("Infer", "Distrib");
            Distrib d = qr.query(root.getVariable());
            timer.stop("Distrib");
            System.out.println(root + " is " + d);
            System.out.println("Cache size: " + cache.size());
        }
        timer.report();
        cache.reportCache();
    }

    /**
     * Check if joint reconstructions are in agreement with marginal reconstructions.
     * Theoretically, they are not the same inference but should be in broad agreement.
     * Each ancestor state inferred by joint is compared (by rank) to the distribution given by marginal reconstruction at the same node.
     * If agreed, the rank is 1, otherwise the joint state ranks 2nd, 3rd, etc by falling probability.
     * The code saves trees which differ to a noticeable extent (average rank > 1.5), and
     * tests the assertion that the average rank is lower than 2 on average.
     * Note that the input columns are random so they display nothing that looks like conservation.
     * This in turn means reconstructions will tend to be more variable than those on natural data.
     */
    @Test
    void compare_Joint_Marg() {
        long SEED = 1L;
        int MIN_SIZE = 2, MAX_SIZE = 100; // number of leaves in trees inspected in this test
        Random rand = new Random(SEED);
        Object[] vals = model.getDomain().getValues();
        double sumavrank = 0;
        for (int NLEAVES = MIN_SIZE; NLEAVES < MAX_SIZE; NLEAVES ++) { // try different tree sizes (number of leaves)
            int cnt_all = 0, cnt_same = 0, rank = 0;
            IdxTree t = getTree(NLEAVES, SEED + NLEAVES);
            PhyloBN pbnm = PhyloBN.create(t, model);
            pbnm.setCache(cache);
            PhyloBN pbnj = PhyloBN.create(t, model);
            pbnj.setCache(cache);
            Object[] values = new Object[t.getSize()];
            Object[][] aln = new Object[t.getNLeaves()][1];
            String[] names = new String[t.getNLeaves()];
            int leafcnt = 0;
            for (int i : t.getLeaves()) { // sequences
                Object sym = vals[rand.nextInt(vals.length)];
                values[i] = sym;
                EnumSeq.Gappy eseq = new EnumSeq.Gappy(model.getDomain());
                BNode leafm = pbnm.getBNode(i);
                BNode leafj = pbnj.getBNode(i);
                leafm.setInstance(sym);
                leafj.setInstance(sym);
                aln[leafcnt][0] = sym;
                names[leafcnt++] = t.getLabel(i).toString();
            }
            VarElim vem = new VarElim();
            vem.instantiate(pbnm.getBN());
            VarElim vej = new VarElim();
            vej.instantiate(pbnj.getBN());
            Query qj = vej.makeMPE();
            CGTable qjt = (CGTable)vej.infer(qj);
            Map<Variable, Object> jrecon = Variable.Assignment.toMap(qjt.getMPE());
            for (int bpidx : t.getAncestors()) {
                BNode nodej = pbnj.getBNode(bpidx);
                t.getBranchPoint(bpidx).setLabel(nodej.getName());
                Object valj = jrecon.get(nodej.getVariable());
                values[bpidx] = valj;
                BNode nodem = pbnm.getBNode(bpidx);
                Query qm = vem.makeQuery(nodem.getVariable());
                CGTable qmt = (CGTable)vem.infer(qm);
                EnumDistrib dm = (EnumDistrib) qmt.query(nodem.getVariable());
                double pm = dm.get(valj);
                Object valm = vals[dm.getMaxIndex()];
                cnt_all += 1;
                for (int i = 0; i < dm.getDomain().size(); i ++) {
                    if (dm.get(i) >= pm)
                        rank += 1;
                }
                if (valj.equals(valm))
                    cnt_same += 1;
                else
                    values[bpidx] = valj.toString() + "_" + valm.toString();
                //System.out.println(nodem.getName() + "=" + valm.toString() + "\t==\t" + nodej.getName() + "=" + valj.toString());
            }
            double avrank = (rank / (double)cnt_all);
            sumavrank += avrank;
            System.out.println("With " + NLEAVES + ":\t" + cnt_same + " are same from " + cnt_all + " with average rank at " + avrank + (avrank > 1.5 ? " (saved)" : ""));
            if (avrank > 1.5) {
                Map<Object, Object> objmap = new HashMap<>();
                for (Map.Entry<Variable, Object> entry : jrecon.entrySet())
                    objmap.put(entry.getKey().getName(), entry.getValue());
                TreeInstance ti = new TreeInstance(t, values);
                try {
                    Newick.save(t, "/Users/mikael/simhome/ASR/infer2022/tst_" + NLEAVES + ".nwk", Newick.MODE_ANCESTOR);
                    FastaWriter writer = new FastaWriter("/Users/mikael/simhome/ASR/infer2022/tst_" + NLEAVES + ".fa");
                    writer.save(names, aln);
                    writer.close();
                    Newick.save(ti, "/Users/mikael/simhome/ASR/infer2022/tstjoint_" + NLEAVES + ".nwk");
                } catch (IOException e) {
                    System.err.println("Failed to save tree instance");
                }
            }
        }
        assertTrue(sumavrank / (MAX_SIZE - MIN_SIZE) < 2.0);
    }

    @Test
    void getProb_Joint2() {
        int NPOS = 100, NLEAVES = 20;
        long SEED = 1L;
        MilliTimer timer = new MilliTimer();
        Random rand = new Random(SEED);
        Object[] vals = model.getDomain().getValues();
        IdxTree t = getTree(NLEAVES, SEED);
        try {
            Newick.save(t, "/Users/mikael/simhome/ASR/infer2022/marg20j.nwk", Newick.MODE_ANCESTOR);
        } catch (IOException e) {
            System.err.println(e);
            System.exit(1);
        }
        PhyloBN pbn2 = PhyloBN.create(t, model);
        System.out.println("Cached nodes: " + pbn2.setCache(cache) + " of " + pbn2.getBN().getNodes().size());
        Object[][] pred = new Object[t.getAncestors().length+t.getNLeaves()][NPOS];
        for (int pos = 0; pos < NPOS; pos ++) { // positions in sequences
            timer.start("Reset");
            pbn2.getBN().resetNodes();
            timer.stopStart("Reset", "setInst");
            for (int i : t.getLeaves()) { // sequences
                BNode leaf = pbn2.getBN().getNode(t.getLabel(i).toString());
                Object sym = vals[rand.nextInt(vals.length)];
                leaf.setInstance(sym);
                pred[i][pos] = sym;
            }
            timer.stopStart("setInst","VE inst");
            VarElim ve = new VarElim();
            ve.instantiate(pbn2.getBN());
            timer.stopStart("VE inst","makeMPE");
            Query q = ve.makeMPE();
            timer.stopStart("makeMPE","Infer");
            CGTable qr = (CGTable)ve.infer(q);
            timer.stopStart("Infer", "Assign");
            Variable.Assignment[] as = qr.getMPE();
            timer.stop("Assign");
            for (Variable.Assignment a : as) {
                int bpidx = t.getIndex(a.var.getName());
                pred[bpidx][pos] = a.val;
            }
//            System.out.println("Cache size: " + cache.size());
        }
        timer.report();
        cache.reportCache();
        for (int i = 0; i < pred.length; i ++) {
            String label = (t.isLeaf(i) ? "" : "N") + t.getLabel(i).toString();
            System.out.print(label + "  \t[" + i + "] \t");
            for (int pos = 0; pos < NPOS; pos ++)
                System.out.print(pred[i][pos]);
            System.out.println();
        }
    }

    @Test
    void fromJSON() {
        Enumerable dom = new Enumerable(new Object[] {1,2,3});
        EnumVariable child = new EnumVariable(dom, "TestChild");
        EnumVariable parent = new EnumVariable(dom, "TestParent");
        EnumVariable child2 = new EnumVariable(aacid, "TestChild2");
        EnumVariable parent2 = new EnumVariable(aacid, "TestParent2");
        SubstNode snode = new SubstNode(child, parent, new JC(1, dom.getValues()), 1);
        SubstNode snode2 = new SubstNode(child2, parent2, SubstModel.createModel("JTT"), 1);
        System.out.println(snode.getStateAsText());
        System.out.println(snode.toJSON());
        System.out.println(snode2.toJSON());
        SubstNode scopy = SubstNode.fromJSON(snode.toJSON(), child, parent);
        assertTrue(scopy.toJSON().toString().equals(snode.toJSON().toString()));
        SubstNode scopy2 = SubstNode.fromJSON(snode2.toJSON(), child2, parent2);
        assertTrue(scopy2.toJSON().toString().equals(snode2.toJSON().toString()));
    }


}