package bn.ctmc;

import bn.BNode;
import bn.Distrib;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import dat.EnumVariable;
import dat.Enumerable;
import dat.phylo.PhyloBN;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.*;
import java.util.concurrent.TimeUnit;

import static org.junit.jupiter.api.Assertions.*;
import util.MilliTimer;

class SubstNodeTest {

    int N = 100;       // number of nodes
    int NDIST = 25;     // number of time distances
    double TIME = 5.0;  // time from root to leaf

    SubstNode[] snarr = new SubstNode[N];
    EnumVariable[] vararr = new EnumVariable[N];
    PhyloBN pbn;
    SubstModel model = SubstModel.createModel("JTT");

    @BeforeEach
    void setUp() {
        long START_TIME = System.currentTimeMillis();
        Random rand = new Random(N);
        double delta_t = TIME / (double)N;
        double sum_t = 0;
        Map<Double, Integer> tmap = new HashMap<>();
        for (int i = 0; i < snarr.length; i ++) {
            int NDIV = NDIST / 2 + 1;
            double adjust = ((rand.nextInt(NDIV - 1) + 1)/(double)NDIV) * delta_t * (rand.nextBoolean() ? -1 : +1);
            double t = delta_t + adjust;
            if (!tmap.containsKey(t))
                tmap.put(t, 0);
            tmap.put(t, tmap.get(t) + 1);
            vararr[i] = new EnumVariable(Enumerable.aacid, "N" + i);
            if (i == 0)
                snarr[i] = new SubstNode(vararr[i], model);
            else
                snarr[i] = new SubstNode(vararr[i], vararr[i - 1], model, t);
            sum_t += t;
        }
        pbn = PhyloBN.createBarebone(snarr);
        System.out.println("Creating BN stretching " + sum_t + " time units, over " + N + " nodes");
        System.out.println("Number of different time distances: " + tmap.size() + " \t");
        for (Map.Entry<Double, Integer> entry : tmap.entrySet()) {
            System.out.print(entry.getKey() + ": " + entry.getValue() + ", ");
        }
        System.out.println();
        long ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
        System.out.println(String.format("Creating %d nodes in %d min, %d s or %d ms", N, TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME)), TimeUnit.MILLISECONDS.toMillis(ELAPSED_TIME)));
        assertTrue(pbn.isValid());
    }

    @Test
    void getProb() {
        BNode root = pbn.getBNode(0);
        root.print();
        Map<Object, Double> selfmap = new HashMap<>();
        ;
        for (Object sym : model.getDomain().getValues()) {
            selfmap.put(sym, 1.0);
        }
        for (int i = 1; i < snarr.length; i ++) {
            BNode node = pbn.getBNode(i);
            for (Object sym : model.getDomain().getValues()) {
                double p = node.get(sym, sym); // P(sym | sym)
                assertTrue(selfmap.get(sym) > p); // make sure the prob of substitution to itself is decreasing with greater distance
                double sum = 0;
                for (Object other : model.getDomain().getValues()) {
                    if (other != sym) {
                        double pp = node.get(other, sym); // P(other | sym)
                        assertTrue(pp <= 1.0 && pp >= 0.0);
                        sum += pp;
                    }
                }
                assertTrue((sum + p) > 0.99 && (sum + p) < 1.01 ); // the probs should sum roughly to 1
                selfmap.put(sym, p);
            }
        }
        BNode leaf = pbn.getBNode(N - 1);
        leaf.print();
    }

    @Test
    void getProb2() {
        MilliTimer timer = new MilliTimer();
        BNode root = pbn.getBNode(0);
        BNode leaf = pbn.getBNode(N - 1);
        for (Object sym : model.getDomain().getValues()) {
            timer.start("Reset");
            pbn.getBN().resetNodes();
            timer.stopStart("Reset", "setInst");
            leaf.setInstance(sym);
            timer.stopStart("setInst","VE inst");
            VarElim ve = new VarElim();
            ve.instantiate(pbn.getBN());
            timer.stopStart("VE inst","makeMPE");
            Query q = ve.makeMPE(root.getVariable());
            timer.stopStart("makeMPE","Infer");
            CGTable qr = (CGTable)ve.infer(q);
            timer.stopStart("Infer", "Distrib");
            Distrib d = qr.query(root.getVariable());
            timer.stop("Distrib");
            System.out.println(root + " is " + d);
        }
        timer.report();
    }

}