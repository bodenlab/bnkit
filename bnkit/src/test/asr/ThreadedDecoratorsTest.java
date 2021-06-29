package asr;

import bn.ctmc.matrix.JC;
import bn.ctmc.matrix.JTT;
import dat.Enumerable;
import dat.phylo.*;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.TimeUnit;

class ThreadedDecoratorsTest {

    class MyDecor implements TreeDecor<Integer> {
        public String name;
        public PhyloBN pbn;
        public int secs = 0;
        public MyDecor(String name, IdxTree tree) {
            this.name = name;
            pbn = PhyloBN.create(tree, new JC(1, Enumerable.aacid.getValues()));

        }
        @Override
        public Integer getDecoration(int idx) {
            return secs;
        }

        @Override
        public void decorate(TreeInstance ti) {
            try {
                this.secs = (int) (Math.random() * 2) + 1;
                //System.out.println("Executing : " + name + " for " + secs + " seconds");
                TimeUnit.SECONDS.sleep(secs);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
    static IdxTree testtree;

    @BeforeAll
    static void setThingsUp() {
        Tree tree;
        try {
            tree = Tree.load("bnkit/src/test/resources/cyp2u1_recon.nwk", "newick");
            testtree = (IdxTree)tree;
        } catch (IOException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        }
    }

    @Test
    void test10xMyDecor()  {
        int N = 10;
        MyDecor[] decors = new MyDecor[N];
        TreeInstance[] tis = new TreeInstance[N];
        for (int i = 0; i < N; i ++) {
            decors[i] = new MyDecor("Noname" + (i+1), testtree);
            tis[i] = testtree.getInstance(new Object[]{},new Object[]{});
        }
        ThreadedDecorators pool = new ThreadedDecorators(decors, tis, 5);
        try {
            Map<Integer, Object> ret = pool.runBatch();
            for (Map.Entry<Integer, Object> entry : ret.entrySet()) {
                Object res = ((MyDecor)entry.getValue()).getDecoration(0);
                //System.out.println("\tCompleted " + entry.getKey() + "@" + res);
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }
    @Test
    void testJoint()  {
        int NPOS = 50;
        int NTHREADS = 5;
        long START_TIME = System.currentTimeMillis();

        Random rand = new Random(START_TIME);
        TreeDecor[] decors = new TreeDecor[NPOS];
        TreeInstance[] tis = new TreeInstance[NPOS];
        for (int i = 0; i < NPOS; i ++) {
            decors[i] = new MaxLhoodJoint(testtree, new JTT());
            int[] leaves = testtree.getLeaves();
            Object[] labels = new Object[leaves.length];
            Object[] values = new Object[leaves.length];
            for (int j = 0; j < leaves.length; j ++) {
                labels[j] = testtree.getLabel(leaves[j]);
                values[j] = Enumerable.aacid.get(rand.nextInt(Enumerable.aacid.size()));
            }
            tis[i] = testtree.getInstance(labels, values);
        }
        ThreadedDecorators pool = new ThreadedDecorators(decors, tis, NTHREADS);
        try {
            START_TIME = System.currentTimeMillis();
            Map<Integer, Object> ret = pool.runBatch();
            long ELAPSED_TIME = (System.currentTimeMillis() - START_TIME);
            System.out.println(String.format("Done in %d min, %d sec", TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME),
                    TimeUnit.MILLISECONDS.toSeconds(ELAPSED_TIME) - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(ELAPSED_TIME))));
            int[] ancs = testtree.getAncestors();
            Object[][] preds = new Object[ancs.length][NPOS];
            for (Map.Entry<Integer, Object> entry : ret.entrySet()) {
                int pos = entry.getKey();
                for (int a = 0; a < ancs.length; a ++)
                    preds[a][pos] = ((TreeDecor)entry.getValue()).getDecoration(ancs[a]);
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }

}