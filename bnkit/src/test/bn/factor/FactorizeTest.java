package bn.factor;

import bn.Predef;
import bn.prob.GaussianDistrib;
import dat.EnumVariable;
import dat.Variable;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class FactorizeTest {

    /** Calculate products linearly, don't consider order */
    protected static final int POOL_OPTION_LINEAR = 0;
    /** Calculate products according to binary tree, optimised to minimise computational cost */
    protected static final int POOL_OPTION_TREE = 1;
    /*************************** Methods to test AbstractFactor **************************/

    protected static Variable[] getVariablePool(long seed, int n) {
        Random random = new Random(seed);
        Variable[] vars = new Variable[n];
        for (int i = 0; i < vars.length; i++) {
            int type = random.nextInt(5);
            Variable var = null;
            switch (type) {
                case 0:
                    var = Predef.Boolean();
                    break;
                case 1:
                    var = Predef.Nominal("a", "b", "c");
                    break;
                case 2:
                    var = Predef.NucleicAcid();
                    break;
                case 3:
                    var = Predef.Real();
                    break;
                case 4:
                    var = Predef.Number(random.nextInt(8) + 2);
                    break;
            }
            vars[i] = var;
        }
        return vars;
    }

    protected static Variable[] getSubset(long seed, Variable[] vars, int n) {
        Random random = new Random(seed);
        Set<Variable> unique = new HashSet<>();
        n = Math.min(vars.length, n);
        Variable[] subset = new Variable[n];
        while (unique.size() < n) {
            unique.add(vars[random.nextInt(vars.length)]);
        }
        unique.toArray(subset);
        return subset;
    }

    protected static AbstractFactor[] getFactorPool(long seed, Variable[] vars, int n) {
        Random random = new Random(seed);
        int M = Math.abs((int) (random.nextGaussian() * n) + 1);
        int N = Math.abs((int) (random.nextGaussian() * n) + 1);
        AbstractFactor[] dfs = new DenseFactor[N];
        for (int i = 0; i < dfs.length; i++) {
            int nvars = random.nextInt(Math.max(1, M));
            if (nvars > 0) {
                Variable[] myvars = getSubset(random.nextInt(), vars, nvars);
                dfs[i] = new DenseFactor(myvars);
            } else {
                dfs[i] = new DenseFactor();
            }
            int npop = dfs[i].getSize(); //random.nextInt(dfs[i].getSize()); // number of entries to populate
            for (int j = 0; j < npop; j++) {
                if (dfs[i].getSize() == 1) {
                    dfs[i].setValue(Math.abs(random.nextGaussian()) / npop);
                    if (dfs[i].isJDF()) {
                        for (Variable nvar : dfs[i].getNonEnumVars()) {
                            dfs[i].setDistrib(nvar, new GaussianDistrib(random.nextGaussian() * random.nextInt(100), Math.abs(random.nextGaussian() * (random.nextInt(10) + 1))));
                        }
                    }
                } else {
                    int index = random.nextInt(npop);
                    dfs[i].setValue(index, Math.abs(random.nextGaussian()) / npop);
                    if (dfs[i].isJDF()) {
                        for (Variable nvar : dfs[i].getNonEnumVars()) {
                            dfs[i].setDistrib(index, nvar, new GaussianDistrib(random.nextGaussian() * random.nextInt(100), Math.abs(random.nextGaussian() * (random.nextInt(10) + 1))));
                        }
                    }
                }
            }
        }
        return dfs;
    }

    protected static AbstractFactor getProductBenchmarked(Factorize.FactorProductTree node) {
        if (node.getFactor() != null) {
            return node.getFactor();
        }
        AbstractFactor X = getProductBenchmarked(node.x);
        AbstractFactor Y = getProductBenchmarked(node.y);
        long startTime = System.nanoTime();
        AbstractFactor f = Factorize.getProduct(X, Y);
        long endTime = System.nanoTime();
        for (int j = 0; j < 20; j++) {
            assertTrue(testProductIntegrity(j, X, Y, f));
/*
            if (testProductIntegrity(j, X, Y, f) == false) {
                System.err.println("Test failed");
                X.display();
                Y.display();
                f.display();
                }
 */
        }
        int overlap = Factorize.getOverlap(X, Y);
        int minevars = Math.min(X.nEVars, Y.nEVars);
        int maxevars = Math.max(Y.nEVars, Y.nEVars);
        System.out.println(maxevars + "\t" + minevars + "\t" + overlap + "\t" + (minevars == 0 ? 0.0 : overlap / (float) minevars) + "\t" + Factorize.getComplexity(X, Y, false) + "\t" + Factorize.getComplexity(X, Y, true) + "\t" + (endTime - startTime) / 100000.0);
        node.setFactor(f);
        return f;
    }

    protected static AbstractFactor getProductBenchmarked(AbstractFactor[] factors) {
        if (factors.length == 0) {
            return null;
        }
        AbstractFactor R = factors[0];
        for (int i = 1; i < factors.length; i++) {
            AbstractFactor X = R;
            AbstractFactor Y = factors[i];
            long startTime = System.nanoTime();
            R = Factorize.getProduct(X, Y);
            long endTime = System.nanoTime();
            for (int j = 0; j < 20; j++) {
                assertTrue(testProductIntegrity(j, X, Y, R));
/*                if (testProductIntegrity(j, X, Y, R) == false) {
                    System.err.println("Test failed");
                    X.display();
                    Y.display();
                    R.display();
                } */
            }
            int overlap = Factorize.getOverlap(X, Y);
            int minevars = Math.min(X.nEVars, Y.nEVars);
            int maxevars = Math.max(Y.nEVars, Y.nEVars);
            System.out.println(maxevars + "\t" + minevars + "\t" + overlap + "\t" + (minevars == 0 ? 0.0 : overlap / (float) minevars) + "\t" + Factorize.getComplexity(X, Y, false) + "\t" + Factorize.getComplexity(X, Y, true) + "\t" + (endTime - startTime) / 100000.0);
        }
        return R;
    }

    protected static AbstractFactor productPool(AbstractFactor[] dfs, int option) {
        if (dfs.length == 0) {
            return null;
        }
        if (dfs.length == 1) {
            return dfs[0];
        }
        AbstractFactor f = null;
        long startTime = System.nanoTime();
        switch (option) {
            case POOL_OPTION_LINEAR:
                f = getProductBenchmarked(dfs);
                break;
            case POOL_OPTION_TREE:
                Factorize.FactorProductTree tree = Factorize.getProductTree(dfs);
                f = getProductBenchmarked(tree);
                break;
        }
        long endTime = System.nanoTime();
        System.out.println("\t\t\t\t\t\t\t" + (endTime - startTime) / 100000.0);
        if (f.nEVars > 0) {
            Random rand = new Random(endTime);
            int n = rand.nextInt(f.nEVars);
            Variable[] sumout = new Variable[n];
            for (int i = 0; i < n; i++) {
                sumout[i] = f.evars[rand.nextInt(n)];
            }
            if (rand.nextBoolean()) {
                Factorize.getMargin(f, sumout);
            } else {
                Factorize.getMaxMargin(f, sumout);
            }
            return f;
        }
        return f;
    }

    protected static boolean testProductIntegrity(long seed, AbstractFactor X, AbstractFactor Y, AbstractFactor PROD) {
        Random rand = new Random(seed);
        int p = rand.nextInt(PROD.getSize());
        Object[] pkey;
        if (PROD.getSize() == 1) {
            pkey = new Object[0];
        } else {
            pkey = PROD.getKey(p);
        }
        EnumVariable[] pevars = PROD.getEnumVars();
        EnumVariable[] xevars = X.getEnumVars();
        EnumVariable[] yevars = Y.getEnumVars();
        Object[] xkey = new Object[X.nEVars];
        Object[] ykey = new Object[Y.nEVars];
        for (int i = 0; i < pkey.length; i++) {
            if (!pevars[i].getDomain().isValid(pkey[i])) {
                return false;
            }
            int xi = -1;
            for (int j = 0; j < X.nEVars; j++) {
                if (xevars[j].equals(pevars[i])) {
                    xi = j;
                }
            }
            int yi = -1;
            for (int j = 0; j < Y.nEVars; j++) {
                if (yevars[j].equals(pevars[i])) {
                    yi = j;
                }
            }
            if (yi != -1) {
                ykey[yi] = pkey[i];
            }
            if (xi != -1) {
                xkey[xi] = pkey[i];
            }
        }
        double xval;
        if (xkey.length != 0) {
            xval = X.getValue(X.getIndex(xkey));
        } else {
            xval = X.getValue();
        }
        double yval;
        if (ykey.length != 0) {
            yval = Y.getValue(Y.getIndex(ykey));
        } else {
            yval = Y.getValue();
        }
        double pval;
        if (pkey.length == 0) {
            pval = PROD.getValue();
        } else {
            pval = PROD.getValue(p);
        }
        return (xval * yval <= pval *1.01 && xval * yval >= pval *0.99);
    }

    protected static void testCrossRef(long seed, Variable[] vars) {
        Random rand = new Random(seed);
        for (int t = 0; t < 20; t++) {
            Variable[] x = new Variable[rand.nextInt(vars.length - 1) + 1];
            Variable[] y = new Variable[rand.nextInt(vars.length - 1) + 1];
            for (int i = 0; i < x.length; i++) {
                x[i] = vars[rand.nextInt(vars.length)];
            }
            for (int i = 0; i < y.length; i++) {
                y[i] = vars[rand.nextInt(vars.length)];
            }
            x = Factorize.getNonredundant(x);
            y = Factorize.getNonredundant(y);
            int[] xcross = new int[x.length];
            int[] ycross = new int[y.length];
            Factorize.getCrossref(x, xcross, y, ycross);
            Variable[] xx = new Variable[x.length];
            Variable[] yy = new Variable[y.length];
            for (int i = 0; i < y.length; i++) {
                if (ycross[i] != -1) {
                    xx[ycross[i]] = y[i];
                }
            }
            for (int i = 0; i < x.length; i++) {
                if (xcross[i] != -1) {
                    yy[xcross[i]] = x[i];
                }
            }
            for (int i = 0; i < x.length; i++) {
                System.out.print(x[i].toString() + ":" + ((xx[i] == null) ? ("-\t") : (xx[i].toString() + "\t")));
            }
            System.out.println();
            for (int i = 0; i < y.length; i++) {
                System.out.print(y[i].toString() + ":" + ((yy[i] == null) ? ("-\t") : (yy[i].toString() + "\t")));
            }
            System.out.println();
            System.out.println();
        }
    }

    @Test
    void getProduct() {
        System.out.println("maxEV\tminEV\tOverlap\tContain\tProduct\tPJoin\tTime (ms)");
        for (long seed = 0; seed < 50; seed++) {
            Variable[] vars = getVariablePool(seed, 10);
            // testCrossRef(seed, vars);

            AbstractFactor[] dfs = getFactorPool(seed, vars, 8);
            AbstractFactor f1 = productPool(dfs, POOL_OPTION_LINEAR);
            AbstractFactor f2 = productPool(dfs, POOL_OPTION_TREE);
            if (f1 == null && f2 == null) {
                continue;
            }
            assertFalse(f1.getSize() != f2.getSize());
            if (f1.getSize() != f2.getSize()) {
                System.err.println("Invalid product size");
            }
            if (f1.getSize() == 1) {
                assertFalse(f1.getValue() < f2.getValue() * 0.999 || f1.getValue() > f2.getValue() * 1.001);
                if (f1.getValue() < f2.getValue() * 0.999 || f1.getValue() > f2.getValue() * 1.001) {
                    System.err.println("Invalid atomic product: " + f1.getValue() + " v " + f2.getValue());
                    System.exit(1);
                }
            } else {
                for (int i = 0; i < f1.getSize(); i++) {
                    assertFalse(f1.getValue(i) < f2.getValue(i) * 0.999 || f1.getValue(i) > f2.getValue(i) * 1.001);
                    if (f1.getValue(i) < f2.getValue(i) * 0.999 || f1.getValue(i) > f2.getValue(i) * 1.001) {
                        System.err.println("Invalid product: " + f1.getValue(i) + " v " + f2.getValue(i));
                        System.exit(1);
                    }
                }
            }

        }
    }

}