package smile.math.special;

import asr.GRASP;
import smile.math.Function;

public class Minimise {

    private static final int MAX = 500;
    private static final double EPS = 10e-4;
    private static final double INVPHI = (1.0 + Math.sqrt(5.0)) / 2.0 - 1.0;
    private static final double RATIO = (3.0 - Math.sqrt(5.0)) / 2.0;
    public final double bestX;
    public final double bestF;
    public final double finalBracketWidth;
    public final boolean hitLowerBound;
    public final boolean hitUpperBound;
    public  static boolean VERBOSE = true;

    public Minimise(double bestX, double bestF, double finalBracketWidth,
                    boolean hitLowerBound, boolean hitUpperBound) {
        this.bestX = bestX;
        this.bestF = bestF;
        this.finalBracketWidth = finalBracketWidth;
        this.hitLowerBound = hitLowerBound;
        this.hitUpperBound = hitUpperBound;
    }

    public static Minimise brent(Function func, double a, double b) {
        double x = b + INVPHI * (a - b);
        double fx = func.apply(x);
        return brentRecursive(func, a, b, x, fx, x, fx, x, fx, 0.0, 0.0, 0);
    }

    /**
     * Recursive loop of Brent's Minimization Method
     *
     * @param func the single variate function
     * @param a left interval
     * @param b right interval
     * @return minimum
     *
     * source: <a href="https://github.com/osveliz/numerical-veliz/blob/master/src/minimization/BrentJarratt.fsx">Brent-Jarratt</a>
     */
//    public static double brent(Function func, double a, double b) {
//        double x = b + INVPHI * (a - b);
//        return brentRecursive(func, a, b, x, x, x, 0, 0, 0);
//    }

    /**
     * Recursive loop of Brent's Minimization Method
     *
     * @param func the single variate function
     * @param a left interval
     * @param b right interval
     * @param v previous iterate value
     * @param w previous iterate value
     * @param x previous iterate value
     * @param dold last delta step
     * @param eold last golden interval size
     * @param i iteration counter
     * @return minimum
     *
     * source: <a href="https://github.com/osveliz/numerical-veliz/blob/master/src/minimization/BrentJarratt.fsx">BrentJarrat.fsx</a>
     */
    public static Minimise brentRecursive(Function func, double a, double b,
                                        double v, double fv,
                                        double w, double fw,
                                        double x, double fx,
                                        double dold, double eold, int i) {

        int newI = i + 1;
        double m = 0.5 * (a + b);
        if (VERBOSE) {
            System.out.println("Iteration: " + newI + " Search interval: [" + a + ", " + b
                    + "] Best point so far: " + x + " f(x)=" + fx);
        }

        if (b - a <= EPS) {
            if (VERBOSE) System.out.println("Converged after " + newI + " iterations. Optimal x=" + x);
            return new Minimise(x, fx, b - a, isNear(x, a, EPS), isNear(x, b, EPS));
        } else if (i > MAX) {
            if (VERBOSE) System.out.println("Exhausted iterations. Approximate x=" + x);
            return new Minimise(x, fx, b - a, isNear(x, a, EPS), isNear(x, b, EPS));
        }

        double r = (x - w) * (fx - fv);
        double tq = (x - v) * (fx - fw);
        double tp = (x - v) * tq - (x - w) * r;
        double tq2 = 2.0 * (tq - r);
        double p = tq2 > 0.0 ? -tp : tp;
        double q = tq2 > 0.0 ? tq2 : -tq2;
        boolean safe = q != 0.0;
        double deltax = safe ? p / q : 0.0;

        boolean parabolic = safe && a < x + deltax && x + deltax < b && Math.abs(deltax) < 0.5 * Math.abs(eold);

        double e;
        if (parabolic) {
            e = dold;
        } else if (x < m) {
            e = b - x;
        } else {
            e = a - x;
        }

        double d = parabolic ? deltax : RATIO * e;
        double u = x + d;
        double fu = func.apply(u);

        if (fu <= fx) {
            double newa = u < x ? a : x;
            double newb = u < x ? x : b;
            return brentRecursive(func, newa, newb, w, fw, x, fx, u, fu, d, e, newI);
        } else {

            double newa = u < x ? u : a;
            double newb = u < x ? b : u;
            if (fu <= fw || w == x) {
                return brentRecursive(func, newa, newb, w, fw, u, fu, x, fx, d, e, newI);
            } else if (fu <= fv || v == x || v == w ) {
                return brentRecursive(func, newa, newb, u, fu, w, fw, x, fx, d, e, newI);
            } else {
                return brentRecursive(func, newa, newb, v, fv, w, fw, x, fx, d, e, newI);
            }
        }
    }

    /**
     * Check if a parameter is hitting the boundary of the search space.
     * @param x
     * @param edge
     * @param eps
     * @return
     */
    private static boolean isNear(double x, double edge, double eps) {
        return Math.abs(x - edge) <= Math.max(EPS, 10 * eps);
    }


    /**
     * Runs brentLogSpace once, and if the result sits on either edge of the
     * bracket (meaning the bracket was too narrow), doubles the log-space
     * window around the same center and retries once.
     */
    public static Minimise brentLogSpaceWithReexpansion(Function realFunc, double minVal, double maxVal) {
        Minimise result = brentLogSpace(realFunc, minVal, maxVal);

        if (result.hitLowerBound || result.hitUpperBound) {
            double logCenter = Math.log(result.bestX);
            double halfWindow = 0.5 * (Math.log(maxVal) - Math.log(minVal));
            double newHalfWindow = halfWindow * 2.0;

            double newMin = Math.exp(logCenter - newHalfWindow);
            double newMax = Math.exp(logCenter + newHalfWindow);

            if (VERBOSE) {
                System.out.println("Bracket [" + minVal + ", " + maxVal
                        + "] was too narrow (best value hit the edge). Re-expanding to ["
                        + newMin + ", " + newMax + "] and retrying.");
            }
            return brentLogSpace(realFunc, newMin, newMax);
        }
        return result;
    }


    public static Minimise brentLogSpace(Function func, double minVal, double maxVal) {

        double logMin = Math.log(minVal);
        double logMax = Math.log(maxVal);

        Function logSpaceFunc = t -> func.apply(Math.exp(t));

        Minimise logResult = brent(logSpaceFunc, logMin, logMax);

        double realBestX = Math.exp(logResult.bestX);
        // finalBracketWidth is still in log units here
        return new Minimise(realBestX, logResult.bestF, logResult.finalBracketWidth,
                logResult.hitLowerBound, logResult.hitUpperBound);
    }


    /**
     * Sizes the next round's log-space bracket around the current best value,
     * using Brent's own final bracket width (in log units) as the curvature
     * signal: a tightly converged (small) width means the objective was sharply
     * peaked there, so next round's window can stay tight; a wide final bracket
     * (e.g. because MAX iterations was hit) means it should stay generous.
     */
    public static double[] nextLogWindowBounds(double bestValue, double lastLogBracketWidth,
                                               double minLogHalfWindow, double maxLogHalfWindow) {
        double halfWindow = Math.max(minLogHalfWindow,
                Math.min(maxLogHalfWindow, 4.0 * Math.max(lastLogBracketWidth, EPS)));
        double logCenter = Math.log(bestValue);
        return new double[] { Math.exp(logCenter - halfWindow), Math.exp(logCenter + halfWindow) };
    }
}
