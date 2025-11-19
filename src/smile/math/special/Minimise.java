package smile.math.special;

import smile.math.Function;

public class Minimise {

    private static final int MAX = 500;
    private static final double EPS = 10e-6;
    private static final double INVPHI = (1.0 + Math.sqrt(5.0)) / 2.0 - 1.0;
    private static final double RATIO = (3.0 - Math.sqrt(5.0)) / 2;



    /**
     * Recursive loop of Brent's Minimization Method
     *
     * @param func the single variate function
     * @param a left interval
     * @param b right interval
     * @return minimum
     *
     * source: <a href="https://github.com/osveliz/numerical-veliz/blob/master/src/minimization/BrentJarratt.fsx">...</a>
     */
    public static double brent(Function func, double a, double b) {
        double x = b + INVPHI * (a - b);
        return brent(func, a, b, x, x, x, 0, 0, 0);
    }

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
     * source: <a href="https://github.com/osveliz/numerical-veliz/blob/master/src/minimization/BrentJarratt.fsx">...</a>
     */
    public static double brent(Function func,  double a, double b, double v,
                               double w, double x, double dold, double eold, int i) {

        double fv = func.apply(v);
        double fw = func.apply(w);
        double fx = func.apply(x);
        int new_i = i + 1;
        double m = 0.5 * (a + b);
        System.out.println("Iteration " + i);
        if (b - a <= EPS) {
            System.out.println("Number of iterations: " + i);
            System.out.println("Solution: " + m);
            return m;
        } else if (i > MAX) {
            System.out.println("Exhausted iterations before convergence");
            System.out.println("Non-optimal solution: " + m);
            return m;
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
            return brent(func, newa, newb, w, x, u, d, e, new_i);
        } else {

            double newa = u < x ? u : a;
            double newb = u < x ? b : u;
            if (fu <= fw || w == x) {
                return brent(func, newa, newb, w, u, x, d, e, new_i);
            } else if (fu <= fv || v == x || v == w ) {
                return brent(func, newa, newb, u, w, x, d, e, new_i);
            } else {
                return brent(func, newa, newb, v, w, x, d, e, new_i);
            }
        }
    }
}
