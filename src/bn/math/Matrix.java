/*
 bnkit -- software for building and using Bayesian networks
 Copyright (C) 2014  M. Boden et al.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package bn.math;

/**
 * This class is based on various maths libraries as indicated for each
 * function--all of which are open source.
 */
public class Matrix {

    /**
     * Matrix exponential of square matrix
     */
    public static class Exp {
        final double[][] evec; // eigenvectors
        final double[] eval_r; // eigenvalues (real)
        final double[] eval_i; // eigenvalues (imaginary)
        final double[][] ievec; // inverse eigenvectors
        final double[][] result;
        /**
         * Determine the exponential of the provided square matrix, storing relevant information in class member variables.
         * @param matrix a square matrix
         */
        public Exp(double[][] matrix) {
            int dim = matrix.length;
            if (dim == 0) 
                throw new IllegalArgumentException("Empty matrix");
            if (dim != matrix[0].length)
                throw new IllegalArgumentException("Not a square matrix");
            result = new double[dim][dim];
            evec = new double[dim][dim];
            for (int i = 0; i < dim; i++) 
                System.arraycopy(matrix[i], 0, result[i], 0, matrix[i].length);
            int[] ordr = new int[dim];
            Matrix.elmhes(result, ordr, dim);
            Matrix.eltran(result, evec, ordr, dim);
            eval_r = new double[dim];
            eval_i = new double[dim];
            Matrix.hqr2(dim, 1, dim, result, evec, eval_r, eval_i);
            ievec = new double[dim][dim];
            Matrix.luinverse(evec, ievec, dim);
        }
        public double[][] getExp() {
            return result;
        }
        public double[][] getEigvec() {
            return evec;
        }
        public double[][] getInvEigvec() {
            return ievec;
        }
        public double[] getEigval() {
            return getEigval(true);
        }
        public double[] getEigval(boolean real) {
            return (real ? eval_r: eval_i);
        }
        
    }

    public static void print(double[][] m) {
        for (int i = 0; i < m.length; i ++) {
            for (int j = 0; j < m[i].length; j ++) {
                System.out.printf("%+5.3f ", m[i][j]);
            }
            System.out.println();
        }
    }

    public static void printLaTeX(double[][] m, Object[] rowsym, Object[] colsym) {
        /*
            \begin{matrix}
    & A       &    R  &    N  &    D  &    C  &    Q  &    E  & ... \\
    A & 0.989 & 0.000 & 0.000 & 0.000 & 0.000 & 0.000 & 0.001 & ...  \\
    R & 0.001 & 0.990 & 0.000 & 0.000 & 0.000 & 0.001 & 0.000 & ...  \\
    N & 0.000 & 0.000 & 0.986 & 0.003 & 0.000 & 0.001 & 0.001 & ...  \\
    D & 0.001 & 0.000 & 0.002 & 0.990 & 0.000 & 0.000 & 0.004 & ...  \\
    C & 0.001 & 0.000 & 0.000 & 0.000 & 0.995 & 0.000 & 0.000 & ...  \\
    Q & 0.001 & 0.001 & 0.001 & 0.000 & 0.000 & 0.986 & 0.003 & ...  \\
    E & 0.001 & 0.000 & 0.000 & 0.004 & 0.000 & 0.002 & 0.988 & ...  \\
    ... & ... \\
    \end{matrix}

         */
        System.out.println("\\begin{matrix}");
        for (int j = 0; j < colsym.length; j ++)
            System.out.printf("& %5s ", colsym[j]);
        System.out.println("\\\\");
        for (int i = 0; i < m.length; i ++) {
            System.out.printf("%s ", rowsym[i]);
            for (int j = 0; j < m[i].length; j ++) {
                System.out.printf("& %5.3f ", m[i][j]);
            }
            System.out.println("\\\\");
        }
        System.out.println("\\end{matrix}");
    }

    /**
     * Based on NETLIB/EISPACK/elmhes.f which is a translation of the Algol 
     * procedure elmhes:
     * http://www.netlib.no/netlib/eispack/elmhes.f
     * Martin and Wilkinson, Num. Math. 12, 349-368 (1968).
     * Given a real general matrix, this subroutine reduces a submatrix to upper hessenberg 
     * form by stabilized elementary similarity transformations.

     * @param a contains the input matrix, and is modified "in place" to on exit contain 
     * the hessenberg matrix.  
     * @param ordr the multipliers which were used in the reduction are stored 
     * in the remaining triangle under the hessenberg matrix.
     * @param n is the order of the matrix.
     */
    public static void elmhes(double[][] a, int[] ordr, int n) {
        int m, j, i;
        double y, x;

        for (i = 0; i < n; i++) {
            ordr[i] = 0;
        }
        for (m = 2; m < n; m++) {
            x = 0.0;
            i = m;
            for (j = m; j <= n; j++) {
                if (Math.abs(a[j - 1][m - 2]) > Math.abs(x)) {
                    x = a[j - 1][m - 2];
                    i = j;
                }
            }
            ordr[m - 1] = i;
            if (i != m) {
                for (j = m - 2; j < n; j++) {
                    y = a[i - 1][j];
                    a[i - 1][j] = a[m - 1][j];
                    a[m - 1][j] = y;
                }
                for (j = 0; j < n; j++) {
                    y = a[j][i - 1];
                    a[j][i - 1] = a[j][m - 1];
                    a[j][m - 1] = y;
                }
            }
            if (x != 0.0) {
                for (i = m; i < n; i++) {
                    y = a[i][m - 2];
                    if (y != 0.0) {
                        y /= x;
                        a[i][m - 2] = y;
                        for (j = m - 1; j < n; j++) {
                            a[i][j] -= y * a[m - 1][j];
                        }
                        for (j = 0; j < n; j++) {
                            a[j][m - 1] += y * a[j][i];
                        }
                    }
                }
            }
        }
    }

    /**
     * Based on NETLIB/EISPACK/eltran.f which is a translation of the Algol 
     * procedure eltran:
     * http://www.netlib.no/netlib/eispack/eltran.f
     * Peters and Wilkinson, Num. Math. 16, 181-204 (1970).
     * This subroutine accumulates the stabilized elementary similarity transformations used 
     * in the reduction of a real general matrix to upper hessenberg form by  elmhes.
     * @param a matrix
     * @param zz contains the transformation matrix produced in the
        reduction by  elmhes.
     * @param ordr contains the multipliers which were used in the
        reduction by  elmhes in its lower triangle
        below the subdiagonal.
     * @param n is the order of the matrix.
     */
    public static void eltran(double[][] a, double[][] zz, int[] ordr, int n) {
        int i, j, m;

        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                zz[i][j] = 0.0;
                zz[j][i] = 0.0;
            }
            zz[i][i] = 1.0;
        }
        if (n <= 2) {
            return;
        }
        for (m = n - 1; m >= 2; m--) {
            for (i = m; i < n; i++) {
                zz[i][m - 1] = a[i][m - 2];
            }
            i = ordr[m - 1];
            if (i != m) {
                for (j = m - 1; j < n; j++) {
                    zz[m - 1][j] = zz[i - 1][j];
                    zz[i - 1][j] = 0.0;
                }
                zz[i - 1][m - 1] = 1.0;
            }
        }
    }

    /**
     * Complex division c = a/b.
     * Based on NETLIB/EISPACK/cdiv.f 
     * http://www.netlib.no/netlib/eispack/cdiv.f

     * @param ar a real part
     * @param ai a imaginary part
     * @param br b real part
     * @param bi b imaginary part
     * @param c an array with two elements: real and imaginary part of the result
     */
    private static void mcdiv(double ar, double ai, double br, double bi, double[] c) {
        double s, ars, ais, brs, bis;

        s = Math.abs(br) + Math.abs(bi);
        ars = ar / s;
        ais = ai / s;
        brs = br / s;
        bis = bi / s;
        s = brs * brs + bis * bis;
        c[0] = (ars * brs + ais * bis) / s;
        c[1] = (ais * brs - ars * bis) / s;
    }

    /**
     * Based on NETLIB/EISPACK/hqr2.f which is a translation of the Algol 
     * procedure hqr2: http://www.netlib.no/netlib/eispack/hqr2.f
     * Peters and Wilkinson, Num. Math. 16, 181-204 (1970).
     * This subroutine finds the eigenvalues and eigenvectors
     * of a real upper hessenberg matrix by the qr method.  the
     * eigenvectors of a real general matrix can also be found
     * if  elmhes  and  eltran have been used to reduce this general matrix 
     * to hessenberg form and to accumulate the similarity transformations.
     * @param n  is the order of the matrix
     * @param low integer determined by the balancing subroutine  balanc. if  balanc  has not been used, set low=1, igh=n.
     * @param hgh integer determined by the balancing subroutine  balanc.  if  balanc  has not been used, set low=1, igh=n.
     * @param h contains the upper hessenberg matrix.
     * @param zz on input contains the transformation matrix produced by eltran after the reduction by  elmhes, or by  ortran  after the
     * reduction by  orthes, if performed. If the eigenvectors of the hessenberg matrix are desired, z must contain the identity matrix.
     * On output, contains the real and imaginary parts of the eigenvectors. If the i-th eigenvalue is real, the i-th column of z 
     * contains its eigenvector. If the i-th eigenvalue is complex with positive imaginary part, the i-th and (i+1)-th columns of z 
     * contain the real and imaginary parts of its eigenvector.  the eigenvectors are unnormalized. If an error exit is made, 
     * none of the eigenvectors has been found.
     * @param wr contain the real parts of the eigenvalues. The eigenvalues are unordered except that complex conjugate pairs
     * of values appear consecutively with the eigenvalue having the positive imaginary part first. If an
     * error exit is made, the eigenvalues should be correct for indices ierr+1,...,n.
     * @param wi contain the imaginary parts of the eigenvalues. The eigenvalues are unordered except that complex conjugate pairs 
     * of values appear consecutively with the eigenvalue having the positive imaginary part first. If an error exit is made, 
     * the eigenvalues should be correct for indices ierr+1,...,n.     
     * @throws ArithmeticException if the limit of 30*n iterations is exhausted while the j-th eigenvalue is being sought.
     */
    public static void hqr2(int n, int low, int hgh, double[][] h, double[][] zz,
            double[] wr, double[] wi) throws ArithmeticException {
        int i, j, k, l = 0, m, en, na, itn, its;
        double p = 0, q = 0, r = 0, s = 0, t, w, x = 0, y, ra, sa, vi, vr, z = 0, norm, tst1, tst2;
        double[] c = new double[2]; // to contain results of complex division
        boolean notLast;

        norm = 0.0;
        k = 1;
        /* store isolated roots and compute matrix norm */
        for (i = 0; i < n; i++) {
            for (j = k - 1; j < n; j++) {
                norm += Math.abs(h[i][j]);
            }
            k = i + 1;
            if (i + 1 < low || i + 1 > hgh) {
                wr[i] = h[i][i];
                wi[i] = 0.0;
            }
        }
        en = hgh;
        t = 0.0;
        itn = n * 30;
        while (en >= low) {	/* search for next eigenvalues */

            its = 0;
            na = en - 1;
            while (en >= 1) {
                /* look for single small sub-diagonal element */
                boolean fullLoop = true;
                for (l = en; l > low; l--) {
                    s = Math.abs(h[l - 2][l - 2]) + Math.abs(h[l - 1][l - 1]);
                    if (s == 0.0) {
                        s = norm;
                    }
                    tst1 = s;
                    tst2 = tst1 + Math.abs(h[l - 1][l - 2]);
                    if (tst2 == tst1) {
                        fullLoop = false;
                        break;
                    }
                }
                if (fullLoop) {
                    l = low;
                }

                x = h[en - 1][en - 1];	/* form shift */

                if (l == en || l == na) {
                    break;
                }
                if (itn == 0) {
                    /* eigenvalues have not converged */
                    System.out.println("Eigenvalues not converged");
                    throw new ArithmeticException();
                }
                y = h[na - 1][na - 1];
                w = h[en - 1][na - 1] * h[na - 1][en - 1];
                /* form exceptional shift */
                if (its == 10 || its == 20) {
                    t += x;
                    for (i = low - 1; i < en; i++) {
                        h[i][i] -= x;
                    }
                    s = Math.abs(h[en - 1][na - 1]) + Math.abs(h[na - 1][en - 3]);
                    x = 0.75 * s;
                    y = x;
                    w = -0.4375 * s * s;
                }
                its++;
                itn--;
                /* look for two consecutive small sub-diagonal elements */
                for (m = en - 2; m >= l; m--) {
                    z = h[m - 1][m - 1];
                    r = x - z;
                    s = y - z;
                    p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
                    q = h[m][m] - z - r - s;
                    r = h[m + 1][m];
                    s = Math.abs(p) + Math.abs(q) + Math.abs(r);
                    p /= s;
                    q /= s;
                    r /= s;
                    if (m == l) {
                        break;
                    }
                    tst1 = Math.abs(p) * (Math.abs(h[m - 2][m - 2]) + Math.abs(z) + Math.abs(h[m][m]));
                    tst2 = tst1 + Math.abs(h[m - 1][m - 2]) * (Math.abs(q) + Math.abs(r));
                    if (tst2 == tst1) {
                        break;
                    }
                }
                for (i = m + 2; i <= en; i++) {
                    h[i - 1][i - 3] = 0.0;
                    if (i != m + 2) {
                        h[i - 1][i - 4] = 0.0;
                    }
                }
                for (k = m; k <= na; k++) {
                    notLast = k != na;
                    if (k != m) {
                        p = h[k - 1][k - 2];
                        q = h[k][k - 2];
                        r = 0.0;
                        if (notLast) {
                            r = h[k + 1][k - 2];
                        }
                        x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                        if (x != 0.0) {
                            p /= x;
                            q /= x;
                            r /= x;
                        }
                    }
                    if (x != 0.0) {
                        if (p < 0.0) {	/* sign */

                            s = -Math.sqrt(p * p + q * q + r * r);
                        } else {
                            s = Math.sqrt(p * p + q * q + r * r);
                        }
                        if (k != m) {
                            h[k - 1][k - 2] = -s * x;
                        } else if (l != m) {
                            h[k - 1][k - 2] = -h[k - 1][k - 2];
                        }
                        p += s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q /= p;
                        r /= p;
                        if (!notLast) {
                            for (j = k - 1; j < n; j++) {	/* row modification */

                                p = h[k - 1][j] + q * h[k][j];
                                h[k - 1][j] -= p * x;
                                h[k][j] -= p * y;
                            }
                            j = (en < (k + 3)) ? en : (k + 3); /* min */

                            for (i = 0; i < j; i++) {	/* column modification */

                                p = x * h[i][k - 1] + y * h[i][k];
                                h[i][k - 1] -= p;
                                h[i][k] -= p * q;
                            }
                            /* accumulate transformations */
                            for (i = low - 1; i < hgh; i++) {
                                p = x * zz[i][k - 1] + y * zz[i][k];
                                zz[i][k - 1] -= p;
                                zz[i][k] -= p * q;
                            }
                        } else {
                            for (j = k - 1; j < n; j++) {	/* row modification */

                                p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
                                h[k - 1][j] -= p * x;
                                h[k][j] -= p * y;
                                h[k + 1][j] -= p * z;
                            }
                            j = (en < (k + 3)) ? en : (k + 3); /* min */

                            for (i = 0; i < j; i++) {	/* column modification */

                                p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
                                h[i][k - 1] -= p;
                                h[i][k] -= p * q;
                                h[i][k + 1] -= p * r;
                            }
                            /* accumulate transformations */
                            for (i = low - 1; i < hgh; i++) {
                                p = x * zz[i][k - 1] + y * zz[i][k]
                                        + z * zz[i][k + 1];
                                zz[i][k - 1] -= p;
                                zz[i][k] -= p * q;
                                zz[i][k + 1] -= p * r;
                            }
                        }
                    }
                }				 /* for k */

            }					 /* while infinite loop */

            if (l == en) {				 /* one root found */

                h[en - 1][en - 1] = x + t;
                wr[en - 1] = h[en - 1][en - 1];
                wi[en - 1] = 0.0;
                en = na;
                continue;
            }
            y = h[na - 1][na - 1];
            w = h[en - 1][na - 1] * h[na - 1][en - 1];
            p = (y - x) / 2.0;
            q = p * p + w;
            z = Math.sqrt(Math.abs(q));
            h[en - 1][en - 1] = x + t;
            x = h[en - 1][en - 1];
            h[na - 1][na - 1] = y + t;
            if (q >= 0.0) {	 /* real pair */

                if (p < 0.0) {	/* sign */

                    z = p - Math.abs(z);
                } else {
                    z = p + Math.abs(z);
                }
                wr[na - 1] = x + z;
                wr[en - 1] = wr[na - 1];
                if (z != 0.0) {
                    wr[en - 1] = x - w / z;
                }
                wi[na - 1] = 0.0;
                wi[en - 1] = 0.0;
                x = h[en - 1][na - 1];
                s = Math.abs(x) + Math.abs(z);
                p = x / s;
                q = z / s;
                r = Math.sqrt(p * p + q * q);
                p /= r;
                q /= r;
                for (j = na - 1; j < n; j++) {	/* row modification */

                    z = h[na - 1][j];
                    h[na - 1][j] = q * z + p * h[en - 1][j];
                    h[en - 1][j] = q * h[en - 1][j] - p * z;
                }
                for (i = 0; i < en; i++) {	/* column modification */

                    z = h[i][na - 1];
                    h[i][na - 1] = q * z + p * h[i][en - 1];
                    h[i][en - 1] = q * h[i][en - 1] - p * z;
                }
                /* accumulate transformations */
                for (i = low - 1; i < hgh; i++) {
                    z = zz[i][na - 1];
                    zz[i][na - 1] = q * z + p * zz[i][en - 1];
                    zz[i][en - 1] = q * zz[i][en - 1] - p * z;
                }
            } else {	/* complex pair */

                wr[na - 1] = x + p;
                wr[en - 1] = x + p;
                wi[na - 1] = z;
                wi[en - 1] = -z;
            }
            en -= 2;
        } /* while en >= low */
        /* backsubstitute to find vectors of upper triangular form */

        if (norm != 0.0) {
            for (en = n; en >= 1; en--) {
                p = wr[en - 1];
                q = wi[en - 1];
                na = en - 1;
                if (q == 0.0) {/* real vector */

                    m = en;
                    h[en - 1][en - 1] = 1.0;
                    if (na != 0) {
                        for (i = en - 2; i >= 0; i--) {
                            w = h[i][i] - p;
                            r = 0.0;
                            for (j = m - 1; j < en; j++) {
                                r += h[i][j] * h[j][en - 1];
                            }
                            if (wi[i] < 0.0) {
                                z = w;
                                s = r;
                            } else {
                                m = i + 1;
                                if (wi[i] == 0.0) {
                                    t = w;
                                    if (t == 0.0) {
                                        tst1 = norm;
                                        t = tst1;
                                        do {
                                            t = 0.01 * t;
                                            tst2 = norm + t;
                                        } while (tst2 > tst1);
                                    }
                                    h[i][en - 1] = -(r / t);
                                } else {	/* solve real equations */

                                    x = h[i][i + 1];
                                    y = h[i + 1][i];
                                    q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
                                    t = (x * s - z * r) / q;
                                    h[i][en - 1] = t;
                                    if (Math.abs(x) > Math.abs(z)) {
                                        h[i + 1][en - 1] = (-r - w * t) / x;
                                    } else {
                                        h[i + 1][en - 1] = (-s - y * t) / z;
                                    }
                                }
                                /* overflow control */
                                t = Math.abs(h[i][en - 1]);
                                if (t != 0.0) {
                                    tst1 = t;
                                    tst2 = tst1 + 1.0 / tst1;
                                    if (tst2 <= tst1) {
                                        for (j = i; j < en; j++) {
                                            h[j][en - 1] /= t;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (q > 0.0) {
                    m = na;
                    if (Math.abs(h[en - 1][na - 1]) > Math.abs(h[na - 1][en - 1])) {
                        h[na - 1][na - 1] = q / h[en - 1][na - 1];
                        h[na - 1][en - 1] = (p - h[en - 1][en - 1]) / h[en - 1][na - 1];
                    } else {
                        mcdiv(0.0, -h[na - 1][en - 1], h[na - 1][na - 1] - p, q, c);
                        h[na - 1][na - 1] = c[0];
                        h[na - 1][en - 1] = c[1];
                    }
                    h[en - 1][na - 1] = 0.0;
                    h[en - 1][en - 1] = 1.0;
                    if (en != 2) {
                        for (i = en - 3; i >= 0; i--) {
                            w = h[i][i] - p;
                            ra = 0.0;
                            sa = 0.0;
                            for (j = m - 1; j < en; j++) {
                                ra += h[i][j] * h[j][na - 1];
                                sa += h[i][j] * h[j][en - 1];
                            }
                            if (wi[i] < 0.0) {
                                z = w;
                                r = ra;
                                s = sa;
                            } else {
                                m = i + 1;
                                if (wi[i] == 0.0) {
                                    mcdiv(-ra, -sa, w, q, c);
                                    h[i][na - 1] = c[0];
                                    h[i][en - 1] = c[1];
                                } else {	/* solve complex equations */

                                    x = h[i][i + 1];
                                    y = h[i + 1][i];
                                    vr = (wr[i] - p) * (wr[i] - p);
                                    vr = vr + wi[i] * wi[i] - q * q;
                                    vi = (wr[i] - p) * 2.0 * q;
                                    if (vr == 0.0 && vi == 0.0) {
                                        tst1 = norm * (Math.abs(w) + Math.abs(q) + Math.abs(x)
                                                + Math.abs(y) + Math.abs(z));
                                        vr = tst1;
                                        do {
                                            vr = 0.01 * vr;
                                            tst2 = tst1 + vr;
                                        } while (tst2 > tst1);
                                    }
                                    mcdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, c);
                                    h[i][na - 1] = c[0];
                                    h[i][en - 1] = c[1];
                                    if (Math.abs(x) > Math.abs(z) + Math.abs(q)) {
                                        h[i + 1][na - 1] = (q * h[i][en - 1]
                                                - w * h[i][na - 1] - ra) / x;
                                        h[i + 1][en - 1] = (-sa - w * h[i][en - 1]
                                                - q * h[i][na - 1]) / x;
                                    } else {
                                        mcdiv(-r - y * h[i][na - 1], -s - y * h[i][en - 1], z, q, c);
                                        h[i + 1][na - 1] = c[0];
                                        h[i + 1][en - 1] = c[1];
                                    }
                                }
                                /* overflow control */
                                t = (Math.abs(h[i][na - 1]) > Math.abs(h[i][en - 1]))
                                        ? Math.abs(h[i][na - 1]) : Math.abs(h[i][en - 1]);
                                if (t != 0.0) {
                                    tst1 = t;
                                    tst2 = tst1 + 1.0 / tst1;
                                    if (tst2 <= tst1) {
                                        for (j = i; j < en; j++) {
                                            h[j][na - 1] /= t;
                                            h[j][en - 1] /= t;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            /* end back substitution. vectors of isolated roots */
            for (i = 0; i < n; i++) {
                if (i + 1 < low || i + 1 > hgh) {
                    for (j = i; j < n; j++) {
                        zz[i][j] = h[i][j];
                    }
                }
            }
            /* multiply by transformation matrix to give vectors of
             * original full matrix. */
            for (j = n - 1; j >= low - 1; j--) {
                m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */

                for (i = low - 1; i < hgh; i++) {
                    z = 0.0;
                    for (k = low - 1; k < m; k++) {
                        z += zz[i][k] * h[k][j];
                    }
                    zz[i][j] = z;
                }
            }
        }
    }

    public static final double EPSILON = 2.220446049250313E-16;

    /**
     * Find the inverse of a matrix using LU decomposition 
     * @param inmat
     * @param imtrx
     * @param size
     * @throws IllegalArgumentException if the matrix is singular
     */
    public static void luinverse(double[][] inmat, double[][] imtrx, int size) throws IllegalArgumentException {
        int i, j, k, l, maxi = 0, idx, ix, jx;
        double sum, tmp, maxb, aw;
        int[] index;
        double[] wk;
        double[][] omtrx;

        index = new int[size];
        omtrx = new double[size][size];

        /* copy inmat to omtrx */
        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                omtrx[i][j] = inmat[i][j];
            }
        }

        wk = new double[size];
        aw = 1.0;
        for (i = 0; i < size; i++) {
            maxb = 0.0;
            for (j = 0; j < size; j++) {
                if (Math.abs(omtrx[i][j]) > maxb) {
                    maxb = Math.abs(omtrx[i][j]);
                }
            }
            if (maxb == 0.0) {
                /* Singular matrix */
                System.out.println("Singular matrix encountered");
                throw new IllegalArgumentException("Singular matrix");
            }
            wk[i] = 1.0 / maxb;
        }
        for (j = 0; j < size; j++) {
            for (i = 0; i < j; i++) {
                sum = omtrx[i][j];
                for (k = 0; k < i; k++) {
                    sum -= omtrx[i][k] * omtrx[k][j];
                }
                omtrx[i][j] = sum;
            }
            maxb = 0.0;
            for (i = j; i < size; i++) {
                sum = omtrx[i][j];
                for (k = 0; k < j; k++) {
                    sum -= omtrx[i][k] * omtrx[k][j];
                }
                omtrx[i][j] = sum;
                tmp = wk[i] * Math.abs(sum);
                if (tmp >= maxb) {
                    maxb = tmp;
                    maxi = i;
                }
            }
            if (j != maxi) {
                for (k = 0; k < size; k++) {
                    tmp = omtrx[maxi][k];
                    omtrx[maxi][k] = omtrx[j][k];
                    omtrx[j][k] = tmp;
                }
                aw = -aw;
                wk[maxi] = wk[j];
            }
            index[j] = maxi;
            if (omtrx[j][j] == 0.0) {
                omtrx[j][j] = EPSILON;
            }
            if (j != size - 1) {
                tmp = 1.0 / omtrx[j][j];
                for (i = j + 1; i < size; i++) {
                    omtrx[i][j] *= tmp;
                }
            }
        }
        for (jx = 0; jx < size; jx++) {
            for (ix = 0; ix < size; ix++) {
                wk[ix] = 0.0;
            }
            wk[jx] = 1.0;
            l = -1;
            for (i = 0; i < size; i++) {
                idx = index[i];
                sum = wk[idx];
                wk[idx] = wk[i];
                if (l != -1) {
                    for (j = l; j < i; j++) {
                        sum -= omtrx[i][j] * wk[j];
                    }
                } else if (sum != 0.0) {
                    l = i;
                }
                wk[i] = sum;
            }
            for (i = size - 1; i >= 0; i--) {
                sum = wk[i];
                for (j = i + 1; j < size; j++) {
                    sum -= omtrx[i][j] * wk[j];
                }
                wk[i] = sum / omtrx[i][i];
            }
            for (ix = 0; ix < size; ix++) {
                imtrx[ix][jx] = wk[ix];
            }
        }
    }

}
