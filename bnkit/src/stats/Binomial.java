/*
 * binfkit -- software for bioinformatics
 * Copyright (C) 2014  M. Boden et al.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package stats;

/**
 *
 * @author mikael
 */
public class Binomial {
    static double MAXIT = 100;
    static double EPS = 3.0e-7;
    static double FPMIN = 1.0e-300;
    static double[] gamma_c = {
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.23173957245,
    0.1208650973866179e-2,
    -0.5395239384953e-5};

//################################################################################
//#       binomial_ncdf
//#
//#       Log of one minus the cumulative distribution function of the binomial 
//#       distribution.  The binomial density gives the probability of
//#       k successes in N independent trials each with probability p of success.
//#
//################################################################################
    public static double log_binomial_ncdf(int N, int k, double p) {
	if (k==0)
            return 0;
	else 
            return log_betai(k, N-k+1, p);
    }
    
    /**
     * The binomial density gives the probability of k successes in N independent trials,
     * each with probability p of success.
     * @param N number of trials 
     * @param k number of successes
     * @param p probability of success
     * @return the (cumulative) probability of (at least) k successes
     */
    public static double binomial_ncdf(int N, int k, double p) {
        return Math.exp(log_binomial_ncdf(N, k, p));
    }
    
//#
//# incomplete beta function
//#
//    public static double betai(a, b, x):
//
//  if (x<0 or x>1): die("Bad x=`" + str(x) + "' in routine betai")
//  if (x==0 or x==1):
//    bt = 0
//  else:
//    bt = exp(gammaln(a+b)-gammaln(a)-gammaln(b)+a*log(x)+b*log(1-x))
//
//  thresh = (a+1)/(a+b+2.0)
//  if (x<thresh): 
//    return(bt*betacf(a,b,x)/a)
//  else: 
//    return(1.0-bt*betacf(b,a,1.0-x)/b)
//#
//# log incomplete beta function
//#
    public static double log_betai(int a, int b, double x) {
        double log_bt = -1e300;
        if (x<0 || x>1)
            throw new RuntimeException("Bad x=" + x + "' in routine betai");
        if (x!=0 && x!=1) 
            log_bt = lgamma(a+b)-lgamma(a)-lgamma(b)+a*Math.log(x)+b*Math.log(1.0-x);
        double thresh = (a+1.0)/(a+b+2.0);
        if (x<thresh) 
            return(log_bt + Math.log(betacf(a,b,x)/a));
        else
            return(Math.log(1.0 - Math.exp(log_bt)*betacf(b,a,1.0-x)/b));
    }
//#
//#       used by betai
//#
    public static double betacf(int a, int b, double x) {

        double qab = a+b;
        double qap = a+1.0;
        double qam = a-1.0;
        double c = 1.0;
        double d = 1.0-qab*x/qap;

        if (Math.abs(d) < FPMIN)
            d = FPMIN;
        d = 1.0/d;
        double h = d;

        for (int m = 1; m < MAXIT+1; m ++) {
            double m2 = 2.0*m;
            double aa = m*(b-m)*x/((qam+m2)*(a+m2));
            d=1.0+aa*d;
            if (Math.abs(d) < FPMIN)
                d=FPMIN;
            c=1.0+aa/c;
            if (Math.abs(c) < FPMIN)
                c=FPMIN;
            d = 1.0/d;
            h *= d*c;
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));

            d=1.0+aa*d;
            if (Math.abs(d) < FPMIN)
                d=FPMIN;
            c=1.0+aa/c;
            if (Math.abs(c) < FPMIN)
                c=FPMIN;
            d = 1.0/d;

            double delta = d*c;
            h *= delta;
            if (Math.abs(delta-1.0) < EPS)
                break;
        }
        return h;
    }
    
    /*
     * Returns an approximation of the Gamma function of x r(x) = integral from 0
     * to infinity of (t^(x-1) * e^(-t) dt) with |error| < 2e-10. Laczos
     * Approximation Reference: Numerical Recipes in C
     * http://www.library.cornell.edu/nr/cbookcpdf.html
     */
    public static double gamma(double x) {
        return Math.exp(lgamma(x));
    }

    /*
     * Returns an approximation of the log of the Gamma function of x. Laczos
     * Approximation Reference: Numerical Recipes in C
     * http://www.library.cornell.edu/nr/cbookcpdf.html
     */
    public static double lgamma(double x) {
        double[] cof = { 76.18009172947146, -86.50532032941677, 24.01409824083091,
                        -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
        double y, z, ser, tmp;
        y = x;
        tmp = x + 5.5;
        tmp -= ((x + 0.5) * Math.log(tmp));
        ser = 1.000000000190015;
        for (int j = 0; j < 6; j += 1) {
            y += 1;
            ser += (cof[j] / y);
        }
        return (-tmp + Math.log(2.5066282746310005 * ser / x));
    }

    public static void main(String[] args) {
        int ntrials = 10;
        int nsuccesses = 7; // or more extreme, so "7" would mean 7, 8, 9 or 10 when ntrials is 10
        double p_success = 0.5; 
        System.out.println("ntrials = " + ntrials);
        System.out.println("nsuccesses = " + nsuccesses);
        System.out.println("p_success = " + p_success);
        System.out.println("Log binomial p = " + log_binomial_ncdf(ntrials, nsuccesses, p_success));
        System.out.println("Binomial p = " + Math.exp(log_binomial_ncdf(ntrials, nsuccesses, p_success)));
        System.out.println("Binomial p = " + binomial_ncdf(ntrials, nsuccesses, p_success));
    }
    
}
