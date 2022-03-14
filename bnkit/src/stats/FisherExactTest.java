package stats;

/**
 * A fast implementation of the Fisher Exact Test based on a java script 
 * by Oyvind Langsrud (Copyright), http://www.matforsk.no/ola/fisher.htm
*/

public class FisherExactTest {

	/**
	 * Reference: "Lanczos, C. 'A precision approximation 
	 * of the gamma function', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
 	 * Translation of  Alan Miller's FORTRAN-implementation
	 * See http://lib.stat.cmu.edu/apstat/245
	 * @param z
	 * @return
	 */
	private static double lngamm(int z) {
		  double x = 0;
		  x += 0.1659470187408462e-06/(z+7);
		  x += 0.9934937113930748e-05/(z+6);
		  x -= 0.1385710331296526    /(z+5);
		  x += 12.50734324009056     /(z+4);
		  x -= 176.6150291498386     /(z+3);
		  x += 771.3234287757674     /(z+2);
		  x -= 1259.139216722289     /(z+1);
		  x += 676.5203681218835     /(z);
		  x += 0.9999999999995183;
		  return (Math.log(x)-5.58106146679532777-z+(z-0.5)*Math.log(z+6.5));
	}



	private static double lnfact(int n) {
		if (n<=1) return 0;
		return lngamm(n+1);
	}

	private static double lnbico(int n, int k) {
		return lnfact(n)-lnfact(k)-lnfact(n-k);
	}

	private static double hyper_323(int n11,int n1_, int n_1, int n) {
		return Math.exp(lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1));
	} 

	private int sn11=0, sn1_=0, sn_1=0, sn=0;
	private double sprob=0;
	
	private double hyper0(int n11i, int n1_i, int n_1i, int ni)	{
		if (!((n1_i|n_1i|ni)!=0)) {
		    if (!(n11i % 10 == 0)) {
			    if (n11i==sn11+1) {
				    sprob *= ((sn1_-sn11)/(double)n11i)*((sn_1-sn11)/(double)(n11i+sn-sn1_-sn_1));  // double conv
				    sn11 = n11i;
				    return sprob;
			    }
			    if (n11i==sn11-1) {
				    sprob *= ((sn11)/(double)(sn1_-n11i))*((sn11+sn-sn1_-sn_1)/(double)(sn_1-n11i)); // double conv
				    sn11 = n11i;
				    return sprob;
			    }
		    }
		    sn11 = n11i;
		} else {
		    sn11 = n11i;
		    sn1_=n1_i;
		    sn_1=n_1i;
		    sn=ni;
		}
		sprob = hyper_323(sn11,sn1_,sn_1,sn);
		return sprob;
	}

	private double hyper(int n11) {
		return hyper0(n11,0,0,0);
	}

	private double sleft,sright,sless,slarg;
	
	private double exact(int n11, int n1_, int n_1, int n) {
		double p,prob;
		int max=n1_;
		if (n_1<max) 
			max=n_1;
		int min = n1_+n_1-n;
		if (min<0) 
			min=0;
		if (min==max) {
		    sless = 1;
		    sright= 1;
		    sleft = 1;
		    slarg = 1;
		    return 1;
		}
		prob=hyper0(n11,n1_,n_1,n);
		sleft=0;
		p=hyper(min);
		int i;
		for (i=min+1; p<0.99999999*prob; i++) {
		    sleft += p;
		    p=hyper(i);
                    if (i > 1000000)
                        break;
		}
		i--;
		if (p<1.00000001*prob) sleft += p;
		else i--;
		sright=0;
		p=hyper(max);
		int j;
		for (j=max-1; p<0.99999999*prob; j--) {
		    sright += p;
		    p=hyper(j);
		}
		j++;
		if (p<1.00000001*prob) sright += p;
		else j++;
		if (Math.abs(i-n11)<Math.abs(j-n11)) {
		    sless = sleft;
		    slarg = 1 - sleft + prob;
		} else {
		    sless = 1 - sright + prob;
		    slarg = sright;
		}
		return prob;
	}

	
	private final double prob;
	
	public FisherExactTest(int[][] matrix) {
		int row1 = matrix[0][0]+matrix[0][1]; // the total number of samples for the null hypothesis (row 1)
		int col1 = matrix[0][0]+matrix[1][0]; // the total number of samples for outcome 1 (col 1)
		int all  = matrix[0][0]+matrix[0][1]+matrix[1][0]+matrix[1][1]; // the total number of samples
		prob=exact(matrix[0][0], row1, col1, all);
	}

	public double getProb() {
		return prob;
	}
	
	public double getPValue(boolean left) {
		if (left)
			return sless;
		else
			return slarg;
	}

	public double get2Tail() {
		return Math.min(1.0, sleft+sright);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("Started "+new java.util.Date());
		FisherExactTest ft=new FisherExactTest(new int[][] {{838, 159}, {78, 29}});
		
		System.out.println("Prob="+ft.getProb());
		System.out.println("p-value="+ft.getPValue(true)+" : "+ft.getPValue(false));
		System.out.println("two-tail="+ft.get2Tail());
		
		System.out.println("Done "+new java.util.Date());

	}

}
