package stats;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The Mann-Whitney U-test (also known as the Mann-Whitney-Wilcoxon test) or
 * the Wilcoxon rank sum test. It is not the signed rank sum.
 * For the calculation of the p-value, an approximation is used. 
 * The approximation is good in cases where the sample sizes are 5 or above (each) and
 * relies on the fact that larger samples will resemble the normal distribution.
 * Alternatively, the test computes the so-called U-value which can be checked against
 * tables and an exact p-value can be determined.
 * @author m.boden
 */
public class MannWhitney {

	private double p;
	private double u;
	private double mean_a;
	private double mean_b;
	
	/** Percent of top sorted data that should be included, set to 0 if not used (default) */  
	public static int SUBLIST = 0;  
	
	private class Data implements Comparable<Data> {
		double measurement;
		boolean group;
		double rank=-1;
		public Data(double measurement, boolean group) {
			this.measurement=measurement;
			this.group=group;
		}
		// sorting is done in descending order (larger values first)
		public int compareTo(Data other) {
			if (this.measurement<other.measurement)
				return +1;
			else if (this.measurement>other.measurement)
				return -1;
			else
				return 0;
		}
	}
	
	/**
	 * Perform a Mann-Whitney test, using two sample-sets.
	 * @param samples the measurements for the two sets
	 * @param group the set membership (true/false) for each measurement
	 */
	public MannWhitney(double[] samples, boolean[] group) {
		if (samples.length!=group.length)
			throw new RuntimeException("Invalid specification");
		ArrayList<Data> all=new ArrayList<Data>(samples.length);
		for (int i=0; i<samples.length; i++) {
			all.add(new Data(samples[i], group[i]));
		}
		run(all);
	}
	
	/**
	 * Perform a Mann-Whitney test, using two sample-sets.
	 * @param sampleA set A
	 * @param sampleB set B
	 */
	public MannWhitney(double[] sampleA, double[] sampleB) {
		ArrayList<Data> all=new ArrayList<Data>(sampleA.length+sampleB.length);
		for (int i=0; i<sampleA.length; i++)
			all.add(new Data(sampleA[i], true));
		for (int i=0; i<sampleB.length; i++)
			all.add(new Data(sampleB[i], false));
		run(all);
	}
	
	private void run(List<Data> all) {
		int rank=0;
		int na=0, nb=0; // the number of points in each set (A & B)
		double sum_a=0, sum_b=0;
		// first sort according to measurement 
		Collections.sort(all);
		if (SUBLIST > 0)
			all = all.subList(0, (int)(all.size()*(SUBLIST/100.0)));
		ArrayList<Data> same=new ArrayList<Data>(); // all entries with the same "measurement"
		double measurement=all.get(0).measurement;
		// assign rank in accordance with sorting order
		for (Data entry:all) {
			if (entry.group) { // positive=true=left
				na++;
				sum_a+=measurement;
			} else { 
				nb++;
				sum_b+=measurement;
			}
			if (measurement!=entry.measurement) { // here's an entry that differed from the previous...
				int firstInGroup=rank+1-same.size(); 
				int lastInGroup=rank; 
				double average=(lastInGroup-firstInGroup)/(double)(2.0); 
				// assign the average of all entries with the same measurement
				for (Data s:same) {
					s.rank=firstInGroup+average;
				}
				same=new ArrayList<Data>();
				measurement=entry.measurement;
			} 
			same.add(entry);
			rank++;
		}
		// the last batch of entries is handled outside the loop...
		int firstInGroup=rank+1-same.size(); 
		int lastInGroup=rank; 
		double average=(lastInGroup-firstInGroup)/(double)(2.0); 
		for (Data s:same) {
			s.rank=firstInGroup+average;
		}
		
		int n=na+nb;		// the total number of measurements
		mean_a=sum_a/na;
		mean_b=sum_b/nb;
		double ta_obs=0; 	// sum of na ranks in group A
		double tb_obs=0; 	// sum of nb ranks in group B
		
		// sum the ranks (replace the measurements)
		for (Data entry:all) {
			if (entry.group)
				ta_obs+=(double)entry.rank;
			else
				tb_obs+=(double)entry.rank;
			//System.out.println(entry.rank+"\t"+entry.measurement+"\t"+(entry.group?"A":"B"));
		}
		
		double tab=ta_obs+tb_obs; 					// sum of n ranks in groups A and B combined 
		double sd=Math.sqrt((na*nb*(n+1.0))/12.0); 	// the standard deviation is the same in both sets
		
		double ta_null=na*(n+1.0)/2.0;				// the sum of the "null" case
		double tb_null=nb*(n+1.0)/2.0;				// the sum of the "null" case
		double ta_max=na*nb+(na*(na+1.0))/2.0;		// the max sum set A can take 
		double tb_max=na*nb+(nb*(nb+1.0))/2.0;		// the max sum set B can take
		double ua=ta_max-ta_obs;					// the "U" value for A which is the mirror of ...
		double ub=tb_max-tb_obs;					// the "U" value for B (we only need one)
		double ua_null=ta_max-ta_null;				// the U value for the null case
		double ub_null=tb_max-tb_null;
		double da=ta_obs>ta_null?-0.5:+0.5;			// a "continuity correction" for A
		double db=tb_obs>tb_null?-0.5:+0.5;			// a "continuity correction" for B
		double za=((ta_obs-ta_null)+da)/sd;			// the z value for A which is the mirror of ...
		double zb=((tb_obs-tb_null)+db)/sd;			// the z value for B (we only need one)
		p=stats.NormalDistribution.f(za);					// figure out the area of the normal distribution
		u=ua;										// remember one of the U values
	}
	
	/**
	 * The Mann-Whitney test computes the U-value, a measure which is equal to the difference between the
	 * maximum possible sum of ranks versus the observed sum.
	 * Use tables to determine p-value from the U-value.
	 * @return the U-value
	 */
	public double getU() {
		return u;
	}
	
	/**
	 * Retrieve the computed p-value or rather an approximation of it.
	 * The approximation is based on the normal distribution and is reliable
	 * when sample sets are of size 5 or larger. 
	 * @param left if true, the area of the left side of the Gaussian, relative the 
	 * estimated z-value, is used
	 * @return the p-value
	 */
	public double getPValue(boolean left) {
		if (left)
			return p;
		else
			return 1-p;
	}

	/**
	 * Retrieve the computed p-value (an approximation of it).
	 * The approximation is based on the normal distribution and is reliable
	 * when sample sets are of size 5 or larger. 
	 * @return the smallest one-tailed p-value
	 */
	public double getPValue() {
		return Math.min(p, 1-p);
	}
	
	/**
	 * Retrieve the two-tailed p-value (which equals double the one-tailed p-value).
	 * @return the two-tailed p-value
	 */
	public double get2Tail() {
		return getPValue()*2;
	}

	public double getMean(boolean left) {
		return (left?mean_a:mean_b);
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MannWhitney mww=new MannWhitney(new double[] {4.6,5.1,5.8,6.5,4.7,5.2,6.1,7.2,4.9,5.5,6.5}, new double[] {5.2,5.6,6.8,8.1,5.3,6.2,7.7,5.4,6.3,8.0});
    	//MannWhitney mww=new MannWhitney(new double[] {0.01,0.02,0.06,1.0,0.05,0.12, 0.13, 0.13, 0.15, 0.051,0.05,0.011}, new double[] {0.7,0.5,0.66, 0.44, 0.13, 0.22, 0.13, 0.20 ,0.10, 0.51,0.52,0.53,0.012});
    	System.out.println("U="+mww.u);
    	System.out.println("Bothtails="+mww.get2Tail()+" Left="+mww.getPValue(true)+" Right="+mww.getPValue(false));
    	System.out.println("Mean Left="+mww.getMean(true)+" Right="+mww.getMean(false));
	}

}
