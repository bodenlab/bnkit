package bn.prob;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Random;

import dat.Enumerable;

import bn.Distrib;


/**
* Mixture Dirichlet distribution is weighted sum of different Dirichlet distribution
* MDir = w1 * Dir1 + w2 * Dir2 + ...
* This class includes learning parameters from data
* However, if you would like to use Dirichlet distribution as prior in BN,
* please look at the method in MixDirichletPrior.java
* @author wangyufei
*
*/

public class MixDirichletDistrib  extends MixtureDistrib implements Serializable{

	/**
	 * The max round the gibbs sampling would run
	 */
	private final int ROUND_LIMITATION = 100;
	/**
	 * The max number of rounds with no update for DL value the algorithm would allow
	 */
	private final int NO_UPDATE = 10;
	private Enumerable domain;

	/**
	 * construct a mixture Dirichlet model from a single component
	 * @param d1
	 * @param weight1
	 */
	public MixDirichletDistrib(DirichletDistrib d1, double weight1) {
		super(d1, weight1);
		domain = (Enumerable) d1.getDomain();
	}

	/**
	 * given the domain of Dirichlet distribution
	 * build an empty Mixture model
	 * @param domain
	 * @param ComponentNum
	 */
	public MixDirichletDistrib(Enumerable domain, int ComponentNum) {
		super();
		Random rand = new Random(System.currentTimeMillis());
		for(int i = 0; i < ComponentNum; i++) {
			super.addDistrib(new DirichletDistrib(EnumDistrib.random(domain,rand.nextInt()), 10), rand.nextDouble());
		}
		this.domain = domain;
	}

	/**
	 * add either a Dirichlet distribution or mixDirichlet distribution
	 */
	public double addDistrib(Distrib d2, double weight2) {

		if(d2 instanceof DirichletDistrib) {
			DirichletDistrib dir = (DirichletDistrib)d2;
			if(!getDomain().equals(dir.getDomain())) {
				throw new RuntimeException("Domain should be the same");
			}
			return super.addDistrib(d2, weight2);
		}

		if(d2 instanceof MixDirichletDistrib) {
			MixDirichletDistrib mixDir = (MixDirichletDistrib)d2;
			if(!getDomain().equals(mixDir.getDomain())) {
				throw new RuntimeException("Domain should be the same");
			}
			return super.addDistrib(d2, weight2);
		}

		throw new RuntimeException("only accept DirichletDistrib or MixDirichletDistrib");

	}

	public Enumerable getDomain() {
		return domain;
	}

	/**
     * The description length used by Ye et al (2011) to score a mixture of Dirichlets.
     * @param data count histograms
     * @return the description length
     */
    private double DL(int[][] data) {
        double outer = 0;
        for (int k = 0; k < data.length; k ++) {
            double inner = Double.MIN_VALUE;
            for (int i = 0; i < distribs.size(); i ++) {
                double log_p_i = ((DirichletDistrib)distribs.get(i)).logLikelihood(data[k]);
                double p_i = Math.exp(log_p_i);
                inner += (this.getWeights(i) * p_i);
            }
            outer += Math.log(inner);
        }
        return -outer;
    }

	/**
	 * learning the parameters, including mixture coefficient and parameters for each component
	 * @param data training data, in this case, this would be a count vector, the weight for all training
	 * point is assumed to be 1
	 */
	public void learnParameters(int[][] data) {

		// necessary parameters
		int nseg = data[0].length;
		int nbins = this.getMixtureSize();
		int dataSize = data.length;
		// init dl value and alpha best
		double dl_best = DL(data);
		double[][] alpha_best = new double[nbins][nseg];
		EnumDistrib m_best = new EnumDistrib(new Enumerable(nbins), this.getAllWeights());
		ArrayList[] bins = new ArrayList[nbins];
		this.getNormalized();

        for (int i = 0; i < nbins; i ++) {
        	DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i);
            System.arraycopy(dirichlet.getAlpha(), 0, alpha_best[i], 0, nseg);
        }
        EnumDistrib p = EnumDistrib.random(new Enumerable(nbins), rand.nextInt()); // probability that sample belongs to bin
        p.setSeed(rand.nextInt());

        int no_update = 0;
        // start iteration
        for (int round = 0; round < ROUND_LIMITATION && no_update < NO_UPDATE; round ++) {
        	// start with empty bins
        	for (int i = 0; i < nbins; i ++) {
                bins[i] = new ArrayList();
        	}

        	// try to put each data points into different bins
        	for(int k = 0; k < dataSize; k++) {
        		try {
                    double[] logprob = new double[nbins];
                    for (int i = 0; i < nbins; i ++) {
                    	DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i);
                        logprob[i] += (Math.log(this.getWeights(i)) + dirichlet.logLikelihood(data[k]));
                    }
                    p.set(EnumDistrib.log2Prob(logprob));
                    Integer index = (Integer)p.sample();
                    bins[index].add(data[k]);
                } catch (RuntimeException ex0) {
                    System.err.println("Problem with data point k = " + k);
                }
        	}
            /*
            Divide the points in the largest bin and place half in the empty bin... (yes, we assume only 1 empty bin)
            What about if there are more than 2 empty bins? We could divide the largest bin multiple times
            or recursively go through the bins and pick the largest and divide until all bins have values?
            This may require further work.
             */
            int largestIndex = 0;
            int largestSize = 0;
            int emptyIndex = -1;
            for (int i = 0; i < nbins; i++){
                if (bins[i].size()> largestSize){
                    largestSize = bins[i].size();
                    largestIndex = i;
                }
                if (bins[i].size() == 0) {
//                    System.out.println("Empty Bin : " + i);
                    emptyIndex = i;
                }
            }
            if (emptyIndex != -1){
                Collections.shuffle(bins[largestIndex]);
                int pointsToSelect = bins[largestIndex].size()/2;
                for (int i = 0; i <pointsToSelect; i++){
                    Object countArray = bins[largestIndex].get(i);
                    bins[emptyIndex].add(countArray);
                    bins[largestIndex].remove(i);
                }
            }

        	// based on the data in each bin, adjust parameters of each dirichlet distribution
        	for(int i = 0; i < nbins; i++) {
        		// update mixture coeffience
        		this.setWeight(i, bins[i].size());
        		// update parameters for model
        		int[][] counts = new int[bins[i].size()][];
        		for(int j = 0; j < bins[i].size(); j++) {
        			counts[j] = (int[])bins[i].get(j);
        		}
        		DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i);
        		dirichlet.setPrior(DirichletDistrib.getAlpha(counts));
        	}
        	this.getNormalized();

        	double dl_cur = DL(data);
//            System.out.println("DL_cur = " + dl_cur + "\tDL_best = " + dl_best);
            if (dl_cur < dl_best) {
                dl_best = dl_cur;
                m_best.set(this.getAllWeights());
                // also save mixing weights and alpha values
                for (int i = 0; i < nbins; i ++) {
                	DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i);
                    System.arraycopy(dirichlet.getAlpha(), 0, alpha_best[i], 0, nseg);
                }
                no_update = 0;
        	} else {
        		no_update ++;
        	}

        }

        // set the best data back
        this.setWeights(m_best.get());
        for(int i = 0; i < nbins; i++) {
        	DirichletDistrib dirichlet = (DirichletDistrib)this.getDistrib(i);
        	dirichlet.setPrior(alpha_best[i]);
        }
	}

	public static void main(String[] args) {

		Enumerable domain = new Enumerable(5);
		/*
		MixDirichletDistrib dis = new MixDirichletDistrib(domain, 9);
		System.out.println(dis.toString());
		int dataNum = 23000;
		int[][] data = new int[dataNum][5];
		for(int i = 0; i < dataNum; i++) {
			EnumDistrib enumDistrib = (EnumDistrib)dis.sample();
			for(int j = 0; j < 5; j++) {
				data[i][j] = (int) (enumDistrib.get(j) * 30);
			}
		}*/

		MixDirichletDistrib dis2 = new MixDirichletDistrib(domain, 9);
		dis2.getNormalized();
		System.out.println(dis2.toXMLString());
		//dis2.learnParameters(loadData("data/mm10_Mixed_NfiX_segmented20_100.out"));
		//System.out.println(dis2.toString());
	}

	public static int[][] loadData(String filename) {
        BufferedReader br = null;
        int[][] data = null;
        try {
            br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();
            List<int[]> alldata = new ArrayList<>();
            while (line != null) {
                String[] tokens = line.split("\t");
                int[] values = new int[tokens.length];
                try {
                    for (int i = 0; i < tokens.length; i ++) {
                        values[i] = Integer.valueOf(tokens[i]);
                    }
                    alldata.add(values);
                } catch (NumberFormatException ex2) {
                    System.err.println("Ignored: " + line);
                }
                line = br.readLine();
            }
            data = new int[alldata.size()][];
            for (int k = 0; k < data.length; k ++) {
                data[k] = new int[alldata.get(k).length];
                for (int j = 0; j < data[k].length; j ++)
                    data[k][j] = alldata.get(k)[j];
            }
            br.close();
        } catch (IOException ex) {
            Logger.getLogger(DirichletDistrib.class.getName()).log(Level.SEVERE, null, ex);
        }
        return data;
    }

}
