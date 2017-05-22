package bn.prob;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Random;

import dat.Enumerable;
import bn.Distrib;

import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * Mixture Dirichlet distribution is weighted sum of different Dirichlet
 * distribution MDir = w1 * Dir1 + w2 * Dir2 + ... This class includes learning
 * parameters from data However, if you would like to use Dirichlet distribution
 * as prior in BN, please look at the method in MixDirichletPrior.java
 *
 * @author wangyufei
 * @author mikael
 * 
*/
public class MixDirichletDistrib extends MixtureDistrib implements Serializable {

    /**
     * The max number of rounds the Gibbs sampler would run
     */
    private final int ROUND_LIMITATION = 200;
    /**
     * The max number of rounds with no update for DL value the algorithm would
     * allow
     */
    private final int NO_UPDATE = 20;
    private Enumerable domain;
    private double components;
    private double letters;
    private double DL_best;
    private EnumDistrib m_best;

    /**
     * Construct a mixture Dirichlet model from a single component
     *
     * @param d1
     * @param weight1
     */
    public MixDirichletDistrib(DirichletDistrib d1, double weight1) {
        super(d1, weight1);
        domain = (Enumerable) d1.getDomain();
    }

    /**
     * Given the domain of Dirichlet distribution build an empty Mixture model
     * Initialises distributions using random points from the data
     *
     * @param domain
     * @param ComponentNum
     */
    public MixDirichletDistrib(Enumerable domain, int ComponentNum, int[][] data) {
        super();
        Random rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < ComponentNum; i++) {
//        	EnumDistrib dataDistrib = new EnumDistrib(domain, data[rand.nextInt(data.length)]);
//        	EnumDistrib randomDistrib = EnumDistrib.random(domain, rand.nextInt());
        	int[] dataPoint = data[rand.nextInt(data.length)]; //Cannot contain any 0's
        	boolean zeros = true;
        	//Ensure dataPoint does not contain 0's
        	while (zeros) {
        		for (int d : dataPoint) {
        			if (d == 0) {
        				dataPoint = data[rand.nextInt(data.length)];
        				continue;
        			}
        		}
        		zeros = false;
        	}
//        	System.out.println(Arrays.toString(dataPoint));
            super.addDistrib(new DirichletDistrib(new EnumDistrib(domain, dataPoint), rand.nextInt(90) + 10), rand.nextDouble());
//            super.addDistrib(new DirichletDistrib(new EnumDistrib(domain, data[rand.nextInt(data.length)]), rand.nextInt(90) + 10), rand.nextDouble());
        }
        this.domain = domain;
        this.components = ComponentNum;
        this.letters = (double)domain.size();
    }
    
    /**
     * Given the domain of Dirichlet distribution build an empty Mixture model
     * Initialises distributions randomly
     *
     * @param domain
     * @param ComponentNum
     */
    public MixDirichletDistrib(Enumerable domain, int ComponentNum) {
        super();
//        Random rand = new Random(System.currentTimeMillis());
//        Random rand = new Random((long)1000); // Mikael: I removed this line to introduce the use of MixtureDistrib's own random number generator. No need for a local one.
        super.setSeed(System.currentTimeMillis());
        for (int i = 0; i < ComponentNum; i++) {
            super.addDistrib(new DirichletDistrib(EnumDistrib.random(domain, nextInt(Integer.MAX_VALUE)), nextInt(90) + 10), nextDouble());
        }
        this.domain = domain;
        this.components = ComponentNum;
        this.letters = (double)domain.size();
    }

    /**
     * Add either a Dirichlet distribution or mixDirichlet distribution
     */
    public double addDistrib(Distrib d2, double weight2) {

        if (d2 instanceof DirichletDistrib) {
            DirichletDistrib dir = (DirichletDistrib) d2;
            if (!getDomain().equals(dir.getDomain())) {
                throw new RuntimeException("Domain should be the same");
            }
            return super.addDistrib(d2, weight2);
        }

        if (d2 instanceof MixDirichletDistrib) {
            MixDirichletDistrib mixDir = (MixDirichletDistrib) d2;
            if (!getDomain().equals(mixDir.getDomain())) {
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
     * Determine the most likely label (index of distribution) to explain the data.
     * 
     * @param data
     * @return index to Dirichlet distribution
     */
    public int getLabel(int[] data) {
        double[] logprob = new double[distribs.size()];
        double[] prob = new double[distribs.size()];
        double[] m = this.getAllWeights();
        int best = 0;
        for (int i = 0; i < logprob.length; i++) {
            logprob[i] = ((DirichletDistrib) distribs.get(i)).logLikelihood(data);
            double p_i = Math.exp(logprob[i]);
            logprob[i] = (Math.log(m[i]) + logprob[i]);
            prob[i] = (m[i] * p_i);
            if (i > 0 && logprob[i] > logprob[best])
                best = i;
        }
        return best;
    }
    
    /**
     * The description length used by Ye et al (2011) to score a mixture of
     * Dirichlets.
     *
     * @param data count histograms
     * @return the description length
     */
    private double DL(int[][] data) {
        double outer = 0;
        double[] m = this.getAllWeights();
        for (int k = 0; k < data.length; k++) {
            double inner = Double.MIN_VALUE;
            for (int i = 0; i < distribs.size(); i++) {
                double log_p_i = ((DirichletDistrib) distribs.get(i)).logLikelihood(data[k]);
                double p_i = Math.exp(log_p_i);
                inner += (m[i] * p_i);
            }
            outer += Math.log(inner);
        }
        return -outer;
    }

    /**
     * Learning the parameters, including mixture coefficient and parameters for
     * each component.
     *
     * @param data training data, in this case, this would be a count vector,
     * the weight for each training point is assumed to be 1
     * @return the description length (smaller means better)
     */
    public double learnParameters(int[][] data) {

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

        for (int i = 0; i < nbins; i++) {
            DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
            System.arraycopy(dirichlet.getAlpha(), 0, alpha_best[i], 0, nseg);
        }
        EnumDistrib p = EnumDistrib.random(new Enumerable(nbins), rand.nextInt()); // probability that sample belongs to bin
        p.setSeed(rand.nextInt());

        int no_update = 0;
        int reintCount = 0;
        boolean anyEmpty = false;
        boolean training = false;
        // start iteration
        likeLoop:
        for (int round = 0; round < ROUND_LIMITATION && no_update < NO_UPDATE; round++) {
            // start with empty bins
            for (int i = 0; i < nbins; i++) {
                bins[i] = new ArrayList();
            }

            // try to put each data points into different bins
            //this.getNormalized();
            double[] m = this.getAllWeights();
            Map<int[], Double> trackPoints = new HashMap();
            for (int k = 0; k < dataSize; k++) {
                try {
                    double[] logprob = new double[nbins];
                    double maxLog = -1000000;
                    for (int i = 0; i < nbins; i++) {
                        DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
                        double log = Math.log(m[i]) + dirichlet.logLikelihood(data[k]);
                        logprob[i] = log;
                        if (log > maxLog)
                        	maxLog = log;
                    }
                    //FIXME some of them would be zero......
                    p.set(EnumDistrib.log2Prob(logprob));
                    Integer index = (Integer) p.sample();
                    bins[index].add(data[k]);
                    //Store each data point and it's lowest log likelihood
                    //FIXME more efficient way to track 'worst' data points
                    //Only need to store nbins-1 data points (worst case scenario)
                    //Could do this with a map and add/remove data points as necessary
                    //Efficiency? Memory?
                    trackPoints.put(data[k], maxLog);
                } catch (RuntimeException ex0) {
                    System.err.println("Problem with data point k = " + k);
                    throw new RuntimeException("Problem with data point k = " + k);
                }
            }
            
            /*
             * Test whether any bins are 'empty' - i.e. data does not comply with what distribution can generate?
             * If there are empty bins -take the data point with the highest log likelihood (lowest probability of 
             * belonging to cluster) and place it in empty bin
             */
            for (int i = 0; i < nbins; i++) {
                if (bins[i].size() == 0) {
                    System.err.println("Empty Bin : " + i);
                    double maxValueInMap=(Collections.max(trackPoints.values()));  // This will return max value in the Hashmap
                    for (Entry<int[], Double> entry : trackPoints.entrySet()) {
                        if (entry.getValue()==maxValueInMap) {
                        	//FIXME better way of removing data point from bin?
                        	for (int k = 0; k < nbins; k++) {
                        		for (int l = 0; l < bins[k].size(); l++) {
                        			if (bins[k].get(l) == entry.getKey()) {
                        				bins[k].remove(l);
                        			}
                        		}
                        	}
                            bins[i].add(entry.getKey());
                            trackPoints.remove(entry.getKey());
                            break;
                        }
                    }
                    anyEmpty = true;
                }
            }
           
            // Report if there were any bins shuffled
            if (anyEmpty) {
            	System.out.println("Empty bins were re-initialised");
            	anyEmpty = false;
            	training = true;
            	round = 0;
            } else {
            	training = true;
            }
            
            // based on the data in each bin, adjust parameters of each dirichlet distribution
//          System.out.println("Updating weights");
            for (int i = 0; i < nbins; i++) {
                // update mixture coefficient
                this.setWeight(i, bins[i].size());
                // update parameters for model
                int[][] counts = new int[bins[i].size()][];
                for (int j = 0; j < bins[i].size(); j++) {
                    counts[j] = (int[]) bins[i].get(j);
                }
                DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
                
                try {
                	dirichlet.setPrior(DirichletDistrib.getAlpha(counts));
                } catch (NullPointerException e) {
                	System.out.println(i);
//                	break mainloop;
                	break likeLoop;
                }
                
//                dirichlet.setPrior(DirichletDistrib.getAlpha(counts));
            }
            this.getNormalized();

            double dl_cur = DL(data);
            System.out.println("DL_cur = " + dl_cur + "\tDL_best = " + dl_best);
            if (dl_cur < dl_best) {
            	setDLBest(dl_cur);
                dl_best = dl_cur;
                m_best.set(this.getAllWeights());
                // also save mixing weights and alpha values1
                for (int i = 0; i < nbins; i++) {
                    DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
                    System.arraycopy(dirichlet.getAlpha(), 0, alpha_best[i], 0, nseg);
                }
                no_update = 0;
            } else {
                no_update++;
            }
        }

        // set the best data back
        this.setWeights(m_best.get());
        
        for (int i = 0; i < nbins; i++) {
            DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
            dirichlet.setPrior(alpha_best[i]);
        }
        
        if (!training) {
        	System.err.print("Failed to remove empty bins\n");
        }
        return dl_best;
    }
    
    /**
     * Learning the parameters, including mixture coefficient and parameters for
     * each component. In the case of ChIP peak data sets, check both the forward
     * and reverse of each vector
     *
     * @param data training data, in this case, this would be a count vector,
     * the weight for each training point is assumed to be 1
     */
    public void learnParametersFlip(int[][] data) {

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

        for (int i = 0; i < nbins; i++) {
            DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
            System.arraycopy(dirichlet.getAlpha(), 0, alpha_best[i], 0, nseg);
        }
        EnumDistrib p = EnumDistrib.random(new Enumerable(nbins), rand.nextInt()); // probability that sample belongs to bin
        p.setSeed(rand.nextInt());

        int no_update = 0;
        int reintCount = 0;
        boolean anyEmpty = false;
        boolean training = false;
        // start iteration
        likeLoop:
        for (int round = 0; round < ROUND_LIMITATION && no_update < NO_UPDATE; round++) {
            // start with empty bins
            for (int i = 0; i < nbins; i++) {
                bins[i] = new ArrayList();
            }

            // try to put each data points into different bins
            //this.getNormalized();
            double[] m = this.getAllWeights();
            Map<int[], Double> trackPoints = new HashMap();
            for (int k = 0; k < dataSize; k++) {
            	//get the reverse count vector
            	int[] reverse = new int[data[k].length];
            	int[] forward = data[k];
            	int stop = data[k].length /2; //assume list is even
            	if (data[k].length % 2 == 1) {
            		stop = (data[k].length/2) +1;
            	}
                for(int r = 0; r < stop; r++) {
                    int temp = data[k][r];
                    reverse[r] = data[k][data[k].length - r - 1];
                    reverse[data[k].length - r - 1] = temp;
                }
                boolean rev = false; //assume forward until proven otherwise
                try {
                    double[] logprob = new double[nbins];
                    double maxLog = -1000000;
                    for (int i = 0; i < nbins; i++) {
                        DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
                        //Here need to check both 'forward' and 'reverse'
                        //Does it need to be tracked? Would it be beneficial to know?
                        //Could be done post hoc if reverse vector is recorded
                        double logForward = Math.log(m[i]) + dirichlet.logLikelihood(data[k]);
                        double logReverse = Math.log(m[i]) + dirichlet.logLikelihood(reverse);
//                        rev = false; //assume forward until proven otherwise
                        double log = logForward;
                        if (logForward < logReverse) {
                        	rev = true;
                        	log = logReverse;
//                        	System.out.println("Reversed");
                        } 
                        logprob[i] = log;
                        if (log > maxLog)
                        	maxLog = log;
                    }
                    p.set(EnumDistrib.log2Prob(logprob));
                    Integer index = (Integer) p.sample();
                    if (rev) {
                    	bins[index].add(reverse);
                    	data[k] = reverse; //data has to reflect current state 
                    					   //of all vectors or DL calculation is incorrect
                    } else {
                    	bins[index].add(data[k]);
                    }
                    
                    //Store each data point and it's lowest log likelihood
                    //FIXME more efficient way to track 'worst' data points?
                    //Only need to store nbins-1 data points (worst case scenario)
                    //Could do this with a map and add/remove data points as necessary
                    //Efficiency? Memory?
                    //What about uniqueness - are all data points unique?
                    //Tracking forward and reverse potentially doubles size of map
                    if (rev) {
                    	trackPoints.put(reverse, maxLog);
                    } else {
                    	trackPoints.put(data[k], maxLog);
                    }
                    
                } catch (RuntimeException ex0) {
                    System.err.println("Problem with data point k = " + k);
                    throw new RuntimeException("Problem with data point k = " + k);
                }
            }
            
            /*
             * Test whether any bins are 'empty' - i.e. data does not comply with what distribution can generate?
             * If there are empty bins -take the data point with the highest log likelihood (lowest probability of 
             * belonging to cluster) and place it in empty bin
             */
            for (int i = 0; i < nbins; i++) {
                if (bins[i].size() == 0) {
                    System.err.println("Empty Bin : " + i);
                    double maxValueInMap=(Collections.max(trackPoints.values()));  // This will return max value in the Hashmap
                    for (Entry<int[], Double> entry : trackPoints.entrySet()) {
                        if (entry.getValue()==maxValueInMap) {
                        	//FIXME better way of removing data point from bin?
                        	for (int k = 0; k < nbins; k++) {
                        		for (int l = 0; l < bins[k].size(); l++) {
                        			if (bins[k].get(l) == entry.getKey()) {
                        				bins[k].remove(l);
                        			}
                        		}
                        	}
                            bins[i].add(entry.getKey());
                            trackPoints.remove(entry.getKey());
                            break;
                        }
                    }
                    anyEmpty = true;
                }
            }
           
            // Report if there were any bins shuffled
            if (anyEmpty) {
            	System.out.println("Empty bins were re-initialised");
            	anyEmpty = false;
            	training = true;
            	round = 0;
            } else {
            	training = true;
            }
            
            // based on the data in each bin, adjust parameters of each dirichlet distribution
//          System.out.println("Updating weights");
            for (int i = 0; i < nbins; i++) {
                // update mixture coefficient
                this.setWeight(i, bins[i].size());
                // update parameters for model
                int[][] counts = new int[bins[i].size()][];
                for (int j = 0; j < bins[i].size(); j++) {
                    counts[j] = (int[]) bins[i].get(j);
                }
                DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
                
                try {
                	dirichlet.setPrior(DirichletDistrib.getAlpha(counts));
                } catch (NullPointerException e) {
                	System.out.println(i);
//                	break mainloop;
                	break likeLoop;
                }
                
//                dirichlet.setPrior(DirichletDistrib.getAlpha(counts));
            }
            this.getNormalized();

            double dl_cur = DL(data);
            System.out.println("DL_cur = " + dl_cur + "\tDL_best = " + dl_best);
            if (dl_cur < dl_best) {
            	setDLBest(dl_cur);
                dl_best = dl_cur;
                m_best.set(this.getAllWeights());
                // also save mixing weights and alpha values1
                for (int i = 0; i < nbins; i++) {
                    DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
                    System.arraycopy(dirichlet.getAlpha(), 0, alpha_best[i], 0, nseg);
                }
                no_update = 0;
            } else {
                no_update++;
            }
        }

        // set the best data back
        this.setWeights(m_best.get());
        
        for (int i = 0; i < nbins; i++) {
            DirichletDistrib dirichlet = (DirichletDistrib) this.getDistrib(i);
            dirichlet.setPrior(alpha_best[i]);
        }
        
        if (!training) {
        	System.err.print("Failed to remove empty bins\n");
        }
    }
    
    public void setDLBest(double dl_best) {
    	this.DL_best = dl_best;
    }
    
    public double getDLBest() {
    	return this.DL_best;
    }
    
    public void setMBest(EnumDistrib m_best) {
    	this.m_best = m_best;
    }
    
    public EnumDistrib getMBest() {
    	return this.m_best;
    }
    
    public static double getEntropy(double[] p) {
        double log_base = Math.log(p.length);
        double ent = 0;
        for (int i = 0; i < p.length; i ++) { 
            double p_i = p[i] + 0.0001;
            ent -= p_i * (Math.log(p_i) / log_base);
        }
        return ent;
    }
    
    /**
     * Approximate the complexity of a Dirichlet Mixture Model.
     * Uses equations 3,4&5 from Ye et al (2011)
     * 
     * @return complexity
     */
    public double getComplexity(int[][]data) {
     	//letters = L
    	//components = M
    	//data.length or nData = n
    	
    	//calculate the average counts across all data
    		//Possible alternative // c = average of (average per bin)
    	Double nData = (double) data.length;
    	Double cCounts = 0.0;
    	for (int d = 0; d < data.length; d++) {
    		for (int col = 0; col < data[d].length; col++) {
    			cCounts += data[d][col];
    		}
    	}
    	Double cAvg = cCounts/(data.length * letters);
    	
    	//delta L,c approaches 0 as c increases. It is a constant. At L = 20, c = 100 the value is 0.057 - a negligible impact
    	// o(1) is a constant with no explanation in either Altschul paper (Ye or Yu)
    	
    	//COMP(Dl, n/M, c), formula 3
    	Double comp3 = ((letters/2)*Math.log(nData/components)) + ((letters-1)/2)*Math.log(cAvg/2) - GammaDistrib.lgamma(letters/2) - ((1/2)*Math.log(letters-1)) ; //excluding 2 constants at end
    	
    	//COMP(Mm, n), formula 4
    	Double comp4 = (((components - 1)/2)*Math.log(nData/2)) + ((1/2)*Math.log(Math.PI)) - GammaDistrib.lgamma(components/2); // + o(1)
    	
    	//calculate M!
    	double mfact = 1.0;
    	for (double i = 1; i <=components; i++) {
    		mfact = mfact*i;
    	}
    	
    	//COMP(DM(m, L), n, c) ~= COMP(Mm, n) + M*COMP(Dl, M/n, c) - log*(M!), formula 5
    	Double complexity = comp4 + components*comp3 - Math.log(mfact);
//    	System.out.println(complexity);
    	
    	return complexity;
    }
    
    public void saveAlphas(String filename) {
    	try {
        	BufferedWriter bd = new BufferedWriter(new FileWriter(filename + "_alpha_"+(int)this.components+".out"));
	        for (int i = 0; i < this.components; i ++) {
	        	DirichletDistrib dd = (DirichletDistrib) this.getDistrib(i);
//	            System.out.println("D"+ i + " = " + dd);
	            bd.write("D"+i+","+dd);
	            bd.newLine();
	        }
	        bd.close();
	    } catch (IOException ex) {
            Logger.getLogger(DirichletDistrib.class.getName()).log(Level.SEVERE, null, ex);
        }   
    }
    
    public void saveClusters(int[][] data, String filename) {
    	
    	int N = data.length;
    	int nbins = (int)this.components;
        List[] bins = new ArrayList[nbins]; // hold components here
        for (int i = 0; i < nbins; i++)
          bins[i] = new ArrayList(); // start with empty bins
        for (int k = 0; k < N; k++) { //N = data points
        	int index = getLabel(data[k]);
        	bins[index].add(k);
        }

        /*
        for (int i = 0; i < dds.length; i ++) {
            System.out.println("D"+ i + " = " + dds[i]);
        }*/

        for (int i = 0; i < nbins; i ++) {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(filename + "_bin_" + i + "_"+nbins+".out"));
                for (int a = 0; a < bins[i].size(); a ++) {
                    int k = (Integer)bins[i].get(a);
                    bw.write(""+k);
//                    bw.write(k + "\t");
//                    for (int j = 0; j < data[k].length; j ++)
//                        bw.write(data[k][j] + "\t");
//                    bw.write(k + "\t");
//                    EnumDistrib d = new EnumDistrib(new Enumerable((int)this.letters), data[k]);
//                    for (int j = 0; j < d.getDomain().size(); j ++)
//                        bw.write(d.get(j) + "\t");
                    bw.newLine();
                }
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(DirichletDistrib.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
       
    public static void main(String[] args) {
    
    	int min = Integer.parseInt(args[1]);
        int max = Integer.parseInt(args[2]);
//        int alphaInit = Integer.parseInt(args[3]);
        
//        String filename = args[0];
//        String filename = "cage_all_expression.out";
        String filename = "wgEncodeH1hescSrf_seg20_500_srf_hg19.out";
//        String filename = "wgEncodeGm12878Max_seg20_500_max_gm_hg19.out";
      
        int[][] data = loadData(filename);
        int N = data.length;
        int nseg = 0;
        for (int i = 0; i < N; i ++) {
            if (nseg == 0) 
                nseg = data[i].length;
            else if (nseg != data[i].length)
                throw new RuntimeException("Error in data: invalid item at data point " + (i + 1));
//            System.out.print("[" + i + "]\t= \t");
//            for (int j = 0; j < data[i].length; j ++) 
//                System.out.print(data[i][j] + "\t");
        }
        
        Enumerable domain = new Enumerable(data[0].length);
        
        double[] dlBestList = new double[max-min];
        double[] complexityList = new double[max-min];
        
        for ( int nbins = min; nbins < max; nbins++) {
        	System.out.println("nbins = " +nbins);
        	MixDirichletDistrib dis = new MixDirichletDistrib(domain, nbins, data);
//        	MixDirichletDistrib dis = new MixDirichletDistrib(domain, nbins);
            dis.learnParametersFlip(data); //For ChIP-seq peak data that has unknown orientation
//            dis.learnParameters(data);
            dis.saveAlphas(filename);
	        dlBestList[nbins - min] = dis.getDLBest();
	        complexityList[nbins - min] = dis.getComplexity(data);
	        //FIXME Track changes here - automate isolation of optimal cluster group
	        dis.saveClusters(data, filename);
        }
        
        //finalComplexity
        double[] fComp = new double[max - min];
        for (int c = 0; c < dlBestList.length; c++ ) {
        	fComp[c] = dlBestList[c] + complexityList[c];
        }
        
        //write list of dl_best scores
//        try {
//        	BufferedWriter bd = new BufferedWriter(new FileWriter(filename + "_ll.out"));
//        	for (double dl : dlBestList) {
//        		bd.write(min + " = " +dl);
//        		min++;
//        	}
//        	bd.close();
//        } catch (IOException ex) {
//        	Logger.getLogger(DirichletDistrib.class.getName()).log(Level.SEVERE, null, ex);
//        } 
        
        //write list of complexity scores
        try {
        	BufferedWriter bd = new BufferedWriter(new FileWriter(filename + "_comp.out"));
        	for (double fc : fComp) {
        		bd.write(min + "\t" +fc);
        		bd.newLine();
        		min++;
        	}
        	bd.close();
        } catch (IOException ex) {
        	Logger.getLogger(DirichletDistrib.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("COMPLETE");
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
                    for (int i = 0; i < tokens.length; i++) {
                        values[i] = Integer.valueOf(tokens[i]);
                    }
                    alldata.add(values);
                } catch (NumberFormatException ex2) {
                    System.err.println("Ignored: " + line);
                }
                line = br.readLine();
            }
            data = new int[alldata.size()][];
            for (int k = 0; k < data.length; k++) {
                data[k] = new int[alldata.get(k).length];
                for (int j = 0; j < data[k].length; j++) {
                    data[k][j] = alldata.get(k)[j];
                }
            }
            br.close();
        } catch (IOException ex) {
            Logger.getLogger(DirichletDistrib.class.getName()).log(Level.SEVERE, null, ex);
        }
        return data;
    }
    
}


