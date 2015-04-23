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

import org.apache.commons.lang3.ArrayUtils;

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
    private final int ROUND_LIMITATION = 100;
    /**
     * The max number of rounds with no update for DL value the algorithm would
     * allow
     */
    private final int NO_UPDATE = 10;
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
        	while (ArrayUtils.contains(dataPoint, 0)) {
        		dataPoint = data[rand.nextInt(data.length)];
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
        Random rand = new Random((long)1000);
//        Random rand = new Random((long)System.currentTimeMillis()/1000000000);
        for (int i = 0; i < ComponentNum; i++) {
            super.addDistrib(new DirichletDistrib(EnumDistrib.random(domain, rand.nextInt()), rand.nextInt(90) + 10), rand.nextDouble());
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
    	//data.length/nData = n
    	
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
    	
    	//COMP(Dl, M/n, c), formula 3
    	Double comp3 = ((letters/2)*Math.log(nData/components)) + ((letters-1)/2)*Math.log(cAvg/2) - GammaDistrib.lgamma(letters/2) - ((1/2)*Math.log(letters-1)) ; //excluding 2 constants at end
    	
    	//COMP(Mm, n), formula 4
    	Double comp4 = (((components - 1)/2)*Math.log(nData)) + ((1/2)*Math.log(Math.PI)) - GammaDistrib.lgamma(components/2); // + o(1)
    	
    	//calculate M!
    	double mfact = 1.0;
    	for (double i = 1; i <=components; i++) {
    		mfact = mfact*i;
    	}
    	
    	//COMP(DM(m, L), n, c) ~= COMP(Mm, n) + M*COMP(Dl, M/n, c) - log*(M!)
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
//        String filename = "wgEncodeRad21_seg20_500_hg19.out";
      
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
            dis.learnParameters(data);
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
    
    public static void main0(String[] args) {
        String filename1 = "wgEncodeRad21_seg20_500_hg19.out";
        String filename2 = "/Users/mikael/simhome/Fantom5/suzy_alltags.csv";
        //String filename2 = "/Users/mikael/simhome/Fantom5/suzy_corr.csv";
        if (args.length > 0)
            filename1 = args[0];
        // this CAGE tag, expresses in condition?
        Object[][] matrix1 = loadObjects(filename1);
        int[][] values1 = new int[matrix1.length][];
        for (int i = 0; i < matrix1.length; i ++) {
            values1[i] = new int[matrix1[i].length - 1];
            for (int j = 0; j < matrix1[i].length - 1; j ++) {
                try {
                    values1[i][j] = (Integer)matrix1[i][j + 1];
                } catch (ClassCastException e) {
                    System.err.println("Invalid integer: " + matrix1[i][j + 1]);
                }
            }
        }
        Object[][] matrix2 = loadObjects(filename2);
        int[][] values2 = new int[matrix2.length][];
        for (int i = 0; i < matrix2.length; i ++) {
            values2[i] = new int[matrix2[i].length - 1];
            for (int j = 0; j < matrix2[i].length - 1; j ++) {
                try {
                    values2[i][j] = (Integer)matrix2[i][j + 1];
                } catch (ClassCastException e) {
                    System.err.println("Invalid integer: " + matrix2[i][j + 1]);
                }
            }
        }
        Enumerable domain = new Enumerable(values1[0].length);
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
        int MAX_CLUSTER = 20;
        int[][] labels = new int[values2.length][MAX_CLUSTER - 1];
        for (int ncluster = 2; ncluster <= MAX_CLUSTER; ncluster ++) {
            MixDirichletDistrib dis = new MixDirichletDistrib(domain, ncluster);
            dis.learnParameters(values1);
            System.out.println("\nDone training with " + ncluster + " clusters");
            String prevgene = "";
            int[] cnt4gene = new int[ncluster]; // counts for a gene
            int tot4gene = 0; // total for a gene
            int[] cnt4all = new int[ncluster]; // counts for a gene
            int tot4all = 0; // total for a gene
            int total = 0; // total number of genes
            double entropy = 0;
            for (int k = 0; k < values2.length; k ++) {
                int currcluster = dis.getLabel(values2[k]);
                labels[k][ncluster - 2] = currcluster;
                cnt4gene[currcluster] ++;
                cnt4all[currcluster] ++;
                tot4gene ++;
                tot4all ++;
                String genename = (String)matrix2[k][0];
                if ((!prevgene.equals(genename) || k == values2.length - 1) && k != 0) { // new gene or last tag, determine the entropy of the previous gene
                    double[] p = new double[ncluster];
                    for (int i = 0; i < ncluster; i ++) {
                        p[i] = cnt4gene[i] / (double)tot4gene;
                        cnt4gene[i] = 0;
                    }
                    tot4gene = 0;
                    entropy += getEntropy(p);
                    total ++; 
                    prevgene = genename;
                }
            }
            System.out.println("Average per-gene entropy is " + entropy / (double)total);
            double[] p = new double[ncluster];
            for (int i = 0; i < ncluster; i ++) {
                p[i] = cnt4all[i] / (double)tot4all;
            }
            System.out.println("Overall entropy is " + getEntropy(p));
            for (int i = 0; i < ncluster; i ++) {
                System.out.println(i + "\t" + cnt4all[i] + "\t" + dis.getDistrib(i));
                
            }
            dis.toXMLString();
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename2 + ".out"));
            for (int k = 0; k < labels.length; k ++) {
                bw.write(matrix2[k][0] + "\t");
                for (int j = 0; j < labels[k].length; j ++)
                    bw.write(labels[k][j] + "\t");
                bw.newLine();
            }
            bw.close();
        } catch (IOException ex) {
            System.err.println("Error writing results");
        }
        
        /*
        MixDirichletDistrib dis2 = new MixDirichletDistrib(domain, 9);
        dis2.getNormalized();
        System.out.println(dis2.toXMLString());
        dis2.learnParameters(values1);
        System.out.println(dis2.toString());
        System.out.println(dis2.toXMLString());
                */
    }

    
    public static Object[][] loadObjects(String filename) {
        BufferedReader br = null;
        Object[][] data = null;
        try {
            br = new BufferedReader(new FileReader(filename));
            String line = br.readLine();
            List<Object[]> alldata = new ArrayList<>();
            while (line != null) {
                String[] tokens = line.split("\t");
                Object[] values = new Object[tokens.length];
                for (int i = 0; i < tokens.length; i++) {
                    try {
                        values[i] = Integer.valueOf(tokens[i]); // value is an int
                    } catch (NumberFormatException e1) {
                        try {
                            values[i] = Double.valueOf(tokens[i]); // value is a double
                        } catch (NumberFormatException e2) {
                            values[i] = tokens[i]; // value is a string
                        }
                    }
                }
                alldata.add(values);
                line = br.readLine();
            }
            data = new Object[alldata.size()][];
            for (int k = 0; k < data.length; k++) {
                data[k] = new Object[alldata.get(k).length];
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


