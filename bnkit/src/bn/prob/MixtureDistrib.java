/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.
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

package bn.prob;

import bn.Distrib;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Mixture of probability distributions.
 * Note that mixtures of mixtures is detected and avoided, leading to a flat mixture.
 * @author mikael
 */
public class MixtureDistrib implements Distrib {

    final Map<Distrib, Double> mixture;
    protected ArrayList<Distrib> distribs;
    protected ArrayList<Double> weights;
    
    private double density;
    Random rand = new Random(1);
    
    /**
     * construct an empty mixture distribution
     * you will have a mixture without any component
     */
    public MixtureDistrib() {
    	mixture = new HashMap<>();
    	distribs = new ArrayList<Distrib>();
    	weights = new ArrayList<Double>();
    	density = 0.0;
    }
    
    /**
     * construct a mixture distribution with a distribution
     * @param d1, could be mixture distribution, or single distribution
     * @param weight1, the corrsponding weight
     */
    public MixtureDistrib(Distrib d1, double weight1) {
    	distribs = new ArrayList<Distrib>();
    	weights = new ArrayList<Double>();
    	density = 0.0;
        mixture = new HashMap<>();
        /*
        try {
            MixtureDistrib packed = (MixtureDistrib) d1;
            for (Map.Entry<Distrib, Double> entry : packed.mixture.entrySet()) {
                double delta = entry.getValue() * weight1;
                mixture.put(entry.getKey(), delta);
                density += delta;
            }
        } catch (ClassCastException e) {
            mixture.put(d1, weight1);
            density = weight1;
        }*/
    	if(d1 instanceof MixtureDistrib) {
    		MixtureDistrib MD = (MixtureDistrib)d1;
    		int size = MD.getMixtureSize();
    		for(int i = 0; i < size; i++) {
    			Distrib distrib = MD.getDistrib(i);
    			double delta = weight1 * MD.getWeights(i);
    			addDistribForced(distrib, delta);
    			mixture.put(distrib, delta);
    		}
    	} else {
    		addDistribForced(d1, weight1);
    		mixture.put(d1, weight1);
    	}
    }
    
    /**
     * Set seed for random number generator.
     * @param seed 
     */
    public void setSeed(long seed) {
        rand = new Random(seed);
    }
    
    /**
     * Get an integer using the random number generator.
     * @param max value
     * @return 
     */
    public int nextInt(int max) {
        return rand.nextInt(max);
    }
    
    /**
     * Get a double using the random number generator.
     * @return 
     */
    public double nextDouble() {
        return rand.nextDouble();
    }
    
    private double addDistribForced(Distrib d2, double weight2) {
    	
        Double prev_weight = mixture.get(d2);
        if (prev_weight == null)
            mixture.put(d2, weight2);
        else
            mixture.put(d2, prev_weight + weight2);
    	if(hasDistrib(d2)) {
    		int index = distribs.indexOf(d2);
    		weights.set(index, weights.get(index) + weight2);
    	} else {
    		distribs.add(d2);
    		weights.add(weight2);
    	}
    	density += weight2;
        return density;
    }
    
    /**
     * add a distribution to this mixture model
     * @param d2 could be mixture distribution, or single distribution
     * @param weight2 the corrsponding weight
     * @return
     */
    public double addDistrib(Distrib d2, double weight2) {
    	/*
        try {
            MixtureDistrib packed = (MixtureDistrib) d2;
            for (Map.Entry<Distrib, Double> entry : packed.mixture.entrySet())
                addDistribForced(entry.getKey(), entry.getValue() * weight2);
        } catch (ClassCastException e) {
            addDistribForced(d2, weight2);
        }*/
    	if(d2 instanceof MixtureDistrib) {
    		MixtureDistrib MD = (MixtureDistrib)d2;
    		int size = MD.getMixtureSize();
    		for(int i = 0; i < size; i++) {
    			Distrib distrib = MD.getDistrib(i);
    			addDistribForced(distrib, weight2 * MD.getWeights(i));
    		}
    	}else {
    		addDistribForced(d2, weight2);
    	}
        return density;
    }
    
    /**
     * check if mixture has this distribution
     * @param d2
     * @return
     */
    public boolean hasDistrib(Distrib d2) {
        return distribs.contains(d2);
    }
    
    /**
     * get the size of mixture
     * @return
     */
    public int getMixtureSize() {
    	return this.distribs.size();
    }
    
    /**
     * get distribution given index
     * @param index
     * @return
     */
    public Distrib getDistrib(int index) {
    	if(index < getMixtureSize()) {
    		return distribs.get(index);
    	}
    	return null;
    }
    
    /**
     * get the weight given index
     * @param index
     * @return
     */
    public Double getWeights(int index) {
    	if(index < getMixtureSize()) {
    		return weights.get(index);
    	}
    	return null;
    }
    
    /**
     * set new weight given the index
     * @param index
     * @param neWeight
     */
    public void setWeight(int index, double neWeight) {
    	this.weights.set(index, neWeight);
    	this.mixture.put(this.distribs.get(index), Double.valueOf(neWeight));
    }
    
    /**
     * set a new set of weights
     * @param neWeight
     */
    public void setWeights(double[] neWeight) {
    	if(neWeight.length != this.weights.size()) {
    		throw new RuntimeException("number of weights invalid");
    	}
    	
    	for(int i = 0; i < neWeight.length; i++) {
    		setWeight(i, neWeight[i]);
    	}
    }
    
    /**
     * get all weights
     * @return
     */
    public double[] getAllWeights() {
    	double[] weight = new double[weights.size()];
        double sum = 0;
    	for (int i = 0; i < weights.size(); i++) {
            weight[i] = weights.get(i).doubleValue();
            sum += weight[i];
    	}
        for (int i = 0; i < weight.length; i ++)
            weight[i] /= sum;
    	return weight;
    }
    
    /**
     * get weights given distribution
     * @param distrib
     * @return
     */
    public Double getWeightsByDistrib(Distrib distrib) {
    	if(mixture.containsKey(distrib)) {
    		return mixture.get(distrib);
    	}
    	return null;
    }
    
    /**
     * Generate a clone which is normalized.
     * @return mixture with weights that add up to 1
     */
    public MixtureDistrib getNormalizedClone() {
        double sum = 0;
        for (Double weight : weights) 
            sum += weight;
        MixtureDistrib md = null;
        //for (Map.Entry<Distrib, Double> entry : mixture.entrySet()) {
        for(int i = 0; i < distribs.size(); i++) {
            if (md == null) 
                md = new MixtureDistrib(distribs.get(i), weights.get(i) / sum);
            else
                md.addDistrib(distribs.get(i), weights.get(i) / sum);
        }
        return md;
    }
    
    /**
     * normalized the mixture weights sum up to be 1
     */
    public void getNormalized() {
    	double sum = 0.0;
    	for(Double weight: weights) {
    		sum += weight;
    	}
    	for(int i = 0; i < weights.size(); i++) {
    		weights.set(i, weights.get(i) / sum);
    		this.mixture.put(distribs.get(i), weights.get(i));
    	}
    }
    
    @Override
    public double get(Object value) {
        double p = 0.0;
        /*
        for (Map.Entry<Distrib, Double> entry : mixture.entrySet()) {
            p += entry.getKey().get(value) * entry.getValue();
        }*/
        for(int i = 0; i < distribs.size(); i++) {
        	p += distribs.get(i).get(value) * weights.get(i);
        }
        return p;
    }
    
    public Distrib componentSample() {
    	double y = rand.nextDouble() * density;
    	Distrib current = null;
        double p = 0.0;
        /*
        for (Map.Entry<Distrib, Double> entry : mixture.entrySet()) {
            current = entry.getKey();
            p += entry.getValue();
            if (p >= y)
                break;
        }*/
        for(int i = 0; i < distribs.size(); i++) {
        	current = distribs.get(i);
        	p += weights.get(i);
        	if (p >= y)
                break;
        }
        return current;
    }
    
    public boolean equals(MixtureDistrib md) {
    	for(int i = 0; i < distribs.size(); i++) {
    		
    		DirichletDistrib dir1 = (DirichletDistrib)getDistrib(i);
    		DirichletDistrib dir2 = (DirichletDistrib)md.getDistrib(i);
    		if(!dir1.equals(dir2)) {
    			return false;
    		}
    		Distrib distrib = getDistrib(i);
    		if(Math.abs(md.getWeightsByDistrib(dir2) - getWeightsByDistrib(dir1)) > 1e-15) {
    			return false;
    		}
    	}
    	
    	return true;
    } 

    /**
     * Sample the mixture distribution.
     * @return the sample
     */
    @Override
    public Object sample() {
        Distrib current = componentSample();
        if (current == null)
            throw new RuntimeException("Invalid MixtureDistrib");
        return current.sample();
    }
    
//    @Override
//    public String toString() {
//        StringBuilder sb = new StringBuilder("mixture size: " + distribs.size() + "\n");
//        //for (Map.Entry<Distrib, Double> entry : mixture.entrySet())
//        for(int i = 0; i < distribs.size(); i++)
//            sb.append("weight: " + String.format("%4.2f", weights.get(i)) + "; distrib: " + distribs.get(i) + "\n");
//        return sb.toString();
//    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("^" + distribs.size());
        //for (Map.Entry<Distrib, Double> entry : mixture.entrySet())
        for(int i = 0; i < distribs.size(); i++)
            sb.append("{" + distribs.get(i) + "*" + String.format("%4.2f", weights.get(i)) + "}");
        return sb.toString();
    }

    public String toXMLString() {
    	StringBuilder sb = new StringBuilder("<MixtureModels>\n");
    	for(int i = 0; i < distribs.size(); i++) {
    		sb.append("<model>\n<weight>" + weights.get(i) + "</weight>\n<distrib>" + distribs.get(i) + "</distrib>\n</model>\n");
    	}
    	sb.append("</MixtureModels>");
    	return sb.toString();
    }
    
    public static void main(String[] args) {
    	
        GaussianDistrib gd1 = new GaussianDistrib(0, 1.0);
        System.out.println("gd1 = " + gd1);
        GaussianDistrib gd2 = new GaussianDistrib(1, 0.5);
        System.out.println("gd2 = " + gd2);
        GaussianDistrib gd3 = new GaussianDistrib(-2, 2.5);
        System.out.println("gd3 = " + gd3);
        
        MixtureDistrib md1 = new MixtureDistrib(gd1, 1.0);
        md1.addDistrib(gd2, 2);
        md1.addDistrib(gd2, 0.5);
        System.out.println("md1 is gd1*1.0 + gd2*2.5 : \n" + md1);
        
        MixtureDistrib md2 = new MixtureDistrib(md1, 1.0);
        System.out.println("mds2 is md1*1.0 : \n" + md2);
        md2.addDistrib(gd1, 0.5);
        System.out.println("md2 += gd1*0.5 : \n" + md2);
        md2.addDistrib(gd3, 2);
        System.out.println("md2 += gd3*2.0 : \n" + md2);
        md2.addDistrib(md1, 2.0);
        System.out.println("md2 += md1*2.0 : \n" + md2);
        md2.addDistrib(gd1, 1.5);
        
        System.out.println("md2 += gd1*1.5 : \n" + md2);
        System.out.println("density = " + md2.density);
    }
}
