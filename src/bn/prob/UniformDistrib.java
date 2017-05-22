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
package bn.prob;

import bn.Distrib;
import dat.Enumerable;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author mikael
 */
public class UniformDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;

    private double lower;
    private double upper;
    private double range;
    private double delta = 0.0;// minimum range quotient
    
    private final Random rand = new Random();

    public UniformDistrib(double lower, double upper) {
        this.lower = lower;
        this.upper = upper;
        if (lower >= upper) 
            throw new RuntimeException("Invalid range for uniform distribution: " + lower + " to " + upper);
        this.range = upper - lower;
    }

    public double getLower() {
        return lower;
    }
    
    public double getUpper() {
        return upper;
    }
    
    public void setLower(double lower) {
        this.lower = lower;
        if (lower >= upper) 
            throw new RuntimeException("Invalid range for uniform distribution: " + lower + " to " + upper);
        this.range = upper - lower;
    }
    
    public void setUpper(double upper) {
        this.upper = upper;
        if (lower >= upper) 
            throw new RuntimeException("Invalid range for uniform distribution: " + lower + " to " + upper);
        this.range = upper - lower;
    }
    
    public void setBounds(double lower, double upper) {
        this.lower = lower;
        this.upper = upper;
        if (lower >= upper) 
            throw new RuntimeException("Invalid range for uniform distribution: " + lower + " to " + upper);
        this.range = upper - lower;
    }
    
    public void setDelta(double delta) {
        this.delta = delta;
    }
    
    @Override
    public String toString() {
        return String.format("%4.2f..%4.2f", lower, upper);
    }
    
    @Override
    public double get(Object value) {
        if ((Double)value >= lower && (Double)value <= upper)
            return 1.0 / (range + delta);
        else
            return 0.0;
    }

    public double getCumulative(Object value) {
        double x = (Double)value;
        if (x < lower)
            return 0;
        else if (x >= upper)
            return 1;
        else
            return (x - lower) / (upper - lower);
    }
    
    public double getMean() {
        return (upper + lower) / 2.0;
    }
    
    public double getVariance() {
        return ((upper - lower) * (upper - lower)) / 12.0;
    }
    
    @Override
    public Double sample() {
        double y = rand.nextDouble() * range + lower;
        return y;
    }
    
    /**
     * Create a density matching specified samples.
     * @param samples samples
     * @return a new distribution which could serve as a starting point for representing a sub-group of those in the samples.
     */
    public static UniformDistrib probe(double[] samples) {
        if (samples.length < 2) 
            throw new RuntimeException("Invalid sample");
        double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY, mean = 0;
        Random myrand = new Random();
        int choose1 = myrand.nextInt(samples.length);
        int choose2 = myrand.nextInt(samples.length);
        while (choose2 == choose1)
            choose2 = myrand.nextInt(samples.length);
        if (choose1 > choose2) {
            int tmp = choose1;
            choose1 = choose2;
            choose2 = tmp;
        }
        for (int i = choose1; i < choose2; i ++) {
            if (samples[i] < min)
                min = samples[i];
            if (samples[i] > max)
                max = samples[i];
            mean += samples[i] / samples.length;
        }
        return new UniformDistrib(min, max);
    }
    

    /**
     * Example that finds a mixture of Uniform densities using Gibbs sampling.
     * @param args 
     */
    public static void main(String[] args) {
        int clusters = 2;
        Enumerable label = new Enumerable(clusters);
        EnumDistrib C = EnumDistrib.random(label); // mixture distrib
        UniformDistrib[] Q = new UniformDistrib[clusters];
        double[] X = {2.9, 3.3, 4.3, 6.1, 6.4, 7.5, 8.1, 8.8};
        
        // initialisation
        Random r = new Random();

        int[] Z_sample = new int[X.length];
        EnumDistrib[] Z = new EnumDistrib[X.length];
        for (int k = 0; k < clusters; k ++) {
            Q[k] = UniformDistrib.probe(X);
            System.out.println("P(Q|Z = " + k + ") = " + Q[k]);
        }
        double[] counts = new double[clusters];
        for (int i = 0; i < X.length; i ++) {
            Z[i] = new EnumDistrib(label);
            Z_sample[i] = r.nextInt(clusters);
            counts[Z_sample[i]] += 1.0;
            System.out.println("X_" + i + " = " + X[i] + ": Z = " + Z_sample[i]);
        }
        C.set(counts);
        for (int k = 0; k < clusters; k ++) {
            System.out.println("P(Z = " + k + ") = " + C.get(k));
        }
            
        // training
        
        System.out.println("Training");
        
        for (int round = 0; round < 10; round ++) {
            int i_z = r.nextInt(X.length);
            System.out.println("Hold-out X_" + i_z + " = " + X[i_z]);
            
            for (int k = 0; k < clusters; k ++) {
                System.out.println("\tSamples labelled Z = " + k + ": ");
                double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
                int cnt = 0;
                for (int i = 0; i < X.length; i ++) {
                    if (i != i_z && k == Z_sample[i]) {
                        if (X[i] < min)
                            min = X[i];
                        if (X[i] > max)
                            max = X[i];
                        System.out.println("\t\tX_" + i + " = " + X[i]);
                        cnt ++;
                    }
                }
                if (min < max)
                    Q[k].setBounds(min, max);
                System.out.println("\t\tNew P(Q|Z = " + k + ") = " + Q[k]);
            }
            
            int i = i_z;
            double[] p = new double[clusters];
            double pseudo = 0.01;
            for (int k = 0; k < clusters; k ++) {
                p[k] = Q[k].get(X[i]) + pseudo;
            }
            Z[i].set(p);
            Z_sample[i] = (int)Z[i].sample();
            System.out.println("\tX_" + i + " = " + X[i] + ": P(Z|Q = " + X[i] + ") = " + Z[i] + " =sample=> " + Z_sample[i]);
            counts = new double[clusters];
            for (i = 0; i < X.length; i ++) {
                counts[Z_sample[i]] += 1.0;
            }
            C.set(counts);
        }
        
        System.out.println("RESULT:");
        for (int i = 0; i < X.length; i ++) {
            System.out.println("X_" + i + " = " + X[i] + ": P(Z|Q = " + X[i] + ") = " + Z[i] + " =sample=> " + Z_sample[i]);
        }
        
    }
    
}
