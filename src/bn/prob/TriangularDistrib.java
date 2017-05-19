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
import java.io.Serializable;
import java.util.Random;

/**
 *
 * @author mikael
 */
public class TriangularDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;

    private double lower;
    private double upper;
    private double mode;
    private double range;
    
    private final Random rand = new Random();

    public TriangularDistrib(double lower, double upper, double mode) {
        this.lower = lower;
        this.upper = upper;
        this.mode = mode;
        if (lower >= upper || mode > upper || mode < lower) 
            throw new RuntimeException("Invalid range for distribution: " + lower + " to " + upper + " with mode " + mode);
        this.range = upper - lower;
    }

    public double getLower() {
        return lower;
    }
    
    public double getUpper() {
        return upper;
    }
    
    public double getMode() {
        return mode;
    }
    
    public void setLower(double lower) {
        this.lower = lower;
        if (lower >= upper || mode < lower) 
            throw new RuntimeException("Invalid range for distribution: " + lower + " to " + upper + " with mode " + mode);
        this.range = upper - lower;
    }
    
    public void setUpper(double upper) {
        this.upper = upper;
        if (lower >= upper || mode > upper) 
            throw new RuntimeException("Invalid range for distribution: " + lower + " to " + upper + " with mode " + mode);
        this.range = upper - lower;
    }
    
    public void setMode(double mode) {
        this.mode = mode;
        if (mode > upper || mode < lower) 
            throw new RuntimeException("Invalid range for distribution: " + lower + " to " + upper + " with mode " + mode);
    }
    
    public void setParams(double lower, double upper, double mode) {
        this.lower = lower;
        this.upper = upper;
        this.mode = mode;
        if (lower >= upper || mode > upper || mode < lower) 
            throw new RuntimeException("Invalid range for distribution: " + lower + " to " + upper + " with mode " + mode);
        this.range = upper - lower;
    }
   
    @Override
    public String toString() {
        return String.format("%4.2f..%4.2f..%4.2f", lower, mode, upper);
    }
    
    @Override
    public double get(Object value) {
        double x = (Double)value;
        if (x >= lower && x <= upper) {
            if (x < mode)
                return 2.0 * (x - lower) / ((upper - lower) * (mode - lower));
            else if (x > mode) 
                return 2.0 * (upper - x) / ((upper - lower) * (upper - mode));
            else // 
                return 2.0 / (upper - lower);
        } else
            return 0.0;
    }

    public double getCumulative(Object value) {
        double x = (Double)value;
        if (x < lower)
            return 0;
        else if (x >= upper)
            return 1;
        else if (x <= mode) {
            return ((x - lower) * (x - lower)) / ((upper - lower) * (mode - lower));
        } else { // x > mode but x < upper
            return 1.0 - ((upper - x) * (upper - x)) / ((upper - lower) * (upper - mode));
        }
    }
    
    public double getMean() {
        return (upper + lower + mode) / 3.0;
    }
    
    public double getVariance() {
        return ((lower * lower) + (upper * upper) + (mode * mode) - (lower * upper) - (lower * mode) - (upper * mode)) / 18.0;
    }
    
    /**
     * Sample from this distribution.
     * https://en.wikipedia.org/wiki/Triangular_distribution
     * @return 
     */
    @Override
    public Double sample() {
        double U = rand.nextDouble();
        double F_c = (mode - lower) / (upper - mode);
        if (U < F_c)
            return lower + Math.sqrt(U * (upper - lower) * (mode - lower));
        else // U >= F_c
            return upper - Math.sqrt((1.0 - U) * (upper - lower) * (upper - mode));
    }
    
    /**
     * Create a density matching specified samples.
     * @param samples samples
     * @return a new distribution matching the samples (not necessarily optimal).
     */
    public static TriangularDistrib probe(double[] samples) {
        if (samples.length < 2) 
            throw new RuntimeException("Invalid sample");
        double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY, mean = 0;
        for (int i = 0; i < samples.length; i ++) {
            if (samples[i] < min)
                min = samples[i];
            if (samples[i] > max)
                max = samples[i];
            mean += samples[i] / samples.length;
        }
        return new TriangularDistrib(min, max, mean);
    }
    
    
}
