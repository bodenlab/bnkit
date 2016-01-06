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
public class UniformDistrib implements Distrib, Serializable {

    private static final long serialVersionUID = 1L;

    private double lower;
    private double upper;
    private double range;
    
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
    
    @Override
    public String toString() {
        return String.format("%4.2f..%4.2f", lower, upper);
    }
    
    @Override
    public double get(Object value) {
        return 1.0 / range;
    }

    @Override
    public Double sample() {
        double y = rand.nextDouble() * range + lower;
        return y;
    }
    
}
