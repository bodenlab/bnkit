/*
 * Copyright (c) 2010-2021 Haifeng Li. All rights reserved.
 *
 * Smile is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Smile is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Smile.  If not, see <https://www.gnu.org/licenses/>.
 */

package smile.math;

import smile.math.random.MersenneTwister;
import smile.math.random.UniversalGenerator;

import java.util.stream.IntStream;

/**
 * This is a high quality random number generator as a replacement of
 * the standard Random class of Java system.
 * 
 * @author Haifeng Li
 */
public class Random {

    private final UniversalGenerator real;
    private final MersenneTwister twister;

    /**
     * Initialize with default random number generator engine.
     */
    public Random() {
        real = new UniversalGenerator();
        twister = new MersenneTwister();
    }

    /**
     * Initialize with given seed for default random number generator engine.
     * @param seed the RNG seed.
     */
    public Random(long seed) {
        real = new UniversalGenerator(seed);
        twister = new MersenneTwister(seed);
    }

    /**
     * Initialize the random generator with a seed.
     * @param seed the RNG seed.
     */
    public void setSeed(long seed) {
        real.setSeed(seed);
        twister.setSeed(seed);
    }

    /**
     * Generator a random number uniformly distributed in [0, 1).
     * @return a pseudo random number
     */
    public double nextDouble() {
        return real.nextDouble();
    }

    /**
     * Generate n uniform random numbers in the range [0, 1)
     * @param d array of random numbers to be generated
     */
    public void nextDoubles(double[] d) {
        real.nextDoubles(d);
    }

    /**
     * Generate a uniform random number in the range [lo, hi)
     * @param lo lower limit of range
     * @param hi upper limit of range
     * @return a uniform random real in the range [lo, hi)
     */
    public double nextDouble(double lo, double hi) {
        return (lo + (hi - lo) * nextDouble());
    }

    /**
     * Generate n uniform random numbers in the range [lo, hi)
     * @param lo lower limit of range
     * @param hi upper limit of range
     * @param d array of random numbers to be generated
     */
    public void nextDoubles(double[] d, double lo, double hi) {
        real.nextDoubles(d);

        double l = hi - lo;        
        int n = d.length;
        for (int i = 0; i < n; i++) {
            d[i] = lo + l * d[i];
        }
    }

    /**
     * Returns a random integer.
     * @return a random integer.
     */
    public int nextInt() {
        return twister.nextInt();
    }
    
    /**
     * Returns a random integer in [0, n).
     * @param n the upper bound of random number.
     * @return a random integer.
     */
    public int nextInt(int n) {
        return twister.nextInt(n);
    }

    /**
     * Returns a random long integer.
     * @return a random long integer.
     */
    public long nextLong() {
        return twister.nextLong();
    }

    /**
     * Returns a permutation of <code>(0, 1, 2, ..., n-1)</code>.
     *
     * @param n the upper bound.
     * @return the permutation of <code>(0, 1, 2, ..., n-1)</code>.
     */
    public int[] permutate(int n) {
        int[] x = IntStream.range(0, n).toArray();
        permutate(x);
        return x;
    }

    /**
     * Permutates an array.
     * @param x the array.
     */
    public void permutate(int[] x) {
        for (int i = 0; i < x.length; i++) {
            int j = i + nextInt(x.length - i);
            MathEx.swap(x, i, j);
        }
    }

    /**
     * Permutates an array.
     * @param x the array.
     */
    public void permutate(float[] x) {
        for (int i = 0; i < x.length; i++) {
            int j = i + nextInt(x.length - i);
            MathEx.swap(x, i, j);
        }
    }

    /**
     * Permutates an array.
     * @param x the array.
     */
    public void permutate(double[] x) {
        for (int i = 0; i < x.length; i++) {
            int j = i + nextInt(x.length - i);
            MathEx.swap(x, i, j);
        }
    }

    /**
     * Permutates an array.
     * @param x the array.
     */
    public void permutate(Object[] x) {
        for (int i = 0; i < x.length; i++) {
            int j = i + nextInt(x.length - i);
            MathEx.swap(x, i, j);
        }
    }
}