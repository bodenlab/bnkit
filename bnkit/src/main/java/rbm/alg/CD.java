/*
    bnkit -- software for building and using probabilistic models
    including Bayesian networks.
    Copyright (C) 2014-2016  M. Boden et al.

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
package rbm.alg;

import dat.EnumVariable;
import rbm.AbstractRBM;
import rbm.BooleanRBM;

import java.util.Random;

/**
 * Contrastive Divergence
 * @author mikael
 */
public class CD<T extends AbstractRBM> {
    
    private final T rbm;
    private final Random rand;

    public CD(T rbm, long seed) {
        this.rbm = rbm;
        this.rand = new Random(seed);
    }
    
    /**
     * Use the probability of state (mean field) to determine gradient of objective function. If false, use state itself.
     * See Hinton G, A Practical Guide to Training Restricted Boltzmann Machines. UTML TR 2010–003, University of Toronto.
     */
    public boolean USE_PROB_GRADIENT = true;
    
    /**
     * Mini-batches are groups of data points for which a weight update is executed.
     * See section 4 in Hinton G, A Practical Guide to Training Restricted Boltzmann Machines. UTML TR 2010–003, University of Toronto.
     */
    public int MINIBATCH_SIZE = 100;
    public double LEARNING_RATE = 0.01;
    public double MOMENTUM = 0.90;
    
    /**
     * Train the RBM using the specified data, when mapped to the specified variables.
     * If variables are not given, all inputs are used.
     * Specific data values that are null represent "absent value", and are not included in training.
     * @param data
     * @param vars 
     */
    public void train(Object[][] data, EnumVariable[] vars) {
        System.out.println("Number of data points: " + data.length);
        if (vars != null)
            rbm.setLinked(vars);
        Double[][] prev = null;
        Object[][] minibatch = new Object[MINIBATCH_SIZE][];
        for (int round = 0; round < data.length / 2; round ++) {
            for (int p = 0; p < MINIBATCH_SIZE; p ++)
                minibatch[p] = data[rand.nextInt(data.length)];
            Double[][] delta = rbm.getCDGradient(minibatch, 1);
            for (int j = 0; j < delta.length; j ++) {
                for (int i = 0; i < delta[j].length; i ++) {
                    if (delta[j][i] == null)
                        continue;
                    if (delta[j][i] != null && prev == null)
                        delta[j][i] *= LEARNING_RATE;
                    else if (delta[j][i] != null && prev[j][i] != null)
                        delta[j][i] = prev[j][i] * MOMENTUM + delta[j][i] * LEARNING_RATE;
                    else
                        delta[j][i] = 0.0;
                }
            }
            rbm.setCDGradient(delta);
            prev = delta;
            if (round % 100 == 0)
                System.out.printf("%05d:\t%10.3f\n", round, rbm.err);
        }

    }

    /**
     * Train the RBM using the specified data.
     * Specific data values that are null represent "absent value", and are not included in training.
     * @param data
     */
    public void train(Object[][] data) {
        train(data, null);
    }
}
