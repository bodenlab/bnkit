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

import dat.Variable;
import rbm.AbstractRBM;

/**
 * Contrastive Divergence
 * @author mikael
 */
public class CD {
    
    private final AbstractRBM rbm;
    
    public CD(AbstractRBM rbm) {
        this.rbm = rbm;
    }
    
    /**
     * Use the probability of state to determine gradient of objective function. If false, use state itself.
     * See Hinton G, A Practical Guide to Training Restricted Boltzmann Machines. UTML TR 2010–003, University of Toronto.
     */
    public boolean USE_PROB_GRADIENT = true;
    
    /**
     * Mini-batches are groups of data points for which a weight update is executed.
     * See section 4 in Hinton G, A Practical Guide to Training Restricted Boltzmann Machines. UTML TR 2010–003, University of Toronto.
     */
    public int MINIBATCH_SIZE = 10; 
    
    
    /**
     * Train the RBM using the specified data, when mapped to the specified variables.
     * @param data
     * @param vars 
     */
    public void train(Object[][] data, Variable[] vars) {
        
    }
}
