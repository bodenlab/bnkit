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
package bn.alg;

import bn.BNet;
import bn.BNode;
import dat.Variable;
import java.util.List;

/**
 * Interface for Bayesian network learning algorithms.
 *
 * @author mikael
 */
public abstract class LearningAlg {

    protected BNet bn;

    public LearningAlg(BNet bn) {
        this.bn = bn;
    }

    /**
     * Randomly set parameters. Useful prior to EM training.
     *
     * @param seed the seed for the random number generator
     */
    public void setRandom(long seed) {
        int cnt = 0;
        for (BNode node : bn.getNodes())
			//node.randomize(seed+(cnt++))
			;
    }

    public abstract void train(Object[][] values, Variable[] vars, long seed);

    public void train(Object[][] values, List<BNode> nodes, long seed) {
        Variable[] vars = new Variable[nodes.size()];
        for (int i = 0; i < vars.length; i++) {
            vars[i] = nodes.get(i).getVariable();
        }
        train(values, vars, seed);
    }

    public void train(Object[][] values, List<BNode> nodes) {
        train(values, nodes, System.currentTimeMillis());
    }

    public void train(Object[][] values) {
        train(values, bn.getOrdered(), System.currentTimeMillis());
    }

}
