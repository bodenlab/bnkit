/**
 *
 */
package bn;

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
