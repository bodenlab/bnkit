package asr;

import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import dat.Variable;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.TimeUnit;

/**
 * "Joint" reconstruction view of Bayesian network.
 */
public class MaxLhoodJoint implements TreeDecor<Object> {

    final private Object[] values;
    final private PhyloBN pbn;

    /**
     * Set-up a joint reconstruction; this creates a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model.
     * @param tree
     * @param model
     */
    public MaxLhoodJoint(IdxTree tree, SubstModel model) {
        pbn = PhyloBN.create(tree, model);
        values = new Object[pbn.getBNodes().length];
    }

    /**
     * Retrieve the inferred state for a specified branch point, as determined by max likelihood
     * @param idx branch point index
     * @return state which jointly with all others assigns the greatest probability to the observed values (at leaves)
     */
    @Override
    public Object getDecoration(int idx) {
        return values[idx];
    }

    /**
     * Determine the joint state that assigns the maximum likelihood to the specified observations
     * @param ti observed tree states
     */
    @Override
    public void decorate(TreeInstance ti) {
        Map<String, Integer> quick = new HashMap<>();
        // instantiate all nodes for which there are values, i.e. leaf nodes most probably but not necessarily.
        // also, leaf nodes do not need to be instantiated.
        BNode[] bnodes = pbn.getBNodes();
        for (int i = 0; i < bnodes.length; i++) {
            quick.put(bnodes[i].getVariable().getName(), i);
            Object y = ti.getInstance(i);
            if (y != null) {
                bnodes[i].setInstance(y);
                values[i] = y;
            } else
                bnodes[i].resetInstance();
        }
        // set-up the inference engine
        VarElim ve = new VarElim();
        ve.instantiate(pbn.getBN());
        Query q_mpe = ve.makeMPE();
        CGTable r1 = (CGTable) ve.infer(q_mpe);
        Variable.Assignment[] assign = r1.getMPE();
        for (Variable.Assignment assign1 : assign) {
            Integer idx = quick.get(assign1.var.getName());
            if (idx != null)
                values[idx] = assign1.val;
        }
    }
}
