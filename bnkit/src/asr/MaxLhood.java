package asr;

import bn.BNet;
import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;
import dat.phylo.BranchPoint;
import dat.phylo.IdxTree;
import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;

import java.util.HashMap;
import java.util.Map;

public class MaxLhood {

    private MaxLhood() {
    }

    /**
     * Class for a Bayesian network that represents branch points,
     * and their relationships as extracted from a phylogenetic tree.
     * Based on deprecated class PhyloBNet.
     *
     * @author mikael
     */
    public static class PhyloBN {

        private final BNet bn;
        private final SubstModel model;
        private double rate = 1.0;
        private BNode[] bnodes = null;

        /**
         * Private constructor so that instances are created for a particular end.
         *
         * @param model substitution model
         */
        private PhyloBN(SubstModel model) {
            bn = new BNet();
            this.model = model;
        }

        /**
         * Get the BN
         *
         * @return Bayesian network instance
         */
        private BNet getBN() {
            return bn;
        }

        /**
         * Construct a BN for specified phylogenetic tree using a single, supplied model.
         *
         * @param tree  indexed phylogenetic tree
         * @param model evolutionary model
         * @return the phylogenetic Bayesian network
         */
        private static PhyloBN create(IdxTree tree, SubstModel model) {
            return create(tree, model, 1.0);
        }

        /**
         * Construct a BN for specified phylogenetic tree using a single, supplied model.
         *
         * @param tree  phylogenetic tree
         * @param model evolutionary model
         * @param rate  the evolutionary rate to be applied
         * @return the phylogenetic Bayesian network
         */
        private static PhyloBN create(IdxTree tree, SubstModel model, double rate) {
            PhyloBN pbn = new PhyloBN(model);
            pbn.rate = rate;
            EnumVariable[] bvars = new EnumVariable[tree.getSize()];
            pbn.bnodes = new BNode[tree.getSize()];
            // first, create variables...
            for (int idx : tree) {
                BranchPoint bp = tree.getBranchPoint(idx);
                EnumVariable rvar = new EnumVariable(model.getDomain(), bp.getID().toString());
                bvars[idx] = rvar;
            }
            // second, create nodes...
            for (int idx : tree) {
                BranchPoint bp = tree.getBranchPoint(idx);
                EnumVariable rvar = bvars[idx];
                int parent = tree.getParent(idx);
                if (parent < 0) { // there's no parent
                    pbn.bnodes[idx] = new SubstNode(rvar, model);
                } else { // there's a parent
                    EnumVariable prvar = bvars[parent];
                    pbn.bnodes[idx] = new SubstNode(rvar, prvar, model, bp.getDistance());
                }
            }
            pbn.bn.add(pbn.bnodes);
            return pbn;
        }
    }
    /**
     * "Joint" reconstruction view of Bayesian network.
     */
    public static class Joint implements TreeDecor<Object> {

        final private Object[] values;
        final private PhyloBN pbn;

        /**
         * Set-up a joint reconstruction; this creates a probabilistic graphical model with the structure and evolutionary
         * distances as defined by the supplied tree and substitution model.
         * @param tree
         * @param model
         */
        public Joint(IdxTree tree, SubstModel model) {
            pbn = PhyloBN.create(tree, model);
            values = new Object[pbn.bnodes.length];
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
            for (int i = 0; i < pbn.bnodes.length; i++) {
                quick.put(pbn.bnodes[i].getVariable().getName(), i);
                Object y = ti.getInstance(i);
                if (y != null) {
                    pbn.bnodes[i].setInstance(y);
                    values[i] = y;
                } else
                    pbn.bnodes[i].resetInstance();
            }
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(pbn.bn);
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

    /**
     * "Marginal" reconstruction view of Bayesian network.
     */
    public static class Marginal implements TreeDecor<EnumDistrib> {

        private EnumDistrib value = null;
        final private PhyloBN pbn;
        final private int ancidx;

        /**
         * Set-up inference of the posterior distribution at the specified ancestor branch point, by creating
         * a probabilistic graphical model with the structure and evolutionary
         * distances as defined by the supplied tree and substitution model.
         * @param bpidx branch point index for the ancestor
         * @param tree phylogenetic tree indexing the branch points and their distances
         * @param model evolutionary model
         */
        public Marginal(int bpidx, IdxTree tree, SubstModel model) {
            this.ancidx = bpidx;
            pbn = PhyloBN.create(tree, model);
        }

        /**
         * Retrieves the already computed distribution
         * @param idx ancestor branch point index (must be the same as when the class instance was created)
         * @return the posterior distribution over states defined by the substitution model
         */
        @Override
        public EnumDistrib getDecoration(int idx) {
            if (idx == ancidx)
                return value;
            else
                throw new ASRRuntimeException("Invalid ancestor index, not defined for marginal inference: " + idx);
        }

        /**
         * Determine the posterior distribution at the ancestor branch point, given the observations at certain branch points in the tree (i.e. leaves)
         * @param ti observed states
         */
        @Override
        public void decorate(TreeInstance ti) {
            // instantiate all nodes for which there are values, i.e. leaf nodes most probably
            for (int i = 0; i < pbn.bnodes.length; i++) {
                Object y = ti.getInstance(i);
                if (y != null)
                    pbn.bnodes[i].setInstance(y);
                else
                    pbn.bnodes[i].resetInstance();
            }
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(pbn.bn);
            Query q = ve.makeQuery(pbn.bnodes[ancidx].getVariable());
            CGTable r1 = (CGTable) ve.infer(q);
            value = (EnumDistrib)r1.query(pbn.bnodes[ancidx].getVariable());
        }
    }

}
