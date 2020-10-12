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

        public Joint(IdxTree tree, SubstModel model) {
            pbn = PhyloBN.create(tree, model);
            values = new Object[pbn.bnodes.length];
        }

        @Override
        public Object getDecoration(int idx) {
            return values[idx];
        }

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
        private int ancidx = 0;

        public Marginal(String ancestorID, IdxTree tree, SubstModel model) {
            pbn = PhyloBN.create(tree, model);
            for (int i = 0; i < pbn.bnodes.length; i++) {
                if (pbn.bnodes[i].getVariable().getName().equals(ancestorID)) {
                    ancidx = i;
                    break;
                }
            }
        }

        public Marginal(int ancidx, IdxTree tree, SubstModel model) {
            this.ancidx = ancidx;
            pbn = PhyloBN.create(tree, model);
        }

        @Override
        public EnumDistrib getDecoration(int idx) {
            if (idx == ancidx)
                return value;
            else
                throw new ASRRuntimeException("Invalid ancestor index, not defined for marginal inference: " + idx);
        }

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
