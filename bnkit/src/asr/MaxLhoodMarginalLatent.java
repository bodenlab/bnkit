package asr;

import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import bn.prob.EnumDistrib;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.TreeDecor;
import dat.phylo.TreeInstance;

/**
 * "Marginal" reconstruction view of Bayesian network.
 */
public class MaxLhoodMarginalLatent implements TreeDecor<EnumDistrib> {

    private EnumDistrib value = null;
    final private IdxTree tree;
    final private SubstModel model;
    final private int bpidx;
    final private PhyloBN pbn;
    final private String[] values; // nominated values for external nodes

    /**
     * Set-up inference of the posterior distribution at the specified ancestor branch point, by creating
     * a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model (which also specifies possible latent states).
     * With this constructor, the tree extends leaves with CPTs, which can be instantiated to the values specified here.
     * @param bpidx branch point index for the ancestor
     * @param tree phylogenetic tree indexing the branch points and their distances
     * @param model evolutionary model (implying an alphabet of latent states)
     * @param values discrete values which are used to instantiate extended nodes (defining an alphabet)
     * @param rate relative evolutionary rate (for the index in POG/alignment)
     */
    public MaxLhoodMarginalLatent(int bpidx, IdxTree tree, SubstModel model, String[] values, double rate) {
        this.tree = tree;
        this.model = model;
        this.bpidx = bpidx;
        this.values = values;
        pbn = PhyloBN.withCPTs(tree, model, values, rate);;
    }

    /**
     * Set-up inference of the posterior distribution at the specified ancestor branch point, by creating
     * a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model (which also specifies possible latent states).
     * With this constructor, the tree extends leaves with CPTs, which can be instantiated to the values specified here.
     * @param bpidx branch point index for the ancestor
     * @param tree phylogenetic tree indexing the branch points and their distances
     * @param model evolutionary model (implying an alphabet of latent states)
     * @param values discrete values which are used to instantiate extended nodes (defining an alphabet)
     */
    public MaxLhoodMarginalLatent(int bpidx, IdxTree tree, SubstModel model, String[] values) {
        this.tree = tree;
        this.model = model;
        this.bpidx = bpidx;
        this.values = values;
        pbn = PhyloBN.withCPTs(tree, model, values, 1);;
    }

    /**
     * Set-up inference of the posterior distribution at the specified ancestor branch point, by creating
     * a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model (which also specifies possible latent states).
     * With this constructor, the tree extends leaves with GDTs, with real values.
     * @param bpidx branch point index for the ancestor
     * @param tree phylogenetic tree indexing the branch points and their distances
     * @param model evolutionary model (implying an alphabet of latent states)
     * @param rate relative evolutionary rate (for the index in POG/alignment)
     */
    public MaxLhoodMarginalLatent(int bpidx, IdxTree tree, SubstModel model, double rate) {
        this.tree = tree;
        this.model = model;
        this.bpidx = bpidx;
        this.values = null; // real
        pbn = PhyloBN.withGDTs(tree, model, rate);;
    }

    /**
     * Set-up inference of the posterior distribution at the specified ancestor branch point, by creating
     * a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model (which also specifies possible latent states).
     * With this constructor, the tree extends leaves with GDTs, with real values.
     * @param bpidx branch point index for the ancestor
     * @param tree phylogenetic tree indexing the branch points and their distances
     * @param model evolutionary model
     */
    public MaxLhoodMarginalLatent(int bpidx, IdxTree tree, SubstModel model) {
        this.tree = tree;
        this.model = model;
        this.bpidx = bpidx;
        this.values = null; // real
        pbn = PhyloBN.withGDTs(tree, model, 1);;
    }

    /**
     * Retrieves the already computed distribution
     * @param idx ancestor branch point index (must be the same as when the class instance was created)
     * @return the posterior distribution over states defined by the substitution model
     */
    @Override
    public EnumDistrib getDecoration(int idx) {
        if (idx == bpidx)
            return value;
        else
            throw new ASRRuntimeException("Invalid ancestor index, not defined for marginal inference: " + idx);
    }

    /**
     * Determine the posterior distribution at the ancestor branch point, given the observations at certain branch points in the tree (i.e. leaves)
     * @param ti observed values; note these are the real values that are attached to leaf nodes in the tree
     */
    @Override
    public void decorate(TreeInstance ti) {
        // instantiate all nodes for which there are values, i.e. leaf nodes most probably
        for (int i = 0; i < ti.getSize(); i++) {
            BNode bnode = pbn.getExtNode(i);
            if (bnode != null) { // not hidden, so can be instantiated and inferred
                Object y = ti.getInstance(i);
                bnode.setInstance(y);
            } // else, this branchpoint is outside of the BN, and will be ignored
        }
        if (pbn.isValid()) {
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(pbn.getBN());
            BNode querynode = pbn.getBNode(bpidx);
            if (querynode == null)
                throw new ASRRuntimeException("Marginal inference of invalid branchpoint: " + bpidx);
            Query q = ve.makeQuery(querynode.getVariable());
            CGTable r1 = (CGTable) ve.infer(q);
            value = (EnumDistrib)r1.query(querynode.getVariable());
        } // else the BN is incapable of performing inference, so just leave values as they are
    }

    public PhyloBN getPhyloBN() {
        return pbn;
    }

}

