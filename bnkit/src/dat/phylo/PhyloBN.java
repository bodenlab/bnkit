package dat.phylo;

import bn.BNet;
import bn.BNode;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import dat.EnumVariable;

/**
 * Class for a Bayesian network that represents branch points,
 * and their relationships as extracted from a phylogenetic tree.
 * Based on deprecated class PhyloBNet.
 *
 * @author mikael
 */
public class PhyloBN {

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
    public BNet getBN() {
        return bn;
    }

    public BNode[] getBNodes() {
        return bnodes;
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     *
     * @param tree  indexed phylogenetic tree
     * @param model evolutionary model
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBN create(IdxTree tree, SubstModel model) {
        return create(tree, model, 1.0);
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     * Note that the tree could in fact be multiple trees, representing separate insertion events for a locus/position in an alignment.
     * There is a single BN regardless but with separate modules, graphically disconnected.
     *
     * @param tree  phylogenetic tree
     * @param model evolutionary model
     * @param rate  the evolutionary rate to be applied
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBN create(IdxTree tree, SubstModel model, double rate) {
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
