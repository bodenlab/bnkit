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
    /** Default relative rate */
    public static double DEFAULT_RATE = 1.0;
    private BNode[] bp2node = null; // branchpoint index to BN node

    /**
     * Private constructor so that instances are created for a particular end.
     *
     */
    private PhyloBN() {
        bn = new BNet();
    }

    /**
     * Get the BN
     *
     * @return Bayesian network instance
     */
    public BNet getBN() {
        return bn;
    }

    /**
     *
     * @return
     * @deprecated
     */
    public BNode[] getBNodes() {
        return bp2node;
    }

    /**
     * Get the Bayesian network node for a given branchpoint index of the phylogenetic tree used to create the network.
     * @param bpidx branchpoint index
     * @return the instance of the BNode if available; note that the network only represents a subset of the branchpoints and
     * if not available null is returned
     */
    public BNode getBNode(int bpidx) {
        return bp2node[bpidx];
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     * The BN will create nodes ONLY for variables that are accessible and valid in the tree.
     * Note that trees are sometimes pruned, which removes the need to include many variables.
     * @param tree  indexed phylogenetic tree
     * @param model evolutionary model
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBN create(IdxTree tree, SubstModel model) {
        return create(tree, model, DEFAULT_RATE);
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     * Note that the tree could in fact be multiple trees, representing separate insertion events for a locus/position in an alignment.
     * There is a single BN regardless but with separate modules, graphically disconnected.
     *
     * @param tree  phylogenetic tree instance
     * @param model evolutionary model
     * @param rate  the evolutionary rate to be applied
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBN create(IdxTree tree, SubstModel model, double rate) {
        PhyloBN pbn = new PhyloBN();
        pbn.bp2node = new BNode[tree.getSize()];
        //  create variables and nodes
        for (int idx : tree) { // iterate through tree depth-first
            if (tree.isConnected(idx)) { // only create variable for branchpoint if it is a parent or a child of one
                BranchPoint bp = tree.getBranchPoint(idx);
                EnumVariable rvar = new EnumVariable(model.getDomain(), bp.getID().toString());
                int parent = tree.getParent(idx);
                if (parent < 0) { // there's no parent
                    pbn.bp2node[idx] = new SubstNode(rvar, model);
                } else { // there's a parent, and because the iterator is "depth-first", the parent must already have been created
                    EnumVariable prvar = (EnumVariable) pbn.bp2node[parent].getVariable();
                    pbn.bp2node[idx] = new SubstNode(rvar, prvar, model, bp.getDistance() * rate); // this is where relative rate is incorporated
                }
                pbn.bn.add(pbn.bp2node[idx]);
            }

        }
        return pbn;
    }

    public static PhyloBN createBarebone(SubstNode[] nodes) {
        PhyloBN pbn = new PhyloBN();
        pbn.bp2node = new BNode[nodes.length];
        for (int i = 0; i < nodes.length; i ++) {
            pbn.bp2node[i] = nodes[i];
        }
        pbn.bn.add(pbn.bp2node);
        return pbn;
    }

    /**
     * Determine if the BN can perform inference, which will only be the case if there is at least one
     * branchpoint which is connected to another.
     * @return true if valid, else false
     */
    public boolean isValid() {
        return this.getBN().getNodes().size() > 0;
    }

}
