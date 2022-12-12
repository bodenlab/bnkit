package dat.phylo;

import bn.BNet;
import bn.BNode;
import bn.Predef;
import bn.alg.EM;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.factor.FactorCache;
import bn.file.BNBuf;
import bn.node.CPT;
import bn.node.GDT;
import bn.prob.EnumDistrib;
import bn.prob.GaussianDistrib;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;
import dat.file.TSVFile;
import json.JSONArray;
import json.JSONObject;

import java.util.*;

/**
 * Class for a Bayesian network that represents branch points,
 * and their relationships as extracted from a phylogenetic tree;
 * branch points can also be extended to represent properties (discrete or continuous).
 * Based on the deprecated class PhyloBNet.
 *
 * @author mikael
 */
public class PhyloBN {

    private final BNet bn;
    private final IdxTree tree;

    /**
     * Default relative rate
     */
    public static double DEFAULT_RATE = 1.0;
    private BNode[] bp2node = null; // branchpoint index to BN node
    private BNode[] bp2ext = null;  // branchpoint index to extra node, e.g. GDT or CPT
    private String[] names = null; // name of node, that could be used to reference external or internal nodes by idx

    /**
     * Private constructor so that instances are created for a particular end.
     */
    private PhyloBN(IdxTree tree) {
        bn = new BNet();
        this.tree = tree;
    }

    private PhyloBN(BNet bn, IdxTree tree) {
        this.bn = bn;
        this.tree = tree;
    }

    /**
     * Get the BN
     *
     * @return Bayesian network instance
     */
    public BNet getBN() {
        return bn;
    }

    public void save(String filename) {
        BNBuf.save(this.getBN(), filename);
    }

    public static PhyloBN load(String filename, IdxTree tree) {
        BNet bn = BNBuf.load(filename);
        return new PhyloBN(bn, tree);
    }

    public IdxTree getTree() {
        return tree;
    }

    /**
     * @return
     * @deprecated
     */
    public BNode[] getBNodes() {
        return bp2node;
    }

    /**
     * Get the Bayesian network node for a given branchpoint index of the phylogenetic tree used to create the network.
     *
     * @param bpidx branchpoint index
     * @return the instance of the BNode if available; note that the network only represents a subset of the branchpoints and
     * if not available null is returned
     */
    public BNode getBNode(int bpidx) {
        return bp2node[bpidx];
    }

    /**
     * Get the Bayesian network GDT node for a given branchpoint index of the phylogenetic tree used to create the network.
     * Note this is NOT the branch point node that takes a discrete state but the GDT that is linked to it.
     *
     * @param bpidx branchpoint index
     * @return the instance of the BNode if available; note that the network only represents a subset of the branch points and
     * if not available null is returned
     */
    public BNode getExtNode(int bpidx) {
        return bp2ext[bpidx];
    }

    /**
     * Checks if this PhyloBN has been set up with external nodes, linked to the internal tree (whose leaves are latent)
     * @return true if external nodes are available
     */
    public boolean isExt() {
        return bp2ext != null;
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     * The BN will create nodes ONLY for variables that are accessible and valid in the tree.
     * Note that trees are sometimes pruned, which removes the need to include many variables.
     *
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
        PhyloBN pbn = new PhyloBN(tree);
        pbn.bp2node = new BNode[tree.getSize()];
        pbn.names = new String[tree.getSize()];
        //  create variables and nodes
        for (int idx : tree) { // iterate through tree depth-first
            if (tree.isConnected(idx)) { // only create variable for branchpoint if it is a parent or a child of one
                BranchPoint bp = tree.getBranchPoint(idx);
                EnumVariable rvar = new EnumVariable(model.getDomain(), bp.getID().toString());
                pbn.names[idx] = bp.getID().toString();
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

    private GDT gdt_master = null;
    private CPT cpt_master = null;

    /**
     * With the specified JSON encoding, override an existing definition of the so-called master node,
     * shared across branch points, extending as a CPT or GDT.
     * @param json
     */
    public void overrideMasterJSON(JSONObject json) {
        if (gdt_master != null) {
            List<EnumVariable> parents = gdt_master.getParents();
            if (parents.size() == 1) {
                EnumVariable par = parents.get(0);
                for (Object state : par.getDomain().getValues()) {
                    JSONObject component = json.getJSONObject(state.toString());
                    GaussianDistrib gd = GaussianDistrib.fromJSON(component);
                    gdt_master.put(new Object[] {state}, gd);
                }
            }
        } else if (cpt_master != null) {
            List<EnumVariable> parents = cpt_master.getParents();
            if (parents.size() == 1) {
                EnumVariable par = parents.get(0);
                for (Object state : par.getDomain().getValues()) {
                    JSONObject component = json.getJSONObject(state.toString());
                    EnumDistrib ed = EnumDistrib.fromJSON(component);
                    cpt_master.put(new Object[] {state}, ed);
                }
            }
        }
    }

    /**
     * Output a text string which describes the state of the master node of this PhyloBN, if available.
     * @return
     */
    public String toString() {
        if (gdt_master != null)
            return gdt_master.getStateAsText();
        else if (cpt_master != null)
            return cpt_master.getStateAsText();
        else
            return "Discrete";
    }

    /**
     * Function to produce a JSON representation of the so-called master node,
     * shared across branch points, extending as a CPT or GDT.
     * @return
     */
    public JSONObject getMasterJSON() {
        if (gdt_master != null) {
            List<EnumVariable> parents = gdt_master.getParents();
            if (parents.size() == 1) {
                EnumVariable par = parents.get(0);
                JSONObject json = new JSONObject();
                for (Object state : par.getDomain().getValues()) {
                    GaussianDistrib gd = (GaussianDistrib)gdt_master.getDistrib(new Object[] {state});
                    json.put(state.toString(), gd.toJSON());
                }
                return json;
            }
        } else if (cpt_master != null) {
            List<EnumVariable> parents = cpt_master.getParents();
            if (parents.size() == 1) {
                EnumVariable par = parents.get(0);
                JSONObject json = new JSONObject();
                for (Object state : par.getDomain().getValues()) {
                    EnumDistrib d = (EnumDistrib) cpt_master.getDistrib(new Object[] {state});
                    json.put(state.toString(), d.toJSON());
                }
                return json;
            }
        }
        return null;
    }

    /**
     * Create a PhyloBN based on the tree, but with GDTs hanging off each of the leaves, forming a BN that can take real
     * values at these extensions and infer latent discrete values at the internal branch point nodes, corresponding to both
     * extants and ancestors in the phylogenetic tree.
     *
     * @param tree
     * @param model
     * @param rate
     * @return
     */
    public static PhyloBN withGDTs(IdxTree tree, SubstModel model, double rate) {
        return withGDTs(tree, model, rate, true, System.currentTimeMillis());
    }

    /**
     * Create a PhyloBN based on the tree, but with GDTs hanging off each of the leaves, forming a BN that can take real
     * values at these extensions and infer latent discrete values at the internal branch point nodes, corresponding to both
     * extants and ancestors in the phylogenetic tree.
     *
     * @param tree
     * @param model
     * @param rate
     * @return
     */
    public static PhyloBN withGDTs(IdxTree tree, SubstModel model, double rate, long SEED) {
        return withGDTs(tree, model, rate, true, SEED);
    }

    /**
     * Create a PhyloBN based on the tree, but with GDTs hanging off each of the leaves, forming a BN that can take real
     * values at these extensions and infer latent discrete values at the internal branch point nodes, corresponding to both
     * extants and ancestors in the phylogenetic tree.
     *
     * @param tree
     * @param model
     * @param rate
     * @return
     */
    public static PhyloBN withGDTs(IdxTree tree, SubstModel model, double rate, boolean leavesOnly, long SEED) {
        PhyloBN pbn = new PhyloBN(tree);
        pbn.bp2node = new BNode[tree.getSize()];
        pbn.names = new String[tree.getSize()];
        // add extra nodes for GDTs
        pbn.bp2ext = new BNode[tree.getSize()];
        //  create variables and nodes
        GDT master = null;
        for (int idx : tree) { // iterate through tree depth-first
            if (tree.isConnected(idx)) { // only create variable for branchpoint if it is a parent or a child of one
                BranchPoint bp = tree.getBranchPoint(idx);
                EnumVariable rvar = new EnumVariable(model.getDomain(), bp.getID().toString());
                pbn.names[idx] = bp.getID().toString();
                int parent = tree.getParent(idx);
                if (parent < 0) { // there's no parent
                    pbn.bp2node[idx] = new SubstNode(rvar, model);
                } else { // there's a parent, and because the iterator is "depth-first", the parent must already have been created
                    EnumVariable prvar = (EnumVariable) pbn.bp2node[parent].getVariable();
                    pbn.bp2node[idx] = new SubstNode(rvar, prvar, model, bp.getDistance() * rate); // this is where relative rate is incorporated
                }
                pbn.bn.add(pbn.bp2node[idx]);
                // now add GDT if the branch point is a leaf...
                if (tree.isLeaf(idx) || !leavesOnly) {
                    Variable real = Predef.Real(bp.getID().toString() + "_Real");
                    GDT gdt = new GDT(real, rvar);
                    if (master == null)
                        master = gdt;
                    else
                        gdt.tieTo(master);
                    gdt.setTrainable(true);
                    gdt.randomize(idx + SEED);
                    pbn.bp2ext[idx] = gdt;
                    pbn.bn.add(pbn.bp2ext[idx]);
                }
            }
        }
        pbn.gdt_master = master;
        return pbn;
    }

    /**
     * Create a PhyloBN based on the tree, but with CPTs hanging off each of the leaves, forming a BN that can take discrete
     * values at these extensions and infer latent discrete values at the internal branch point nodes, corresponding to both
     * extants and ancestors in the phylogenetic tree.
     *
     * @param tree
     * @param model
     * @param rate
     * @return
     */
    public static PhyloBN withCPTs(IdxTree tree, SubstModel model, String[] alphas, double rate) {
        return withCPTs(tree, model, alphas, rate, true, System.currentTimeMillis());
    }
    /**
     * Create a PhyloBN based on the tree, but with CPTs hanging off each of the leaves, forming a BN that can take discrete
     * values at these extensions and infer latent discrete values at the internal branch point nodes, corresponding to both
     * extants and ancestors in the phylogenetic tree.
     *
     * @param tree
     * @param model
     * @param rate
     * @return
     */
    public static PhyloBN withCPTs(IdxTree tree, SubstModel model, String[] alphas, double rate, long SEED) {
        return withCPTs(tree, model, alphas, rate, true, SEED);
    }

    /**
     * Create a PhyloBN based on the tree, but with CPTs hanging off each of the leaves, forming a BN that can take discrete
     * values at these extensions and infer latent discrete values at the internal branch point nodes, corresponding to both
     * extants and ancestors in the phylogenetic tree.
     *
     * @param tree
     * @param model
     * @param rate
     * @return
     */
    public static PhyloBN withCPTs(IdxTree tree, SubstModel model, String[] alphas, double rate, boolean leavesOnly, long SEED) {
        PhyloBN pbn = new PhyloBN(tree);
        pbn.bp2node = new BNode[tree.getSize()];
        pbn.names = new String[tree.getSize()];
        // add extra nodes for GDTs
        pbn.bp2ext = new BNode[tree.getSize()];
        //  create variables and nodes
        CPT master = null;
        for (int idx : tree) { // iterate through tree depth-first
            if (tree.isConnected(idx)) { // only create variable for branchpoint if it is a parent or a child of one
                BranchPoint bp = tree.getBranchPoint(idx);
                EnumVariable rvar = new EnumVariable(model.getDomain(), bp.getID().toString());
                pbn.names[idx] = bp.getID().toString();
                int parent = tree.getParent(idx);
                if (parent < 0) { // there's no parent
                    pbn.bp2node[idx] = new SubstNode(rvar, model);
                } else { // there's a parent, and because the iterator is "depth-first", the parent must already have been created
                    EnumVariable prvar = (EnumVariable) pbn.bp2node[parent].getVariable();
                    pbn.bp2node[idx] = new SubstNode(rvar, prvar, model, bp.getDistance() * rate); // this is where relative rate is incorporated
                }
                pbn.bn.add(pbn.bp2node[idx]);
                // now add GDT if the branch point is a leaf...
                if (tree.isLeaf(idx) || !leavesOnly) {
                    EnumVariable dvar = Predef.Nominal(alphas, bp.getID().toString() + "_Discrete");
                    CPT cpt = new CPT(dvar, rvar);
                    if (master == null)
                        master = cpt;
                    else
                        cpt.tieTo(master);
                    cpt.setTrainable(true);
                    cpt.randomize(idx + SEED);
                    pbn.bp2ext[idx] = cpt;
                    pbn.bn.add(pbn.bp2ext[idx]);
                }
            }
        }
        pbn.cpt_master = master;
        return pbn;
    }

    /**
     * Train the so-called master node of the BN.
     * @param datasetWithHeaders representation straight from TSV file
     * @param seed random seed to reproduce stochastic training decisions
     */
    public void trainEM(TSVFile datasetWithHeaders, long seed) {
        trainEM(datasetWithHeaders.getHeaders(), datasetWithHeaders.getRows(), seed);
    }
    /**
     * Train the so-called master node of the BN.
     * @param labels names of features, which must match variables in the BN
     * @param data data matrix (rows represent samples, columns represent features)
     * @param seed random seed to reproduce stochastic training decisions
     */

    public void trainEM(String[] labels, Object[][] data, long seed) {
        Variable[] vars = new Variable[labels.length];
        for (int i = 0; i < labels.length; i++) {
            boolean found = false;
            for (int j = 0; j < names.length; j ++) {
                if (names[j].equals(labels[i])) {
                    vars[i] = bp2ext == null ? bp2node[j].getVariable() : bp2ext[j].getVariable();
                    found = true;
                    break;
                }
            }
            if (!found)
                throw new RuntimeException("Could not find label in BN");
        }
        EM em = new EM(bn);
        em.setEMOption(1);
        em.train(data, vars, seed);
    }

    /**
     * Construct a BN from a list of substitution nodes, assuming that they have been organised into a tree-like structure.
     * @param nodes
     * @return
     */
    public static PhyloBN createBarebone(SubstNode[] nodes) {
        PhyloBN pbn = new PhyloBN(null); // mark this BN as one that is not created from a tree
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

    /**
     * Enable cache for all nodes that can be cached.
     * @return number of nodes cache-enabled
     */
    public int setCache(FactorCache cache) {
        int cnt = 0;
        for (int i = 0; i < bp2node.length; i ++) {
            try { // presently, only SubstNode can be cached
                ((SubstNode)bp2node[i]).setCache(cache);
                if (((SubstNode)bp2node[i]).isCache())
                    cnt += 1;
            } catch (ClassCastException e) {
                // a node that is NOT cached
            }
        }
        return cnt;
    }
}
