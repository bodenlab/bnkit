package asr;

import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
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
 * Note that the BN represents a single position (one rate).
 */
public class MaxLhoodJoint implements TreeDecor<Object> {

    private SubstModel model = null;
    private SubstModel.ModelCache modelcache = null;
    final private IdxTree tree;
    private PhyloBN pbn = null;
    private Inference inf = null;


    /**
     * Set-up a joint reconstruction; this creates a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model.
     * @param tree
     * @param model
     * @param rate relative evolutionary rate
     */
    public MaxLhoodJoint(IdxTree tree, SubstModel model, double rate) {
        this.tree = tree;
        this.model = model;
        this.pbn = PhyloBN.create(tree, model, rate);
    }

    /**
     * Set-up a joint reconstruction; this creates a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and substitution model. Uses default rate.
     * @param tree
     * @param model
     */
    public MaxLhoodJoint(IdxTree tree, SubstModel model) {
        this.tree = tree;
        this.model = model;
        this.pbn = PhyloBN.create(tree, model);
    }

    /**
     * Set-up a joint reconstruction but hold-off with introducing what values are being optimised;
     * this prepares for a probabilistic graphical model with the structure and evolutionary
     * distances as defined by the supplied tree and (later supplied) substitution model.
     * @param tree
     * //@param rate relative evolutionary rate
     */
    public MaxLhoodJoint(IdxTree tree, SubstModel.ModelCache modelcache) {
        this.tree = tree;
        this.modelcache = modelcache;
    }

    /**
     * Set-up joint inference based directly on a PhyloBN instance (with defined/pre-trained nodes).
     * The BN always has nodes for the nodes of the phylogenetic tree, and if isExt() is true it has nodes for accessory variables.
     * Accessory variables are "tied" meaning that they share the parameters. Moreover, these variables may be confined to leaves on the
     * tree, or apply to all nodes of the tree.
     * @param pbn
     */
    public MaxLhoodJoint(PhyloBN pbn) {
        this.tree = pbn.getTree();
        this.pbn = pbn;
    }


    /**
     * Retrieve the inferred state for a specified branch point, as determined by max likelihood
     * @param idx branch point index
     * @return state which jointly with all others assigns the greatest probability to the observed values (at leaves)
     */
    @Override
    public Object getDecoration(int idx) {
        return inf.values[idx];
    }

    public long ELAPSED_TIME = 0;

    /**
     * Determine the joint state that assigns the maximum likelihood to the specified observations.
     * This is the bit that will be multi-threaded usually.
     * @param ti observed tree states
     */
    @Override
    public void decorate(TreeInstance ti) {
        long START_TIME = System.currentTimeMillis();
        if (model == null && !pbn.isExt())
            inf = infer(ti, true);
        else
            inf = infer(ti);
        ELAPSED_TIME += (System.currentTimeMillis() - START_TIME);
    }

    public String toElapsedTime() {
        return String.format("Elapsed %d ms",
                TimeUnit.MILLISECONDS.toMillis(ELAPSED_TIME));
    }

    /**
     * Determine the joint state that assigns the maximum likelihood to the specified observations.
     * This is the bit that will be multi-threaded usually.
     * @param ti observed tree states; these correspond to either the nodes of the phylogenetic tree, or to variables attached below them
     */
    public Inference infer(TreeInstance ti) {
        Inference myinf = new Inference(ti);            // setting up a structure to keep the values
        Map<Variable, Integer> quick = new HashMap<>();   // creating a map so that variables (whatever they may represent) can be linked to a specific index in the tree
        // instantiate all nodes for which there are values, i.e. leaf nodes most probably but not necessarily or exclusively
        for (int i = 0; i < ti.getSize(); i++) {
            myinf.values[i] = ti.getInstance(i);  // input value is always output value, so transferred to result here
            // we may infer two different values for the same index in the phylogenetic tree!
            // inferred values are in two places; in the main tree (always), ...
            if (pbn.isExt()) { // will have accessory (ext) variables
                BNode bnode = pbn.getExtNode(i); // but they are not always there (e.g. accessory vars for ancestors)
                if (bnode != null) {             // not hidden, so can be instantiated and inferred
                    quick.put(bnode.getVariable(), i);
                    if (myinf.values[i] != null)
                        bnode.setInstance(myinf.values[i]);
                    BNode parent_bnode = pbn.getBNode(i);
                    quick.put(parent_bnode.getVariable(), i);
                } else {       // additional complication is that some indices do NOT have accessory variables, so we rely on the main tree instead
                    bnode = pbn.getBNode(i);
                    if (bnode != null) {
                        quick.put(bnode.getVariable(), i);
                        if (myinf.values[i] != null)
                            bnode.setInstance(myinf.values[i]);
                    } // else, this branchpoint is outside of the BN
                }
            } else { // will NOT have accessory variables, so much more straightforward...
                BNode bnode = pbn.getBNode(i);
                if (bnode != null) {
                    quick.put(bnode.getVariable(), i);
                    bnode.setInstance(myinf.values[i]);
                } // else, this branchpoint is outside of the BN
            }
        }
        if (pbn.isValid()) {
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(pbn.getBN());
            Query q_mpe = ve.makeMPE();
            CGTable r1 = (CGTable) ve.infer(q_mpe);
            Variable.Assignment[] assign = r1.getMPE();
            for (Variable.Assignment assign1 : assign) {
                Integer idx = quick.get(assign1.var);
                if (idx != null) {
                    if (ti.getInstance(idx) == null) { // was not instantiated
                        myinf.values[idx] = assign1.val;
                    }
                } else {
                    System.out.println("Failed to fetch variable from recon: " + assign1.var);
                }
            }
        } // else the BN is incapable of performing inference, so just leave values as they are
        return myinf;
    }

    /**
     * Get the model for the specified values (size N).
     * Only to be used internally, if the model is NOT assigned through the constructor.
     * TODO: pool models so that they can be re-used; also consider optional model types (currently uniform JC-like for any N)
     * @param values
     * @return
     */
    private SubstModel getModel(Object[] values) {
        if (modelcache != null) {
            Integer tag = values.length;
            synchronized (modelcache) {
                SubstModel m = modelcache.get(tag);
                if (m == null) {
                    m = new JC(1, values);
                    modelcache.add(m, tag);
                    m.CACHE_SIZE = (int) (this.tree.getSize());
                }
                return m;
            }
        }
        return new JC(1, values);
    }

    /**
     * Determine the joint state that assigns the maximum likelihood to the specified observations.
     * This implementation of inference will convert values to a standard list and create a model for it internally.
     * This will be multi-threaded usually.
     * @param ti observed tree states
     * @param recodeNull if true, convert null at leaves to a distinct value to be optimised,
     *                   null at ancestors are kept since they may be queried;
     *                   if false, all null values are treated as uninstantiated/query nodes
     * @throws ASRRuntimeException if a model has been instantiated already
     */
    public Inference infer(TreeInstance ti, boolean recodeNull) {
        if (model != null) {
            boolean all_valid = true;
            for (Object sym : ti.getPossible()) {
                if (!model.getDomain().isValid(sym)) {
                    all_valid = false;
                    break;
                }
            }
            if (!all_valid)
                throw new ASRRuntimeException("Invalid inference setting: model is given, still recoding of values required");
            return infer(ti);
        }
        // construct a model suited to the values in the tree instance
        Object[] key = ti.encode(recodeNull); // convert instances to standardised list
        model = getModel(ti.getPossible());
        this.pbn = PhyloBN.create(tree, model);
        Inference myinf = new Inference(ti);
        Map<String, Integer> quick = new HashMap<>();
        // instantiate all nodes for which there are values, i.e. leaf nodes most probably but not necessarily.
        // also, leaf nodes do not need to be instantiated.
        for (int i = 0; i < ti.getSize(); i++) {
            myinf.values[i] = ti.getInstance(i);
            BNode bnode = pbn.getBNode(i);
            if (bnode != null) { // not hidden, so can be instantiated and inferred
                quick.put(bnode.getVariable().getName(), i);
                bnode.setInstance(myinf.values[i]);
            } // else, this branchpoint is outside of the BN, and will be ignored
        }
        if (pbn.isValid()) {
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(pbn.getBN());
            Query q_mpe = ve.makeMPE();
            CGTable r1 = (CGTable) ve.infer(q_mpe);
            Variable.Assignment[] assign = r1.getMPE();
            for (Variable.Assignment assign1 : assign) {
                Integer idx = quick.get(assign1.var.getName());
                if (idx != null) { // need to convert back to original alphabet
                    try {
                        myinf.values[idx] = key[(Integer)assign1.val];
                    } catch (ClassCastException e) {
                        throw new ASRRuntimeException("Invalid symbol inferred by ML: " + assign1.val);
                    }
                }
            }
        } // else the BN is incapable of performing inference, so just leave values as they are
        return myinf;
    }

    class Inference {
        final private Object[] values;
        final private IdxTree tree;
        public Inference(TreeInstance ti) {
            this.values = new Object[ti.getSize()];
            this.tree = ti.getTree();
        }
        @Override
        public String toString() {
            return getTreeInstance().toString();
        }

        public TreeInstance getTreeInstance() {
            return new TreeInstance(tree, values);
        }
    }

}
