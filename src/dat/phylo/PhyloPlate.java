package dat.phylo;

import api.JSONUtils;
import bn.BNet;
import bn.BNode;
import bn.TiedNode;
import bn.alg.EM;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.ctmc.matrix.JC;
import dat.EnumVariable;
import dat.Enumerable;
import dat.Variable;
import dat.file.TSVFile;
import json.JSONArray;
import json.JSONObject;

import java.util.ArrayList;
import java.util.List;


/**
 * This is a class representing a Bayesian network that centres on branch points in a phylogenetic tree.
 * Each branch point is essentially represented by a "plate" (a shared/repeated sub-BN),
 * which contains (a) nodes that reflect properties of the branch point (observed or latent states), and
 * links-to (b) nodes that are part of a tree with evolutionary distances determining their (latent) states.
 * Importantly, the parameters of the plates are shared between all nodes, enabling the (EM) learning of
 * universal principles between them.
 *
 * @author mikael
 *
 * TODO: in to- and fromJSON methods, enable "empty" fields to be left-out, allow non-listed fields to default to empty
 */
public class PhyloPlate {

    private final BNet bn; // the Bayesian network, sharing (potentially repeated) topology with tree, and including the plates
    private final IdxTree tree; // the tree used as input, deciding the topology of the backbone BN representing branchpoints in tree, and the reference points for the "plate"
    /**
     * Default relative rate
     */
    public static double DEFAULT_RATE = 1.0; // multiplier of evolutionary time
    private BNode[][] bp2node = null;   // branch point index to latent nodes (one set for each branch point; that set has one for each "tree")
    private Plate[] bp2plate = null;    // branch point index to plate instances (

    private SubstModel[] models = null; // substitution models, one for each tree/enumerable

    private String[] names = null; // name of node, that could be used to reference phylo- or plate nodes by idx

    public int EM_ROUNDS = 25;

    /**
     * Create a Bayesian network (BN) from a tree, preparing plates with branch point connectors for each enumerable data type (as identified by a template),
     * @param tree the tree with which to determine topology of BN
     * @param template the type of each substitution node, and prepares for plate nodes connected to them
     */
    public PhyloPlate(IdxTree tree, Modes template) {
        this.tree = tree;
        bn = new BNet();
        Enumerable[] latents = template.mode_types;
        bp2node = new BNode[tree.getSize()][latents.length]; // array with all the branchpoint nodes
        bp2plate = new Plate[tree.getSize()];
        names = new String[tree.getSize()];
        models = new SubstModel[template.mode_types.length];
        for (int i = 0; i < models.length; i ++)
            models[i] = new JC(1, template.mode_types[i]);
        //  create variables and nodes
        for (int idx : tree) { // iterate through tree depth-first
            if (tree.isConnected(idx)) { // only create variable for branchpoint if it is a parent or a child of one (should be all)
                BranchPoint bp = tree.getBranchPoint(idx);
                names[idx] = bp.getLabel().toString(); // get name of branchpoint
                Plate plate = new Plate(names[idx], template); // create a plate from modes spec, with branchpoint name as a prefix
                EnumVariable[] evars = plate.getParentVariables(); // retrieve plate's branchpoint variables, so we can connect them
                int parent_idx = tree.getParent(idx); //
                if (parent_idx < 0) { // there's no parent
                    for (int i = 0; i < evars.length; i ++)
                        bp2node[idx][i] = new SubstNode(evars[i], models[i]); // the branchpoint is a root, so the state is not conditioned
                } else { // there's a parent, and because the iterator is "depth-first", the parent must already have been created
                    for (int i = 0; i < evars.length; i ++) {
                        EnumVariable prvar = (EnumVariable) bp2node[parent_idx][i].getVariable(); // the branchpoint has a parent that provides a condition
                        bp2node[idx][i] = new SubstNode(evars[i], prvar, models[i], bp.getDistance()); // this is where relative rate needs to be incorporated (TODO)
                    }
                }
                for (int i = 0; i < evars.length; i ++) // add the branchpoint nodes to the actual BN
                    bn.add(bp2node[idx][i]);
                bp2plate[idx] = plate; // add the plate itself, but note that nodes (within it) need to be added separately to the plate...
            }
        }
    }

    public void addNode(BNode node) {
        this.bn.add(node);
    }

    /**
     * FIXME: this constructor is not fully implemented
     * Create a Bayesian network from a tree, based on cloning plates
     * @param tree
     * @param plate the nodes defining the plate including variables that specify the type of each substitution node
     */
    public PhyloPlate(IdxTree tree, Plate plate) {
        this.tree = tree;
        bn = new BNet();
        bp2node = new BNode[tree.getSize()][plate.bpcvars.length];
        bp2plate = new Plate[tree.getSize()];
        names = new String[tree.getSize()];
        models = new SubstModel[plate.bpcvars.length];
        for (int i = 0; i < plate.bpcvars.length; i ++)
            models[i] = new JC(1, plate.modetypes[i]);
        //  create variables and nodes
        for (int idx : tree) { // iterate through tree depth-first
            if (tree.isConnected(idx)) { // only create variable for branchpoint if it is a parent or a child of one
                BranchPoint bp = tree.getBranchPoint(idx);
                names[idx] = bp.getID().toString();
                Plate plate_local = plate.createAnother(names[idx]); // <--------
                EnumVariable[] evars = plate_local.getParentVariables();
                int parent_idx = tree.getParent(idx);
                if (parent_idx < 0) { // there's no parent
                    for (int i = 0; i < evars.length; i ++)
                        bp2node[idx][i] = new SubstNode(evars[i], models[i]);
                } else { // there's a parent, and because the iterator is "depth-first", the parent must already have been created
                    for (int i = 0; i < evars.length; i ++) {
                        EnumVariable prvar = (EnumVariable) bp2node[parent_idx][i].getVariable();
                        bp2node[idx][i] = new SubstNode(evars[i], prvar, models[i], bp.getDistance()); // this is where relative rate is incorporated
                    }
                }
                for (int i = 0; i < evars.length; i ++)
                    bn.add(bp2node[idx][i]);
                // add plate
                bp2plate[idx] = plate_local;
            }
        }
    }

    /**
     * The "Modes" class specifies the latent connection between "plates" and the backbone that is governed by evolutionary models.
     * It essentially represents the type of discrete states in preparation for creating plates.
     */
    public static class Modes {
        Enumerable[] mode_types;
        public Modes(Enumerable[] mode_types) {
            this.mode_types = mode_types;
        }

        public JSONObject toJSON() {
            JSONObject json = new JSONObject();
            JSONArray jtypes = new JSONArray();
            for (Enumerable e : mode_types)
                jtypes.put(e.toJSON());
            json.put("Modetypes", jtypes);
            return json;
        }
        public static Modes fromJSON(JSONObject json) {
            JSONArray jtypes = json.getJSONArray("Modetypes");
            Enumerable[] mode_types = new Enumerable[jtypes.length()];
            for (int i = 0; i < mode_types.length; i ++)
                mode_types[i] = Enumerable.fromJSON(jtypes.getJSONObject(i));
            return new Modes(mode_types);
        }
    }

    /**
     * The "Plate" class defines an individual plate; this is first created and attached to each branchpoint
     */
    public static class Plate {

        String name;

        // The types of modes--connecting the plate with the evolutionary inference; the order determines the order of their branch point connector variables (below)
        Enumerable[] modetypes = null;
        // The branch point connector (BPC) variables; the order determines the order in which they connect with BNodes below
        EnumVariable[] bpcvars = null;

        // The nodes within the plate, some of which will have parent variables that are BPC variables
        BNode[] bnodes = null;
        // branch point connector indices [bpc][bnode]
        int[][] bpcidxs = null;
        // All variables that can be assigned values (via learning say)
        Variable[] vars = null;

        private Plate(String name, Enumerable[] modetypes) {
            this(name, modetypes, extractNames(name, modetypes));
        }

        private Plate(String name, Enumerable[] modetypes, String[] names) {
            this.name = name;
            this.modetypes = modetypes;
            this.bpcvars = new EnumVariable[modetypes.length];
            for (int i = 0; i < modetypes.length; i ++)
                this.bpcvars[i] = new EnumVariable(modetypes[i], names[i]);
            this.bnodes = new BNode[0];
        }

        /**
         * Create an empty plate but allocate the branch point connector variables, which
         * prepares for substitution nodes to be created outside of the plate and for plate nodes
         * to be created inside.
         * @param name
         * @param template
         */
        public Plate(String name, Modes template) {
            this(name, template.mode_types);
        }

        @Override
        public String toString() {
            StringBuffer sb = new StringBuffer();
            for (int i = 0; i < modetypes.length; i ++)
                sb.append(this.bpcvars[i].getName() + ((i < modetypes.length - 1) ? ";" : ""));
            return name + "[" +sb.toString() + "]";
        }

        public EnumVariable[] getParents(int[] idxs) {
            EnumVariable[] ret = new EnumVariable[idxs.length];
            for (int i = 0; i < idxs.length; i ++)
                ret[i] = bpcvars[idxs[i]];
            return ret;
        }

        public BNode[] getNodes() {
            return bnodes;
        }

        private static String prefix = "__";
        private String replacePrefix(String nameWithPrefix) {
            String platename = this.name;
            int shift = nameWithPrefix.indexOf(prefix);
            if (shift >= 0)
                return platename + "__" + nameWithPrefix.substring(shift + prefix.length());
            else
                return platename + "__" + nameWithPrefix;
        }

        public String expungePrefix(String nameWithPrefix) {
            String myprefix = this.name + prefix;
            int shift = nameWithPrefix.indexOf(myprefix);
            if (shift >= 0)
                return nameWithPrefix.substring(shift + myprefix.length());
            else
                return nameWithPrefix;
        }

        /**
         * Add a node to extend the plate.
         * Presently the node needs to have one or more branchpoint variables as parents.
         * @param node the node to add
         */
        public synchronized void addNode(BNode node, PhyloPlate bn) {
            addNode(node);
            bn.addNode(node);
        }
        /**
         * Add a node to extend the plate.
         * Presently the node needs to have one or more branchpoint variables as parents.
         * @param node the node to add
         */
        public synchronized void addNode(BNode node) {
            Variable var = node.getVariable();
            var.setName(replacePrefix(var.getName()));
            List<EnumVariable> parent_vars = node.getParents();
            int[] idxs = new int[parent_vars.size()];
            for (int i = 0; i < parent_vars.size(); i ++) {
                for (int j = 0; j < bpcvars.length; j++) {
                    if (bpcvars[j].equals(parent_vars.get(i))) {
                        idxs[i] = j;
                        break;
                    }
                }
            }
            // record keeping: update the arrays holding pointers to nodes and the indices to their parents
            BNode[] copy_nodes = new BNode[bnodes.length + 1];
            int[][] copy_idxs = new int[bnodes.length + 1][];
            for (int i = 0; i < bnodes.length; i ++) {
                copy_nodes[i] = bnodes[i];
                copy_idxs[i] = bpcidxs[i];
            } // old ones done, ... next the new one node and parents of it
            bnodes = copy_nodes;
            bnodes[bnodes.length - 1] = node;
            bpcidxs = copy_idxs;
            bpcidxs[bnodes.length - 1] = idxs;
        }

        public synchronized Plate createAnother(String name) {
            Plate clone = new Plate(name, this.modetypes);
            BNode[] newnodes = new BNode[this.bnodes.length];
            for (int i = 0; i < this.bnodes.length; i ++) {
                /*
                if (this.bnodes[i] instanceof bn.node.CPT) {
                    Variable var = bnodes[i].getVariable();
                    newnodes[i] = new CPT(var, );

                } else if (this.bnodes[i] instanceof bn.node.GDT) {

                }
                 */
            }
            return clone;
        }

        private static String[] extractNames(String name, Enumerable[] modes) {
            String[] names = new String[modes.length];
            for (int i = 0; i < modes.length; i ++)
                names[i] = name + "_" + (i + 1);
            return names;
        }
        public static Enumerable[] extractEnumerable(EnumVariable[] vars) {
            Enumerable[] modes = new Enumerable[vars.length];
            for (int i = 0; i < vars.length; i ++)
                modes[i] = vars[i].getDomain();
            return modes;
        }

        public EnumVariable[] getParentVariables() {
            return bpcvars;
        }

        public Variable[] getVariables() {
            if (vars == null && bnodes != null) {
                vars = new Variable[bnodes.length];
                for (int i = 0; i < bnodes.length; i ++)
                    vars[i] = bnodes[i].getVariable();
            }
            return vars;
        }

        public void setInstances(Object[] values) {
            if (bnodes != null && values.length == bnodes.length) {
                for (int i = 0; i < bnodes.length; i++)
                    bnodes[i].setInstance(values[i]);
            } else
                throw new RuntimeException("Invalid instance array for plate");
        }

        public void setMaster(Plate master) {
            BNode[] masternodes = master.getNodes();
            if (masternodes.length != bnodes.length)
                throw new RuntimeException("Master plate not compatible with current plate");
            for (int i = 0; i < bnodes.length; i ++) {
                try {
                    if (!((TiedNode)bnodes[i]).tieTo(masternodes[i]))
                        throw new RuntimeException("Current plate cannot be mastered");
                } catch (ClassCastException e) {
                    throw new RuntimeException("Current plate cannot be mastered");
                }
            }
        }

        /**
         * Generate JSON for this plate.
         * @return JSON representation of this plate, including all the nodes with parameter values
         */
        public JSONObject toJSON() {
            JSONObject json = new JSONObject();
            json.put("Name", name);
            JSONArray jtypes = new JSONArray();
            for (Enumerable e : this.modetypes)
                jtypes.put(e.toJSON());
            json.put("Modetypes", jtypes);
            JSONArray jconn = new JSONArray();
            for (int[] targetnodes : this.bpcidxs) {
                JSONArray jtarg = new JSONArray();
                for (int tnode : targetnodes)
                    jtarg.put(tnode);
                jconn.put(jtarg);
            }
            json.put("Targets", jconn);
            JSONArray jarray = new JSONArray();
            for (BNode bnode : bnodes)
                jarray.put(bnode.toJSON());
            json.put("Nodes", jarray);
            return json;
        }

        /**
         * Create a plate from a JSON representation, including nodes with parameter values
         * @param json
         * @return new instance of this class
         */
        public static Plate fromJSON(JSONObject json) {
            // json could describe Modes (with variables) or a Plate (with parameterised nodes)
            String name = json.getString("Name");
            JSONArray jtypes = json.getJSONArray("Modetypes");
            JSONArray jconn = json.getJSONArray("Targets");
            Enumerable[] tree_types = new Enumerable[jtypes.length()];
            for (int i = 0; i < tree_types.length; i ++)
                tree_types[i] = Enumerable.fromJSON(jtypes.getJSONObject(i));
            int[][] bpcidxs = new int[jconn.length()][];
            for (int i = 0; i < bpcidxs.length; i ++) {
                JSONArray jtarg = jconn.getJSONArray(i);
                bpcidxs[i] = new int[jtarg.length()];
                for (int j = 0; j < jtarg.length(); j ++)
                    bpcidxs[i][j] = jtarg.getInt(j);
            }
            Plate plate = new Plate(name, tree_types);
            JSONArray jnodes = json.optJSONArray("Nodes");
            if (jnodes != null) {
                BNode[] bnodes = new BNode[jnodes.length()];
                for (int i = 0; i < jnodes.length(); i++) {
                    bnodes[i] = BNode.fromJSON(jnodes.getJSONObject(i), plate.getParents(bpcidxs[i]));
                    plate.addNode(bnodes[i]);
                }
            }
            return plate;
        }

    }

    /**
     * Get the BN
     *
     * @return Bayesian network instance
     */
    public BNet getBN() {
        return bn;
    }

    public int getNModels() {
        return models.length;
    }
    
    public IdxTree getTree() {
        return tree;
    }

    /**
     * Get the plate (sub-BN) for a given branchpoint index of the phylogenetic tree used to create the network.
     *
     * @param bpidx branchpoint index
     * @return the instance of the plate if available; note that the network may only represent a subset of the branch points and
     * if not available null is returned
     */
    public Plate getPlate(int bpidx) {
        return bp2plate[bpidx];
    }

    /**
     * Get the nodes at a given branchpoint index of the phylogenetic tree used to create the network.
     *
     * @param bpidx branchpoint index
     * @return the latent nodes at the branchpoint, if not available, null is returned
     */
    public BNode[] getBNodes(int bpidx) {
        if (bp2node[bpidx][0] == null)
            return null;
        return bp2node[bpidx];
    }

    private Plate master = null;

    public void setMaster(Plate master) {
        this.master = master;
    }

    public Plate getMaster() {
        return master;
    }

    public void overrideMaster(JSONObject json) {
        if (master != null) {
            /*
            GDT replace = GDT.fromJSON(json, gdt_master.getVariable(), gdt_master.getParents());
            List<EnumVariable> parents = gdt_master.getParents();
            if (parents.size() == 1) {
                EnumVariable par = parents.get(0);
                for (Object state : par.getDomain().getValues()) {
                    Object[] cond = new Object[] {state};
                    gdt_master.put(cond, replace.getDistrib(cond));
                }
            }
             */
        }
    }

    /**
     * Function to produce a JSON representation of the so-called master node,
     * shared across branch points, extending as a CPT or GDT.
     * @return
     */
    public JSONObject getMasterJSON() {
        if (master != null) {
            return master.toJSON();
        }
        return null;
    }




    /**
     * Train the so-called master plate of the BN.
     * @param labels names of features, which correspond to variables in the BN
     * @param data data tensor (first dim represents samples, second represents items, third represents features)
     * @param seed random seed to reproduce stochastic training decisions
     */
    // note: not the same signature as in PhyloBN; only ONE sample possible here, with rows corresponding to entry (ONE plate)
    public void trainEM(String[] labels, Object[][] data, long seed) {
        List<Variable> varlist = new ArrayList<>();
        List<Object> datlist = new ArrayList<>();
        for (int i = 0; i < labels.length; i ++) {
            for (int j = 0; j < names.length; j ++) {
                if (names[j].equals(labels[i])) {                       // found plate!
                    if (data[i].length == bp2plate[j].bnodes.length) {  // same number of values as nodes in the plate
                        for (int n = 0; n < bp2plate[j].bnodes.length; n++) {
                            varlist.add(bp2plate[j].bnodes[n].getVariable());
                            datlist.add(data[i][n]);
                        }
                        break;
                    }
                }
            }
        }
        Object[][] dats = new Object[1][datlist.size()];
        datlist.toArray(dats[0]);
        Variable[] vars = new Variable[varlist.size()];
        varlist.toArray(vars);
        EM em = new EM(bn);
        em.setMaxRounds(EM_ROUNDS);
        em.setEMOption(1);
        em.train(dats, vars, seed);
    }
    /**
     * Train the so-called master plate of the BN.
     * @param dataset the dataset that is used to train the BN
     * @param seed random seed to reproduce stochastic training decisions
     */
    public void trainEM(JSONUtils.DataSet dataset, long seed) {
        List<Variable> varlist = new ArrayList<>();
        // dataset contains: names of items, which correspond to branch points in the tree
        // dataset contains: names of features, which expand to variables in the BN
        // dataset contains: data tensor (first dim represents samples, second represents items, third represents features)
        Object[][] data = dataset.getFlattenedData();
        String[] headers = dataset.getFlattenedHeaders();
        Variable[] vars = new Variable[headers.length];
        for (int i = 0; i < headers.length; i ++) {
            BNode node = bn.getNode(headers[i]);
            if (node != null)
                vars[i] = node.getVariable();
            else
                System.err.println("Missing variable for " + headers[i]);
        }
        EM em = new EM(bn);
        em.setMaxRounds(EM_ROUNDS);
        em.setEMOption(1);
        em.train(data, vars, seed);
    }
    /**
     * Train the so-called master plate of the BN.
     * @param dataset the dataset that instantiates the BN
     * @param sampleidx the sample in the dataset that is used
     */
    public void instantiate(JSONUtils.DataSet dataset, int sampleidx) {
        // dataset contains: names of items, which correspond to branch points in the tree
        // dataset contains: names of features, which expand to variables in the BN
        // dataset contains: data tensor (first dim represents samples, second represents items, third represents features)
        Object[][] data = dataset.getFlattenedData();
        String[] headers = dataset.getFlattenedHeaders();
        for (int i = 0; i < headers.length; i ++) {
            BNode node = bn.getNode(headers[i]);
            if (node != null)
                node.setInstance(data[sampleidx][i]);
        }
    }


}
