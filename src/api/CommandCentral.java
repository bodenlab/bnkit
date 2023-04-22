package api;

import asr.*;
import bn.BNode;
import bn.Distrib;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.Variable;
import dat.file.TSVFile;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.PhyloPlate;
import dat.phylo.TreeInstance;
import dat.pog.POGTree;
import dat.pog.POGraph;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.util.*;

public class CommandCentral {

    final private GServer server;

    /**
     * Create a command central, which parses requests formatted as JSON, and creates a GRequest.
     * The GRequest extracts parameters for the job in its constructor, and has a "run" method, which
     * will execute as a thread, once the job has been scheduled by the server.
     * @param server
     */
    public CommandCentral(GServer server) {
        this.server = server;
    }

    /**
     * The server passes on a request to this method, which in turn decides what specific
     * GRequest class that should be instantiated as a result.
     * If you write a new command, this is where you need to add a hook; that hook is
     * usually a new definition of a GRequest further below.
     * @param json
     * @return
     */
    public GRequest createRequest(JSONObject json) {
        GRequest request = null;
        try {
            String command = json.getString("Command");
            String authtoken = json.getString("Auth");
            JSONObject jparams = json.optJSONObject("Params");

            if (command.equals("Recon")) {
                request = new GRequest_Recon(command, authtoken, jparams);

            } else if (command.equals("Pogit")) {
                request = new GRequest_Pogit(command, authtoken, jparams);

            } else if (command.equals("Label-tree")) {
                request = new GRequest_LabelTree(command, authtoken, jparams);

            } else if (command.equals("Train")) {
                request = new GRequest_Train(command, authtoken, jparams);

            } else if (command.equals("TrainModes")) {
                request = new GRequest_TrainModes(command, authtoken, jparams);

            } else if (command.equals("InferModes")) {
                request = new GRequest_InferModes(command, authtoken, jparams);

            } else if (command.equals("Infer")) {
                request = new GRequest_InferBP(command, authtoken, jparams);

            } else if (command.equals("Fake")) {
                request = new GRequest_Fake(command, authtoken, jparams);

            } else {

            }
        } catch (GRequestRuntimeException syntax) {
            throw new GRequestRuntimeException("Invalid request: " + syntax.getMessage());
        } catch (JSONException e) {
            e.printStackTrace();
            throw new GRequestRuntimeException("Invalid request: missing Command or Auth");
        }
        return request;
    }


    /**
     * Definition of the "Pogit" command.
     *
     * "Pogit" combines a tree and an alignment to construct a POGTree, i.e.
     * POGs representing sequences in alignment, attached to leaves of the tree (with congruent labels)
     */
    public static class GRequest_Pogit extends GRequest {
        private IdxTree idxTree = null;
        private EnumSeq.Alignment aln = null;

        public GRequest_Pogit(String command, String auth, JSONObject params) {
            super(command, auth);
            // extract necessary detail from command params etc
            try {
                JSONObject tree = params.getJSONObject("Tree");
                idxTree = IdxTree.fromJSON(tree);
                JSONObject jaln = params.getJSONObject("Alignment");
                aln = EnumSeq.Alignment.fromJSON(jaln);
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command : " + command + "; " + e.getMessage());
            }
        }

        @Override
        public boolean isQueued() {
            return false; // run immediately
        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public void run() {
            POGTree pogtree = new POGTree(aln, idxTree);
            this.setResult(pogtree.toJSON());
        }
    }
    /**
     * Definition of the "Recon" command.
     *
     * "Recon" uses a tree and alignment to perform a reconstruction, under a range of conditions
     * specified by a set of parameters.
     */
    public static class GRequest_Recon extends GRequest {
        private IdxTree idxTree = null;
        private EnumSeq.Alignment aln = null;
        private GRASP.Inference MODE = null;
        private int[] ancestors = null; // if marginal inference, a list of ancestors of at least one
        private Enumerable alpha = null;
        private SubstModel MODEL = null;
        String[] INDELS = new String[] {"BEP", "BEML", "SICP", "SICML", "PSP", "PSML"};
        int INDEL_IDX = 0; // default indel approach is that above indexed 0
        private POGTree pogTree = null;
        private Prediction indelpred = null;
        private double[] RATES = null;
        private long START_TIME;

        public GRequest_Recon(String command, String auth, JSONObject params) {
            super(command, auth);
            // extract necessary detail from command params etc
            try {
                pogTree = POGTree.fromJSON(params);
                idxTree = pogTree.getTree();
            } catch (ASRRuntimeException e1) { // incomplete pogtree
                try {
                    JSONObject tree = params.getJSONObject("Tree");
                    idxTree = IdxTree.fromJSON(tree);
                    JSONObject jaln = params.getJSONObject("Alignment");
                    aln = EnumSeq.Alignment.fromJSON(jaln);
                    pogTree = new POGTree(aln, idxTree);
                } catch (JSONException e2) {
                    throw new GRequestRuntimeException("Invalid JSON in command: " + command + "; " + e2.getMessage());
                } catch (RuntimeException e3) {
                    throw new GRequestRuntimeException("Invalid parameters in command: " + command + "; " + e3.getMessage());
                }
            }
            try {
                String infmode = params.optString("Inference", "Joint");
                MODE = infmode.equals("Joint") ? GRASP.Inference.JOINT : (infmode.equals("Marginal") ? GRASP.Inference.MARGINAL : null);
                if (MODE == GRASP.Inference.MARGINAL) {
                    try {
                        Integer ancspec1 = params.optInt("Ancestor");
                        if (ancspec1 == null) {
                            JSONArray ancspec2 = params.optJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i++)
                                ancestors[i] = ancspec2.getInt(i);
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = ancspec1;
                        }
                    } catch (ClassCastException e) { // ancestor IDs specified with "N" prefix
                        String ancspec1 = params.optString("Ancestor");
                        if (ancspec1 == null) {
                            JSONArray ancspec2 = params.optJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i++)
                                ancestors[i] = Integer.parseInt(ancspec2.getString(i).substring(1));
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = Integer.parseInt(ancspec1.substring(1));
                            ;
                        }
                    }
                }
                String indels = params.optString("Indels", "BEP");
                for (int i = 0; i < INDELS.length; i++) {
                    if (INDELS[i].equalsIgnoreCase(indels)) {
                        INDEL_IDX = i;
                        break;
                    }
                }
                // other params:
                // model (default is JTT)
                String modelname = params.optString("Model", "JTT");
                MODEL = SubstModel.createModel(modelname);
                if (MODEL == null)
                    throw new ASRRuntimeException("Invalid model");
//              TODO: more user data error checking
//                if (MODEL.getDomain().)
                // rates (optional)
                JSONArray jrates = params.optJSONArray("Rates");
                if (jrates != null) {
                    RATES = new double[jrates.length()];
                    for (int i = 0; i < RATES.length; i++)
                        RATES[i] = jrates.getDouble(i);
                }
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command : " + command + "; " + e.getMessage());
            }
            //
            // TODO: figure out compute resources, threads, memory, priority, queueing strategy etc.
            //
        }

        @Override
        public boolean isQueued() {
            return true;
        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public void run() {
            START_TIME = System.currentTimeMillis();
            switch (INDEL_IDX) {
                case 0:
                    indelpred = Prediction.PredictByBidirEdgeParsimony(pogTree);
                    break;
                case 1:
                    indelpred = Prediction.PredictByBidirEdgeMaxLhood(pogTree);
                    break;
                case 2:
                    indelpred = Prediction.PredictBySICP(pogTree);
                    break;
                case 3:
                    indelpred = Prediction.PredictBySICML(pogTree);
                    break;
                case 4:
                    indelpred = Prediction.PredictByParsimony(pogTree);
                    break;
                case 5:
                    indelpred = Prediction.PredictByMaxLhood(pogTree);
                    break;
                default:
                    break;
            }
            if (MODE == GRASP.Inference.JOINT)
                indelpred.getJoint(MODEL, RATES);
            else if (MODE == GRASP.Inference.MARGINAL) {
                for (int i = 0; i < ancestors.length; i ++) {
                    if (indelpred.getTree().getIndex(ancestors[i]) < 0)
                        ; // System.out.println(ancestors[i] + " is not a valid ancestor number");
                    else
                        indelpred.getMarginal(ancestors[i], MODEL, RATES);
                }
            }
            Map<Object, POGraph> pogs = indelpred.getAncestors(MODE);
            this.setResult(indelpred.toJSONJustAncestors());
        }
    }

    public static class GRequest_LabelTree extends GRequest {
        private final IdxTree idxTree;

        public GRequest_LabelTree(String command, String auth, JSONObject params) {
            super(command, auth);
            try {
                JSONObject tree = params.getJSONObject("Tree");
                idxTree = IdxTree.fromJSON(tree);
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command: " + command + "; " + e.getMessage());
            }
        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public boolean isQueued() {
            return false;
        }

        @Override
        public void run() {
            idxTree.getBranchPoint(0).setInternalLabels(0);
            JSONObject myres = new JSONObject();
            myres.put("Tree", idxTree.toJSON());
            this.setResult(myres);
        }
    }

    public static class GRequest_Train extends GRequest {

        private final IdxTree idxTree;
        private final JSONUtils.DataSet dataset;
        private enum DATATYPE {CONTINUOUS, ENUMERABLE};
        private final DATATYPE TIP_TYPE;
        private PhyloBN pbn = null;
        private SubstModel MODEL = null;
        private int OBSERVED_NSTATES;
        private String[] OBSERVED_ALPHA = null;
        private Object[] LATENT_STATES = null;
        private Double GAMMA = 1.0;
        private Double RATE = 1.0;
        private Long SEED = 1L;
        private Boolean LEAVES_ONLY = true;

        public GRequest_Train(String command, String auth, JSONObject params) {
            super(command, auth);
            try {
                JSONObject tree = params.getJSONObject("Tree");
                idxTree = IdxTree.fromJSON(tree);
                JSONObject jdataset = params.getJSONObject("Dataset");
                dataset = JSONUtils.DataSet.fromJSON(jdataset);
                JSONObject TIP_DISTRIB = params.optJSONObject("Distrib");
                JSONArray JSTATES = params.optJSONArray("States");
                if (JSTATES != null) {
                    LATENT_STATES = new Object[JSTATES.length()];
                    for (int i = 0; i < JSTATES.length(); i ++)
                        LATENT_STATES[i] = JSTATES.get(i);
                }
                LEAVES_ONLY = params.optBoolean("Leaves-only", LEAVES_ONLY);
                SEED = params.optLong("Seed", SEED);
                RATE = params.optDouble("Rate", RATE);
                GAMMA = params.optDouble("Gamma", GAMMA);
                if (TSVFile.isDoubleOrInt(dataset.getNonitemisedData())) { // real values for tip nodes
                    TIP_TYPE = DATATYPE.CONTINUOUS;
                    if (LATENT_STATES != null) {
                        MODEL = new JC(GAMMA, LATENT_STATES);
                    } else
                        throw new GRequestRuntimeException("Latent states are invalid for " + command + "; States are " + (JSTATES != null ? JSTATES : "not given"));
                } else { // the column has discrete values (not real)
                    TIP_TYPE = DATATYPE.ENUMERABLE;
                    Set observed = new HashSet(new TSVFile(dataset.getNonitemisedData()).getValues());
                    OBSERVED_NSTATES = observed.size();
                    OBSERVED_ALPHA = new String[observed.size()];
                    observed.toArray(OBSERVED_ALPHA);
                    if (TIP_DISTRIB == null && LATENT_STATES == null) { // no extra nodes
                        MODEL = new JC(GAMMA, OBSERVED_ALPHA);
                    } else if (LATENT_STATES != null) {
                        MODEL = new JC(GAMMA, LATENT_STATES);
                    } else
                        throw new GRequestRuntimeException("Latent states are invalid for " + command + "; States are " + (JSTATES != null ? JSTATES : "not given"));
                }

                if (TIP_TYPE == DATATYPE.CONTINUOUS) {
                    pbn = PhyloBN.withGDTs(idxTree, MODEL, RATE, LEAVES_ONLY, SEED);
                    if (TIP_DISTRIB != null) // GDT given?
                        pbn.overrideMasterJSON(TIP_DISTRIB);
                } else { // discrete states only
                    if (TIP_DISTRIB == null && LATENT_STATES == null) { // no extra nodes
                        pbn = PhyloBN.create(idxTree, MODEL, RATE);
                    } else {
                        pbn = PhyloBN.withCPTs(idxTree, MODEL, OBSERVED_ALPHA, RATE, LEAVES_ONLY, SEED);
                        if (TIP_DISTRIB != null) // CPT given?
                            pbn.overrideMasterJSON(TIP_DISTRIB);
                    }
                }
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command: " + command + "; " + e.getMessage());
            } catch (JSONUtils.JSONUtilsException e) {
                throw new GRequestRuntimeException("Invalid format: " + command + "; " + e.getMessage());
            }
            //
            // TODO: figure out compute resources, threads, memory, priority, queueing strategy etc.
            //

        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public boolean isQueued() {
            return true;
        }

        @Override
        public void run() {
            pbn.trainEM(dataset.getFeatures(), dataset.getNonitemisedData(), SEED);
            JSONObject myres = new JSONObject();
            myres.put("Distrib", pbn.getMasterJSON());
            this.setResult(myres);
        }
    }
    public static class GRequest_TrainModes extends GRequest {

        private final IdxTree idxTree;
        private final JSONUtils.DataSet dataset;
        private PhyloPlate pbn = null;
        private int EM_ROUNDS = 50;
        private Double GAMMA = 1.0;
        private Double RATE = 1.0;
        private Long SEED = 1L;
        private Boolean LEAVES_ONLY = true;

        public GRequest_TrainModes(String command, String auth, JSONObject params) {
            super(command, auth);
            try {
                JSONObject tree = params.getJSONObject("Tree");
                idxTree = IdxTree.fromJSON(tree);
                JSONObject jdataset = params.getJSONObject("Dataset");
                dataset = JSONUtils.DataSet.fromJSON(jdataset);
                JSONObject jdistrib = params.getJSONObject("Distrib");
                JSONArray modetypes = jdistrib.getJSONArray("Modetypes");
                Enumerable[] modes = new Enumerable[modetypes.length()];
                for (int i = 0; i < modes.length; i ++)
                    modes[i] = Enumerable.fromJSON(modetypes.getJSONObject(i));
                PhyloPlate.Modes template = new PhyloPlate.Modes(modes);
                //System.out.println("Created template: " + template.toJSON());
                LEAVES_ONLY = params.optBoolean("Leaves-only", LEAVES_ONLY);
                SEED = params.optLong("Seed", SEED);
                RATE = params.optDouble("Rate", RATE);
                GAMMA = params.optDouble("Gamma", GAMMA);
                EM_ROUNDS = params.optInt("Rounds", EM_ROUNDS);
                pbn = new PhyloPlate(idxTree, template);
                PhyloPlate.Plate master = null;
                for (int idx : idxTree) {
                    if (LEAVES_ONLY && !idxTree.isLeaf(idx))
                        continue;
                    //System.out.println("Adding plate to branchpoint " + idx);
                    PhyloPlate.Plate plate = pbn.getPlate(idx);
                    JSONArray jnodes = jdistrib.getJSONArray("Nodes");
                    JSONArray jtargets = jdistrib.getJSONArray("Targets");
                    BNode[] nodes = new BNode[jnodes.length()];
                    for (int j = 0; j < nodes.length; j ++) {
                        JSONArray jtargnode = jtargets.getJSONArray(j);
                        int[] targets = new int[jtargnode.length()];
                        for (int k = 0; k < targets.length; k ++)
                            targets[k] = jtargnode.getInt(k);
                        plate.addNode(BNode.fromJSON(jnodes.getJSONObject(j), plate.getParents(targets)), pbn);
                    }
                    if (master == null) {
                        master = plate;
                        pbn.setMaster(master);
                    } else {
                        plate.setMaster(master);
                    }
                    //System.out.println("\tPlate: " + plate.toJSON());
                }
                dataset.curateFeatures(pbn.getBN());
                pbn.EM_ROUNDS = EM_ROUNDS;
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command: " + command + "; " + e.getMessage());
            } catch (JSONUtils.JSONUtilsException e) {
                throw new GRequestRuntimeException("Invalid format: " + command + "; " + e.getMessage());
            } catch (RuntimeException e) {
                throw new GRequestRuntimeException("Technical problem: " + command + "; " + e.getMessage());
            }
            //
            // TODO: figure out compute resources, threads, memory, priority, queueing strategy etc.
            //

        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public boolean isQueued() {
            return true;
        }

        @Override
        public void run() {
            pbn.trainEM(dataset, SEED);
            JSONObject myres = new JSONObject();
            myres.put("Distrib", pbn.getMasterJSON());
            this.setResult(myres);
        }
    }
    public static class GRequest_InferModes extends GRequest {

        private final IdxTree idxTree;
        private final JSONUtils.DataSet dataset;
        private PhyloPlate pbn = null;
        private Object[] LATENT_STATES = null;
        private Double GAMMA = 1.0;
        private Double RATE = 1.0;
        private Long SEED = 1L;
        private Boolean LEAVES_ONLY = true;
        private Boolean INFER_LATENT = true; // set to false for inferring the plate nodes
        private GRASP.Inference MODE = null;
        private int[] ancestors = null; // if marginal inference, a list of ancestors of at least one

        private Map<String, Integer> latent2idx = new HashMap<>();
        private Map<String, Integer> feat2idx = new HashMap<>();

        public GRequest_InferModes(String command, String auth, JSONObject params) {
            super(command, auth);
            try {
                JSONObject tree = params.getJSONObject("Tree");
                idxTree = IdxTree.fromJSON(tree);
                JSONObject jdataset = params.getJSONObject("Dataset");
                dataset = JSONUtils.DataSet.fromJSON(jdataset);
                JSONObject jdistrib = params.getJSONObject("Distrib");
                JSONArray modetypes = jdistrib.getJSONArray("Modetypes");
                Enumerable[] modes = new Enumerable[modetypes.length()];
                for (int i = 0; i < modes.length; i ++)
                    modes[i] = Enumerable.fromJSON(modetypes.getJSONObject(i));
                PhyloPlate.Modes template = new PhyloPlate.Modes(modes);
                LEAVES_ONLY = params.optBoolean("Leaves-only", LEAVES_ONLY);
                INFER_LATENT = params.optBoolean("Latent", INFER_LATENT);
                SEED = params.optLong("Seed", SEED);
                RATE = params.optDouble("Rate", RATE);
                GAMMA = params.optDouble("Gamma", GAMMA);
                String infmode = params.optString("Inference", "Marginal");
                MODE = infmode.equals("Joint") ? GRASP.Inference.JOINT : (infmode.equals("Marginal") ? GRASP.Inference.MARGINAL : null);
                if (MODE == GRASP.Inference.MARGINAL) {
                    try {
                        Integer ancspec1 = params.optInt("Ancestor", -1);
                        if (ancspec1 == -1) {
                            JSONArray ancspec2 = params.getJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i++)
                                ancestors[i] = ancspec2.getInt(i);
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = ancspec1;
                        }
                    } catch (ClassCastException e) { // ancestor IDs specified with "N" prefix
                        String ancspec1 = params.optString("Ancestor");
                        if (ancspec1 == null) {
                            JSONArray ancspec2 = params.getJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i++)
                                ancestors[i] = Integer.parseInt(ancspec2.getString(i).substring(1));
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = Integer.parseInt(ancspec1.substring(1));
                            ;
                        }
                    }
                } else { // joint inference
                    ancestors = new int[idxTree.getSize()];
//                    for (int idx : idxTree)
                    //                        ancestors[idx] = idxTree.getLabel(idx);
                }
                pbn = new PhyloPlate(idxTree, template);
                PhyloPlate.Plate master = null;
                for (int idx : idxTree) {
                    if (LEAVES_ONLY && !idxTree.isLeaf(idx))
                        continue;
                    //System.out.println("Adding plate to branchpoint " + idx);
                    PhyloPlate.Plate plate = pbn.getPlate(idx);
                    JSONArray jnodes = jdistrib.getJSONArray("Nodes");
                    JSONArray jtargets = jdistrib.getJSONArray("Targets");
                    BNode[] nodes = new BNode[jnodes.length()];
                    for (int j = 0; j < nodes.length; j ++) {
                        JSONArray jtargnode = jtargets.getJSONArray(j);
                        int[] targets = new int[jtargnode.length()];
                        for (int k = 0; k < targets.length; k ++)
                            targets[k] = jtargnode.getInt(k);
                        plate.addNode(BNode.fromJSON(jnodes.getJSONObject(j), plate.getParents(targets)), pbn);
                    }
                    if (master == null) {
                        master = plate;
                        pbn.setMaster(master);
                    } else {
                        plate.setMaster(master);
                    }
                    //System.out.println("\tPlate: " + plate.toJSON());
                }
                dataset.curateFeatures(pbn.getBN());

                // create name-to-idx map
                for (int idx : idxTree) {
                    BNode[] nodes_on_plate = pbn.getBNodes(idx);
                    for (int j = 0; j < nodes_on_plate.length; j++)
                        latent2idx.put(nodes_on_plate[j].getVariable().toString(), idx);
                    BNode[] nodes_in_plate = pbn.getPlate(idx).getNodes();
                    for (int j = 0; j < nodes_in_plate.length; j++)
                        feat2idx.put(nodes_in_plate[j].getVariable().toString(), idx);
                }

            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command: " + command + "; " + e.getMessage());
            } catch (JSONUtils.JSONUtilsException e) {
                throw new GRequestRuntimeException("Invalid format: " + command + "; " + e.getMessage());
            } catch (RuntimeException e) {
                throw new GRequestRuntimeException("Technical problem: " + command + "; " + e.getMessage());
            }
            //
            // TODO: figure out compute resources, threads, memory, priority, queueing strategy etc.
            //

        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public boolean isQueued() {
            return true;
        }

        @Override
        public void run() {
            int[] predictidxs = new int[ancestors.length];
            String[] predictitems = new String[ancestors.length];
            for (int i = 0; i < ancestors.length; i ++) {
                if (MODE == GRASP.Inference.MARGINAL) {
                    predictidxs[i] = idxTree.getIndex(ancestors[i]); // retrieve the branchpoint index for (each of) the nominated ancestors
                    predictitems[i] = "N" + ancestors[i];
                } else {
                    predictidxs[i] = i; // retrieve the branchpoint index for (each of) the nominated ancestors
                    predictitems[i] = (!idxTree.isLeaf(i) ? "N" : "") + idxTree.getLabel(i).toString();
                }
            }
            String[] predictfeats = null;
            if (INFER_LATENT) {
                predictfeats = new String[pbn.getNModels()];
                for (int i = 0; i < predictfeats.length; i++)
                    predictfeats[i] = "State " + i;
            } else {
                Variable[] vars = pbn.getMaster().getVariables();
                predictfeats = new String[vars.length];
                for (int i = 0; i < vars.length; i++) {
                    predictfeats[i] = pbn.getMaster().expungePrefix(vars[i].getName());
                }
            }
            JSONObject jpred = new JSONObject();
            JSONArray jsamples = new JSONArray();
            JSONUtils.DataSet predict = new JSONUtils.DataSet(predictitems, predictfeats);
            for (int i = 0; i < dataset.getNSamples(); i ++) {
                JSONObject jmarg = new JSONObject();
                Object[][] output = new Object[predictitems.length][predictfeats.length];
                if (MODE == GRASP.Inference.MARGINAL) {
                    for (int k = 0; k < predictidxs.length; k ++) {
                        int bpidx = predictidxs[k];
                        if (bpidx >= 0) { // did find it...
                            JSONArray jarr = new JSONArray();
                            // inference below; first create the inference instance
                            // perform marginal inference
                            VarElim ve = new VarElim();
                            ve.instantiate(pbn.getBN());
                            BNode[] querynodes = INFER_LATENT ? pbn.getBNodes(bpidx) : pbn.getPlate(bpidx).getNodes();
                            for (int j = 0; j < querynodes.length; j++) {
                                pbn.instantiate(dataset, i);
                                Query q = ve.makeQuery(querynodes[j].getVariable());
                                CGTable r1 = (CGTable) ve.infer(q);
                                Object anydistrib = r1.query(querynodes[j].getVariable());
                                if (anydistrib instanceof EnumDistrib) {
                                    jarr.put(((EnumDistrib) anydistrib).toJSON());
                                    output[k][j] = ((EnumDistrib) anydistrib).toJSON();
                                } else {
                                    jarr.put(anydistrib.toString());
                                    output[k][j] = anydistrib;
                                }
                            }
                            jmarg.put(predictitems[k], jarr);
                        } else // did not find ancestor
                            ;
                    }
                    predict.addItemisedSample(output);
                    jsamples.put(jmarg);
                } else { // JOINT
                    VarElim ve = new VarElim();
                    ve.instantiate(pbn.getBN());
                    pbn.instantiate(dataset, i);
                    Query q = ve.makeMPE();
                    CGTable r1 = (CGTable) ve.infer(q);
                    Variable.Assignment[] assign = r1.getMPE();
                    for (Variable.Assignment assign1 : assign) {
                        String name = assign1.var.toString();
                        if (INFER_LATENT) {
                            Integer idx = latent2idx.get(name);
                            if (idx != null) { // inferring state and this variable is a match
                                BNode[] nodes_on_plate = pbn.getBNodes(idx);
                                for (int j = 0; j < nodes_on_plate.length; j++)
                                    if (nodes_on_plate[j].getVariable().equals(assign1.var)) {
                                        output[idx][j] = assign1.val;
                                        break;
                                    }
                            }
                        } else { // collect inferred features
                            Integer idx = feat2idx.get(name);
                            if (idx != null) { // inferring state and this variable is a match
                                BNode[] nodes_in_plate = pbn.getPlate(idx).getNodes();
                                for (int j = 0; j < nodes_in_plate.length; j++)
                                    if (nodes_in_plate[j].getVariable().equals(assign1.var)) {
                                        output[idx][j] = assign1.val;
                                        break;
                                    }
                            }
                        }
                    }
                    predict.addItemisedSample(output);
                }
                jpred.put("Predict", predict.toJSON());
                this.setResult(jpred);
            }
        }
    }

    public static class GRequest_InferBP extends GRequest {

        private final IdxTree idxTree;
        private final JSONUtils.DataSet dataset;
        private enum DATATYPE {CONTINUOUS, ENUMERABLE};
        private final DATATYPE TIP_TYPE;
        private PhyloBN pbn = null;
        private SubstModel MODEL = null;
        private int OBSERVED_NSTATES;
        private String[] OBSERVED_ALPHA = null;
        private Object[] LATENT_STATES = null;
        private Double GAMMA = 1.0;
        private Double RATE = 1.0;
        private Long SEED = 1L;
        private Boolean LEAVES_ONLY = true;
        private Integer MARG_LABEL = 0;
        private GRASP.Inference MODE = null;
        private int[] ancestors = null; // if marginal inference, a list of ancestors of at least one

        public GRequest_InferBP(String command, String auth, JSONObject params) {
            super(command, auth);
            try {
                JSONObject tree = params.getJSONObject("Tree");
                idxTree = IdxTree.fromJSON(tree);
                JSONObject jdataset = params.getJSONObject("Dataset");
                dataset = JSONUtils.DataSet.fromJSON(jdataset);
                JSONObject TIP_DISTRIB = params.getJSONObject("Distrib");
                JSONArray JSTATES = params.optJSONArray("States");
                if (JSTATES != null) {
                    LATENT_STATES = new Object[JSTATES.length()];
                    for (int i = 0; i < JSTATES.length(); i ++)
                        LATENT_STATES[i] = JSTATES.get(i);
                }
                LEAVES_ONLY = params.optBoolean("Leaves-only", LEAVES_ONLY);
                SEED = params.optLong("Seed", SEED);
                String infmode = params.optString("Inference", "Marginal");
                MODE = infmode.equals("Joint") ? GRASP.Inference.JOINT : (infmode.equals("Marginal") ? GRASP.Inference.MARGINAL : null);
                if (MODE == GRASP.Inference.MARGINAL) {
                    try {
                        Integer ancspec1 = params.optInt("Ancestor");
                        if (ancspec1 == null) {
                            JSONArray ancspec2 = params.optJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i++)
                                ancestors[i] = ancspec2.getInt(i);
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = ancspec1;
                        }
                    } catch (ClassCastException e) { // ancestor IDs specified with "N" prefix
                        String ancspec1 = params.optString("Ancestor");
                        if (ancspec1 == null) {
                            JSONArray ancspec2 = params.optJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i++)
                                ancestors[i] = Integer.parseInt(ancspec2.getString(i).substring(1));
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = Integer.parseInt(ancspec1.substring(1));
                            ;
                        }
                    }
                }
                RATE = params.optDouble("Rate", RATE);
                GAMMA = params.optDouble("Gamma", GAMMA);
                if (TSVFile.isDoubleOrInt(dataset.getNonitemisedData())) { // real values for tip nodes
                    TIP_TYPE = DATATYPE.CONTINUOUS;
                    if (LATENT_STATES != null) {
                        MODEL = new JC(GAMMA, LATENT_STATES);
                    } else
                        throw new GRequestRuntimeException("Latent states are invalid for " + command + "; States are " + (JSTATES != null ? JSTATES : "not given"));
                } else { // the column has discrete values (not real)
                    TIP_TYPE = DATATYPE.ENUMERABLE;
                    Set observed = new HashSet(new TSVFile(dataset.getNonitemisedData()).getValues());
                    OBSERVED_NSTATES = observed.size();
                    OBSERVED_ALPHA = new String[observed.size()];
                    observed.toArray(OBSERVED_ALPHA);
                    if (TIP_DISTRIB == null && LATENT_STATES == null) { // no extra nodes
                        MODEL = new JC(GAMMA, OBSERVED_ALPHA);
                    } else if (LATENT_STATES != null) {
                        MODEL = new JC(GAMMA, LATENT_STATES);
                    } else
                        throw new GRequestRuntimeException("Latent states are invalid for " + command + "; States are " + (JSTATES != null ? JSTATES : "not given"));
                }
                if (TIP_TYPE == DATATYPE.CONTINUOUS) {
                    pbn = PhyloBN.withGDTs(idxTree, MODEL, RATE, LEAVES_ONLY, SEED);
                    if (TIP_DISTRIB != null) // GDT given?
                        pbn.overrideMasterJSON(TIP_DISTRIB);
                } else { // discrete states only
                    if (TIP_DISTRIB == null && LATENT_STATES == null) { // no extra nodes
                        pbn = PhyloBN.create(idxTree, MODEL, RATE);
                    } else {
                        pbn = PhyloBN.withCPTs(idxTree, MODEL, OBSERVED_ALPHA, RATE, LEAVES_ONLY, SEED);
                        if (TIP_DISTRIB != null) // CPT given?
                            pbn.overrideMasterJSON(TIP_DISTRIB);
                    }
                }
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command : " + command + "; " + e.getMessage());
            }
            //
            // TODO: figure out compute resources, threads, memory, priority, queueing strategy etc.
            //

        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public boolean isQueued() {
            return true;
        }

        @Override
        public void run() {
            if (MODE == GRASP.Inference.MARGINAL) {
                JSONObject jpred = new JSONObject();
                for (int MARG_LABEL : ancestors) {
                    int bpidx = idxTree.getIndex(MARG_LABEL); // retrieve the branchpoint index for (each of) the nominated ancestors
                    if (bpidx < 0) // did not find it...
                        continue;
                    JSONArray jarr = new JSONArray();
                    // inference below; first create the inference instance
                    MaxLhoodMarginal<EnumDistrib> inf = new MaxLhoodMarginal(bpidx, pbn);
                    // perform marginal inference
                    String[] headers = dataset.getFeatures();
                    Object[][] rows = dataset.getNonitemisedData();
                    for (int i = 0; i < rows.length; i ++) {
                        TreeInstance ti = idxTree.getInstance(headers, rows[i]);
                        inf.decorate(ti);
                        // retrieve the distribution at the node previously nominated
                        Distrib anydistrib = inf.getDecoration(bpidx);
                        if (anydistrib instanceof EnumDistrib) {
                            jarr.put(((EnumDistrib) anydistrib).toJSON());
                        } else {
                            jarr.put(anydistrib.toString());
                        }
                    }
                    jpred.put("N" + MARG_LABEL, jarr);
                }
                this.setResult(jpred);
            } else { // joint:
                JSONObject jpred = new JSONObject();
                String[] headers = dataset.getFeatures();
                Object[][] rows = dataset.getNonitemisedData();
                Object[][] preds = new Object[idxTree.getSize()][rows.length];
                MaxLhoodJoint inf = new MaxLhoodJoint(pbn);
                Object[][] values = new Object[rows.length][idxTree.getSize()];
                String[] features = new String[idxTree.getSize()];
                for (int i = 0; i < rows.length; i ++) {
                    TreeInstance ti = idxTree.getInstance(headers, rows[i]);
                    inf.decorate(ti);
                    // inference will include all nodes (a subset instantiated before inference)
                    for (int j : idxTree) {
                        preds[j][i] = inf.getDecoration(j); //ti.getInstance(j);
                        values[i][j] = preds[j][i];
                    }
                }
                for (int NODE_LABEL : idxTree) {
                    JSONArray jarr = new JSONArray();
                    jarr.put(preds[NODE_LABEL]);
                    if (idxTree.isLeaf(NODE_LABEL)) {
                        features[NODE_LABEL] = idxTree.getLabel(NODE_LABEL).toString();
                    } else {
                        features[NODE_LABEL] = "N" + idxTree.getLabel(NODE_LABEL).toString();
                    }
                }
                JSONUtils.DataSet result = new JSONUtils.DataSet(features, values);
                jpred.put("Predict", JSONUtils.toJSON(result));
                this.setResult(jpred);
            }
        }
    }

    public static class GRequest_Fake extends GRequest {
        int sleepfor;

        public GRequest_Fake(String command, String auth, JSONObject params) {
            super(command, auth);
            try {
                sleepfor = params.getInt("Sleep");
            } catch (JSONException e) {
                throw new GRequestRuntimeException("Invalid JSON in command : " + command + "; " + e.getMessage());
            }
        }

        /**
         * Execute the job
         * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
         */
        @Override
        public boolean isQueued() {
            return true;
        }

        @Override
        public void run() {
            try {
                Thread.sleep(sleepfor);
            } catch (InterruptedException e) {
                System.err.println("Server was interrupted");
            }
            JSONObject myres = new JSONObject();
            myres.put("Fake", "Completed in " + sleepfor + " ms");
            this.setResult(myres);
        }
    }


    public static class GRequestRuntimeException extends RuntimeException {
        public GRequestRuntimeException(String msg) {
            super(msg);
        }
    }
}
