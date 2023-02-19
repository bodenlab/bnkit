package api;

import asr.*;
import bn.Distrib;
import bn.ctmc.SubstModel;
import bn.ctmc.matrix.JC;
import bn.prob.EnumDistrib;
import dat.EnumSeq;
import dat.Enumerable;
import dat.file.TSVFile;
import dat.phylo.IdxTree;
import dat.phylo.PhyloBN;
import dat.phylo.TreeInstance;
import dat.pog.POGTree;
import dat.pog.POGraph;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.io.IOException;
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
                if (TSVFile.isDoubleOrInt(dataset.values)) { // real values for tip nodes
                    TIP_TYPE = DATATYPE.CONTINUOUS;
                    if (LATENT_STATES != null) {
                        MODEL = new JC(GAMMA, LATENT_STATES);
                    } else
                        throw new GRequestRuntimeException("Latent states are invalid for " + command + "; States are " + (JSTATES != null ? JSTATES : "not given"));
                } else { // the column has discrete values (not real)
                    TIP_TYPE = DATATYPE.ENUMERABLE;
                    Set observed = new HashSet(new TSVFile(dataset.values).getValues());
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
            pbn.trainEM(dataset.headers, dataset.values, SEED);
            JSONObject myres = new JSONObject();
            myres.put("Distrib", pbn.getMasterJSON());
            this.setResult(myres);
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
                if (TSVFile.isDoubleOrInt(dataset.values)) { // real values for tip nodes
                    TIP_TYPE = DATATYPE.CONTINUOUS;
                    if (LATENT_STATES != null) {
                        MODEL = new JC(GAMMA, LATENT_STATES);
                    } else
                        throw new GRequestRuntimeException("Latent states are invalid for " + command + "; States are " + (JSTATES != null ? JSTATES : "not given"));
                } else { // the column has discrete values (not real)
                    TIP_TYPE = DATATYPE.ENUMERABLE;
                    Set observed = new HashSet(new TSVFile(dataset.values).getValues());
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
                    String[] headers = dataset.headers;
                    Object[][] rows = dataset.values;
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
                String[] headers = dataset.headers;
                Object[][] rows = dataset.values;
                Object[][] preds = new Object[idxTree.getSize()][rows.length];
                MaxLhoodJoint inf = new MaxLhoodJoint(pbn);
                JSONUtils.DataSet result = new JSONUtils.DataSet();
                result.values = new Object[rows.length][idxTree.getSize()];
                result.headers = new String[idxTree.getSize()];
                for (int i = 0; i < rows.length; i ++) {
                    TreeInstance ti = idxTree.getInstance(headers, rows[i]);
                    inf.decorate(ti);
                    // inference will include all nodes (a subset instantiated before inference)
                    for (int j : idxTree) {
                        preds[j][i] = inf.getDecoration(j); //ti.getInstance(j);
                        result.values[i][j] = preds[j][i];
                    }
                }
                for (int NODE_LABEL : idxTree) {
                    JSONArray jarr = new JSONArray();
                    jarr.put(preds[NODE_LABEL]);
                    if (idxTree.isLeaf(NODE_LABEL)) {
                        result.headers[NODE_LABEL] = idxTree.getLabel(NODE_LABEL).toString();
                    } else {
                        result.headers[NODE_LABEL] = "N" + idxTree.getLabel(NODE_LABEL).toString();
                    }
                }
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
