package api;

import asr.ASRRuntimeException;
import asr.GRASP;
import asr.GServer;
import asr.Prediction;
import bn.ctmc.SubstModel;
import dat.EnumSeq;
import dat.Enumerable;
import dat.phylo.IdxTree;
import dat.pog.POGTree;
import dat.pog.POGraph;
import json.JSONArray;
import json.JSONException;
import json.JSONObject;

import java.util.HashMap;
import java.util.Map;

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

            } else if (command.equals("Label-tree")) {
                request = new GRequest_LabelTree(command, authtoken, jparams);

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
        private Prediction indelpred = null;
        private double[] RATES = null;
        private long START_TIME;

        public GRequest_Recon(String command, String auth, JSONObject params) {
            super(command, auth);
            // extract necessary detail from command params etc
            try {
                JSONObject jinput = params.optJSONObject("POGTree");
                if (jinput != null) {
                    POGTree pogtree = POGTree.fromJSON(jinput);
                    JSONArray jancs = params.optJSONArray("Ancestors");
                    if (jancs != null) {
                        if (jancs.length() == pogtree.getTree().getNParents()) {
                            Map<Object, POGraph> ancestors = new HashMap<>();
                            for (int i = 0; i < jancs.length(); i++) {
                                JSONObject obj = jancs.getJSONObject(i);
                                POGraph pog = POGraph.fromJSON(obj);
                                ancestors.put(pog.getName(), pog);
                            }
                            indelpred = new Prediction(pogtree, ancestors);
                        }
                    }
                    idxTree = pogtree.getTree();
                } else {
                    JSONObject tree = params.getJSONObject("Tree");
                    idxTree = IdxTree.fromJSON(tree);
                    JSONObject jaln = params.getJSONObject("Alignment");
                    aln = EnumSeq.Alignment.fromJSON(jaln);
                }
                String infmode = params.optString("Inference", "Joint");
                MODE = infmode.equals("Joint") ? GRASP.Inference.JOINT : (infmode.equals("Marginal") ? GRASP.Inference.MARGINAL : null);
                if (MODE == GRASP.Inference.MARGINAL) {
                    try {
                        Integer ancspec1 = params.optInt("Ancestor");
                        if (ancspec1 == null) {
                            JSONArray ancspec2 = params.optJSONArray("Ancestors");
                            ancestors = new int[ancspec2.length()];
                            for (int i = 0; i < ancestors.length; i ++)
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
                            for (int i = 0; i < ancestors.length; i ++)
                                ancestors[i] = Integer.parseInt(ancspec2.getString(i).substring(1));
                        } else {
                            ancestors = new int[1];
                            ancestors[0] = Integer.parseInt(ancspec1.substring(1));;
                        }
                    }
                }
                if (indelpred == null) { // no POGTree AND no ancestor POGs provided
                    String indels = params.optString("Indels", "BEP");
                    for (int i = 0; i < INDELS.length; i ++) {
                        if (INDELS[i].equalsIgnoreCase(indels)) {
                            INDEL_IDX = i;
                            break;
                        }
                    }
                }
                // other params:
                // model (default is JTT)
                String modelname = params.optString("Model", "JTT");
                MODEL = SubstModel.createModel(modelname);
                if (MODEL == null)
                    throw new ASRRuntimeException("Invalid model");
                // rates (optional)
                JSONArray jrates = params.optJSONArray("Rates");
                if (jrates != null) {
                    RATES = new double[jrates.length()];
                    for (int i = 0; i < RATES.length; i ++)
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
            if (indelpred == null) {
                POGTree pogtree = new POGTree(aln, idxTree);
                switch (INDEL_IDX) {
                    case 0:
                        indelpred = Prediction.PredictByBidirEdgeParsimony(pogtree);
                        break;
                    case 1:
                        indelpred = Prediction.PredictByBidirEdgeMaxLhood(pogtree);
                        break;
                    case 2:
                        indelpred = Prediction.PredictBySICP(pogtree);
                        break;
                    case 3:
                        indelpred = Prediction.PredictBySICML(pogtree);
                        break;
                    case 4:
                        indelpred = Prediction.PredictByParsimony(pogtree);
                        break;
                    case 5:
                        indelpred = Prediction.PredictByMaxLhood(pogtree);
                        break;
                    default:
                        break;
                }
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
            // TODO: collate results
            JSONObject myres = new JSONObject();
            myres.put("Prediction", indelpred.toJSON());
            this.setResult(myres);
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
                throw new GRequestRuntimeException("Invalid JSON in command : " + command + "; " + e.getMessage());
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
