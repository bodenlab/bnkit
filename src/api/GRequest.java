package api;

import dat.file.Newick;
import dat.phylo.IdxTree;
import json.JSONException;
import json.JSONObject;

import java.util.*;
import java.util.concurrent.PriorityBlockingQueue;

/**
 * Class to represent a request to GRASP server
 */
public class GRequest extends Thread implements Comparable<GRequest> {

    private static int JOBCOUNTER = 0;
    private final int JOB;
    private final String command;
    private final String authtoken;

    private Map<String, Object> params;
    private Object param;
    private JSONObject result; // container for result of job, when complete

    private GRequest(String command, String authtoken) {
        this.command = command;
        this.authtoken = authtoken;
        this.JOB = JOBCOUNTER;
        GRequest.JOBCOUNTER += 1;
    }

    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        json.put("Command", command);
        json.put("Auth", authtoken);
        if (params != null) {
            JSONObject jparams = new JSONObject();
            for (Map.Entry<String, Object> entry : params.entrySet())
                jparams.put(entry.getKey(), entry.getValue());
            json.put("Params", jparams);
        } else if (param != null)
            json.put("Param", param);
        return json;
    }

    public void setParams(Map<String, Object> params) {
        this.param = null;
        this.params = params;
    }

    public void setParam(Object param) {
        this.params = null;
        this.param = param;
    }

    @Override
    public String toString() {
        return "Command \"" + this.command + "\" [" + authtoken + "]";
    }

    public static GRequest fromJSON(JSONObject json) {
        try {
            String command = json.getString("Command");
            String authtoken = json.getString("Auth");
            GRequest greq = new GRequest(command, authtoken);
            JSONObject jparams = json.optJSONObject("Params");
            if (jparams != null) {
                Map<String, Object> params = new HashMap<>();
                for (String key : jparams.keySet())
                    params.put(key, jparams.opt(key));
                greq.setParams(params);
            } else {
                Object jparam = json.opt("Param");
                if (jparam != null)
                    greq.setParam(jparam);
            }
            return greq;
        } catch (JSONException e) {
            e.printStackTrace();
            throw new RuntimeException("Invalid request: missing Command or Auth");
        }
    }

    /**
     * Execute the job
     * and place results in the instance of this class for later retrieval with {@see GRequest#getResult()}
     */
    @Override
    public void run() {
        if (command.equals("Label-tree")) {
            JSONObject tree = (JSONObject) params.get("Tree");
            IdxTree idxTree = IdxTree.fromJSON(tree);
            System.out.println("Labelling: " + idxTree);
        } else if (command.equals("Retrieve")) {
            int job = (int) param;


        } else {
            System.out.println("Server started job " + JOB + " using thread: " + toString());
            try {
                Thread.sleep(10000);
            } catch (InterruptedException e) {
                System.err.println("Server was interrupted");
            }
        }
        System.out.println("Job done.");
    }

    public JSONObject getResult() {
        return result;
    }

    public int getJob() {
        return JOB;
    }

    public boolean isQueued() {
        if (command.startsWith("Recon-")) {
            return true;
        } else if (command.startsWith("Label-")) {
            return false;
        } else
            return false;
    }

    @Override
    public int compareTo(GRequest o) {
        return this.JOB - o.JOB;
    }

    /**
     * Basic job queue for requests
     */
    public static class JobQueue {
        List<GRequest> queue = new ArrayList<>();
        public JobQueue() {
        }
        public void add(GRequest request) {
            queue.add(request);
        }
        public synchronized GRequest poll() {
            if (queue.size() > 0) {
                GRequest req = queue.get(0);
                queue.remove(0);
                return req;
            }
            return null;
        }

        /**
         * Find job request entry
         * @param job job number
         * @return request
         */
        public GRequest getRequest(int job) {
            for (GRequest req : queue) {
                if (req.getJob() == job)
                    return req;
            }
            return null;
        }

        /**
         * Retrieve place in queue of job.
         * @param request job request
         * @return place 1-size of queue; 0 if the job request does not exist
         */
        public int getPlace(GRequest request) {
            return queue.indexOf(request) + 1;
        }
        /**
         * Retrieve place in queue of job.
         * @param job job number
         * @return place 1-size of queue; 0 if the job request does not exist
         */
        public int getPlace(int job) {
            int place = 1;
            for (GRequest req : queue) {
                if (req.getJob() == job)
                    return place;
                place += 1;
            }
            return 0;
        }
    }
}
