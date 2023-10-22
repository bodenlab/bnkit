package api;

import asr.Job;
import dat.EnumSeq;
import dat.phylo.IdxTree;
import json.JSONException;
import json.JSONObject;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Class to represent a request to GRASP server
 */
public class GRequest extends Thread implements Comparable<GRequest>, Job {

    private static int JOBCOUNTER = 1;
    private final int JOB;
    enum STATUS {NONE, WAITING, RUNNING, COMPLETED};
    STATUS status;
    private String ERRMSG = null;
    private final String command;
    private final String authtoken;

    private Map<String, Object> params;
    private int PRIORITY = 0; // default -10 to +10 with -10 being lowest priority and +10 being the highest
    private int NTHREADS = 1; // default number of threads that will be allowed to spawn from this
    private int MEMORY   = 1; // expected GB requirement of the request

    private JSONObject result = null; // container for result of job, when complete
    private String savedAsFile = null;  // filename, if result has been written

    public GRequest(String command, String authtoken) {
        this.command = command;
        this.authtoken = authtoken;
        this.JOB = JOBCOUNTER;
        this.status = STATUS.NONE;
        GRequest.JOBCOUNTER += 1;
    }

    /**
     * Retrieve a JSON object containing meta info about request
     * @return
     */
    public JSONObject job2JSON() {
        JSONObject json = new JSONObject();
        json.put("Command", command);
        json.put("Auth", authtoken);
        json.put("Job", getJob());
        json.put("Status", getStatus());
        json.put("Priority", PRIORITY);
        json.put("Threads", NTHREADS);
        json.put("Memory", MEMORY);
        return json;
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
        }
        return json;
    }

    public void setParams(Map<String, Object> params) {
        this.params = params;
    }

    @Override
    public String toString() {
        return "Command \"" + this.command + "\" [" + authtoken + "]";
    }

    protected void setError(String msg) {
        ERRMSG = msg;
    }

    public boolean isFailed() {
        return ERRMSG != null;
    }

    public String getError() {
        return ERRMSG == null ? "" : ERRMSG;
    }

    @Override
    public void run() {
    }

    public void runnow() {
        status = STATUS.RUNNING;
        run();
        status = STATUS.COMPLETED;
        System.out.println("Server completed job " + JOB);
    }

    public STATUS getStatus() {
        return status;
    }

    public boolean isWaiting() {
        return status == STATUS.WAITING;
    }

    public JSONObject getResult() {
        if (status == STATUS.COMPLETED) {
            if (this.result != null) // in memory still
                return result;
            else if (this.savedAsFile != null) { // not in memory, but saved
                try {
                    // read from a temporary file
                    BufferedReader br = new BufferedReader(new FileReader(savedAsFile));
                    StringBuilder sb = new StringBuilder();
                    String line = br.readLine();
                    while (line != null) {
                        sb.append(line);
                        line = br.readLine();
                    }
                    JSONObject json = new JSONObject(sb.toString());
                    return json;
                } catch (JSONException e) {
                    throw new RuntimeException(e);
                } catch (FileNotFoundException e) {
                    throw new RuntimeException(e);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            } else { // job has been removed entirely
                // TODO: respond?
            }
        } else { // job not complete
            // TODO: respond?
        }
        return null;
    }

    protected void setResult(JSONObject json) {
        result = json;
    }

    public int getJob() {
        return JOB;
    }

    public boolean isQueued() {
        if (command.startsWith("Recon")) {
            return true;
        } else if (command.startsWith("Fake")) {
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
    public static class JobQueue extends Thread {
        private final File directory;
        private final List<GRequest> currentjobs = new ArrayList<>();

        /**
         * Open a job queue, with storage (of completed jobs) located in specified directory
         * @param storage
         */
        public JobQueue(File storage) {
            this.directory = storage;
        }
        public synchronized void add(GRequest request) {
            request.status = STATUS.WAITING;
            currentjobs.add(request);
        }

        /**
         * Retrieve next job
         * TODO: possibly use priority instead of just polling the next in queue, then later when selected, ensure that CPU/memory requirements are managed
         *
         * @return next job to execute
         */
        public synchronized GRequest poll() {
            if (currentjobs.size() > 0) {
                for (GRequest req : currentjobs) {
                    if (req.status == STATUS.WAITING) {
                        return req;
                    }
                }
            }
            return null;
        }

        public List<GRequest> getJobs() {
            return currentjobs;
        }

        /**
         * Find job request entry in queue
         * @param job job number
         * @return request
         */
        public GRequest getRequest(int job) {
            for (GRequest req : currentjobs) {
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
            return getPlace(request.JOB);
        }
        /**
         * Retrieve place in queue of job.
         * @param job job number
         * @return place 1-size of queue; 0 if the job request does not exist
         */
        public int getPlace(int job) {
            int place = 1;
            for (GRequest req : currentjobs) {
                if (req.status == STATUS.WAITING) {
                    if (req.getJob() == job)
                        return place;
                    else
                        place += 1;
                }
            }
            return 0;
        }

        /**
         * Cancel job in queue (when WAITING).
         * @param request job request
         * @return true if cancelled, else false
         */
        public boolean cancel(GRequest request) {
            return cancel(request.JOB);
        }
        /**
         * Cancel job in queue (when WAITING).
         * @param job job number
         * @return true if cancelled, else false
         */
        public boolean cancel(int job) {
            for (GRequest req : currentjobs) {
                if (req.status == STATUS.WAITING) {
                    if (req.getJob() == job) {
                        return currentjobs.remove(req);
                    }
                }
            }
            return false;
        }

        @Override
        public void run() {
            // job queue thread
            // TODO: additional tasks not yet implemented include
            //  1a. cleaning-up the queue from client-retrieved jobs,
            //  1b. keep a pointer for queue to avoid re-scanning from top,
            //  2. cleaning-up filed jobs once lapsed
            try {
                while (true) {
                    GRequest req = this.poll();
                    if (req == null) {
                        // System.out.println("Waiting (sleeping) for jobs to be submitted");
                        Thread.sleep(500);
                    } else {
                        System.out.println("Found job; " + currentjobs.size() + " requests in queue");
                        // TODO: more advanced scheduling could be done; here, 1 job at a time
                        req.status = STATUS.RUNNING;
                        System.out.println("Server started job " + req.JOB);
                        req.start(); // could just be "run", so that current call stack is used
                        req.join();  // always wait for the request to finish
                        req.status = STATUS.COMPLETED;
                        if (req.isFailed())
                            System.out.println("Server failed to complete job " + req.JOB + " because: " + req.ERRMSG);
                        else
                            System.out.println("Server completed job " + req.JOB);
                        JSONObject result = req.getResult();
                        result.put("Error", req.ERRMSG);
                        if (result != null) {
                            // create a temporary file
                            try {
                                String filename = directory + File.separator + "GRequest_" + req.JOB + ".json";
                                BufferedWriter br = new BufferedWriter(new FileWriter(filename));
                                // Writes a string to the above temporary file
                                br.write(result.toString());
                                System.out.println("Wrote job " + req.JOB + " to file \"" + filename + "\"");
                                br.close();
                                req.savedAsFile = filename;
                                req.result = null;
                            } catch (IOException e2) {
                                System.err.println("Failed to create temp file for completed job " + req.JOB);
                            }
                        } else {
                            // TODO: how to act if the job completed without result?
                        }
                    }
                }
            } catch (InterruptedException e) {
                System.err.println("Job queue interrupted: " + e.getMessage());
                // TODO: save state of queue?

            }
        }
    }
}
