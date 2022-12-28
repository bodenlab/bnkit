package api;

import json.JSONObject;

import java.util.HashMap;
import java.util.Map;

public class GMessage {

    static final Map<Integer, String> errmsg = new HashMap<>() {
        {
            put(1, "Request is not in JSON");
            put(2, "Job is not available");
            put(3, "Unknown command");
            put(4, "Request is invalid");
            put(5, "Connection denied");
        }
    };

    public static JSONObject errorToJSON(int error) {
        return errorToJSON(error, null);
    }
    public static JSONObject errorToJSON(int error, String comment) {
        JSONObject json = new JSONObject();
        json.put("Error",   error);
        String msg = errmsg.get(error);
        if (msg == null)
            msg = "Unknown error";
        json.put("Message", msg);
        if (comment != null)
            json.put("Comment", comment);
        return json;
    }

    public static JSONObject server2clientReJob(int job, String msg) {
        JSONObject json = new JSONObject();
        json.put("Job",   job);
        if (msg == null)
            msg = "Unknown";
        json.put("Message", msg);
        return json;
    }

    public static JSONObject server2clientReJob(int job, String key, Object value) {
        JSONObject json = new JSONObject();
        json.put("Job",   job);
        json.put(key, value);
        return json;
    }
    public static JSONObject server2clientReJob(int job, String attachment, JSONObject jsonObject) {
        JSONObject json = new JSONObject();
        json.put("Job",   job);
        json.put(attachment, jsonObject);
        return json;
    }

    /**
     * Extract job number from JSON message
     * @param json
     * @return
     */
    public static int fromJSON2Job(JSONObject json) {
        Integer job = json.optInt("Job");
        return job;
    }

    /**
     * Extract result from JSON message
     * @param json
     * @return
     */
    public static JSONObject fromJSON2Result(JSONObject json) {
        return json.optJSONObject("Result");
    }

    public static JSONObject client2serverReJob(int job, String command) {
        JSONObject json = new JSONObject();
        json.put("Job",   job);
        json.put("Command", command);
        return json;
    }

}
