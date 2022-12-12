package dat.pog;

import bn.prob.EnumDistrib;
import dat.Enumerable;
import dat.colourschemes.Clustal;
import json.JSONObject;

import java.io.Serializable;

public class EnumNode extends Node implements Serializable {

    private Enumerable domain;
    private EnumDistrib distrib;
    private int[] counts;
    private int tot = 0;


    private boolean notUpdated = true;

    public EnumNode() {

    }

    public EnumNode(Enumerable domain) {
        this.domain = domain;
        this.distrib = new EnumDistrib(domain);
        counts = new int[domain.size()];
    }

    public EnumNode(EnumDistrib distrib) {
        this.distrib = distrib;
        this.domain = distrib.getDomain();
        this.counts = null;
        this.notUpdated = false;
    }

    public void add(Object sym) {
        if (counts == null)
            return;
        this.notUpdated = true;
        counts[domain.getIndex(sym)] += 1;
        tot += 1;
    }

    public EnumDistrib getDistrib() {
        if (notUpdated)
            distrib = new EnumDistrib(domain, counts);
        notUpdated = false;
        return distrib;
    }

    public int get(Object sym) {
        if (counts == null)
            return 0;
        return counts[domain.getIndex(sym)];
    }

    public String getFillcolor() {
        Object consensus = domain.get(getDistrib().getMaxIndex());
        return Clustal.getColour((Character)consensus);
    }

    public String getLabel() {
        if (notUpdated)
            distrib = new EnumDistrib(domain, counts);
        notUpdated = false;
        if (tot == 0 && counts != null)
            return null;
        Object consensus = distrib.getMax();
        return consensus.toString();
    }

    public String getStyle() {
        if (getDistrib().get(getDistrib().getMaxIndex()) > 0.5)
            return "bold,rounded,filled";
        return "rounded,filled";
    }
    public Enumerable getDomain() {
        return domain;
    }

    public void setDomain(Enumerable domain) {
        this.domain = domain;
    }

    public void setDistrib(EnumDistrib distrib) {
        this.distrib = distrib;
    }

    public int[] getCounts() {
        return counts;
    }

    public void setCounts(int[] counts) {
        this.counts = counts;
    }

    public int getTot() {
        return tot;
    }

    public void setTot(int tot) {
        this.tot = tot;
    }

    public boolean isNotUpdated() {
        return notUpdated;
    }

    public void setNotUpdated(boolean notUpdated) {
        this.notUpdated = notUpdated;
    }

    public static Node[] toArray(EnumDistrib[] distribs) {
        EnumNode[] arr = new EnumNode[distribs.length];
        for (int i = 0; i < arr.length; i ++) {
            if (distribs[i] != null) {
                arr[i] = new EnumNode(distribs[i]);
                arr[i].xlabel = i + 1;
            }
        }
        return arr;
    }
    /**
     * Save the instance as a JSON object
     * @return
     */
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        if (label != null)
            json.put("Label", label);
        EnumDistrib d = getDistrib();
        if (d != null) {
            for (String key : d.toJSON().keySet())
                json.put(key, d.toJSON().getJSONObject(key));
        }
        return json;
    }

    /**
     * Convert a JSONObject into an instance of this class
     * @param json
     * @return
     */
    public static EnumNode fromJSON(JSONObject json) {
        String label = json.optString("Label", null);
        EnumDistrib d = EnumDistrib.fromJSON(json);
        EnumNode n = new EnumNode(d);
        if (label != null)
            n.setLabel(label);
        return n;
    }

}
