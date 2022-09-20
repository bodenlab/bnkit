package dat.pog;

import dat.colourschemes.Clustal;
import json.JSONObject;

import java.io.Serializable;

public class SymNode extends Node implements Serializable {

    private Object value = null;

    public SymNode() {
        this.value = null;
    }

    public SymNode(Object value) {
        this.value = value;
    }

    public Object getValue() {
        return value;
    }
    public void setValue(Object sym) {
        this.value = sym;
    }

    public Object get() {
        return this.value;
    }
    @Override
    public String getFillcolor() {
        return Clustal.getColour((char)value);
    }
    @Override
    public String getStyle() {
        return "rounded,filled";
    }

    /**
     * Save the instance as a JSON object
     * @return
     */
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        if (label != null)
            json.put("Label", label);
        if (value != null)
            json.put("Value", value);
        return json;
    }

    /**
     * Convert a JSONObject into an instance of this class
     * @param json
     * @return
     */
    public static SymNode fromJSON(JSONObject json) {
        String label = json.optString("Label", null);
        Object value = json.optString("Value", null);
        SymNode n = new SymNode(value);
        if (label != null)
            n.setLabel(label);
        return n;
    }

    public static Node[] toArray(Object[] values) {
        SymNode[] arr = new SymNode[values.length];
        for (int i = 0; i < arr.length; i ++) {
            if (values[i] != null) {
                arr[i] = new SymNode(values[i]);
                arr[i].xlabel = i + 1;
            }
        }
        return arr;
    }

    @Override
    public Object[] getObjectArray() {
        if (label != null)
            return new Object[] {get(), label};
        else
            return new Object[] {get()};
    }

    public static SymNode fromObjectArray(Object[] input) {
        SymNode n = new SymNode(input[0]);
        if (input.length >= 2)
            n.setLabel(input[1].toString());
        return n;
    }

}
