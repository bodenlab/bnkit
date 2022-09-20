package dat.pog;

import json.JSONObject;

/**
 * Default, place-holder implementation of an edge
 */
public class Edge {

    public Edge() {

    }

    /**
     *      *   edge  [style=bold, fontname="Arial", weight=100]
     *      *   Empty     -> Announced [ label="announce"   ];
     * @return
     */
    public String toDOT() {
        return  ((label == null) ? "" : ("label=\"" + label + "\",")) +
//                ((getWeight() == null) ? "" : ("penwidth=" + getWeight() + ",")) +
                ((fontname == null) ? "" : ("fontname=\"" + fontname + "\",")) +
                ((style == null) ? "" : ("style=\"" + style + "\""));
    }

    public String getLabel() {
        return label;
    }

    public void setLabel(String label) {
        this.label = label;
    }
    public String label = null;
    public String fontname = null; // "Arial";
    public String style = null; // "bold", "dotted";

    public JSONObject toJSON() {
        JSONObject edge = new JSONObject();
        if (getLabel() != null)
            edge.put("Label", this.getLabel());
        return edge;
    }

    public static Edge fromJSON(JSONObject json) {
        Edge e = new Edge();
        String label = json.optString("Label", null);
        if (label != null)
            e.setLabel(label);
        return e;
    }

    public Object[] toObjectArray() {
        if (label != null)
            return new Object[] {label};
        else
            return new Object[] {};
    }
    public static Edge fromObjectArray(Object[] input) {
        Edge e = new Edge();
        if (input != null)
            if (input.length > 0)
                if (input[0] != null)
                    e.setLabel(input[0].toString());
        return e;
    }
}
