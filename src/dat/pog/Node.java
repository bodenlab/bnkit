package dat.pog;

import json.JSONObject;

import java.io.Serializable;

/**
 * Default, place-holder implementation of a node
 */
public class Node implements Serializable {

    /**
     *      *   node  [style="rounded,filled,bold", shape=box, fixedsize=true, width=1.3, fontname="Arial"];
     *      *   Created   [fillcolor=black, shape=circle, label="", width=0.25];
     *      *   Destroyed [fillcolor=black, shape=doublecircle, label="", width=0.3];
     *      *   Empty     [fillcolor="#a0ffa0"];
     * @return
     */

    protected String label = null;

    /** Many features of Node class are specific to DOT */
    private static String shape = null; // box, circle, doublecircle
    protected Object xlabel = null;
    private Boolean fixedsize = null;
    private Double width = null;
    private static String fontname = null; // "Arial";

    public Node() {
    }

    public String getLabel() {
        return label;
    };
    public void setLabel(String label) { this.label = label; }

    public Object getXlabel() {
        return xlabel;
    }

    public void setXlabel(Object xlabel) {
        this.xlabel = xlabel;
    }

    public Boolean getFixedsize() {
        return fixedsize;
    }

    public void setFixedsize(Boolean fixedsize) {
        this.fixedsize = fixedsize;
    }

    public Double getWidth() {
        return width;
    }

    public void setWidth(Double width) {
        this.width = width;
    }

    public static String getFontname() {
        return fontname;
    }

    public static void setFontname(String fontname) {
        Node.fontname = fontname;
    }

    public static String getShape() {
        return shape;
    }

    public static void setShape(String shape) {
        Node.shape = shape;
    }

    @Override
    public boolean equals(Object other) {
        if (! (other instanceof Node))
            return false;
        String otherlabel = ((Node) other).getLabel();
        if (label == null)
            return otherlabel == null;
        else if (otherlabel != null)
            return (label.equalsIgnoreCase(((Node) other).getLabel()));
        return false;
    }

    public String toDOT() {
        return  ((getFillcolor() == null) ? "" : ("fillcolor=\"" + getFillcolor() + "\",")) +
                ((shape == null) ? "" : ("shape=" + shape + ",")) +
                ((getLabel() == null) ? "" : ("label=\"" + getLabel() + "\",")) +
                ((xlabel == null) ? "" : ("xlabel=\"" + xlabel + "\",")) +
                ((fixedsize == null) ? "" : ("fixedsize=" + fixedsize + ",")) +
                ((width == null) ? "" : ("width=" + width + ",")) +
                ((fontname == null) ? "" : ("fontname=\"" + fontname + "\",")) +
                ((getStyle() == null) ? "" : ("style=\"" + getStyle() + "\""));
    }

    public Object[] getObjectArray() {
        if (label != null)
            return new Object[] {label};
        else
            return new Object[0];
    }

    public static Node fromObjectArray(Object[] input) {
        Node n = new Node();
        if (input.length >= 1)
            n.setLabel(input[0].toString());
        return n;
    }

    /**
     * Save the instance as a JSON object
     * @return
     */
    public JSONObject toJSON() {
        JSONObject json = new JSONObject();
        if (label != null)
            json.put("Label", label);
        return json;
    }

    /**
     * Convert a JSONObject into an instance of this class
     * @param json
     * @return
     */
    public static Node fromJSON(JSONObject json) {
        String label = json.optString("Label", null);
        Node n = new Node();
        if (label != null)
            n.setLabel(label);
        return n;
    }

    public void setXLabel(Object xlabel) {
        this.xlabel = xlabel;
    }

    public String getFillcolor() {
        return null;
    };
    public String getStyle() {
        return null; // "rounded,filled,bold";
    }

}
