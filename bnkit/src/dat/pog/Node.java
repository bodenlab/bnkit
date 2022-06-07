package dat.pog;

/**
 * Default, place-holder implementation of a node
 */
public class Node {
    /**
     *      *   node  [style="rounded,filled,bold", shape=box, fixedsize=true, width=1.3, fontname="Arial"];
     *      *   Created   [fillcolor=black, shape=circle, label="", width=0.25];
     *      *   Destroyed [fillcolor=black, shape=doublecircle, label="", width=0.3];
     *      *   Empty     [fillcolor="#a0ffa0"];
     * @return
     */

    public Node() {

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

    public void setXLabel(Object xlabel) {
        this.xlabel = xlabel;
    }

    public String getFillcolor() {
        return null;
    };
    public static String shape = null; // box, circle, doublecircle
    public String label = null;
    public Object xlabel = null;
    public String getLabel() {
        return label;
    };
    public void setLabel(String label) { this.label = label; }
    public Boolean fixedsize = null;
    public Double width = null;
    public static String fontname = null; // "Arial";
    public String getStyle() {
        return null; // "rounded,filled,bold";
    }

}
