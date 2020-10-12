package dat.pog;

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

}
