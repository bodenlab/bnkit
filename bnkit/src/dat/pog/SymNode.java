package dat.pog;

import dat.colourschemes.Clustal;

public class SymNode extends Node {
    Object value = null;

    public SymNode(Object value) {
        this.value = value;
    }

    public void set(Object sym) {
        this.value = sym;
    }

    public Object get() {
        return this.value;
    }

    public String getFillcolor() {
        return Clustal.getColour((char)value);
    }

    public String getLabel() {
        return value.toString();
    }

    public String getStyle() {
        return "rounded,filled";
    }

    public static Node[] toArray(Object[] values) {
        SymNode[] arr = new SymNode[values.length];
        for (int i = 0; i < arr.length; i ++) {
            arr[i] = new SymNode(values[i]);
        }
        return arr;
    }
}
