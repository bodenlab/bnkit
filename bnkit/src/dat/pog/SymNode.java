package dat.pog;

import dat.Enumerable;
import dat.colourschemes.Clustal;

public class SymNode<E extends Enumerable> extends Node {

    Enumerable domain;
    Object value = null;

    public SymNode(E domain, Object value) {
        this.domain = domain;
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
}
