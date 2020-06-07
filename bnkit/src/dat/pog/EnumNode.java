package dat.pog;

import bn.prob.EnumDistrib;
import dat.Enumerable;
import dat.colourschemes.Clustal;

public class EnumNode extends Node {

    Enumerable domain;
    EnumDistrib distrib;
    int[] counts;
    boolean needs_update = true;

    public EnumNode(Enumerable domain) {
        this.domain = domain;
        this.distrib = new EnumDistrib(domain);
        counts = new int[domain.size()];
    }

    public void add(Object sym) {
        this.needs_update = true;
        counts[domain.getIndex(sym)] += 1;
    }

    public EnumDistrib getDistrib() {
        if (needs_update)
            distrib = new EnumDistrib(domain, counts);
        needs_update = false;
        return distrib;
    }

    public int get(Object sym) {
        return counts[domain.getIndex(sym)];
    }

    public String getFillcolor() {
        Object consensus = domain.get(getDistrib().getMaxIndex());
        return Clustal.getColour((Character)consensus);
    }

    public String getLabel() {
        Object consensus = distrib.getMax();
        return consensus.toString();
    }

    public String getStyle() {
        if (getDistrib().get(getDistrib().getMaxIndex()) > 0.5)
            return "bold,rounded,filled";
        return "rounded,filled";
    }
}
