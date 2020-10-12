package dat.pog;

import bn.prob.EnumDistrib;
import dat.Enumerable;
import dat.colourschemes.Clustal;

public class EnumNode extends Node {

    Enumerable domain;
    EnumDistrib distrib;
    int[] counts;
    int tot = 0;
    boolean needs_update = true;

    public EnumNode(Enumerable domain) {
        this.domain = domain;
        this.distrib = new EnumDistrib(domain);
        counts = new int[domain.size()];
    }

    public EnumNode(EnumDistrib distrib) {
        this.distrib = distrib;
        this.domain = distrib.getDomain();
        this.counts = null;
        this.needs_update = false;
    }

    public void add(Object sym) {
        if (counts == null)
            return;
        this.needs_update = true;
        counts[domain.getIndex(sym)] += 1;
        tot += 1;
    }

    public EnumDistrib getDistrib() {
        if (needs_update)
            distrib = new EnumDistrib(domain, counts);
        needs_update = false;
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
        if (needs_update)
            distrib = new EnumDistrib(domain, counts);
        needs_update = false;
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

    public static Node[] toArray(EnumDistrib[] distribs) {
        EnumNode[] arr = new EnumNode[distribs.length];
        for (int i = 0; i < arr.length; i ++) {
            if (distribs[i] != null)
                arr[i] = new EnumNode(distribs[i]);
        }
        return arr;
    }

}
