package dat.pog;

import dat.EnumSeq;
import dat.Enumerable;

import java.util.HashSet;
import java.util.Set;

/**
 * Edge that forms part of a given sequence
 */
public class SeqEdge extends Edge {

    private Set<EnumSeq<Enumerable>> seqs;
    private int[] counts = null;
    private int count_total = 0;
    private Enumerable domain;
    private int normalise_by = 1;

    public SeqEdge(Enumerable domain) {
        super();
        this.domain = domain;
        seqs = new HashSet<>();
    }

    /**
     * Set the normalisation constant with which sequence counts are normalised
     * @param normalise_by
     */
    public void setTotal(int normalise_by) {
        this.normalise_by = normalise_by;
    }

    public void add(EnumSeq<Enumerable> seq) {
        seqs.add(seq);
    }
    public void add(EnumSeq<Enumerable> seq, Object from_sym, Object to_sym) {
        add(seq);
        if (counts == null) {
            counts = new int[domain.size() * domain.size()];
            counts[domain.getIndex(from_sym) * domain.size() + domain.getIndex(to_sym)] += 1;
            count_total += 1;
        }
    }

    public double getProb(Enumerable from_sym, Enumerable to_sym) {
        return count_total > 0 ? ((double)getCount(from_sym, to_sym)) / ((double) count_total) : 0;
    }

    public int getCount(Enumerable from_sym, Enumerable to_sym) {
        return counts[domain.getIndex(from_sym) * domain.size() + domain.getIndex(to_sym)];
    }

    public Double getWeight() {
        if (seqs.size() > 0)
            return (2 * Double.valueOf(seqs.size() / (double)normalise_by));
        else return null;
    }

}