package dat.pog;

import java.util.Objects;

public class POGEdge {
    final int idx1, idx2; // direction idx1 --> idx2

    /**
     * Create an edge from idx1 to idx2
     * @param idx1 first index
     * @param idx2 second index
     */
    public POGEdge(int idx1, int idx2) {
        this.idx1 = idx1;
        this.idx2 = idx2;
    }
    @Override
    public String toString() {
        return idx1 + "-->" + idx2;
    }
    @Override
    public boolean equals(Object o) {
        try {
            POGEdge e = (POGEdge)o;
            return (e.idx1 == this.idx1 && e.idx2 == this.idx2);
        } catch (ClassCastException e) {
            e.printStackTrace();
            return false;
        }
    }
    @Override
    public int hashCode() {
        return Objects.hash(idx1, idx2);
    }
    public int[] getIndices() {
        return new int[] {idx1, idx2};
    }
}
