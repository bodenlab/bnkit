package dat.pog;

import java.util.Objects;

public class EnumEdge {
    final EnumNode node1, node2; // direction node1 --> node2

    /**
     * Create an edge from node1 to node2
     * @param node2 first node
     * @param node2 second node
     */
    public EnumEdge(EnumNode node1, EnumNode node2) {
        this.node1 = node1;
        this.node2 = node2;
    }
    @Override
    public String toString() {
        return node1 + "-->" + node2;
    }
    @Override
    public boolean equals(Object o) {
        try {
            EnumEdge e = (EnumEdge)o;
            return (e.node1 == this.node1 && e.node2 == this.node2);
        } catch (ClassCastException e) {
            e.printStackTrace();
            return false;
        }
    }
    @Override
    public int hashCode() {
        return Objects.hash(node1, node2);
    }
    public EnumNode[] getPair() {
        return new EnumNode[] {node1, node2};
    }

}
