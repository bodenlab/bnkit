package dat.pog;

import java.util.ArrayList;
import java.util.List;

/**
 * A "combined nodes" node.
 * Useful for graphs which are the result of merging one or more graphs.
 */
public class AggregateNode extends Node {

    private List<Node> unpacked;

    public List<Node> unpack() {
        return unpacked;
    }

    public AggregateNode(Node... nodes) {
        this.unpacked = new ArrayList<>(nodes.length);
        for (Node n : nodes) {
            if (n instanceof AggregateNode)
                unpacked.addAll(((AggregateNode)n).unpack());
            else
                unpacked.add(n);
        }
    }

}
