package dat.pog;

public class POG extends IdxGraph {

    /**
     * Create partial order graph, which is a directed and terminated IdxGraph with some specific properties.
     * @param nNodes     number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
     */
    public POG(int nNodes) {
        super(nNodes, false, true);
    }

    
}
