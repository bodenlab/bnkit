package dat.pog;

public class POGraph<E extends dat.Enumerable> extends IdxEdgeGraph {

    /**
     * Create partial order graph, which is a directed and terminated IdxGraph with some specific properties.
     * @param nNodes     number of nodes excluding (optional) virtual nodes to have multiple entry and exit points in the graph
     */
    public POGraph(int nNodes) {
        super(nNodes, false, true);
    }

    public String toString() {
        return super.toString();
    }
    
}
