/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNet;
import bn.BNode;
import bn.EnumVariable;
import bn.Predef;
import bn.Variable;
import bn.alg.CGTable;
import bn.alg.CGVarElim;
import bn.alg.Query;
import com.mxgraph.layout.hierarchical.mxHierarchicalLayout;
import com.mxgraph.layout.mxGraphLayout;
import com.mxgraph.layout.mxParallelEdgeLayout;
import com.mxgraph.layout.mxPartitionLayout;
import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.util.mxConstants;
import com.mxgraph.util.mxRectangle;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxStylesheet;
import java.awt.BorderLayout;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Random;
import javax.swing.JFileChooser;
import javax.swing.JPanel;

/**
 *
 * @author jun View class GraphPanel visually represents the graph structure.
 *
 * Mouse and keyboard input are captured in this class to allow node selection
 * and movement, edge insertion, key events eg. zoom, save/load.
 */
public final class GraphPanel extends JPanel implements Serializable, Observer {

    private mxGraph graph;
    final mxGraphComponent graphComponent;
    private final Map<String, Object> allVertices = new HashMap<String, Object>();
    public final Object lastPressedVertex = null;
    private List<Object> selectedCells = new ArrayList<>();
    private final List<NodeModel> nodeModels = new ArrayList<>();
    private BNModel model;
    private final Map<String, Integer> nodeCounts = new HashMap<>(); // Track number of each type of node to ensure
                                                                     // new nodes are assigned unique names.
    private NodeModel queryNode; // For now, there is only one query node

    public mxGraphComponent getGraphComponent() {
        return graphComponent;
    }

    public mxGraph getGraph() {
        return graph;
    }

    public void setGraph(mxGraph graph) {
        this.graph = graph;
    }

    public BNContainer getBNContainer() {
        return model.getBNC();
    }

    public NodeModel getQueryNode() {
        return queryNode;
    }

    public void setQueryNode(NodeModel node) {
        queryNode = node;
    }

    public Map getNodeCounts() {
        return nodeCounts;
    }

    @Deprecated
    /**
     * Not presently used. Pre-defines visual
     * styles for nodes. 
     */
    public void defStyleSheets(mxGraph graph) {
        String STRING_STYLE = "STRING_STYLE";
        String BOOL_STYLE = "BOOL_STYLE";
        String EDGE_ORTH = "EDGE_ORTH";

        Hashtable<String, Object> style;
        mxStylesheet stylesheet = graph.getStylesheet();

        // String nodes are elliptical
        style = new Hashtable<String, Object>();
        style.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_ELLIPSE);
        stylesheet.putCellStyle(STRING_STYLE, style);

        // bool nodes are hexagonal
        style = new Hashtable<String, Object>();
        style.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_HEXAGON);
        stylesheet.putCellStyle(BOOL_STYLE, style);

        // Edges are orthogonal
        style = new Hashtable<String, Object>();
        style.put(mxConstants.STYLE_EDGE, mxConstants.EDGESTYLE_ORTHOGONAL);
        stylesheet.putCellStyle(EDGE_ORTH, style);
    }

    public Object[] getAllVertices() {
        Object[] all = new Object[allVertices.size()];
        int i = 0;
        for (Object vertex : allVertices.values()) {
            all[i++] = vertex;
        }
        return all;
    }

    public void addVertex(String name, Object vertex) {
        allVertices.put(name, vertex);
    }

    public void removeVertex(Object vertex) {
        for (Map.Entry<String, Object> entry : allVertices.entrySet()) {
            if (entry.getValue().equals(vertex)) {
                allVertices.remove(entry.getKey());
                return;
            }
        }
    }

    public List<Object> getSelectedCells() {
        return selectedCells;
    }

    public void addCellSelection(Object o) {
        selectedCells.add(o);
    }

    public void removeVertex(String name) {
        allVertices.remove(name);
    }

    public Object getVertex(String name) {
        return allVertices.get(name);
    }

    public Object[] getAllCells() {
        Object[] arrEdges = graph.getAllEdges(getAllVertices());
        Object[] all = new Object[allVertices.size() + arrEdges.length];
        int i = 0;
        for (Object vertex : allVertices.values()) {
            all[i++] = vertex;
        }
        for (Object edge : arrEdges) {
            all[i++] = edge;
        }
        return all;
    }

    public String getSelectedCellType() {
        mxIGraphModel model = graph.getModel();
        if (model.isVertex(selectedCells)) {
            return "cell";
        } else if (model.isEdge(selectedCells)) {
            return "edge";
        }
        return null;
    }

    /**
     * Automatically arrange graph into layout.
     * Passed index specifies layout configuration. 
     * Presently, only mxParallelEdgeLayout is used.
     * @param layout_index 
     */
    public void executeLayout(int layout_index) {
        final mxGraphLayout layout;
        graph.getModel().beginUpdate();
        try {
            if (layout_index == 1) {
                layout = new mxParallelEdgeLayout(graph);
                layout.execute(graph.getDefaultParent());
            } else if (layout_index == 2) {
                layout = new mxPartitionLayout(graph);
                layout.execute(graph.getDefaultParent());
            } else {
                layout = new mxHierarchicalLayout(graph);
                ((mxHierarchicalLayout) layout).setInterHierarchySpacing(50);
                ((mxHierarchicalLayout) layout).setIntraCellSpacing(50);
                layout.execute(graph.getDefaultParent());
            }
        } finally {
            graph.getModel().endUpdate();
        }
    }

    /**
     * Displays a save dialog for saving network to file.
     * File structure is XML.
     */
    public void saveNetwork() {
        JFileChooser c = new JFileChooser();
        int rVal = c.showSaveDialog(c);
        if (rVal == JFileChooser.APPROVE_OPTION) {
            File file = c.getSelectedFile();
            System.out.println(c.getSelectedFile().toString());
            if (file != null) {
                try {
                    System.out.println("file saved!");
                    model.getBNC().save(file.getCanonicalPath());
                } catch (IOException ex) {

                }
            }
        }
        if (rVal == JFileChooser.CANCEL_OPTION) {
        }
    }

    /** 
     * Displays a load dialog to load network from XML-format file.
     * 
     */
    public void loadNetwork() {
        JFileChooser c = new JFileChooser();
        int rVal = c.showOpenDialog(c);
        if (rVal == JFileChooser.APPROVE_OPTION) {
            File file = c.getSelectedFile();
            if (file != null) {
                try {
                    model.getBNC().loadnm(file.getCanonicalPath(), true);
                    renderNetwork(model.getBNC());
                    executeLayout(1);
                } catch (IOException ex) {

                }
            }
        }
        if (rVal == JFileChooser.CANCEL_OPTION) {

        }
    }

    GraphPanel(BNModel mod) {
        this();
        model = mod;
//        bnc = model.getBNC();
    }

    public GraphPanel() {
        super();
        this.graph = new mxGraph();
        Object parent = this.graph.getDefaultParent();
        this.graph.setMinimumGraphSize(new mxRectangle(0, 0, 800, 600));
        this.graph.setAllowDanglingEdges(false);
        this.graph.setAllowLoops(false);
        this.graph.setCellsEditable(false);
        this.graph.setCellsResizable(false);
        this.graph.setKeepEdgesInBackground(true);
        this.graph.setAllowNegativeCoordinates(true);
        this.setAutoscrolls(true);

        clearNodeCounts();

        Map<String, Object> stil = new HashMap<String, Object>();
        stil.put(mxConstants.STYLE_ROUNDED, false);
//        stil.put(mxConstants.STYLE_EDGE, mxConstants.EDGESTYLE_ORTHOGONAL);
        stil.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_CONNECTOR);
        stil.put(mxConstants.STYLE_ENDARROW, mxConstants.ARROW_CLASSIC);
        stil.put(mxConstants.STYLE_VERTICAL_ALIGN, mxConstants.ALIGN_MIDDLE);
        stil.put(mxConstants.STYLE_ALIGN, mxConstants.ALIGN_CENTER);
        stil.put(mxConstants.STYLE_STROKECOLOR, "#6482B9");
        mxStylesheet foo = new mxStylesheet();
        foo.setDefaultEdgeStyle(stil);
        graph.setStylesheet(foo);
        graphComponent = new mxGraphComponent(this.graph);
        this.add(graphComponent, BorderLayout.CENTER);

    }

    /**
     * Perform Inference. Selected node is set as Query node, and NodeModels
     * which are so specified set as Evidence. Presently only variable
     * elimination is implemented so this is used by default for inference.
     */
    public void doInference() {

        if (queryNode == null) {
            System.err.println("queryNode is null. Returning...");
            return;
        }
        BNet bn = model.getBNC().getBNet();
        CGVarElim ve = new CGVarElim();
        ve.instantiate(bn);
        Query q = ve.makeQuery(queryNode.getVariable());
        CGTable res = (CGTable) ve.infer(q);
        res.display();
    }
    /**
     * Reset nodeCounts count values.
     */
    public void clearNodeCounts(){
        for (String predef : Predef.getVariableTypes()) {
            nodeCounts.put(predef, 0);
        }
    } 
    
    /**
     * Re-renders the network. This involves removing all vertices in View and
     * generating them again from underlying Model.
     *
     * @param bnc
     */
    public void renderNetwork(BNContainer bnc) {
        // Clear graph then repopulate.
        graph.removeCells(graph.getChildVertices(graph.getDefaultParent()));
        
        // Iterate through list of nodes stored in Model and draw them.
        for (BNode node : bnc.getBNet().getNodes()) {
            NodeModel nm = new NodeModel(node);
            Variable var = nm.getVariable();
            if (var != null) {
                String predef = var.getPredef();
                String name = var.getName();
                String params = var.getParams();
                graph.getModel().beginUpdate();
                System.out.println("<<predef: " + predef
                        + ", name: " + name
                        + ", params: " + params + ">>");
                try {
                    String color = (Predef.isEnumerable(predef) ? "yellow" : "orange");
                    Object newvertex = graph.insertVertex(graph.getDefaultParent(),
                            null, name, 50, 50, 100, 50, "ROUNDED;strokeColor=black;fillColor=" + color);
                    this.addVertex(name, newvertex);
                } finally {
                    graph.getModel().endUpdate();
                }
            } else {
                System.err.println("In renderNetwork, node var is null");
            }
        } 
        // Insert edges.
        for (BNode node : bnc.getBNet().getNodes()) {
            NodeModel nm = new NodeModel(node);
            String child_name = nm.getVariable().getName();
            Object child_vertex = this.getVertex(child_name);
            if (nm.getParents() != null) {
                for (EnumVariable parent : nm.getParents()) {
                    String parent_name = parent.getName();
                    Object parent_vertex = this.getVertex(parent_name);
                    graph.getModel().beginUpdate();
                    System.out.println("Inserting edge between " + parent_name + " and " + child_name);
                    try {
                        Object newedge = graph.insertEdge(graph.getDefaultParent(), null, "", parent_vertex, child_vertex);
                    } finally {
                        graph.getModel().endUpdate();
                    }
                }
            } else {
            }
        }
        this.executeLayout(1);
    }

    /**
     * Adds node to model. 
     *
     * @param name
     * @param predef
     * @param params
     */
    public void addNodetoBNC(String name, String predef, String params) {
        if (name == null) {
            name = predef + " node-" + nodeCounts.get(predef);
        }
        // Set default parameters.
        if (params == null) {
            params = predef.equalsIgnoreCase("String") ? "a;b"
                    : predef.equalsIgnoreCase("Number") ? "5"
                    : null;
        }

        String type = Predef.getBNodeType(predef);
        nodeCounts.put(predef, nodeCounts.get(predef) + 1);
        Variable var = Predef.getVariable(name, predef, params);
        NodeModel nm = Predef.getNodeModel(var, new ArrayList<Variable>(), type);

        if (predef.equalsIgnoreCase("Number")) {
            nm.getVariable().setPredef("Number");
            nm.getVariable().setName("Number node-" + (nodeCounts.get("Number") - 1));
            System.out.println("name is:" + "Number node-" + nodeCounts.get("Number"));
        }
        model.getBNC().addNode(nm);
    }
    
    @Override
    // Method called by Observable mxModel NodeModel.
    public void update() {

    }

    @Override
    public void setSubject(Observable sub) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
