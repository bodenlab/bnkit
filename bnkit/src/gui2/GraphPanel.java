/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNode;
import bn.EnumVariable;
import bn.Predef;
import bn.Variable;
import com.mxgraph.layout.hierarchical.mxHierarchicalLayout;
import com.mxgraph.layout.mxGraphLayout;
import com.mxgraph.layout.mxParallelEdgeLayout;
import com.mxgraph.layout.mxPartitionLayout;
import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxGeometry;
import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.swing.handler.mxRubberband;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.swing.util.mxGraphTransferable;
import com.mxgraph.util.mxConstants;
import com.mxgraph.util.mxEvent;
import com.mxgraph.util.mxEventObject;
import com.mxgraph.util.mxEventSource.mxIEventListener;
import com.mxgraph.util.mxRectangle;
import com.mxgraph.util.mxResources;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxGraphSelectionModel;
import com.mxgraph.view.mxStylesheet;
import java.awt.BorderLayout;
import java.awt.datatransfer.Transferable;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.TransferHandler;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 *
 * @author jun
 */
public class GraphPanel extends JPanel implements Serializable, Observer {

    private mxGraph graph;
    final mxGraphComponent graphComponent;
    private Map<String, Object> allVertices = new HashMap<String, Object>();
    private BNContainer bnc;
    public final Object lastPressedVertex = null;
    public List<Object> selectedCells = new ArrayList<>();
    private List<NodeModel> nodeModels = new ArrayList<>();
    private BNModel model;
    private Map<String, Integer> nodeCounts = new HashMap<>();

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
        return bnc;
    }

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

        // bool nodes are elliptical
        style = new Hashtable<String, Object>();
        style.put(mxConstants.STYLE_SHAPE, mxConstants.SHAPE_HEXAGON);
        stylesheet.putCellStyle(BOOL_STYLE, style);

        style = new Hashtable<String, Object>();
        style.put(mxConstants.STYLE_EDGE, mxConstants.EDGESTYLE_ORTHOGONAL);
        stylesheet.putCellStyle(EDGE_ORTH, style);

    }

    public void setBNContainer(BNContainer bnc) {
        this.bnc = bnc;
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

    public void setLayout(String lay) {
        mxGraph graph = getGraph();
        final mxGraphLayout layout;
        graph.getModel().beginUpdate();
        layout = new mxHierarchicalLayout(graph);
        ((mxHierarchicalLayout) layout).setInterHierarchySpacing(50);
        ((mxHierarchicalLayout) layout).setIntraCellSpacing(50);
        layout.execute(graph.getDefaultParent());

        graph.getModel().endUpdate();

    }

    /**
     * Creates a new Vertex 'node' in the View.
     *
     * @param name
     * @param predef
     * @param params
     */
    public void createNode(String name, String predef, String params) {
        // Set defaults for initialising node...
        if (name == null) {
            name = predef + " node-" + nodeCounts.get(predef);
        }
        // Set default parameters.
        if (params == null) {
            params = predef.equalsIgnoreCase("String") ? "a;b"
                    : predef.equalsIgnoreCase("Number") ? "5" : // this bugs out when real int provided
                    null;
        }

        String type = Predef.getBNodeType(predef);
        try {
            String color = (Predef.isEnumerable(predef) ? "yellow" : "orange"); // yellow for enumerable nodes, orange for continuous
//                String cellStyle = (Predef.parameterName(predef).equals("String") ? "STRING_STYLE" : "BOOL_STYLE");
            graph.getModel().beginUpdate();
            Variable var = Predef.getVariable(name, predef, params);

            try {
                // Create visible node
                Random rand = new Random();
                this.defStyleSheets(graph); // custom vertex and edge styles
                Object newvertex = graph.insertVertex(graph.getDefaultParent(), null, name, 10 + rand.nextInt(50),
                        10 + rand.nextInt(50), 100, 50, "ROUNDED;strokeColor=black;fillColor=" + color);
                
                // TODO: find a new home for this
                this.addVertex(name, newvertex);
                
                // 'Select' the new node.
                selectedCells.clear();
                this.addCellSelection(newvertex);
                System.out.println("bnode is: " + Predef.getBNode(var, new ArrayList<Variable>(), type));
                System.out.println("type is: " + type);
                
                // Add Node to BNodeMap.

//                nm.register(mainFrame);
            } finally {
                graph.getModel().endUpdate();
            }
//            }
        } catch (RuntimeException e) {
//            error_msg = e.getLocalizedMessage();
        }
    }

    public Transferable createTransferableNode(String predef) {
        // Changes in GraphPanel model
//        this.addVertex(name, newvertex);
//        selectedCells.clear();
//        this.addCellSelection(newvertex);

        String name = predef + " node-" + nodeCounts.get(predef);
        
        // Changes in GraphPanel View
        mxRectangle bounds = new mxRectangle(0, 0, 100, 50);
        mxGeometry geometry = new mxGeometry(0, 0, 100, 50);
        geometry.setRelative(false);
        String color = (Predef.isEnumerable(predef) ? "yellow" : "orange");
        mxCell vertex = new mxCell(name, geometry, "ROUNDED;strokeColor=black;fillColor=" + color);
//                vertex.setId(null);
        vertex.setVertex(true);
        vertex.setConnectable(true);
        return new mxGraphTransferable(new Object[]{vertex}, bounds);
    }

    public Object generateNode() {

        return null;
    }

    public void addNodetoBNC(String name, String predef, String params) {
        System.out.println("node added to bnc!!");
        if (name == null) {
            name = predef + " node-" + nodeCounts.get(predef);
        }
        // Set default parameters.
        if (params == null) {
            params = predef.equalsIgnoreCase("String") ? "a;b"
                    : predef.equalsIgnoreCase("Number") ? "5" : // this bugs out when real int provided
                    null;
        }

        String type = Predef.getBNodeType(predef);
        nodeCounts.put(predef, nodeCounts.get(predef) + 1);
        Variable var = Predef.getVariable(name, predef, params);
        if (Predef.getBNodeType(predef).equalsIgnoreCase("CPT")) {
            BNode newBNode = Predef.getBNode(var, new ArrayList<Variable>(), type);
            NodeModel nm = Predef.getNodeModel(var, new ArrayList<Variable>(), type);

            //TODO: investigate this...
            // This should not be necessary. Predef.getBNode is not correctly storing
            // number node predef and name
            if (predef.equalsIgnoreCase("Number")) {
                nm.getVariable().setPredef("Number");
                nm.getVariable().setName("Number node");
            }

            model.getBNC().addNode(newBNode);
            model.getBNC().addNode(nm);
        } else {
            System.out.println("GPT case");
            BNode newBNode = Predef.getBNode(var, new ArrayList<Variable>(), type);
            NodeModel nm = Predef.getNodeModel(var, new ArrayList<Variable>(), type);
        }

    }

    private GraphPanel getGraphPanel() {
        return this;
    }

    public void deleteSelected() {
        if (selectedCells.isEmpty()) {
            return;
        }
        System.out.println("Delete selection ");
        mxIGraphModel model = graph.getModel();
        for (Object cell : selectedCells) {

            // Edges must be deleted first, otherwise removeParent fails.
            if (graph.getModel().isEdge(cell)) {
                Object child = model.getTerminal(cell, false);
                Object parent = model.getTerminal(cell, true);
                System.out.println("deleting edge between " + ((mxCell) parent).getValue()
                        + " and " + ((mxCell) child).getValue());
                BNode childnode = bnc.getNode(graph.getLabel(child));

                Variable parentvar = bnc.getVariable(graph.getLabel(parent));
//                graph.removeCells(new Object[]{cell});
                System.out.println("delete vert:: childnode: " + childnode + " parentvar: " + parentvar);
                bnc.removeParent(childnode, parentvar);
            }

        }
        for (Object cell : selectedCells) {
            if (model.isVertex(cell)) {
                bnc.removeNode(bnc.getNodeModel(graph.getLabel(cell)));
                String nodename = graph.getLabel(cell);
                BNode node = bnc.getNode(nodename);
                nodeCounts.put(node.getVariable().getPredef(), nodeCounts.get(node.getVariable().getPredef()) - 1);
                bnc.removeNode(node);
                removeVertex(cell);
            }
        }
        graph.removeCells(selectedCells.toArray(new Object[selectedCells.size()]));
    }

    public void deleteAll() {
        for (Object cell : graph.getChildVertices(graph.getDefaultParent())) {
            addCellSelection(cell);
            removeVertex(cell);
        }
        deleteSelected();
        graph.removeCells(graph.getChildVertices(graph.getDefaultParent()));
        bnc.clear();
    }

    public void executeLayout(int layout_index) {
        final mxGraph graph = getGraph();
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

    public void saveNetwork() {
        JFileChooser c = new JFileChooser();

        // For now only allow .xml files
        FileFilter filter = new FileNameExtensionFilter("XML file", "xml");
        c.setFileFilter(filter);
        int rVal = c.showSaveDialog(c);
        if (rVal == JFileChooser.APPROVE_OPTION) {
            File file = c.getSelectedFile();
            System.out.println(c.getSelectedFile().toString());
            if (file != null) {
                try {
                    System.out.println("file saved!");
                    bnc.save(file.getCanonicalPath());
                } catch (IOException ex) {

                }
            }
        }
        if (rVal == JFileChooser.CANCEL_OPTION) {
            // do nothing for now
        }
    }

    public void loadNetwork() {
        JFileChooser c = new JFileChooser();

        // For now only allow .xml files
        FileFilter filter = new FileNameExtensionFilter("XML file", "xml");
        c.setFileFilter(filter);
        int rVal = c.showOpenDialog(c);
        if (rVal == JFileChooser.APPROVE_OPTION) {
            File file = c.getSelectedFile();
            if (file != null) {
                try {
//                    bnc.load(file.getCanonicalPath());
                    bnc.loadnm(file.getCanonicalPath(), true);
                    renderNetwork(bnc);
                    this.setLayout("");
                    // also need to add nodes to graph.
                } catch (IOException ex) {

                }
            }
            // open
        }
        if (rVal == JFileChooser.CANCEL_OPTION) {

        }
        // Later replace this with
        // graphPanel.loadNetwork();
    }

    GraphPanel(BNModel mod) {
        this();
        model = mod;
        bnc = model.getBNC();
    }

    public GraphPanel() {
        super();
        this.graph = new mxGraph();
        Object parent = this.graph.getDefaultParent();
        // To have scroll bars by default, set min graph size smaller than frame container
        this.graph.setMinimumGraphSize(new mxRectangle(0, 0, 400, 500));
        this.graph.setMaximumGraphBounds(new mxRectangle(0, 0, 2000, 1600));
        this.graph.setAllowDanglingEdges(false);
        this.graph.setAllowLoops(false);
        this.graph.setCellsEditable(false);
        this.graph.setCellsResizable(false);
        this.graph.setKeepEdgesInBackground(true);
        this.graph.setAllowNegativeCoordinates(true);
        this.setAutoscrolls(true);

        for (String predef : Predef.getVariableTypes()) {
            nodeCounts.put(predef, 0);
        }

        Map<String, Object> stil = new HashMap<String, Object>();
        stil.put(mxConstants.STYLE_ROUNDED, false);
        stil.put(mxConstants.STYLE_EDGE, mxConstants.EDGESTYLE_ORTHOGONAL);
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

        // Enable rubberband (multiple) selection
        mxRubberband rubberband = new mxRubberband(graphComponent);
        graphComponent.setPanning(true);

        // Set tooltips for nodes.
        graphComponent.setToolTips(true);

//        graphComponent.getToolTipText();
        // Listener for edge creation
        graphComponent.getConnectionHandler().addListener(mxEvent.CONNECT, new mxIEventListener() {
            public void invoke(Object sender, mxEventObject evt) {
                Object edge = evt.getProperty("cell");
                Object parentNode = ((mxCell) edge).getTerminal(true);
                Object childNode = ((mxCell) edge).getTerminal(false);

                System.out.println("Created edge: " + edge);
                System.out.println("Between " + parentNode + " and " + childNode);

                selectedCells.clear();
                selectedCells.add(edge);

                // update parent-child relationships
                NodeModel nm = bnc.getNodeModel(graph.getLabel(childNode));
                Variable parentvar = bnc.getVariable(graph.getLabel(parentNode));
                bnc.addParent(nm, parentvar);

//                ((mxCell) childNode).setValue("test remove me");
            }
        });

        // Mouse scroll handler for zooming
        graphComponent.getGraphControl().addMouseWheelListener(new MouseAdapter() {
            public void mouseWheelMoved(MouseWheelEvent e) {
                System.out.println("Mouse scrolled");
                if (e.getWheelRotation() < 0) {
                    graphComponent.zoomIn();
                } else {
                    graphComponent.zoomOut();
                }
                System.out.println(mxResources.get("scale") + ": "
                        + (int) (100 * graphComponent.getGraph().getView().getScale())
                        + "%");
            }
        });

        // Register listeners for key presses.
        graphComponent.getGraphControl().addKeyListener(new KeyAdapter() {
            public void keyDown(KeyEvent e) {
                // required
            }

            public void keyTyped(KeyEvent e) {
                //required
            }

            public void keyPressed(KeyEvent e) {

                int key = e.getKeyCode();
                if (key == KeyEvent.VK_DELETE) {
                    deleteSelected();
                } else if ((key == KeyEvent.VK_S) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    saveNetwork();
                } else if ((key == KeyEvent.VK_O) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    loadNetwork();
                } else if ((key == KeyEvent.VK_EQUALS) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    graphComponent.zoomIn();
                } else if ((key == KeyEvent.VK_MINUS) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    graphComponent.zoomOut();
                }
            }

            public void keyReleased(KeyEvent e) {
            }
        });

        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                // Clear previous selections
                if (selectedCells != null) {
                    selectedCells.clear();
                }
                selectedCells.add(cell);

                if (cell != null) {

                    // If right-click, open properties dialog
                    if (SwingUtilities.isRightMouseButton(e)) {
                        BNode node = bnc.getNode(graph.getLabel(cell));

                        NodeModel nm = bnc.getNodeModel(graph.getLabel(cell));

                        if (nm == null) {
                            System.out.println("getLabel is: " + graph.getLabel(cell));
                        }

                        // register listener here.
                        NodePropertiesDialog npp = new NodePropertiesDialog(null, true, nm, getGraphPanel());
                        npp.setVisible(true);
                    } else {
                        // Single-clicked single cell
                        graphComponent.getGraphControl().requestFocus();
                        mxIGraphModel model = graph.getModel();
                        // select whatever is pointed at
                        if (model.isVertex(cell)) {
                            System.out.println("Cell is: " + cell
                                    + ", Selected node is:" + bnc.getNode(graph.getLabel(cell)));

                            // now add this to MainJFrame's propertiesPanel
//                        nodeProps.setBNode(bnc.getNode(graph.getLabel(cell)));
//                        nodeProps.updateDisplay();
                            if (e.getClickCount() == 2) {
                                System.out.println("Double-clicked-to-select Vertex=" + graph.getLabel(cell));
                                BNode node = bnc.getNode(graph.getLabel(cell));
                                NodeParamsDialog dialog = new NodeParamsDialog(null, true);
                                NodeModel nm = bnc.getNodeModel(graph.getLabel(cell));

                                //TODO: find a way to correctly setmodel
                                dialog.setModel(nm);
                                dialog.setVisible(true);
                            }

                        } else if (graph.getModel().isEdge(cell)) {
                            System.out.println("Clicked-to-select Edge=" + graph.getLabel(cell));
                        }
                    }
                } else {
                    selectedCells.clear();
                    System.out.println("deselected");
                    graphComponent.getGraphControl().requestFocus();
                }
            }
        });

        // Multiple selection
        // This handles general click events (by mouse selection or Ctrl+A)
        graph.getSelectionModel().addListener(mxEvent.CHANGE, new mxIEventListener() {
            @Override
            public void invoke(Object sender, mxEventObject evt) {
                if (sender instanceof mxGraphSelectionModel) {
                    selectedCells.clear();
                    for (Object cell : ((mxGraphSelectionModel) sender).getCells()) {
                        System.out.println("selected cell=" + graph.getLabel(cell));
                        selectedCells.add(cell);
                        graphComponent.getGraphControl().requestFocus();
                    }
                }
            }
        });

        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            public void mouseReleased(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                if (cell != null) {
                }
            }

        });

    }

    /**
     * Re-renders the network. This involves removing all vertices in View and
     * generating them again from underlying Model.
     *
     * @param bnc
     */
    public void renderNetwork(BNContainer bnc) {
        // Clear graph then repopulate.
        graph.removeCells(this.getAllCells());

        for (BNode node : bnc.getBNetnm().getNodes()) {
            // Just iterate through elements in nodems??

            NodeModel nm = new NodeModel(node);
            // check if parents is null beforehand
//            NodeModel nm = new NodeModel( node.getVariable(), node.getParents());

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
        } // variables done... now connect them
        for (BNode node : bnc.getBNetnm().getNodes()) {
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
//                System.out.println("No parents :(");
            }
        }
        this.executeLayout(1);
    }

    @Override
    // Method called by Observable model NodeModel.
    public void update() {

    }

    @Override
    public void setSubject(Observable sub) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
