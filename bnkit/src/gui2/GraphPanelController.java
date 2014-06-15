/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui2;

import bn.BNode;
import bn.Predef;
import bn.Variable;
import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.swing.handler.mxRubberband;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.util.mxEvent;
import com.mxgraph.util.mxEventObject;
import com.mxgraph.util.mxEventSource;
import com.mxgraph.util.mxResources;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxGraphSelectionModel;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import javax.swing.SwingUtilities;

/**
 *
 * @author jun
 * GraphPanelController implements handlers for GraphPanel,.
 * 
 * Mouse and keyboard input are captured to allow node selection
 * and movement, edge insertion as well as keyboard and mouse events.
 */
public class GraphPanelController {
    final private GraphPanel graphPanel;
    final private BNModel model;
    
    public GraphPanelController(GraphPanel gp, BNModel mdl){
        this.graphPanel = gp;
        this.model = mdl;
    }
    
    public void initHandlers(){
        final mxGraphComponent graphComponent = graphPanel.getGraphComponent();
        final BNContainer bnc = graphPanel.getBNContainer();
        final mxGraph graph = graphPanel.getGraph();
        
        // Enable rubberband (multiple) selection
        mxRubberband rubberband = new mxRubberband(graphComponent);
        graphComponent.setPanning(true);

        // Set tooltips for nodes.
        graphComponent.setToolTips(true);

        // Listener for edge creation
        graphComponent.getConnectionHandler().addListener(mxEvent.CONNECT, new mxEventSource.mxIEventListener() {
            @Override
            public void invoke(Object sender, mxEventObject evt) {
                Object edge = evt.getProperty("cell");
                Object parentNode = ((mxCell) edge).getTerminal(true);
                Object childNode = ((mxCell) edge).getTerminal(false);
                
                NodeModel nmp = bnc.getNodeModel(graph.getLabel(parentNode));
                NodeModel nm = bnc.getNodeModel(graph.getLabel(childNode));
                Variable parentvar = bnc.getVariable(graph.getLabel(parentNode));
                
                // The following pre-conditions exist:
                // - Graph may not have cycles.
                // - Continuous nodes may not have children.
                if (!Predef.isEnumerable(nmp.getVariable().getPredef())){
                    // Number nodes can not have children
                    System.err.println("Continuous nodes may not have children.");
                    graph.removeCells(new Object[]{edge});
                    return;
                }
                
                // Constraints not violated, create edge.
                System.out.println("Created edge: " + edge);
                System.out.println("Between " + parentNode + " and " + childNode);
                
                graphPanel.getSelectedCells().clear();
                graphPanel.getSelectedCells().add(edge);

                // update parent-child relationships
                bnc.addParent(nm, parentvar);
            }
        });

        // Mouse scroll handler for zooming
        graphComponent.getGraphControl().addMouseWheelListener(new MouseAdapter() {
            @Override
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
            @Override
            public void keyPressed(KeyEvent e) {
                int key = e.getKeyCode();
                if (key == KeyEvent.VK_DELETE) {
                    deleteSelected();
                } else if ((key == KeyEvent.VK_S) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    graphPanel.saveNetwork();
                } else if ((key == KeyEvent.VK_O) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    graphPanel.loadNetwork();
                } else if ((key == KeyEvent.VK_EQUALS) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    graphComponent.zoomIn();
                } else if ((key == KeyEvent.VK_MINUS) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    graphComponent.zoomOut();
                }
            }
        });

        // Add a listener for mouse-clicks. 
        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                // If a single node is clicked, clear previous selections and select node.
                if (graphPanel.getSelectedCells() != null) {
                    graphPanel.getSelectedCells().clear();
                }
                graphPanel.getSelectedCells().add(cell);               
                if (cell != null) {

                    // If right-click, open properties dialog
                    if (SwingUtilities.isRightMouseButton(e)) {
                        System.out.println("graph children are: ");
                        for (Object node: graph.getChildVertices(graph.getDefaultParent())){
                            System.out.println(" >"  + ((mxCell) node).getValue() );
                        }
                        
                        NodeModel nm = bnc.getNodeModel(graph.getLabel(cell));

                        // register listener here.
                        NodePropertiesDialog npp = new NodePropertiesDialog(null, true, nm, graphPanel);
                        npp.setVisible(true);
                    } else {
                        // Single-clicked single cell
                        graphComponent.getGraphControl().requestFocus();
                        mxIGraphModel model = graph.getModel();
                        // select whatever is pointed at
                        if (model.isVertex(cell)) {
                            if (e.getClickCount() == 2) {
                                System.out.println("Double-clicked-to-select Vertex=" + graph.getLabel(cell));
                                NodeParamsDialog dialog = new NodeParamsDialog(null, true);
                                NodeModel nm = bnc.getNodeModel(graph.getLabel(cell));
                                dialog.setModel(nm);
                                dialog.setVisible(true);
                            }

                        } else if (graph.getModel().isEdge(cell)) {
                            System.out.println("Clicked-to-select Edge=" + graph.getLabel(cell));
                        }
                    }
                } else {
                    graphPanel.getSelectedCells().clear();
                    System.out.println("deselected");
                    graphComponent.getGraphControl().requestFocus();
                }
            }
        });

        // Add listener for multiple cell selection.
        graph.getSelectionModel().addListener(mxEvent.CHANGE, new mxEventSource.mxIEventListener() {
            @Override
            public void invoke(Object sender, mxEventObject evt) {
                if (sender instanceof mxGraphSelectionModel) {
                    graphPanel.getSelectedCells().clear();
                    for (Object cell : ((mxGraphSelectionModel) sender).getCells()) {
                        System.out.println("selected cell=" + graph.getLabel(cell));
                        graphPanel.getSelectedCells().add(cell);
                        graphComponent.getGraphControl().requestFocus();
                    }
                }
            }
        });
    }
    
    /**
     * Add node to GraphPanel View.
     *
     * @param name
     * @param predef
     * @param params
     */
    public void createNode(String name, String predef, String params) {
        // Set defaults for initialising node...
        mxGraph graph = graphPanel.getGraph();
        if (name == null) {
            name = predef + " node-" + graphPanel.getNodeCounts().get(predef);
        }
        // Set default parameters.
        if (params == null) {
            params = predef.equalsIgnoreCase("String") ? "a;b"
                    : predef.equalsIgnoreCase("Number") ? "5"
                    : null;
        }

        String type = Predef.getBNodeType(predef);
        try {
            String color = (Predef.isEnumerable(predef) ? "yellow" : "orange"); // yellow for enumerable nodes, orange for continuous
            graph.getModel().beginUpdate();
            Variable var = Predef.getVariable(name, predef, params);

            try {
                // Create visible node
                Random rand = new Random();
//                graphPanel.defStyleSheets(graph); // custom vertex and edge styles, presently unused.
                Object newvertex = graph.insertVertex(graph.getDefaultParent(), null, name, 10 + rand.nextInt(50),
                        10 + rand.nextInt(50), 100, 50, "ROUNDED;strokeColor=black;fillColor=" + color);

                graphPanel.addVertex(name, newvertex);

                // 'Select' the new node.
                graphPanel.getSelectedCells().clear();
                graphPanel.addCellSelection(newvertex);
                System.out.println("bnode is: " + Predef.getBNode(var, new ArrayList<Variable>(), type));
                System.out.println("type is: " + type);

                // Add Node to BNodeMap.
//                nm.register(mainFrame);
            } finally {
                graph.getModel().endUpdate();
            }
//            }
        } catch (RuntimeException e) {
        }
    }

    /**
     * Adds node to Model. 
     *
     * @param name
     * @param predef
     * @param params
     */
    public void addNodetoBNC(String name, String predef, String params) {
        Map<String, Integer> nodeCounts = graphPanel.getNodeCounts();
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

    /**
     * Delete the selected nodes from the View and Model.
     */
    public void deleteSelected() {
        mxGraph graph = graphPanel.getGraph();
        Map<String, Integer> nodeCounts = graphPanel.getNodeCounts();
        BNContainer bnc = model.getBNC();
        if (graphPanel.getSelectedCells().isEmpty()) {
            return;
        }
        System.out.println("Delete selection ");
        mxIGraphModel mxModel = graph.getModel();
        for (Object cell : graphPanel.getSelectedCells()) {
            // Edges must be deleted first, otherwise removeParent fails.
            if (graph.getModel().isEdge(cell)) {
                Object child = mxModel.getTerminal(cell, false);
                Object parent = mxModel.getTerminal(cell, true);
                NodeModel childnode = bnc.getNodeModel(graph.getLabel(child));
                Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                bnc.removeParent(childnode, parentvar);
            }

        }
        for (Object cell : graphPanel.getSelectedCells()) {
            if (mxModel.isVertex(cell)) {
                String nodename = graph.getLabel(cell);
                NodeModel node = bnc.getNodeModel(nodename);
                bnc.removeNode(node); 
                nodeCounts.put(node.getVariable().getPredef(), nodeCounts.get(node.getVariable().getPredef()) - 1);
                graphPanel.removeVertex(cell);
            }
        }
        graph.removeCells(graphPanel.getSelectedCells().toArray(new Object[graphPanel.getSelectedCells().size()]));
    }

    /**
     * Delete all nodes from View and Model.
     */
    public void deleteAll() {
        mxGraph graph = graphPanel.getGraph();
        for (Object cell : graph.getChildVertices(graph.getDefaultParent())) {
            graphPanel.addCellSelection(cell);
            graphPanel.removeVertex(cell);
        }
        graph.removeCells(graph.getChildVertices(graph.getDefaultParent()));
        model.getBNC().clear();
        graphPanel.clearNodeCounts();
    }
    
}
