/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNode;
import bn.Variable;
import com.mxgraph.layout.hierarchical.mxHierarchicalLayout;
import com.mxgraph.layout.mxGraphLayout;
import com.mxgraph.layout.mxParallelEdgeLayout;
import com.mxgraph.layout.mxPartitionLayout;
import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.swing.handler.mxRubberband;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.util.mxConstants;
import com.mxgraph.util.mxEvent;
import com.mxgraph.util.mxEventObject;
import com.mxgraph.util.mxEventSource.mxIEventListener;
import com.mxgraph.util.mxRectangle;
import com.mxgraph.util.mxResources;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxGraphSelectionModel;
import com.mxgraph.view.mxStylesheet;
import static gui2.StyleExample.NEW_CUSTOM_EDGE;
import java.awt.BorderLayout;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import javax.swing.JPanel;

/**
 *
 * @author mikael
 */
public class MyGraphPanel extends JPanel implements Serializable {

    private mxGraph graph;
    final mxGraphComponent graphComponent;
    private Map<String, Object> allVertices = new HashMap<String, Object>();
    private BNContainer bnc;
    public final Object lastPressedVertex = null;
    public List<Object> selectedCell = new ArrayList<Object>();
    private int mousex;
    private int mousey;

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

    public void removeAll() {
        graph.removeCells(graph.getChildVertices(graph.getDefaultParent()));
    }

    public void setSelected(Object o) {
        selectedCell.add(o);
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
        if (model.isVertex(selectedCell)) {
            return "cell";
        } else if (model.isEdge(selectedCell)) {
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

    public void deleteSelected() {
        if (selectedCell.isEmpty()) return;
        System.out.println("Delete selection ");
        for (Object cell : selectedCell) {
            mxIGraphModel model = graph.getModel();

            if (model.isVertex(cell)) {
                String nodename = graph.getLabel(cell);
                graph.removeCells(new Object[]{cell});
                BNode node = bnc.getNode(nodename);
                bnc.removeNode(node);
                removeVertex(cell);
            } else if (graph.getModel().isEdge(cell)) {
                Object child = model.getTerminal(cell, false);
                Object parent = model.getTerminal(cell, true);
                BNode childnode = bnc.getNode(graph.getLabel(child));
                Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                graph.removeCells(new Object[]{cell});
                bnc.removeParent(childnode, parentvar);
            }
        }
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

    public MyGraphPanel() {
        super();
        this.graph = new mxGraph();
        Object parent = this.graph.getDefaultParent();
        // To have scroll bars by default, set min graph size smaller than frame container
        this.graph.setMinimumGraphSize(new mxRectangle(0, 0, 485, 565));
        this.graph.setMaximumGraphBounds(new mxRectangle(0, 0, 1000, 800));
        this.graph.setAllowDanglingEdges(false);
        this.graph.setAllowLoops(false);
        this.graph.setCellsEditable(false);
        this.graph.setCellsResizable(false);
        this.graph.setKeepEdgesInBackground(true);
        this.graph.setAllowNegativeCoordinates(true);

        //this.graph.getView().setTranslate(new mxPoint(500,400)); // middle?
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

        // enable rubberband selection
        mxRubberband rubberband = new mxRubberband(graphComponent);
        // Selection Listener
//        graphComponent.getGraphControl().addGraphSelectionListener(new MyListener()
//        {
//            
//        })
        graphComponent.setPanning(true);

        // handler for edge creation
        graphComponent.getConnectionHandler().addListener(mxEvent.CONNECT, new mxIEventListener() {
            public void invoke(Object sender, mxEventObject evt) {
                System.out.println("Selected edge" + evt.getProperty("cell"));
                setSelected(evt.getProperty("cell"));
            }
        });
        
       

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

        graphComponent.getGraphControl().addKeyListener(new KeyAdapter() {

            public void keyDown(KeyEvent e) {
                // required
            }

            public void keyTyped(KeyEvent e) {
                //required
            }

            public void keyPressed(KeyEvent e) {
                System.out.println("a key was pressed");
                int key = e.getKeyCode();
                if (key == KeyEvent.VK_DELETE) {
                    Object cell = graphComponent.getCellAt(mousex, mousey);
//                   System.out.println(graphComponent.);
                    mxIGraphModel model = graph.getModel();

                    System.out.println("Delete pressed ");
                    if (model.isVertex(cell)) {
                        String nodename = graph.getLabel(cell);
                        System.out.println("Delete-pressed-to-delete Vertex: " + nodename);
                        graph.removeCells(new Object[]{cell});
                        BNode node = bnc.getNode(nodename);
                        bnc.removeNode(node);
                        removeVertex(cell);
                    } else if (graph.getModel().isEdge(cell)) {
                        System.out.println("Delete-pressed-to-delete Edge: " + graph.getLabel(cell));
                        Object child = model.getTerminal(cell, false);
                        Object parent = model.getTerminal(cell, true);
                        BNode childnode = bnc.getNode(graph.getLabel(child));
                        Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                        graph.removeCells(new Object[]{cell});
                        bnc.removeParent(childnode, parentvar);
                    }
                }
            }
        });

        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {

                // TODO: Need to make changes in MainJFrame
                // show side pane, etc.
                //testing
                nodePropertiesPane nn = new nodePropertiesPane();
//                mxStylesheet sty = graph.getStylesheet().getDefaultVertexStyle();

                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                // Clear previous selections
                if (selectedCell != null) {
                    selectedCell.clear();
                }
                selectedCell.add(cell);
                for (Object cells : selectedCell) {
                    System.out.println("selected cells:" + cells);
                }

                if (cell != null) {
                    System.out.println("Single clicked: " + cell.toString());
                    mxIGraphModel model = graph.getModel();
                    // select whatever is pointed at
                    if (model.isVertex(cell)) {
                        if (e.getClickCount() == 2) {
                            System.out.println("Double-clicked-to-select Vertex=" + graph.getLabel(cell));
                            BNode node = bnc.getNode(graph.getLabel(cell));
                            NodeParamsDialog dialog = new NodeParamsDialog(null, true);
                            dialog.setModel(node);
                            dialog.setVisible(true);
                        }
                    } else if (graph.getModel().isEdge(cell)) {
                        System.out.println("Clicked-to-select Edge=" + graph.getLabel(cell));
                    }
                } else {
                    selectedCell.clear();
                    System.out.println("deselected");
                }
            }
        });

        //Multiple selection
        // This handles general click events (by mouse selection or Ctrl+A)
        graph.getSelectionModel().addListener(mxEvent.CHANGE, new mxIEventListener() {
            @Override
            public void invoke(Object sender, mxEventObject evt) {
                System.out.println("evt.toString() = " + evt.toString());
                System.out.println("Selection in graph component");
                if (sender instanceof mxGraphSelectionModel) {
                    for (Object cell : ((mxGraphSelectionModel) sender).getCells()) {
                        System.out.println("cell=" + graph.getLabel(cell));
                        //TODO: SELECT MULTI CELLS
                        selectedCell.add(cell);

                    }
                }

                for (Object o : selectedCell) {
                    System.out.println("Obj: " + o);
                }
            }
        });

        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            public void mouseReleased(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                if (cell != null) {
                    mxIGraphModel model = graph.getModel();
                    // select whatever is pointed at
                    if (model.isVertex(cell)) {
                        System.out.println("Released-to-select Vertex=" + graph.getLabel(cell));
                        // could arrive here because the user wants to...
                        // 1. connect TO this node (so add parent to it)
                        // 2. just selected it, maybe moved it (so ignore)
                        for (int i = 0; i < model.getEdgeCount(cell); i++) {
                            Object edge = model.getEdgeAt(cell, i);
                            Object child = model.getTerminal(edge, false);
                            Object parent = model.getTerminal(edge, true);
                            if (parent != cell) {
                                System.out.println("\tParent (#" + i + ")=" + graph.getLabel(parent));
                                BNode childnode = bnc.getNode(graph.getLabel(child));
                                Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                                bnc.addParent(childnode, parentvar);
                            }
                        }
                    } else if (graph.getModel().isEdge(cell)) {
                        System.out.println("Released-to-select Edge=" + graph.getLabel(cell));
                        // could arrive here because the user wants to...
                        // 1. disconnect the edge (so need to update child)
                        // 2. just selected it (so ignore)
                        Object child = model.getTerminal(cell, false);
                        Object parent = model.getTerminal(cell, true);
                        BNode childnode = bnc.getNode(graph.getLabel(child));
                        Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                        bnc.removeParent(childnode, parentvar);
                    }
                }
            }

        });
        graphComponent.getGraphControl().addMouseMotionListener(new MouseAdapter() {
            public void mouseDragged(MouseEvent e) {
                // Track mouse location...
                mousex = e.getX();
                mousey = e.getY();
            }
        });

    }

}
