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
import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxICell;
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
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import javax.swing.TransferHandler;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

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
    private final nodePropertiesPane nodeProps = new nodePropertiesPane();

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

    public void addCellSelection(Object o) {
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
        if (selectedCell.isEmpty()) {
            return;
        }
        System.out.println("Delete selection ");
        mxIGraphModel model = graph.getModel();
        for (Object cell : selectedCell) {

            // Edges must be deleted first, otherwise removeParent fails.
            if (graph.getModel().isEdge(cell)) {
                Object child = model.getTerminal(cell, false);
                Object parent = model.getTerminal(cell, true);
                System.out.println("deleting edge between " + parent + " and " + child);
                BNode childnode = bnc.getNode(graph.getLabel(child));
                Variable parentvar = bnc.getVariable(graph.getLabel(parent));
//                graph.removeCells(new Object[]{cell});
                System.out.println("delete vert:: childnode: " + childnode + " parentvar: " + parentvar);
                bnc.removeParent(childnode, parentvar);
            }

        }
        for (Object cell : selectedCell) {
            if (model.isVertex(cell)) {
                String nodename = graph.getLabel(cell);
//                graph.removeCells(new Object[]{cell});
                BNode node = bnc.getNode(nodename);
                bnc.removeNode(node);
                removeVertex(cell);
            }
        }
        graph.removeCells(selectedCell.toArray(new Object[selectedCell.size()]));
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

    public void saveNetwork(){
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
    
    public void loadNetwork(){
         JFileChooser c = new JFileChooser();

        // For now only allow .xml files
      FileFilter filter = new FileNameExtensionFilter("XML file", "xml");
      c.setFileFilter(filter);
        int rVal = c.showOpenDialog(c);
        if (rVal == JFileChooser.APPROVE_OPTION) {
            File file = c.getSelectedFile();
            if (file != null) {
                try {
                    bnc.load(file.getCanonicalPath());
//                    renderNetwork(bnc);
                    setLayout("");
                } catch (IOException ex) {

                }
            }
            // open
        }
        if (rVal == JFileChooser.CANCEL_OPTION) {
            // do nothing for now
        }
    }
    
    public MyGraphPanel() {
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

        // refactor this later
        // add nodeProps to parent frame's propertiespane
//        MainJFrame mainFrame = (MainJFrame) SwingUtilities.getUnwrappedParent(this);
//        javax.swing.JPanel panel = mainFrame.getPanelContainer();
//        panel.add(nodeProps);
//        panel.validate();
//        this.getParent().getComponents();
        nodeProps.setVisible(false);

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

        // Enable rubberband (multiple) selection
        mxRubberband rubberband = new mxRubberband(graphComponent);
        graphComponent.setPanning(true);

        // Listener for edge creation
        graphComponent.getConnectionHandler().addListener(mxEvent.CONNECT, new mxIEventListener() {
            public void invoke(Object sender, mxEventObject evt) {
                Object edge = evt.getProperty("cell");
                Object parentNode = ((mxCell) edge).getTerminal(true);
                Object childNode = ((mxCell) edge).getTerminal(false);

                System.out.println("Created edge: " + edge);
                System.out.println("Between " + parentNode + " and " + childNode);

                selectedCell.clear();
                selectedCell.add(edge);

                // update parent-child relationships
                BNode childnode = bnc.getNode(graph.getLabel(childNode));
                Variable parentvar = bnc.getVariable(graph.getLabel(parentNode));
                bnc.addParent(childnode, parentvar);

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

        // Key press handler for graph
        // TODO: if global solution required try this
        // http://stackoverflow.com/a/1379517
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
                    System.err.println("opening file using ctrl+o, must manually render layout.");
                } else if ((key == KeyEvent.VK_C) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    System.out.println("ctrl+c");
                    TransferHandler.getCopyAction();
                } else if ((key == KeyEvent.VK_V) && ((e.getModifiers() & KeyEvent.CTRL_MASK) != 0)) {
                    System.out.println("ctrl+v");
                }
            }

            public void keyReleased(KeyEvent e) {
            }
        });

        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                // TODO: Need to make changes in MainJFrame
                // show side pane, etc.

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
                    // Single-clicked single cell
                    graphComponent.getGraphControl().requestFocus();
                    mxIGraphModel model = graph.getModel();
                    // select whatever is pointed at
                    if (model.isVertex(cell)) {
                        System.out.println("Selected node is:" + bnc.getNode(graph.getLabel(cell)));

                        // now add this to MainJFrame's propertiesPanel
//                        nodeProps.setBNode(bnc.getNode(graph.getLabel(cell)));
                        nodeProps.updateDisplay();

                        if (e.getClickCount() == 2) {
                            System.out.println("Double-clicked-to-select Vertex=" + graph.getLabel(cell));
                            System.out.print(bnc.getBNet());
                            BNode node = bnc.getNode(graph.getLabel(cell));
                            NodeParamsDialog dialog = new NodeParamsDialog(null, true);

                            //TODO: find a way to correctly setmodel
                            dialog.setModel(node);
                            dialog.setVisible(true);
                        }
                    } else if (graph.getModel().isEdge(cell)) {
                        System.out.println("Clicked-to-select Edge=" + graph.getLabel(cell));
                    }
                } else {
                    selectedCell.clear();
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
                    selectedCell.clear();
                    for (Object cell : ((mxGraphSelectionModel) sender).getCells()) {
                        System.out.println("selected cell=" + graph.getLabel(cell));
                        selectedCell.add(cell);
                        graphComponent.getGraphControl().requestFocus();
                    }
                }
            }
        });

        graphComponent.getGraphControl().addMouseListener(new MouseAdapter() {
            public void mouseReleased(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                if (cell != null) {
                    mxIGraphModel model = graph.getModel();

                }
            }

        });

    }

}
