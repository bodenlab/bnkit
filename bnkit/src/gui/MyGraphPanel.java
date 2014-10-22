/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui;

import bn.BNode;
import bn.Predef;
import dat.Variable;
import com.mxgraph.layout.hierarchical.mxHierarchicalLayout;
import com.mxgraph.layout.mxCompactTreeLayout;
import com.mxgraph.layout.mxGraphLayout;
import com.mxgraph.layout.mxParallelEdgeLayout;
import com.mxgraph.layout.mxPartitionLayout;
import com.mxgraph.layout.mxStackLayout;
import com.mxgraph.model.mxGraphModel;
import com.mxgraph.model.mxIGraphModel;
import java.beans.*;
import java.io.Serializable;
import javax.swing.JPanel;

import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.util.mxConstants;
import com.mxgraph.util.mxPoint;
import com.mxgraph.util.mxRectangle;
import com.mxgraph.view.mxGraph;
import com.mxgraph.view.mxStylesheet;
import java.awt.BorderLayout;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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
    
    public void setBNContainer(BNContainer bnc) {
        this.bnc = bnc;
    }
    
    public Object[] getAllVertices() {
        Object[] all = new Object[allVertices.size()];
        int i = 0;
        for (Object vertex : allVertices.values())
            all[i ++] = vertex;
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
        for (Object vertex : allVertices.values())
            all[i ++] = vertex;
        for (Object edge : arrEdges)
            all[i ++] = edge;
        return all;
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
                ((mxHierarchicalLayout)layout).setInterHierarchySpacing(50);
                ((mxHierarchicalLayout)layout).setIntraCellSpacing(50);
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
        this.graph.setMinimumGraphSize(new mxRectangle(0,0,900,700));
        this.graph.setMaximumGraphBounds(new mxRectangle(0,0,1000,800));
        this.graph.setAllowDanglingEdges(false);
        this.graph.setAllowLoops(false);
        this.graph.setCellsEditable(false);
        this.graph.setCellsResizable(false);
        this.graph.setKeepEdgesInBackground(true);
        this.graph.setAllowNegativeCoordinates(true);
        
        //this.graph.getView().setTranslate(new mxPoint(500,400)); // middle?

        Map<String, Object> stil = new HashMap<String, Object>();
        stil.put(mxConstants.STYLE_ROUNDED, true);
        stil.put(mxConstants.STYLE_EDGE, mxConstants.EDGESTYLE_ELBOW);
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
        
        graphComponent.getGraphControl().addMouseListener(new MouseAdapter()
        {
            public void mouseClicked(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                if (cell != null) {
                    mxIGraphModel model = graph.getModel();
                    if (MainJFrame.isDeleteMode()) { // delete, but we only do that when released
                        ;
                    } else { // select whatever is pointed at
                        if (model.isVertex(cell)) {
                            if (e.getClickCount() == 2) {
                                System.out.println("Double-clicked-to-select Vertex="+graph.getLabel(cell));
                                BNode node = bnc.getNode(graph.getLabel(cell));
                                NodeParamsDialog dialog = new NodeParamsDialog(null, true);
                                dialog.setModel(node);
                                dialog.setVisible(true);
                            }
                        } else if (graph.getModel().isEdge(cell)) {
                            System.out.println("Clicked-to-select Edge="+graph.getLabel(cell));
                        }
                    }
                }
            }
        });
        graphComponent.getGraphControl().addMouseListener(new MouseAdapter()
        {
            public void mouseReleased(MouseEvent e) {
                Object cell = graphComponent.getCellAt(e.getX(), e.getY());
                if (cell != null) {
                    mxIGraphModel model = graph.getModel();
                    if (MainJFrame.isDeleteMode()) { // delete whatever is pointed at
                        if (model.isVertex(cell)) {
                            String nodename = graph.getLabel(cell);
                            System.out.println("Released-to-delete Vertex=" + nodename);
                            graph.removeCells(new Object[] {cell});
                            BNode node = bnc.getNode(nodename);
                            bnc.removeNode(node);
                            removeVertex(cell);
                        } else if (graph.getModel().isEdge(cell)) {
                            System.out.println("Released-to-delete Edge="+graph.getLabel(cell));
                            Object child = model.getTerminal(cell, false);
                            Object parent = model.getTerminal(cell, true);
                            BNode childnode = bnc.getNode(graph.getLabel(child));
                            Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                            graph.removeCells(new Object[] {cell});
                            bnc.removeParent(childnode, parentvar);
                        }                    
                    } else { // select whatever is pointed at
                        if (model.isVertex(cell)) {
                            System.out.println("Released-to-select Vertex="+graph.getLabel(cell));
                            // could arrive here because the user wants to...
                            // 1. connect TO this node (so add parent to it)
                            // 2. just selected it, maybe moved it (so ignore)
                            for (int i = 0; i < model.getEdgeCount(cell); i ++) {
                                Object edge = model.getEdgeAt(cell, i);
                                Object child = model.getTerminal(edge, false);
                                Object parent = model.getTerminal(edge, true);
                                if (parent != cell) {
                                    System.out.println("\tParent (#"+i+")="+graph.getLabel(parent));
                                    BNode childnode = bnc.getNode(graph.getLabel(child));
                                    Variable parentvar = bnc.getVariable(graph.getLabel(parent));
                                    bnc.addParent(childnode, parentvar);
                                }
                            }
                        } else if (graph.getModel().isEdge(cell)) {
                            System.out.println("Released-to-select Edge="+graph.getLabel(cell));
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
            }
        });
        graphComponent.getGraphControl().addMouseMotionListener(new MouseAdapter()
        {
            public void mouseDragged(MouseEvent e) {
            }
        });
        

    }

}
