/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.Predef;
import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxGeometry;
import com.mxgraph.swing.util.mxGraphTransferable;
import com.mxgraph.util.mxRectangle;
import java.awt.Color;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.border.Border;

/**
 *
 * @author jun 
 * Extension of JLabel which can be dragged and dropped onto GraphPanel to
 * add nodes to the network.
 */
public class NodeLabel extends JLabel implements Transferable,
        DragSourceListener, DragGestureListener {

    //marks this NodeLabel as the source of the Drag

    private DragSource source;
    private GraphPanel graphPanel;
    private Transferable transferable;

    public NodeLabel(gui2.GraphPanel graph, String msg) {
        this(msg);
        graphPanel = graph;
    }

    private NodeLabel(String message) {
        super(message);

        //The Drag will copy the DnDButton rather than moving it
        source = new DragSource();
        source.createDefaultDragGestureRecognizer(this, DnDConstants.ACTION_COPY, this);
        Border border = BorderFactory.createDashedBorder(Color.BLACK, 2, 2);
        getThis().setBorder(border);
    }

    private NodeLabel getThis() {
        return this;
    }

     //The DataFlavor is a marker to let the DropTarget know how to
    //handle the Transferable
    public DataFlavor[] getTransferDataFlavors() {
        return new DataFlavor[]{new DataFlavor(mxCell.class, "mxCell")};

    }

    public boolean isDataFlavorSupported(DataFlavor flavor) {
        return true;
    }

    public Object getTransferData(DataFlavor flavor) {
        return this;
    }

    public void dragEnter(DragSourceDragEvent dsde) {
    }

    public void dragOver(DragSourceDragEvent dsde) {
    }

    public void dropActionchanged(DragSourceDragEvent dsde) {
    }

    public void dragExit(DragSourceEvent dse) {
    }

    /**
     * When NodeLabel is dropped, checks whether drop is a success and handles
     * accordingly.
     *
     * @param dsde
     */
    public void dragDropEnd(DragSourceDropEvent dsde) {
        if (dsde.getDropSuccess()) {
            // Drop success, add Node to Model.
            String type = getThis().getText();
            graphPanel.addNodetoBNC(null, type, null);
            repaint();
        } else {
        }
    }

    /**
     * Triggers when a NodeLabel is clicked and dragged. Generates the visual
     * Node to be added to GraphPanel.
     * @param dge 
     */
    public void dragGestureRecognized(DragGestureEvent dge) {
        transferable = createTransferableNode(getThis().getText());
        source.startDrag(dge, DragSource.DefaultCopyDrop, transferable, this);
    }

    @Override
    public void dropActionChanged(DragSourceDragEvent dsde) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Renders dragged label as a Node.
     * @param predef
     * @return 
     */
    public Transferable createTransferableNode(String predef) {

        String name = predef + " node-" + graphPanel.getNodeCounts().get(predef);

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

}
