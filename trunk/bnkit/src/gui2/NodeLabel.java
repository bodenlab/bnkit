/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui2;

import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxGeometry;
import com.mxgraph.swing.util.mxGraphTransferable;
import com.mxgraph.util.mxRectangle;
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
import javax.swing.JComponent;
import javax.swing.JLabel;

/**
 *
 * @author Jun
 * JLabel which can be dragged and dropped to return an mxCell.
 */
public class NodeLabel extends JLabel implements Transferable,
DragSourceListener, DragGestureListener{
    //marks this JLabel as the source of the Drag
    private DragSource source;
    private gui2.GraphPanel graphPanel;
    private Transferable transferable;
    
    public NodeLabel(gui2.GraphPanel graph, String msg){
        this(msg);
        graphPanel = graph;
    }

    private NodeLabel(String message){
        super(message);

        //The Drag will copy the DnDButton rather than moving it
        source = new DragSource();
        source.createDefaultDragGestureRecognizer(this, DnDConstants.ACTION_COPY, this);
    }

    private NodeLabel getThis(){
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

    public void dragEnter(DragSourceDragEvent dsde) {}
    public void dragOver(DragSourceDragEvent dsde) {}
    public void dropActionchanged(DragSourceDragEvent dsde) {}
    public void dragExit(DragSourceEvent dse) {}
    
    public void dragDropEnd(DragSourceDropEvent dsde) {
        if (dsde.getDropSuccess()){
            String type = getThis().getText();
            graphPanel.addNodetoBNC(null, type, null);
            repaint();
        } else {
        }
    }

    //when a DragGesture is recognized, initiate the Drag
    public void dragGestureRecognized(DragGestureEvent dge) {
        transferable = graphPanel.createTransferableNode(getThis().getText());
        source.startDrag(dge, DragSource.DefaultCopyDrop, transferable, this);
    }

    protected Transferable createTransferable(JComponent comp) {
                NodeLabel parentComponent = (NodeLabel) comp;
                String type = parentComponent.getText();
                mxRectangle bounds = new mxRectangle(0, 0, 100, 50);
                mxGeometry geometry = new mxGeometry(0, 0, 100, 50);
                geometry.setRelative(false);
                mxCell vertex = new mxCell(type, new mxGeometry(0, 0, 100, 50), null);
//                vertex.setId(null);
                vertex.setVertex(true);
                vertex.setConnectable(true);

                // TODO is it neccessary to set the parent here??
//                Object parent = graphPanel.getGraphComponent().getGraph().getDefaultParent();
//                vertex.setParent((mxICell) parent);
                return new mxGraphTransferable(new Object[]{vertex}, bounds);
            
        }
    @Override
    public void dropActionChanged(DragSourceDragEvent dsde) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}