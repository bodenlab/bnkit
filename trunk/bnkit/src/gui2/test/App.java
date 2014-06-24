package gui2.test;



import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Point;
import java.awt.datatransfer.Transferable;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.TransferHandler;
import javax.swing.border.EmptyBorder;

import com.mxgraph.model.mxCell;
import com.mxgraph.model.mxGeometry;
import com.mxgraph.model.mxICell;
import com.mxgraph.swing.mxGraphComponent;
import com.mxgraph.swing.util.mxGraphTransferable;
import com.mxgraph.util.mxRectangle;
import com.mxgraph.view.mxGraph;

@SuppressWarnings("serial")
public class App extends JFrame {

    private JPanel contentPane;
    private JTextField txtDragHere;
    private mxGraphComponent mxGraphComponent;

    private mxICell parentCell;

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    App frame = new App();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }

    /**
     * Create the frame.
     */
    public App() {
        initComponents();
    }

    private void initComponents() {
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(100, 100, 450, 300);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        contentPane.setLayout(new BorderLayout(0, 0));
        setContentPane(contentPane);
        contentPane.add(getTxtDragHere(), BorderLayout.NORTH);
        contentPane.add(getGraphComponent(), BorderLayout.CENTER);
    }

    private JTextField getTxtDragHere() {
        if (txtDragHere == null) {
            txtDragHere = new JTextField();
            txtDragHere.setEditable(false);
            txtDragHere.setText("Mark & Drag Here");
            txtDragHere.setColumns(10);
            txtDragHere.setDragEnabled(true);
            txtDragHere.setTransferHandler(new TextFieldTransferHandler());
        }
        return txtDragHere;
    }

    private mxGraphComponent getGraphComponent() {
        if (mxGraphComponent == null) {
            mxGraphComponent = new mxGraphComponent(new mxGraph()) {

                public Object[] importCells(Object[] cells, double dx, double dy, Object target, Point location) {

                    // super code
                    Object[] objList = super.importCells( cells, dx, dy, target, location);

                    try {

                        mxICell dragCell = (mxICell) objList[0];
                        Object parent = getGraphComponent().getGraph().getDefaultParent();
                        getGraphComponent().getGraph().insertEdge(parent, null, "Edge", parentCell, dragCell);

                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                    return objList;
                  }

            };

            Object parent = mxGraphComponent.getGraph().getDefaultParent();

            parentCell = (mxICell) mxGraphComponent.getGraph().
                    insertVertex(parent, null, "Parent", 10, 10, 50, 50);
        }
        return mxGraphComponent;
    }

    protected class TextFieldTransferHandler extends TransferHandler {

        @Override
        public int getSourceActions(JComponent comp) {
            return COPY;
        }

        @Override
        protected Transferable createTransferable(JComponent comp) {

            if (comp instanceof JTextField) {
                JTextField sourceComp = (JTextField) comp;

                mxRectangle bounds = new mxRectangle(0, 0, 100, 50);

                mxGeometry geometry = new mxGeometry(0, 0, 100, 50);
                geometry.setRelative(false);

                mxCell vertex = new mxCell(sourceComp.getText(), geometry, null);
                vertex.setId(null);
                vertex.setVertex(true);
                vertex.setConnectable(true);

                // TODO is it neccessary to set the parent here??
                Object parent = getGraphComponent().getGraph().getDefaultParent();
                vertex.setParent((mxICell) parent);

                return new mxGraphTransferable(new Object[] { vertex }, bounds);
            }
            return null;
        }

    }
}