/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.EnumDistrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.Enumerable;
import java.util.List;
import javax.swing.table.AbstractTableModel;

/**
 *
 * @author mikael
 * modified by jun
 */
public class NodeParamsDialog extends javax.swing.JDialog {

    MyTableModel myTableModel;

    /**
     * Creates new form NodeParamsDialog
     */
    public NodeParamsDialog(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        initComponents();
    }

    public void setModel(NodeModel node) {
        if (node == null) {
            System.out.println("node null, do nothing");
            return;
        }
        myTableModel = new MyTableModel(node);
        myParamsTable.setModel(myTableModel);
        this.setTitle(node.getBNode().getVariable().getPredef() + "(" + node.getType() + ") " +
                "node " + "\"" + node.getName() + "\""); // add information about Q/E/I state
    }

    class MyTableModel extends AbstractTableModel {

        private final NodeModel node;
        private EnumTable table;

        public MyTableModel(NodeModel nd) {
            this.node = nd;
            this.table = node.getTable();
            System.out.println("@@Printing node");
            node.print();
            if (table == null) {
                System.out.println("node table is null");
            } else {
                System.out.println("In MyTableModel, print table");
                table.display();
            }
        }

        public int getColumnCount() {
            if (table != null) {
                return node.getParents().size() + 1;
            } else {
                return 1;
            }
        }

        public int getRowCount() {
            if (table != null) {
                return table.getSize();
            } else {
                return 1;
            }
        }

        public String getColumnName(int col) {
            if (table != null) {
                List<EnumVariable> parents = node.getParents();
                if (col < parents.size()) {
                    return parents.get(col).getName();
                }
            }
            return node.getVariable().getName();
        }

        public Object getValueAt(int row, int col) {
            if (table != null) {
                if (col < node.getParents().size()) {
                    Object[] key = table.getKey(row);
                    return key[col];
                } else {
                    return table.getValue(row);
                }
            } else {
                return node.getDistrib();
            }
        }

        public Class getColumnClass(int c) {
            Object value = this.getValueAt(0, c);
            // By default, booleans will be rendered as checkbox. Override this
            // behaviour
//            if (value.getClass().equals(boolean.class)){
//                return Object.class;
//            }
            return (value == null ? Object.class : value.getClass());
        }

        /*
         * Don't need to implement this method unless your table's
         * editable.
         */
        public boolean isCellEditable(int row, int col) {
            if (table != null) {
                if (col < node.getParents().size()) {
                    return false;
                }
            }
            return true;
        }

        /*
         * Don't need to implement this method unless your table's
         * data can change.
         */
        public void setValueAt(Object value, int row, int col) {
            if (table == null) {
                return;
            }
            if (node.getBNode().getType().equalsIgnoreCase("CPT")){
                // Generate an EnumVariable.
                Enumerable domain = (Enumerable) node.getBNode().getVariable().getDomain();
                table.setValue(myParamsTable.getSelectedRow(), 
                        getEnumVar(value, domain));
            }
            table.display();
            
            fireTableCellUpdated(row, col);
        }
    }

    private EnumDistrib getEnumVar(Object input, Enumerable domain){
        // EnumVariable v1 = Predef.Boolean();
        // cpt1.put(new Object[]{true, false}, new EnumDistrib(v1.getDomain(), new double[]{1, 0}));
        String[] strArr = ((String) input).split(",");
        double[] dblArr = new double[strArr.length];
        for (int i = 0; i < strArr.length; i++) {
            dblArr[i] = Double.parseDouble(strArr[i]);
        }
        return new EnumDistrib(domain, dblArr);
    }
    
    private double[] parseInput(Object input){
        String[] strArr = ((String) input).split(",");
        double[] dblArr = new double[strArr.length];
        for (int i = 0; i < strArr.length; i++) {
            dblArr[i] = Double.parseDouble(strArr[i]);
        }
        return dblArr;
    }
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        myScrollPane = new javax.swing.JScrollPane();
        myParamsTable = new javax.swing.JTable();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        myParamsTable.setModel(new javax.swing.table.DefaultTableModel(
                new Object[][]{
                    {null, null, null, null},
                    {null, null, null, null},
                    {null, null, null, null},
                    {null, null, null, null}
                },
                new String[]{
                    "Title 1", "Title 2", "Title 3", "Title 4"
                }
        ));
        myScrollPane.setViewportView(myParamsTable);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(myScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addComponent(myScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 300, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(NodeParamsDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(NodeParamsDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(NodeParamsDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(NodeParamsDialog.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the dialog */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                NodeParamsDialog dialog = new NodeParamsDialog(new javax.swing.JFrame(), true);
                dialog.addWindowListener(new java.awt.event.WindowAdapter() {
                    @Override
                    public void windowClosing(java.awt.event.WindowEvent e) {
                        System.exit(0);
                    }
                });
                dialog.setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTable myParamsTable;
    private javax.swing.JScrollPane myScrollPane;
    // End of variables declaration//GEN-END:variables
}