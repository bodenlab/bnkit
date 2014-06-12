/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui2;

/**
 *
 * @author Jun
 * Instantiate the Model, Views, Controllers.
 */
public class ProjectMain {
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
            java.util.logging.Logger.getLogger(MainJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(MainJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(MainJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(MainJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                
                // Initialise the views.
                MainJFrame mainFrame = new MainJFrame(false);
                mainFrame.setVisible(true);
                
                // Initialise the model.
                BNModel model = new BNModel();
                
                // Initialise the controllers.
                GraphPanel graphPanel = new GraphPanel(model);
                GraphPanelController gpController = new GraphPanelController(graphPanel, model);
                graphPanel.getGraphComponent().requestFocus();
                BNController bnController = new BNController(mainFrame, graphPanel, model, gpController);

                mainFrame.setGraphPanel(graphPanel);
                mainFrame.initDrawPanel();
                mainFrame.initNodeControls();
                gpController.initHandlers();
                bnController.InitButtonHandlers();
                
                
            }
        });
    }
}
