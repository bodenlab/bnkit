/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;

/**
 *
 * @author Jun All the controller logic lives here. Receives calls from View
 * classes MainJFrame and GraphPanel.
 */
public class BNController implements Observer {

    private final MainJFrame mainFrame;
    private final BNModel model;
    private final GraphPanel graphPanel;

    public BNController(MainJFrame mf, GraphPanel gp, BNModel mod) {
        this.mainFrame = mf;
        this.graphPanel = gp;
        this.model = mod;
    }

    /**
     * Initialises handler events for View class MainJFrame.
     */
    public void InitButtonHandlers() {
        mainFrame.initNodeControls();
        // THis msut be called after InitNodeButtons has been calle din MainJFrame
        
        // Only used if MainJFrame using Button node implementation
        if (mainFrame.usingButtons) {
            // Configure action handlers for add node buttons
            for (JButton button : mainFrame.getNodeButtons()) {
                button.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        JButton thisbtn = (JButton) e.getSource();
                        String type = thisbtn.getText();
                        graphPanel.createNode(null, type, null);
                    }
                });
            }
        } else {
            // Configure action handlers for drag-drop labels
            
        }
        mainFrame.getDeleteButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.deleteSelected();
            }
        });

        mainFrame.getDeleteAllButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.deleteAll();
            }
        });

        mainFrame.getRefreshButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.renderNetwork(model.getBNC());
                graphPanel.setLayout("");
            }
        });

        mainFrame.getMenuOpen().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.loadNetwork();
            }
        });

        mainFrame.getMenuSave().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.saveNetwork();
            }
        });

        mainFrame.getLayoutButton().addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                graphPanel.setLayout("");
            }
        });
    }

    @Override
    public void update() {
        // When update occurs, update View
        // and update the tables in BNode list.
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setSubject(Observable sub) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
