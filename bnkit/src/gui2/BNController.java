/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNode;
import bn.Predef;
import bn.Variable;
import com.mxgraph.model.mxIGraphModel;
import com.mxgraph.view.mxGraph;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Random;
import javax.swing.JButton;

/**
 *
 * @author Jun All the controller logic lives here. Receives calls from View
 * classes MainJFrame and GraphPanel.
 */
public class BNController {

    private MainJFrame mainFrame;
    private BNModel model;
    private GraphPanel graphPanel;

    public BNController(MainJFrame mf, GraphPanel gp, BNModel mod) {
        this.mainFrame = mf;
        this.graphPanel = gp;
        this.model = mod;
    }

    /**
     * Initialises handler events for View class MainJFrame.
     */
    public void InitButtonHandlers() {
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
            public void actionPerformed(ActionEvent e) {
                graphPanel.saveNetwork();
            }
        });
    }

}
