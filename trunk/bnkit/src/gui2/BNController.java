/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import gui2.test.Observable;
import bn.BNet;
import bn.JPT;
import bn.alg.CGTable;
import bn.alg.Query;
import com.mxgraph.model.mxCell;
import com.mxgraph.view.mxGraph;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;

/**
 *
 * @author jun 
 * Controller logic for MainJFrame lives here in BNController. 
 * Receives calls from View
 * classes MainJFrame and GraphPanel.
 */
public class BNController{ //implements Observer 

    private final MainJFrame mainFrame;
    private final BNModel model;
    private final GraphPanel graphPanel;
    private final GraphPanelController graphPanelController;

    public BNController(MainJFrame mf, GraphPanel gp, BNModel mod, GraphPanelController gpc) {
        this.mainFrame = mf;
        this.graphPanel = gp;
        this.model = mod;
        this.graphPanelController = gpc;
    }

    /**
     * Initialises handler events for View class MainJFrame.
     */
    public void InitButtonHandlers() {

        // Only used if MainJFrame using Button node implementation
        if (mainFrame.usingButtons) {
            // Configure action handlers for add node buttons
            for (JButton button : mainFrame.getNodeButtons()) {
                button.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        JButton thisbtn = (JButton) e.getSource();
                        String type = thisbtn.getText();
                        graphPanelController.createNode(null, type, null);
                        graphPanelController.addNodetoBNC(null, type, null);
                    }
                });
            }
        } else {
            // Configure action handlers for drag-drop labels

        }
        mainFrame.getDeleteButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanelController.deleteSelected();
            }
        });

        mainFrame.getDeleteAllButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanelController.deleteAll();
            }
        });

        mainFrame.getRefreshButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.renderNetwork(model.getBNC());
                graphPanel.executeLayout(0);
                for (NodeModel node : model.getBNC().getNodeModelArr().values()) {
                    System.out.println(" >Node is: " + node.getName());
                }
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
            @Override
            public void actionPerformed(ActionEvent e) {
                graphPanel.executeLayout(0);
            }
        });

        mainFrame.getSetQueryButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                java.util.List<Object> cells = graphPanel.getSelectedCells();
                if (cells.size() == 1) {
                    mainFrame.setQueryLbl("<html>Query Node: " + ((mxCell) cells.get(0)).getValue());
                    
                    BNContainer bnc = graphPanel.getBNContainer();
                    NodeModel node = bnc.getNodeModel((String)((mxCell) cells.get(0)).getValue());
                    
                    if (node == null) {
                        System.err.println("Oops, done soemthing wrong");
                    } else {
                        node.setInstance(null);
                        node.setInferenceModel("QUERY");
                        graphPanel.setQueryNode(node);
                    }
                    
                } else {
                    mainFrame.setQueryLbl("<html>Query Node: error, select one query cell.");
                }
            }

        });

        mainFrame.getSetInferButton().addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
//                graphPanel.doInference();
                doInference();
//                graphPanel.testInfer();
            }
        });
    }

//    @Override
    public void update() {
        // When update occurs, update View
        // and update the tables in BNode list.
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

//    @Override
    public void setSubject(Observable sub) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    public void doInference(){
        if (graphPanel.getQueryNode() == null) {
            System.err.println("queryNode is null. Returning...");
            return;
        }
        
        BNet bn = graphPanel.getBNContainer().getBNet();
        //aa
        bn.alg.CGVarElim ve = new bn.alg.CGVarElim();
        ve.instantiate(bn);
        Query q = ve.makeQuery(graphPanel.getQueryNode().getVariable());
        CGTable res = (CGTable)ve.infer(q);
//        res.display();
        
//        System.out.println("Variable elimination--------------");
//		bn.alg.CGVarElim ve = new bn.alg.CGVarElim();
//		ve.instantiate(bn);
//		Query q = ve.makeQuery(B);
//		JPT jpt = ve.infer(q).getJPT();
//		jpt.display();
        if (graphPanel.getQueryNode().getType().equalsIgnoreCase("CPT")){
            JPT jpt = ve.infer(q).getJPT();
            jpt.display();
            mainFrame.getResultLabel().setText("<html>" + jpt.toString());
        } else {
            mainFrame.getResultLabel().setText("<html>" + res.toString());
        }
        
        
        //set mainFrame's panel??
        
        
        mainFrame.validate();
    }
}
