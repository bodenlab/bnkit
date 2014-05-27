/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNet;
import bn.Predef;
import bn.Variable;
import com.mxgraph.model.mxCell;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.AbstractButton;
import javax.swing.ButtonGroup;
import javax.swing.GroupLayout.ParallelGroup;
import javax.swing.GroupLayout.SequentialGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JRootPane;
import javax.swing.SwingUtilities;

/**
 *
 * @author Jun
 */
public class NodePropertiesDialog extends javax.swing.JDialog {

    private javax.swing.JLabel nodeDescriptionLabel;
    private javax.swing.JTextField nodeNameField;
    private javax.swing.JTextField nodeParametersField;
    private javax.swing.JLabel lbl1;
    private javax.swing.JLabel lbl2;
    private javax.swing.JButton applyBtn;
    private javax.swing.JButton cancelBtn;
    private GraphPanel graphPanel;
    private NodeModel nodeModel;
    
    private BNContainer bnc;
    private Object selectedRadioButton;

    public NodePropertiesDialog(java.awt.Frame parent, boolean modal, NodeModel nm, GraphPanel gp) {
        super(parent, modal);
        graphPanel = gp;
        bnc = graphPanel.getBNContainer();
        this.setTitle("Node Properties");
        
        // Set default button
//        JRootPane rootpane = SwingUtilities.getRootPane(applyBtn);
//        rootPane.setDefaultButton(applyBtn);
        
        nodeModel = nm;
        init();
        System.out.println("NodePropertiesDialog Constructor");
        if (nodeModel == null) {
            System.err.println("constructor of NodePropertiesDialog, bnode is null now, good going.");
        }
        System.out.println(nodeModel + " name: " + nodeModel.getName() + ", predef: "
                + nodeModel.getVariable().getPredef() + ", params: "
                + nodeModel.getVariable().getParams());

        nodeNameField.setText(nodeModel.getName());
        nodeParametersField.setText(nodeModel.getVariable().getParams());
        nodeDescriptionLabel.setText(getTypeDescription(nodeModel.getVariable().getPredef()));
    }

    private String getTypeDescription(String predef) {
        String header = "<html>" + predef + " (" + nodeModel.getType() + ")<br>";
        if (Predef.isParameterised(predef)) {
            if (predef.equalsIgnoreCase("String")) {
                return header + "Parameters takes a list of strings, delimited by semicolon.<br>"
                        + "e.g. A;B;C";
            } else if (predef.equalsIgnoreCase("Number")) {
                return header + "Parameters takes a single integer value which "
                        + "represents the maximum value.";
            }
        } else {
            nodeParametersField.setEnabled(false);
            return header + "No parameters.";
        }

        return null;
    }

    private void init() {

        // Set up the 'node actions' panel
        // Set up the 'node properties' panel
        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        nodeDescriptionLabel = new javax.swing.JLabel();
        nodeNameField = new javax.swing.JTextField();
        nodeParametersField = new javax.swing.JTextField();
        lbl1 = new javax.swing.JLabel("Node name");
        lbl2 = new javax.swing.JLabel("Parameters");
        applyBtn = new javax.swing.JButton("Apply");
        cancelBtn = new javax.swing.JButton("Cancel");

        // Loop through model.getModelNames() to generate JRadioButtons...
        JRadioButton dummybtn;
        ArrayList<JRadioButton> buttonArr = new ArrayList<>();
        ButtonGroup bgroup = new ButtonGroup();
        
        for (String s: nodeModel.getModelNames()){
            dummybtn = new JRadioButton(s);
            dummybtn.addActionListener(new ActionListener(){
                public void actionPerformed(ActionEvent e){
                    AbstractButton thisbtn = (AbstractButton) e.getSource();
                    selectedRadioButton = thisbtn;
                    System.out.println("selected radiobutton: " + 
                            ((JRadioButton)selectedRadioButton).getText());
                }
            });
            buttonArr.add(dummybtn);
            bgroup.add(dummybtn);
        }
        
        applyBtn.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                applyBtnPressed(evt);
            }
        });

        cancelBtn.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                cancelBtnPressed(evt);
            }
        });

        this.setBackground(new java.awt.Color(204, 204, 204));
//        this.setBorder(javax.swing.BorderFactory.createTitledBorder("Node Properties"));

        
        
        javax.swing.GroupLayout propertiesPanelLayout = new javax.swing.GroupLayout(getContentPane());
        
//        ParallelGroup hgroup =  propertiesPanelLayout.createParallelGroup();
//        SequentialGroup vgroup = propertiesPanelLayout.createSequentialGroup();
//        for (JRadioButton b: buttonArr){
//            hgroup.addComponent(b);
//            vgroup.addComponent(b);
//        }
        
        this.setLayout(propertiesPanelLayout);
        propertiesPanelLayout.setHorizontalGroup(
                propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(propertiesPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                .addComponent(nodeDescriptionLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(nodeNameField)
                                .addComponent(nodeParametersField)
                                .addGroup(propertiesPanelLayout.createSequentialGroup()
                                        .addGroup(propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                .addComponent(lbl1)
                                                .addComponent(lbl2))
                                        .addGap(0, 0, Short.MAX_VALUE)))
                        .addContainerGap())
//                .addGroup(hgroup)
                .addGroup(propertiesPanelLayout.createSequentialGroup()
                        .addGap(36, 36, 36)
                        .addComponent(applyBtn)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 48, Short.MAX_VALUE)
                        .addComponent(cancelBtn)
                        .addGap(34, 34, 34))
        );
        propertiesPanelLayout.setVerticalGroup(
                propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(propertiesPanelLayout.createSequentialGroup()
                        .addGap(6, 6, 6)
                        .addComponent(lbl1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nodeNameField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(lbl2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nodeParametersField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(nodeDescriptionLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 61, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 7, Short.MAX_VALUE)
//                        .addGroup(vgroup)
                        .addGroup(propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(applyBtn)
                                .addComponent(cancelBtn)))
        // Pack the two frames into the parent
        );
        pack();
    }

    private void applyBtnPressed(java.awt.event.MouseEvent evt) {
        // If old name is not the same as the new name, need to remove
        // the old <key,value> pair from the map and re-enter it.
        String oldNodeName = nodeModel.getVariable().getName();
        String newNodeName = nodeNameField.getText();
        if (!oldNodeName.equals(newNodeName)) {
            
            // Bayesian Network requires at least one parent-child relationship
            // check in BNet constructor failing
//            BNet bn = bnc.getBNetnm(); 
            
            
            // The node will have no children, because the 'child' will point to old name.
            
            // Make children point to new parent name
//            System.out.println("==bn is: " + Arrays.toString(bn.getNodes().toArray()));
//            if (bn.getChildren(nodeModel) != null) { // always null.... nodemodel has changed.
//                for (String nodename : bn.getChildren(nodeModel)) {
//                    System.out.println("!!Child is: " + nodename);
//                    bnc.getNodeModel(nodename).setName(newNodeName);
//                }
//            } else {
//                System.out.println(nodeModel.getName() + " has no children");
//            }
            
            
            // Alternate brute-force method for getting chjildren to point to new parent name
            for (NodeModel nm: bnc.getNodeModelArr().values()){
                System.out.println("!!nm is:" + nm.getName());
                if (nm.getParents() != null && !nm.getParents().isEmpty()){ // null object pattern...
                    for (Variable v: nm.getParents()){
                        System.out.println("parent node is: " + v.getName());
                        if (v.getName().equals(oldNodeName)){
                            v.setName(newNodeName);
                        }
                    }
                }
            }
            
            
            //TODO: use a placeholder node for search?
            NodeModel tempnode = new NodeModel(nodeModel);
            
            // Since get by ID is not exposed in Java implementation of JGraphX,
            // have to bruteforce search for matching cell...
            // TODO: replace this with listener, and have View handle the name update
            for (Object cell: graphPanel.getAllVertices()){
                if (((mxCell) cell).getValue().equals(oldNodeName)){
                    graphPanel.getGraph().getModel().setValue(
                        cell, newNodeName);
                }
            }
            
            bnc.removeNode(nodeModel);
            nodeModel.getVariable().setName(newNodeName);
            bnc.addNode(nodeModel);
            
            
            
            
        }
        // Update variables
        nodeModel.getVariable().setParams(nodeParametersField.getText());
//        nodeModel.setModel(((JRadioButton) selectedRadioButton).getText());
        dispose();
    }

    private void cancelBtnPressed(java.awt.event.MouseEvent evt) {
        dispose();
    }

    public void setCellName() {

    }
}
