/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNode;
import bn.CPT;
import bn.GDT;
import java.awt.Color;
import java.awt.event.ActionEvent;

/**
 *
 * @author Jun
 */
public class nodePropertiesPane extends javax.swing.JPanel {

    
    // TODO: change this to a JDialog
    private javax.swing.JLabel nodeDescriptionLabel;
    private javax.swing.JTextField nodeNameField;
    private javax.swing.JTextField nodeParametersField;
    private javax.swing.JLabel lbl1;
    private javax.swing.JLabel lbl2;
    private javax.swing.JButton applyBtn;
    private javax.swing.JButton cancelBtn;
    private BNode bnode;

    public nodePropertiesPane(BNode node) {
        this();
        bnode = node;
    }
    
    public nodePropertiesPane() {
        super();
        init();
        this.setForeground(Color.red); //testing

    }

    private void init() {

        // Set up the 'node actions' panel
        // Set up the 'node properties' panel
        nodeDescriptionLabel = new javax.swing.JLabel();
        nodeNameField = new javax.swing.JTextField();
        nodeParametersField = new javax.swing.JTextField();
        lbl1 = new javax.swing.JLabel("Node name");
        lbl2 = new javax.swing.JLabel("Parameters");
        applyBtn = new javax.swing.JButton("Apply");
        cancelBtn = new javax.swing.JButton("Cancel");

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
        this.setBorder(javax.swing.BorderFactory.createTitledBorder("Node Properties"));

        javax.swing.GroupLayout propertiesPanelLayout = new javax.swing.GroupLayout(this);
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
                        .addGroup(propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(applyBtn)
                                .addComponent(cancelBtn)))
        // Pack the two frames into the parent

        );
    }

    private void applyBtnPressed(java.awt.event.MouseEvent evt) {
        // update the node
    }

    private void cancelBtnPressed(java.awt.event.MouseEvent evt) {
        nodeNameField.setText("");
        nodeParametersField.setText("");
        this.setVisible(false);
        // dispose when this is JDialog
    }

    // Update display based on BNode contents.
    public void updateDisplay() {
        if (bnode != null) {
            // try casting to GDT and CPT and retrieve values.
            try {
                try {
                    nodeNameField.setText(((CPT) bnode).getName());
                    nodeDescriptionLabel.setText(((CPT) bnode).getType());
                } finally {
                    nodeNameField.setText(((GDT) bnode).getName());
                }
            } catch (ClassCastException ex) {
                
            }
        } else {
            // 19/05/14 -- This will always display because setBNode has been removed.
            // Later when this is JDialog, BNode will be passed to constructor.
            System.err.println("in updateDisplay(): bnode is null");
        }
    }


    public void setCellName() {

    }
}
