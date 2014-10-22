/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNet;
import dat.Enumerable;
import bn.Predef;
import dat.Variable;
import com.mxgraph.model.mxCell;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import javax.swing.AbstractButton;
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;

/**
 *
 * @author jun Dialog presented when node is right-clicked on.
 *
 * Allows nodes to be parameterised and evidenced. Model and View are updated
 * accordingly.
 */
public class NodePropertiesDialog extends javax.swing.JDialog {

    private javax.swing.JLabel nodeDescriptionLabel;
    private javax.swing.JTextField nodeNameField;
    private javax.swing.JTextField nodeParametersField;
    private javax.swing.JLabel lbl1;
    private javax.swing.JLabel lbl2;
    private javax.swing.JButton applyBtn;
    private javax.swing.JButton cancelBtn;
    private javax.swing.JPanel radioPanel;
    private javax.swing.JTextField optField;
    private ButtonGroup checkButtonGroup;
    private final GraphPanel graphPanel;
    private final NodeModel nodeModel;

    private final BNContainer bnc;
    private JRadioButton selectedRadioButton;

    public NodePropertiesDialog(java.awt.Frame parent, boolean modal, NodeModel nm, GraphPanel gp) {
        super(parent, modal);
        graphPanel = gp;
        bnc = graphPanel.getBNContainer();
        this.setTitle("Node Properties");

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
//        nodeNameField.setEditable(false); // Temporary
        nodeParametersField.setText(nodeModel.getVariable().getParams());
        nodeDescriptionLabel.setText(getTypeDescription(nodeModel.getVariable().getPredef()));
    }

    /**
     * Returns a short description about the queried predef. This needs to be
     * maintained manually.
     *
     * @param predef
     * @return predef-description
     */
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

    /**
     * Initialises GUI elements.
     */
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

        radioPanel = new javax.swing.JPanel();
        radioPanel.setVisible(true);

        // Generate radioButtons for the parameters.
        JRadioButton dummybtn;
        ArrayList<JRadioButton> buttonArr = new ArrayList<>();
        checkButtonGroup = new ButtonGroup();
        javax.swing.JLabel evidenceLbl = new javax.swing.JLabel("Evidenced value");
        ArrayList<String> paramsList = new ArrayList<>();

        if (nodeModel == null) {
            System.out.println("NodePropertiesDialog constructed with null NodeModel");
        }
        
        String params = nodeModel.getVariable().getParams();

        // Populate the parameter checkboxes.
        if (Predef.isEnumerable(nodeModel.getVariable().getPredef())) {
            paramsList.add("None");
            if (nodeModel.getVariable().getPredef().equalsIgnoreCase("String")) {
                paramsList.addAll(Arrays.asList(params.split(";")));
            } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Boolean")) {
                paramsList.addAll(Arrays.asList("True", "False"));
            } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Number")) {
                for (int i = 1; i < Integer.parseInt(params); i++) {
                    paramsList.add(String.valueOf(i));
                }
            } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Amino acid")) {
                for (int i = 0; i < Predef.AminoAcid().size(); i++) {
                    paramsList.add(String.valueOf(
                            Predef.AminoAcid().getDomain().get(i)
                    ));
                }
            } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Nucleic acid")) {
                for (int i = 0; i < Predef.NucleicAcid().size(); i++) {
                    paramsList.add(String.valueOf(
                            Predef.NucleicAcid().getDomain().get(i)));
                }
            } else {
                paramsList.add("In development");
            }

            // Generate JRadioButtons for parameters.
            for (String s : paramsList) {
                dummybtn = new JRadioButton(s);
                dummybtn.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        AbstractButton thisbtn = (AbstractButton) e.getSource();
                        selectedRadioButton = (JRadioButton) thisbtn;
                        System.out.println("selected radiobutton: "
                                + ((JRadioButton) selectedRadioButton).getText());
                    }
                });
                buttonArr.add(dummybtn);
                checkButtonGroup.add(dummybtn);
                radioPanel.add(dummybtn);
            }
        } else { //GDT
            optField = new javax.swing.JTextField();
            optField.setColumns(10);
            radioPanel.add(optField);
            if (nodeModel.getInstance() != null) {
                optField.setText(String.valueOf(nodeModel.getInstance()));
            }
        }

        // Select checkboxes if node is evidenced.
        if (Predef.isEnumerable(nodeModel.getVariable().getPredef())) {
            setParameterBoxes();
        }

        // Listener for 'Apply' button click
        applyBtn.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                applyBtnPressed(evt);
            }
        });

        // Listener for 'Cancel' button click
        cancelBtn.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                cancelBtnPressed(evt);
            }
        });

        this.setBackground(new java.awt.Color(204, 204, 204));

        javax.swing.GroupLayout propertiesPanelLayout = new javax.swing.GroupLayout(getContentPane());

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
                .addComponent(evidenceLbl)
                .addComponent(radioPanel)
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
                        .addComponent(evidenceLbl)
                        .addComponent(radioPanel)
                        .addGroup(propertiesPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(applyBtn)
                                .addComponent(cancelBtn)))
        );
        pack();
    }

    /**
     * Handler for 'Apply' button press. Performs updates to Model and View.
     *
     * @param evt
     */
    private void applyBtnPressed(java.awt.event.MouseEvent evt) {

        // Apply name change
        Variable oldVar = Predef.getVariable(nodeModel.getName(),
                nodeModel.getVariable().getPredef(), nodeModel.getVariable().getParams());
        NodeModel nmCopy = new NodeModel(oldVar);
        
        String newNodeName = nodeNameField.getText();
        if (!oldVar.getName().equals(newNodeName)) {
            nodeModel.getVariable().setName(newNodeName);

            // Iterate through all nodes and check if their parentvars match the old var
            for (NodeModel nm : bnc.getNodeModelArr().values()) {
                if (!(nm.getParents() == null) ) {
                    for (Variable pvars : nm.getParents()) {
                        if (pvars.getName().equals(oldVar.getName())) { // this is a child
                            bnc.removeParent(nm, oldVar);
                            bnc.addParent(nm, nodeModel.getVariable());
                        }
                    }
                } 
            }
            
            //Rename cells in View.
            for (Object cell : graphPanel.getGraph().getChildVertices(graphPanel.getGraph().getDefaultParent())) {
                if (((mxCell) cell).getValue().equals(oldVar.getName())) {
                    graphPanel.getGraph().getModel().setValue(
                            cell, newNodeName);
                }
            }

        }

        // Update done, remove old node from bnc, add new node.
        bnc.removeNode(nmCopy);
        bnc.addNode(nodeModel);
        
        // Modify oldvar to track change in params
        oldVar = Predef.getVariable(nodeModel.getName(),
                nodeModel.getVariable().getPredef(), nodeModel.getVariable().getParams());
        
        // Update node variables
        nodeModel.getVariable().setParams(nodeParametersField.getText());
        // now update parent params of child nodes.
        // Update parentVars of children of nodeModel
//            for (NodeModel nm : bnc.getNodeModelArr().values()) {
//                if (!(nm.getParents() == null) ) {
//                    for (Variable pvars : nm.getParents()) {
//                        if (pvars.getName().equals(oldVar.getName()) &&
//                                !(pvars.getParams().equals(oldVar.getParams()))) {
//                            System.err.println(nodeModel.getVariable().getParams());
//                            System.err.println(pvars.getParams());
//                            bnc.addParent(nm, nodeModel.getVariable());
//                            bnc.removeParent(nm, oldVar);
//                            System.err.println(pvars.getParams());
//                        }
//                    }
//                } else{
//                } 
//            }
//        
         

        // Update state
        if (Predef.isEnumerable(nodeModel.getVariable().getPredef())) {
            if (checkParameterBox().equals("null")) {
                nodeModel.setInferenceModel("IGNORE");
                nodeModel.setInstance(null);
            } else {
                nodeModel.setInferenceModel("Evidence");
                nodeModel.setInstance(checkParameterBox());
            }
        } else {
            if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Real")) {
                if (optField.getText().isEmpty()) {
                    nodeModel.setInferenceModel("IGNORE");
                    nodeModel.setInstance(null);
                }
                nodeModel.setInferenceModel("Evidence");
                try {
                    nodeModel.setInstance(Double.parseDouble(optField.getText()));
                } catch (NumberFormatException e) {

                }
            }
        }
        
 
        
        dispose();
    }

    /**
     *
     */
    private void setParameterBoxes() {
        String instanceString = "";
        java.util.Enumeration<AbstractButton> radioButtons = checkButtonGroup.getElements();

        if (nodeModel.getInstance() == null) {
            while (radioButtons.hasMoreElements()) {
                AbstractButton b = radioButtons.nextElement();
                if (((JRadioButton) b).getText().equalsIgnoreCase("None")) {
                    checkButtonGroup.setSelected(b.getModel(), true);
                    selectedRadioButton = (JRadioButton) b;
                }
            }
            return;
        }

        instanceString = nodeModel.getInstance().toString();

        // TODO: wrap casts in try catches
        if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Boolean")) {
            if ((Boolean) nodeModel.getInstance()) {
                instanceString = "True";
            } else {
                instanceString = "False";
            }
        } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Number")) {
            instanceString = String.valueOf((Integer) nodeModel.getInstance());
        } else {
            instanceString = (String) nodeModel.getInstance();
        }
        while (radioButtons.hasMoreElements()) {
            AbstractButton b = radioButtons.nextElement();
            if (((JRadioButton) b).getText().equalsIgnoreCase(instanceString)) {
                checkButtonGroup.setSelected(b.getModel(), true);
                selectedRadioButton = (JRadioButton) b;
            }
        }

    }

    private Object checkParameterBox() {
        String checkTxt;
        if (Predef.isEnumerable(nodeModel.getVariable().getPredef())) {
            if (selectedRadioButton.getText().isEmpty()) {
                return "null";
            } else {
                checkTxt = selectedRadioButton.getText();
            }
            if (nodeModel.getVariable().getPredef().equalsIgnoreCase("String")) {
                return checkTxt;
            } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Boolean")) {
                if (checkTxt.equalsIgnoreCase("true")) {
                    return true;
                } else if (checkTxt.equalsIgnoreCase("false")) {
                    return false;
                } else {
                    return "null";
                }
            } else if (nodeModel.getVariable().getPredef().equalsIgnoreCase("Number")) {
                try {
                    return Integer.parseInt(checkTxt);
                } catch (NumberFormatException e) {
                    return "null";
                }
            }
        } else {
            System.err.println("Calling checkParameterBox() for a continuous variable.");
        }
        return "null";
    }

    /**
     * Handler for 'Cancel' button press. Discards changes and closes dialog.
     *
     * @param evt
     */
    private void cancelBtnPressed(java.awt.event.MouseEvent evt) {
        dispose();
    }

    public void setCellName() {

    }
}
