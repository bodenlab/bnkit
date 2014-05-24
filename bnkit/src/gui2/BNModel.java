/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package gui2;

import bn.BNode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Study Room
 */
public class BNModel {
    private BNContainer bnc;
    private HashMap<String, BNode> BNodeMap = new HashMap<>();
    private List<Object> selectedCells = new ArrayList<>();
    
    public BNContainer getBNC(){
        return bnc;
    }

    public BNModel(){
        this.bnc = new BNContainer();
    }
    /**
     * 
     * @return BNodeMap
     */
    public HashMap<String, BNode> getBNodeMap() {
        return BNodeMap;
    }

    public void setBNodeMap(HashMap<String, BNode> BNodeMap) {
        this.BNodeMap = BNodeMap;
    }

    public List<Object> getSelectedCells() {
        return selectedCells;
    }

    
    public void setSelectedCells(List<Object> selectedCells) {
        this.selectedCells = selectedCells;
    }
    
    public void setSelectedCells(Object selectedCells) {
        this.selectedCells.clear();
        this.selectedCells.add(selectedCells);
    }
 
    
}

