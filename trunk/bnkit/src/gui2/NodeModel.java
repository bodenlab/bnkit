/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import gui2.test.Observer;
import bn.BNet;
import bn.BNode;
import bn.node.CPT;
import dat.Continuous;
import bn.Distrib;
import dat.EnumTable;
import dat.EnumVariable;
import bn.Factor;
import bn.node.GDT;
import bn.SampleTable.*;
import dat.Variable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 *
 * @author jun 
 * Wrapper around the BNode class.
 * NodeModel implements Observable in the Observer pattern. 
 */

public class NodeModel {

    public void countInstance(Object[] key, Object value) {
        bnode.countInstance(key, value);
    }

    public Distrib getDistrib(Object[] key) {
        return bnode.getDistrib(key);
    }

    private enum DistributionType {
        GDT, CPT
    };
    private DistributionType DistType;

    private enum inferenceModel {
        QUERY, EVIDENCE, IGNORE
    };
    private BNode bnode;
    private List<Observer> observers;
    private boolean changed;
    private final Object MUTEX = new Object();
    private String cellName = null;
    private inferenceModel state = inferenceModel.QUERY;

    private NodeModel(String type) {
        if (type.equalsIgnoreCase("GDT")) {
            DistType = DistributionType.GDT;
        } else if (type.equalsIgnoreCase("CPT")) {
            DistType = DistributionType.CPT;
        }
        this.observers = new ArrayList<>();
    }
    
    public BNode getBNode() {
        return bnode;
    }

    public NodeModel(BNode bnod) {
        this.bnode = bnod;
        this.observers = new ArrayList<>();
    }

    public NodeModel(Variable var) {
        try {
            var = (EnumVariable) var;
            this.bnode = new CPT((EnumVariable) var);
            // change the 'String' constructor to a method.
        } catch(ClassCastException e) {
            var = (Variable<Continuous>) var;
            this.bnode = new GDT(var);
        }
        this.observers = new ArrayList<>();
    }
    

    public NodeModel(Variable var, List<EnumVariable> parents) {
        try {
            var = (EnumVariable) var;
            this.bnode = new CPT((EnumVariable) var, parents);
            // change the 'String' constructor to a method.
        } catch(ClassCastException e) {
            var = (Variable<Continuous>) var;
            this.bnode = new GDT(var, parents);
        }
//        this("CPT");
        this.observers = new ArrayList<>();
    }
    
    public NodeModel(Variable var, EnumVariable... parents) {
        try {
            var = (EnumVariable) var;
            this.bnode = new CPT((EnumVariable) var, parents);
            // change the 'String' constructor to a method.
        } catch(ClassCastException e) {
            var = (Variable<Continuous>) var;
            this.bnode = new GDT(var, parents);
        }
//        this("CPT");
        this.observers = new ArrayList<>();
    }


    public String getCellName(){
        return cellName;
    }
    
    public void setCellName(String name){
        cellName = name;
    }
    
    /**
     * Set inference model of the node: 'Query', 'Evidence', 'Ignore'
     * @param state 
     */
    public void setInferenceModel(String state){
        this.state = inferenceModel.valueOf(state.toUpperCase());
    }
    
    public String getInferenceModel(){
        return state.name();
    }
    
    public ArrayList<String> getModelNames(){
        ArrayList<String> names = new ArrayList<>();
        for (inferenceModel s: inferenceModel.values()){
            names.add(s.name());
        }
        return names;
    }
  

    public void setName(String name) {
        this.bnode.getVariable().setName(name);
        changed = true;
//        notifyObservers();
    }

    public void setParams(String name) {
        this.bnode.getVariable();
        changed = true;
//        notifyObservers();
    }
//    @Override
//    public void register(Observer obj) {
//        if (obj == null) {
//            throw new NullPointerException("Null Observer");
//        }
//        synchronized (MUTEX) {
//            if (!observers.contains(obj)) {
//                observers.add(obj);
//            }
//        }
//    }
//
//    @Override
//    public void unregister(Observer obj) {
//        synchronized (MUTEX) {
//            observers.remove(obj);
//        }
//    }
//
//    @Override
//    public void notifyObservers() {
//        List<Observer> observersLocal = null;
//        // This needs to be thread-safe!
//        synchronized (MUTEX) {
//            if (!changed) {
//                return;
//            }
//            observersLocal = new ArrayList<>(this.observers);
//            this.changed = false;
//        }
//        for (Observer obj : observersLocal) {
//            obj.update();
//        }
//    }
//    
//    @Override
//    public Object getUpdate(Observer obj) {
//        return this.getVariable();
//    }

    public String getName() {
        return bnode.getName();
    }

    public Double get(Object[] key, Object value) {
        return bnode.get(key, value);
    }

    public Double get(Object value, Object... key) {
        return bnode.get(value, key);
    }

    public Double get(Object value) {
        return bnode.get(value);
    }

    public Variable getVariable() {
        changed = true;
        return bnode.getVariable();
    }

    public List<EnumVariable> getParents() {
        return bnode.getParents();
    }

    public EnumTable getTable() {
        changed = true;
//        notifyObservers();
        return bnode.getTable();
    }

    public Distrib getDistrib() {
        return bnode.getDistrib();
    }

    public void print() {
        bnode.print();
    }

    public String getType() {
        return bnode.getType();
    }

    public String getStateAsText() {
        return bnode.getStateAsText();
    }

    public boolean setState(String dump) {
        changed = true;
//        notifyObservers();
        boolean success = bnode.setState(dump);
        if (success){
            System.out.println("set state success");
        } else {
            System.out.println("set state failure");
        }
        return success;
    }
    
    public boolean isRoot() {
        return bnode.isRoot();
    }

    public void setInstance(Object value) {
        changed = true;
//        notifyObservers();
        bnode.setInstance(value);
    }

    public void resetInstance() {
        bnode.resetInstance();
    }

    public Object getInstance() {
        return bnode.getInstance();
    }

    public void countInstance(Object[] key, Object value, Double prob) {
        bnode.countInstance(key, value, prob);
    }

    public void maximizeInstance() {
        bnode.maximizeInstance();
    }

    public boolean isTrainable() {
        return bnode.isTrainable();
    }

    public void randomize(long seed) {
        bnode.randomize(seed);
    }
    
    // update parent variable
    public void update(NodeModel oldparent, NodeModel newparent){
        
    }
    
}
