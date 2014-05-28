/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gui2;

import bn.BNet;
import bn.BNode;
import bn.CPT;
import bn.Continuous;
import bn.Distrib;
import bn.EnumTable;
import bn.EnumVariable;
import bn.Factor;
import bn.FactorTable;
import bn.GDT;
import bn.Variable;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Jun NodeModel implements Observable in the Observer pattern. Changes
 * in EnumVariable var's members are observed.
 */
//http://stackoverflow.com/questions/1658702/how-do-i-make-a-class-extend-observable-when-it-has-extended-another-class-too
public class NodeModel implements Observable, BNode {

    @Override
    public Factor makeFactor(BNet bn) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void countInstance(Object[] key, Object value) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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

//    /**
//     * Copy constructor for NodeModel
//     * @param oldNode 
//     */
//    public NodeModel (NodeModel oldNode){
//        this(oldNode.getVariable(), oldNode.getParents());
//    }
    
    public BNode getBNode() {
        return bnode;
    }

    public NodeModel(BNode bnode) {
        try {
            EnumVariable evar = (EnumVariable) bnode.getVariable();
            this.bnode = new CPT(evar, bnode.getParents());
            // change the 'String' constructor to a method.
        } catch(ClassCastException e) {
            Variable var = (Variable<Continuous>) bnode.getVariable();
            this.bnode = new GDT(var, bnode.getParents());
        }
        this.observers = new ArrayList<>();
    }

    public NodeModel(EnumVariable var) {
        this("CPT");
        this.bnode = new CPT(var);
    }

    public NodeModel(EnumVariable var, List<EnumVariable> parents) {
        this("CPT");
        this.bnode = new CPT(var, parents);
    }

    public NodeModel(EnumVariable var, EnumVariable... parents) {
        this("CPT");
        this.bnode = new CPT(var, parents);
    }

    public NodeModel(Variable<Continuous> var) {
        this("GDT");
        this.bnode = new GDT(var);
    }

    public NodeModel(Variable var, List<EnumVariable> parents) {
        try {
            var = (EnumVariable) var;
            this.bnode = new CPT((EnumVariable) var, bnode.getParents());
            // change the 'String' constructor to a method.
        } catch(ClassCastException e) {
            var = (Variable<Continuous>) var;
            this.bnode = new GDT(var, bnode.getParents());
        }
        this.observers = new ArrayList<>();
    }

    public NodeModel(Variable<Continuous> var, EnumVariable... parents) {
        this("GDT");
        this.bnode = new GDT(var, parents);
    }

    public String getCellName(){
        return cellName;
    }
    
    public void setCellName(String name){
        cellName = name;
    }
    
    public void setModel(String state){
        // TODO: throw exception if string doesn't match?
        this.state = inferenceModel.valueOf(state.toUpperCase());
    }
    
    public String getModel(){
        return state.name();
    }
    
    public ArrayList<String> getModelNames(){
        ArrayList<String> names = new ArrayList<>();
        for (inferenceModel s: inferenceModel.values()){
            names.add(s.name());
        }
        return names;
    }
    
    @Override
    public void register(Observer obj) {
        if (obj == null) {
            throw new NullPointerException("Null Observer");
        }
        synchronized (MUTEX) {
            if (!observers.contains(obj)) {
                observers.add(obj);
            }
        }
    }

    @Override
    public void unregister(Observer obj) {
        synchronized (MUTEX) {
            observers.remove(obj);
        }
    }

    @Override
    public void notifyObservers() {
        List<Observer> observersLocal = null;
        // This needs to be thread-safe!
        synchronized (MUTEX) {
            if (!changed) {
                return;
            }
            observersLocal = new ArrayList<>(this.observers);
            this.changed = false;
        }
        for (Observer obj : observersLocal) {
            obj.update();
        }
    }

    public void setName(String name) {
        this.bnode.getVariable().setName(name);
        changed = true;
        notifyObservers();
    }

    public void setParams(String name) {
        this.bnode.getVariable();
        changed = true;
        notifyObservers();
    }

    @Override
    public Object getUpdate(Observer obj) {
        return this.getVariable();
    }

    @Override
    public String getName() {
        return bnode.getName();
    }

    @Override
    public Double get(Object[] key, Object value) {
        return bnode.get(key, value);
    }

    @Override
    public Double get(Object value, Object... key) {
        return bnode.get(value, key);
    }

    @Override
    public Double get(Object value) {
        return bnode.get(value);
    }

    @Override
    public Variable getVariable() {
        // need to do this more intelligently...
        // should watch for a change in stored BNode
        changed = true;
        notifyObservers();
        return bnode.getVariable();
    }

    @Override
    public List<EnumVariable> getParents() {
        return bnode.getParents();
    }

    @Override
    public EnumTable getTable() {
        changed = true;
        notifyObservers();
        return bnode.getTable();
    }

    @Override
    public Distrib getDistrib() {
        return bnode.getDistrib();
    }

    @Override
    public void print() {
        bnode.print();
    }

    @Override
    public String getType() {
        return bnode.getType();
    }

    @Override
    public String getStateAsText() {
        return bnode.getStateAsText();
    }

    @Override
    public boolean setState(String dump) {
        changed = true;
        notifyObservers();
        return bnode.setState(dump);
    }

    @Override
    public boolean isRoot() {
        return bnode.isRoot();
    }

    @Override
    public void setInstance(Object value) {
        changed = true;
        notifyObservers();
        bnode.setInstance(value);
    }

    @Override
    public void resetInstance() {
        bnode.resetInstance();
    }

    @Override
    public Object getInstance() {
        return bnode.getInstance();
    }

//    @Override
//    public FactorTable makeFactor(BNet bn) {
//        return bnode.makeFactor(bn);
//    }

    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        bnode.countInstance(key, value, prob);
    }

    @Override
    public void maximizeInstance() {
        bnode.maximizeInstance();
    }

    @Override
    public boolean isTrainable() {
        return bnode.isTrainable();
    }

    @Override
    public void randomize(long seed) {
        bnode.randomize(seed);
    }

}