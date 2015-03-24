/*
    bnkit -- software for building and using Bayesian networks
    Copyright (C) 2014  M. Boden et al.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package bn.ctmc;

import bn.prob.EnumDistrib;
import bn.factor.Factor;
import dat.EnumVariable;
import dat.Variable;
import dat.EnumTable;
import bn.*;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.matrix.*;
import bn.factor.AbstractFactor;
import bn.factor.DenseFactor;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 *
 * @author mikael
 */
public class SubstNode implements BNode, TiedNode {
    
    private final EnumVariable var;
    private final EnumVariable parent;
    private final List<EnumVariable> parentAsList;
    private SubstModel model;
    private final Object[] values;
    
    private double time = 0.0;
    private Object instance = null;
    private boolean relevant = false; //for inference, track whether the node is relevant to the query

    /**
     * Create a conditional "substitution node" for an enumerable variable. 
     * The variable is conditioned on a set of Enumerable variables.
     *
     * @param var variable
     * @param parent parent variable
     * @param model the substitution model to use to generate conditional probabilities
     * @param t time for substitution (typically greater time, means greater chance for substitution)
     */
    public SubstNode(EnumVariable var, EnumVariable parent, SubstModel model, double t) {
        this.var = var;
        this.values = var.getDomain().getValues();
        this.parent = parent;
        this.parentAsList = Collections.singletonList(parent);
        this.time = t;
        this.model = model;
    }

    /**
     * Create a root "substitution node", ie a node modeling the ancestor at which the tree starts.
     * The substitution model is only used to provide the zero-order statistics.
     * @param var the variable
     * @param model substitution model but only the priors are used
     */
    public SubstNode(EnumVariable var, SubstModel model) {
        this.var = var;
        this.values = var.getDomain().getValues();
        this.model = model;
        this.parent = null; // this is a root node
        this.parentAsList = null;
    }
    
    @Override
    public String getName() {
        return getVariable().toString();
    }

    public String toString() {
        if (parent == null) {
            return "SN(" + var.getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            return "SN(" + var.getName() + "|" + parent.getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        }
    }

    /**
     * Determine the conditional probability P(X=value|Y=key).
     * @param key the condition Y, one value only
     * @param value the value assigned to the variable X represented by this node
     * @return the probability
     */
    @Override
    public Double get(Object[] key, Object value) {
        if (key != null) {
            if (key.length == 1) {
                Object Y = key[0];
                Object X = value;
                return model.getProb(X, Y, time);
            }
        }
        throw new RuntimeException("Invalid key: " + key);
    }

    @Override
    public Double get(Object value, Object... key) {
        return get(key, value);
    }

    @Override
    public Double get(Object value) {
        return model.getProb(value);
    }

    @Override
    public Variable getVariable() {
        return var;
    }

    @Override
    public List<EnumVariable> getParents() {
        return this.parentAsList;
    }

    public double getTime() {
        return time;
    }
    
    public SubstModel getModel() {
        return model;
    }
    
    @Override
    public EnumTable getTable() {
        throw new UnsupportedOperationException("Not supported."); 
    }

    @Override
    public Distrib getDistrib(Object[] key) {
        return model.getDistrib(key[0], time);
    }

    @Override
    public Distrib getDistrib() {
        return new EnumDistrib(var.getDomain(), model.F); 
    }

    @Override
    public void print() {
        if (parent == null) { // root
            for (Object value : this.values)
                System.out.print(String.format("[%5s]", value.toString()));
            System.out.println();
            for (Object X : this.values) {
                System.out.print(String.format(" %-5.3f ", get(X)));
            }
            System.out.println();
        } else {
            System.out.print("Idx ");
            System.out.print("Parent ");
            for (Object value : this.values)
                System.out.print(String.format("[%5s]", value.toString()));
            System.out.println();
            for (int i = 0; i < this.values.length; i++) {
                System.out.print(String.format("%3d ", i));
                Object Y = this.values[i];
                System.out.print(String.format("%-5s  ", Y.toString()));
                for (Object X : this.values) {
                    System.out.print(String.format(" %-5.3f ", get(X, Y)));
                }
                System.out.println();
            }
        }
    }

    @Override
    public String getType() {
        return "SubstNode";
    }

    @Override
    public String getStateAsText() {
        StringBuilder sbuf = new StringBuilder("\n");
        sbuf.append(time + "; (time)\n");
        sbuf.append(model.getName() + "; (model)\n");
        return sbuf.toString();
    }

    @Override
    public boolean setState(String dump) {
        int row = 1;
        boolean error = false;
        for (String line : dump.split("\n")) {
            String[] specline = line.split(";");
            switch (row) {
                case 2:
                    try {
                        this.time = Double.parseDouble(specline[0]);
                    } catch (NumberFormatException e) {
                        System.err.println("Time incorrect in state for SubstNode: " + this + " ==> " + specline[0]);
                        error = true;
                    }
                    break;
                case 3:
                    String modelName = specline[0].trim();
                    if (modelName.equalsIgnoreCase("JTT"))
                        this.model = new JTT();
                    else if (modelName.equalsIgnoreCase("WAG"))
                        this.model = new WAG();
                    else if (modelName.equalsIgnoreCase("LG"))
                        this.model = new LG();
                    else if (modelName.equalsIgnoreCase("Dayhoff"))
                        this.model = new Dayhoff();
                    else {
                        System.err.println("Model incorrect in state for SubstNode: " + this + " ==> " + specline[0]);
                        error = true;
                    }
                    break;
            }
            row ++;
        }
        return !error;
    }

    @Override
    public boolean isRoot() {
        return parent == null; 
    }

    @Override
    public void setInstance(Object value) {
        this.instance = value;
    }

    @Override
    public void resetInstance() {
        this.instance = null;
    }

    @Override
    public Object getInstance() {
        return instance;
    }

    @Override
    public Distrib makeDistrib(Collection<Sample> samples) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Make a Factor out of this node. If a variable is instantiated it will
     * be factored out.
     * If a parent is not relevant, it will not be included in the factor
     *
     * @param relevant only include relevant nodes, with instantiations if available
     * @return factor of node considering if parents are relevant
     */
    @Override
    public Factor makeFactor(Map<Variable, Object> relevant) {
        EnumVariable myvar = var;
        EnumVariable mypar = parent; // parent variable, if any
        // get value of this node, if any assigned
        Object varinstance = relevant.get(myvar);
        if (mypar != null) { // there is a parent variable
            // three possibilities:
            // 1. neither node variable or parent variable are instantiated, in which case, the parent is either
            //      a. relevant to the query, or
            //      b. irrelevant and needs to be marginalised subsequently
            // 2. node variable is instantiated, not parent, in which case, the parent is either
            //      a. relevant to the query, or
            //      b. irrelevant and needs to be marginalised subsequently
            // 3. parent variable is instantiated, not node variable

            // Check if parent variable is relevant
            boolean parent_is_relevant = relevant.containsKey(mypar);
            // get value of the parent node, if any assigned
            Object parinstance = relevant.get(mypar);
            
            // option 1
            if (varinstance == null && parinstance == null) { 
                List<Variable> vars_new = new ArrayList<>();
                vars_new.add(myvar);
                if (parent_is_relevant) {
                    vars_new.add(mypar);
                    Factor ft = new Factor(vars_new);
                    // F(X, Y) from P(X | Y)
                    for (Object X : values) {
                        for (Object Y : values) {
                            ft.setFactor(new Object[] {X, Y}, model.getProb(X, Y, time));
                        }
                    }
                    return ft;
                } else { // parent is NOT relevant
                    Factor ft = new Factor(vars_new);
                    // F(X) from P(X | Y)
                    for (Object X : values) {
                        double sum = 0;
                        for (Object Y : values) {
                            sum += model.getProb(X, Y, time);
                        }
                        ft.setFactor(new Object[] {X}, sum);
                    }
                    return ft;
                }
            } 
            // option 2
            else if (varinstance != null && parinstance == null) {
                if (parent_is_relevant) {
                    List<Variable> vars_new = new ArrayList<>();
                    vars_new.add(mypar);
                    Factor ft = new Factor(vars_new);
                    ft.evidenced = true;
                    // F(Y) from P(X = x| Y)
                    Object X = varinstance;
                    for (Object Y : values) {
                        ft.setFactor(new Object[] {Y}, model.getProb(X, Y, time));
                    }
                    return ft;
                } else { // parent is NOT relevant
                    Factor ft = new Factor();
                    ft.evidenced = true;
                    // F() from P(X = x| Y) and Y is irrelevant to query
                    Object X = varinstance;
                    double sum = 0;
                    for (Object Y : values) {
                        sum += model.getProb(X, Y, time);
                    }
                    ft.setFactor(sum);
                    return ft;
                }
            } 
            // option 3
            else if (varinstance == null && parinstance != null) {
                List<Variable> vars_new = new ArrayList<>();
                vars_new.add(myvar);
                Factor ft = new Factor(vars_new);
                ft.evidenced = true;
                // F(X) from P(X | Y)
                Object Y = parinstance;
                for (Object X : values) {
                    ft.setFactor(new Object[] {X}, model.getProb(X, Y, time));
                }
                return ft;
            }
            throw new RuntimeException("Invalid setting");
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                Factor ft = new Factor();
                ft.setFactor(model.getProb(varinstance));
                return ft;
            } else { // not instantiated
                List<Variable> vars_new = new ArrayList<>(1);
                vars_new.add(myvar);
                Factor ft = new Factor(vars_new);
                for (Object X : values) {
                    ft.setFactor(new Object[] {X}, model.getProb(X));
                }
                return ft;
            }
        }
    }

    @Override
    public AbstractFactor makeDenseFactor(Map<Variable, Object> relevant) {
        EnumVariable myvar = var;
        EnumVariable mypar = parent; // parent variable, if any
        // get value of this node, if any assigned
        Object varinstance = relevant.get(myvar);
        if (mypar != null) { // there is a parent variable
            // three possibilities:
            // 1. neither node variable or parent variable are instantiated, in which case, the parent is either
            //      a. relevant to the query, or
            //      b. irrelevant and needs to be marginalised subsequently
            // 2. node variable is instantiated, not parent, in which case, the parent is either
            //      a. relevant to the query, or
            //      b. irrelevant and needs to be marginalised subsequently
            // 3. parent variable is instantiated, not node variable

            // Check if parent variable is relevant
            boolean parent_is_relevant = relevant.containsKey(mypar);
            // get value of the parent node, if any assigned
            Object parinstance = relevant.get(mypar);
            
            // option 1
            if (varinstance == null && parinstance == null) { 
                if (parent_is_relevant) {
                    AbstractFactor ft = new DenseFactor(new Variable[] {mypar, myvar});
                    EnumVariable[] evars = ft.getEnumVars();
                    boolean cond_first = mypar.equals(evars[0]); // check if the condition is listed as first variable in factor
                    // F(X, Y) from P(X | Y)
                    for (Object X : values) {
                        for (Object Y : values) {
                            if (cond_first)
                                ft.setValue(new Object[] {Y, X}, model.getProb(X, Y, time));
                            else
                                ft.setValue(new Object[] {X, Y}, model.getProb(X, Y, time));
                        }
                    }
                    return ft;
                } else { // parent is NOT relevant
                    AbstractFactor ft = new DenseFactor(new Variable[] {myvar});
                    // F(X) from P(X | Y)
                    for (Object X : values) {
                        double sum = 0;
                        for (Object Y : values) {
                            sum += model.getProb(X, Y, time);
                        }
                        ft.setValue(new Object[] {X}, sum);
                    }
                    return ft;
                }
            } 
            // option 2
            else if (varinstance != null && parinstance == null) {
                if (parent_is_relevant) {
                    AbstractFactor ft = new DenseFactor(new Variable[] {mypar});
                    ft.evidenced = true;
                    // F(Y) from P(X = x| Y)
                    Object X = varinstance;
                    for (Object Y : values) {
                        ft.setValue(new Object[] {Y}, model.getProb(X, Y, time));
                    }
                    return ft;
                } else { // parent is NOT relevant
                    AbstractFactor ft = new DenseFactor();
                    ft.evidenced = true;
                    // F() from P(X = x| Y) and Y is irrelevant to query
                    Object X = varinstance;
                    double sum = 0;
                    for (Object Y : values) {
                        sum += model.getProb(X, Y, time);
                    }
                    ft.setValue(sum);
                    return ft;
                }
            } 
            // option 3
            else if (varinstance == null && parinstance != null) {
                AbstractFactor ft = new DenseFactor(new Variable[] {myvar});
                ft.evidenced = true;
                // F(X) from P(X | Y)
                Object Y = parinstance;
                for (Object X : values) {
                    ft.setValue(new Object[] {X}, model.getProb(X, Y, time));
                }
                return ft;
            }
            throw new RuntimeException("Invalid setting");
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                AbstractFactor ft = new DenseFactor();
                ft.setValue(model.getProb(varinstance));
                return ft;
            } else { // not instantiated
                AbstractFactor ft = new DenseFactor(new Variable[] {myvar});
                for (Object X : values) {
                    ft.setValue(new Object[] {X}, model.getProb(X));
                }
                return ft;
            }
        }
    }

    @Override
    public void countInstance(Object[] key, Object value, Double prob) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void countInstance(Object[] key, Object value) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void maximizeInstance() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isTrainable() {
        return false; // perhaps in the future "time" and aspects of model will be trainable
    }

    @Override
    public void randomize(long seed) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isRelevant() {
        return relevant;
    }

    @Override
    public void setRelevant(boolean relevant) {
        this.relevant = relevant;
    }

    @Override
    public void setTrainable(boolean trainable) {
        // no effect
    }

    @Override
    public void tieTo(BNode source) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public BNode getTieSource() {
        return null;
    }

    public static void main(String[] args) {
        // here are all the variables, representing the sequence at a position
        EnumVariable r1 = Predef.AminoAcid("R1");
        EnumVariable x1 = Predef.AminoAcid("X1");
        EnumVariable x2 = Predef.AminoAcid("X2");
        EnumVariable x3 = Predef.AminoAcid("X3");
        EnumVariable x4 = Predef.AminoAcid("X4");
        EnumVariable x5 = Predef.AminoAcid("X5");
        EnumVariable x6 = Predef.AminoAcid("X6");
        EnumVariable x7 = Predef.AminoAcid("X7");
        EnumVariable y1 = Predef.AminoAcid("Y1");
        EnumVariable y2 = Predef.AminoAcid("Y2");
        EnumVariable y3 = Predef.AminoAcid("Y3");
        EnumVariable y4 = Predef.AminoAcid("Y4");
        EnumVariable y5 = Predef.AminoAcid("Y5");
        EnumVariable y6 = Predef.AminoAcid("Y6");
        EnumVariable y7 = Predef.AminoAcid("Y7");
        // here's the phylogenetic tree represented by a BN
        SubstNode sn1r = new SubstNode(r1, new JTT());
        SubstNode sn1x = new SubstNode(x1, r1, new JTT(), 1.0);
        SubstNode sn2x = new SubstNode(x2, x1, new JTT(), 1.0);
        SubstNode sn3x = new SubstNode(x3, x2, new JTT(), 1.0);
        SubstNode sn4x = new SubstNode(x4, x2, new JTT(), 1.0);
        SubstNode sn5x = new SubstNode(x5, x1, new JTT(), 1.0);
        SubstNode sn6x = new SubstNode(x6, x5, new JTT(), 1.0);
        SubstNode sn7x = new SubstNode(x7, x5, new JTT(), 1.0);
        SubstNode sn1y = new SubstNode(y1, r1, new JTT(), 1.0);
        SubstNode sn2y = new SubstNode(y2, y1, new JTT(), 1.0);
        SubstNode sn3y = new SubstNode(y3, y2, new JTT(), 1.0);
        SubstNode sn4y = new SubstNode(y4, y2, new JTT(), 1.0);
        SubstNode sn5y = new SubstNode(y5, y1, new JTT(), 1.0);
        SubstNode sn6y = new SubstNode(y6, y5, new JTT(), 1.0);
        SubstNode sn7y = new SubstNode(y7, y5, new JTT(), 1.0);
        BNet bn = new BNet();
        bn.add(sn1x, sn2x, sn3x, sn4x, sn5x, sn6x, sn7x, sn1y, sn2y, sn3y, sn4y, sn5y, sn6y, sn7y, sn1r);
        // assign values to some variables, what we observe at the leaves of the tree
        sn3x.setInstance('Q');
        sn4x.setInstance('A');
        sn6x.setInstance('R');
        sn7x.setInstance('Q');
        sn3y.setInstance('L');
        sn4y.setInstance('K');
        sn6y.setInstance('K');
        sn7y.setInstance('R');
        
        //sn1.print();
        //sn2.print();
        //sn3.print();
        
        // Use variable elimination to perform inference in the tree
        VarElim ve = new VarElim();
        ve.instantiate(bn);
        
        // construct a query, this one to determine the "most probable explanation" (MPE) to the observations above
        // where each variable is assigned a value
//        Query q = ve.makeMPE(r1, x1, y1); // we only look at three
        Query q = ve.makeMPE(); // we only look at three
        CGTable qr = (CGTable) ve.infer(q);
        Variable.Assignment[] assigned = qr.getMPE();
        for (Variable.Assignment a : assigned) {
            System.out.println(a.var + " = " + a.val.toString());
        }
        
        // construct other queries, this time to determine the normal joint probability
        q = ve.makeQuery(r1);
        qr = (CGTable) ve.infer(q);
        qr.display();
        q = ve.makeQuery(x1);
        qr = (CGTable) ve.infer(q);
        qr.display();
        q = ve.makeQuery(y1);
        qr = (CGTable) ve.infer(q);
        qr.display();
    }

	@Override
	public List<Sample> getConditionDataset(int conditionIndex) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Distrib getlikelihoodDistrib() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void put(Object[] key, Distrib distr) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void put(Distrib prob) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void put(Distrib prob, Object... key) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void put(int index, Distrib distr) {
		// TODO Auto-generated method stub
		
	}
}
