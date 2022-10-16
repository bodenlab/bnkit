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

import bn.factor.Factorize;
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
import dat.Enumerable;

import java.util.*;

/**
 * This class makes a branch point in a phylogenetic tree into a node in a Bayesian network.
 *
 * It has been rewritten so that an instance does not make explicit reference to the substitution model as
 * accesses for probs are dependent on time, and since BN nodes use distinct time distances, it is better to
 * keep probs here (and leave the model separate).
 *
 * @author mikael
 */
public class SubstNode implements BNode, TiedNode {
    
    private final EnumVariable var;
    private final EnumVariable parent;
    private final List<EnumVariable> parentAsList;
    protected final EnumTable<EnumDistrib> table; // table of (enumerable) probability distributions
    private final double[][] probs;
    private final EnumDistrib prior; // one (enumerable) probability distribution that is used if this variable is NOT conditioned

    //private SubstModel model;
    private String modelname = null;
    private final Enumerable alpha;
    private final Object[] values;
    private boolean gap;

    private double time;
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
        this.parent = parent;
        this.parentAsList = Collections.singletonList(parent);
        this.time = t;
        //this.model = model;
        this.modelname = model.getName();
        this.alpha = model.getDomain();
        this.values = this.alpha.getValues();
        this.table = new EnumTable<>(parent);
        this.probs = model.getProbs(t);
        for (int i = 0; i < this.values.length; i ++) {
            table.setValue(i, new EnumDistrib(this.alpha, probs[i]));
        }
        this.prior = null; // this is a cond prob
    }

    /**
     * Create a root "substitution node", ie a node modeling the ancestor at which the tree starts.
     * The substitution model is only used to provide the zero-order statistics.
     * @param var the variable
     * @param model substitution model but only the priors are used
     */
    public SubstNode(EnumVariable var, SubstModel model) {
        this.var = var;
        //this.model = model;
        this.alpha = model.getDomain();
        this.values = this.alpha.getValues();
        this.table = null;
        this.probs = null;
        this.prior = new EnumDistrib(var.getDomain(), model.F);
        this.parent = null; // this is a root node
        this.parentAsList = null;
        this.time = 0;
    }
    
    @Override
    public String getName() {
        return getVariable().toString();
    }

    public boolean getGap() {
        return gap;
    }
    public void setGap(boolean gap) {
        this.gap = gap;
    }

    public String toString() {
        if (parent == null) {
            return "SN(" + var.getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        } else {
            return "SN(" + var.getName() + "|" + parent.getName() + ")" + (getInstance() == null ? "" : "=" + getInstance());
        }
    }

    /**
     * Get conditional probability P(X=x|Y=y,time)
     * @param X
     * @param Y
     * @return
     */
    public double getProb(Object X, Object Y) {
        int index_X = alpha.getIndex(X);
        int index_Y = alpha.getIndex(Y);
        return probs[index_Y][index_X];
        //return this.table.getValue(index_X).get(index_Y);
    }

    /**
     * Get conditional probability P(X=x|Y=y,time) by index to value
     * @param index_X
     * @param index_Y
     * @return
     */
    private double getProb(int index_X, int index_Y) {
        return probs[index_Y][index_X];
    }

    /**
     * Get probability P(X=x)
     * @param X
     * @return
     */
    public double getProb(Object X) {
        int index_X = alpha.getIndex(X);
        return this.prior.get(index_X);
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
                //return table.getValue(key).get(value);
                Object Y = key[0];
                Object X = value;
                return getProb(X, Y);
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
        return this.prior.get(value);
        //return getProb(value);
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

    /**
     *
     * @deprecated pog2020 -- remove all refs to model to make thread safe
     */
    public SubstModel getModel() {
        return SubstModel.createModel(modelname);
    }
    
    @Override
    public EnumTable getTable() {
        throw new UnsupportedOperationException("Not supported."); 
    }

    /**
     * Retrieve the distribution for this node that applies GIVEN the parents' instantiations.
     * Requires all parent nodes to be instantiated.
     * @param key the parent values
     * @return the distribution of the variable for this node
     */
    @Override
    public Distrib getDistrib(Object[] key) {
        if (this.table == null || key == null)
            return this.getDistrib();
        try {
            return this.table.getValue(key);
        } catch (RuntimeException e) {
            throw new RuntimeException("Evaluation of SubstNode " + this.toString() + " failed since condition was not fully specified: " + e.getMessage());
        }
    }

    @Override
    public EnumDistrib getDistrib() {
        return prior;
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
        sbuf.append(modelname + "; (model)\n");
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
                    modelname = specline[0].trim();
                    SubstModel model = SubstModel.createModel(modelname);
                    if (model == null) {
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
                            ft.setFactor(new Object[] {X, Y}, getProb(X, Y));
                        }
                    }
                    return ft;
                } else { // parent is NOT relevant
                    Factor ft = new Factor(vars_new);
                    // F(X) from P(X | Y)
                    for (Object X : values) {
                        double sum = 0;
                        for (Object Y : values) {
                            sum += getProb(X, Y);
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
                        ft.setFactor(new Object[] {Y}, getProb(X, Y));
                    }
                    return ft;
                } else { // parent is NOT relevant
                    Factor ft = new Factor();
                    ft.evidenced = true;
                    // F() from P(X = x| Y) and Y is irrelevant to query
                    Object X = varinstance;
                    double sum = 0;
                    for (Object Y : values) {
                        sum += getProb(X, Y);
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
                    ft.setFactor(new Object[] {X}, getProb(X, Y));
                }
                return ft;
            }
            throw new RuntimeException("Invalid setting");
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                Factor ft = new Factor();
                ft.setFactor(getProb(varinstance));
                return ft;
            } else { // not instantiated
                List<Variable> vars_new = new ArrayList<>(1);
                vars_new.add(myvar);
                Factor ft = new Factor(vars_new);
                for (Object X : values) {
                    ft.setFactor(new Object[] {X}, getProb(X));
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
            // Consider what defines this factor:
            // - parinstance value (null or value from mypar.getDomain()) AND the domain (but NOT the variable name, so use domain hash-code)
            // - varinstance value (null or value from myvar.getDomain()) AND the domain (but NOT the variable name, so use domain hash-code)
            // - the probabilities that are used to set factor
            // option 1
            if (varinstance == null && parinstance == null) { 
                if (parent_is_relevant) {
                    AbstractFactor ft = new DenseFactor(new Variable[] {mypar, myvar});
                    EnumVariable[] evars = ft.getEnumVars();
                    boolean cond_first = mypar.equals(evars[0]); // check if the condition is listed as first variable in factor
                    for (int index_X = 0; index_X < values.length; index_X ++) {
                        for (int index_Y = 0; index_Y < values.length; index_Y ++) {
                            if (cond_first)         // F(Y, X) from P(X | Y)
                                ft.setValue(new Object[] {values[index_Y], values[index_X]}, getProb(index_X, index_Y));
                            else                    // F(X, Y) from P(X | Y)
                                ft.setValue(new Object[] {values[index_X], values[index_Y]}, getProb(index_X, index_Y));
                        }
                    }
                    return ft;
                } else { // parent is NOT relevant
                    AbstractFactor ft = new DenseFactor(new Variable[] {myvar});
                    // F(X) from P(X | Y)
                    for (Object X : values) {
                        double sum = 0;
                        for (Object Y : values) {
                            sum += getProb(X, Y);
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
                        ft.setValue(new Object[] {Y}, getProb(X, Y));
                    }
                    return ft;
                } else { // parent is NOT relevant
                    AbstractFactor ft = new DenseFactor();
                    ft.evidenced = true;
                    // F() from P(X = x| Y) and Y is irrelevant to query
                    Object X = varinstance;
                    double sum = 0;
                    for (Object Y : values) {
                        sum += getProb(X, Y);
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
                    ft.setValue(new Object[] {X}, getProb(X, Y));
                }
                return ft;
            }
            throw new RuntimeException("Invalid setting");
        } else { // no parents, just a prior
            if (varinstance != null) { // instantiated prior
                AbstractFactor ft = new DenseFactor();
                ft.setValue(getProb(varinstance));
                return ft;
            } else { // not instantiated
                AbstractFactor ft = new DenseFactor(new Variable[] {myvar});
                for (Object X : values) {
                    ft.setValue(new Object[] {X}, getProb(X));
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

    public static void main0(String[] args) {
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

    public static void myPrintAADistrib(Distrib d) {
        Character[] S = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
        for (Character c : S)
            System.out.print(c + "    ");
        System.out.println();
        for (Character c : S)
            System.out.print(String.format("%4.2f", ((int)((d.get(c)*100))/100.0)) + " ");
        System.out.println();
    }
    
    public static void main1(String[] args) {
        // all the variables, representing the sequence at a position
        EnumVariable x1 = Predef.AminoAcid("X1");
        EnumVariable x2 = Predef.AminoAcid("X2");
        EnumVariable x3 = Predef.AminoAcid("X3");
        // next link variables to nodes forming a tree
        SubstNode snx1 = new SubstNode(x1, new WAG());
        SubstNode snx2 = new SubstNode(x2, x1, new WAG(), 1.0);
        SubstNode snx3 = new SubstNode(x3, x1, new WAG(), 1.0);
        BNet bn = new BNet();
        bn.add(snx1, snx2, snx3);
        // assign values to some variables, what we observe at the leaves of the tree
        snx2.setInstance('K');
        snx3.setInstance('Q');

        myPrintAADistrib(snx1.getDistrib());
        Character[] S = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
        double[] myd = new double[S.length];
        Arrays.fill(myd, 1.0);
        
        for (int i = 0; i < S.length; i ++)
            myd[i] *= snx1.getDistrib().get(S[i]);
        
        for (int i = 0; i < S.length; i ++) 
            System.out.print(S[i] + "    ");
        System.out.println();
        for (int i = 0; i < S.length; i ++) {
            Distrib d = snx2.getDistrib(new Object[] {S[i]});
            System.out.print(String.format("%4.2f", ((int)((d.get('K')*100))/100.0)) + " ");
            myd[i] *= d.get('K');
        }
        System.out.println();
                
        for (int i = 0; i < S.length; i ++) 
            System.out.print(S[i] + "    ");
        System.out.println();
        for (int i = 0; i < S.length; i ++) {
            Distrib d = snx3.getDistrib(new Object[] {S[i]});
            System.out.print(String.format("%4.2f", ((int)((d.get('Q')*100))/100.0)) + " ");
            myd[i] *= d.get('Q');
        }
        System.out.println();
        
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
        q = ve.makeQuery(x1);
        qr = (CGTable) ve.infer(q);
        myPrintAADistrib(qr.query(x1));
        
        
        System.out.println();
        EnumDistrib mydistr = new EnumDistrib(Enumerable.aacid_alt, myd);
        myPrintAADistrib(mydistr);
        for (int i = 0; i < S.length; i ++) {
            System.out.println(S[i] + "\t" + myd[i]);
        }
    }

    /**
     * Example in tutorial
     * @param args
     */
    public static void main(String[] args) {
        // all the variables, representing the sequence at a position
        // each variable can be instantiated to one of twenty amino acids
        // we use a predefined alphabet
        EnumVariable n0 = Predef.AminoAcid("N0");
        EnumVariable n1 = Predef.AminoAcid("N1");
        EnumVariable n2 = Predef.AminoAcid("N2");
        EnumVariable a18 = Predef.AminoAcid("A18");
        EnumVariable a57 = Predef.AminoAcid("A57");
        EnumVariable a25 = Predef.AminoAcid("A25");
        EnumVariable a59 = Predef.AminoAcid("A59");
        // next link variables to nodes forming a tree;
        // each node will represent a conditional probability table;
        // the structure they form collectively is specific to a position,
        // below is index 6 for the example in the tutorial
        SubstNode sn0 = new SubstNode(n0, new LG());
        SubstNode sa18 = new SubstNode(a18, n0, new LG(), 0.07);
        SubstNode sn1 = new SubstNode(n1, n0, new LG(), 0.01);
        SubstNode sa57 = new SubstNode(a57, n1, new LG(), 0.07);
        // Option 0: be naive and loop through all possible instantiations of ancestral variables n0 and n1
        Object[] aa = Predef.AminoAcid().getDomain().getValues();
        EnumTable<Double> table = new EnumTable<>(n0, n1); // table for storing all products
        double sum = 0; // keep total so we can normalise later
        for (int i = 0; i < aa.length; i ++) {
            for (int j = 0; j < aa.length; j ++) {
                double p = // here's the product
                        sn0.get(aa[i]) *             // P(N0)
                        sn1.get(aa[j], aa[i]) *      // P(N1 | N0)
                        sa18.get('E', aa[i]) * // P(a18=E | N0)
                        sa57.get('D', aa[j]);  // P(a57=D | N1)
                table.setValue(new Object[] {aa[i], aa[j]}, p); // we save the product in a table keyed by n0 and n1
                sum += p;
            }
        }
        // table.display(); // prints out all entries we've saved, NOT normalised
        // below we normalise the table using the sum above
        for (int i = 0; i < aa.length; i ++) {
            for (int j = 0; j < aa.length; j ++) {
                int index = table.getIndex(new Object[] {aa[i], aa[j]});
                // we can access the entry in the table by a separate index, defined by the n0|n1 key
                table.setValue(index, table.getValue(index) / sum);
            }
        }
        // table.display(); // prints out all entries again but now normalised
        // below we access only a sub-set of the entres of specific interest
        Object[] select = new Object[] {'D','E'};
        for (int i = 0; i < select.length; i ++) {
            for (int j = 0; j < select.length; j ++)
                System.out.println("P(N0=" + select[i] + ", N1=" + select[j] + " | a18=E, a57=D) = " + table.getValue(new Object[] {select[i], select[j]}));
        }
        // Option 1: use variable elimination as implemented in bnkit, so make a BN and perform inference much much more efficiently
        BNet bn = new BNet();
        bn.add(sn0, sa18, sn1, sa57);
        // assign values to some variables, what we observe at the leaves of the tree
        sa18.setInstance('E');
        sa57.setInstance('D');
        // instantiate inference engine on BN
        VarElim ve = new VarElim();
        ve.instantiate(bn);
        // construct a query, this one to determine the "most probable explanation" (MPE) to the observations above
        // where each variable is assigned the value which jointly is most probable
        Query q = ve.makeMPE();  // the inference engine will figure out what variables that have been instantiated to form the condition
        CGTable qr = (CGTable) ve.infer(q);
        Variable.Assignment[] assigned = qr.getMPE(); // extract the answer; for MPE this means all un-instantiated variable assignments
        for (Variable.Assignment a : assigned) {
            System.out.println(a.var + " = " + a.val.toString());
        }

        // lets work on index 17; we re-use the variables above but we need a few more nodes
        SubstNode sn2 = new SubstNode(n2, n1, new LG(), 0.02);
        SubstNode sa25 = new SubstNode(a25, n2, new LG(), 0.01);
        SubstNode sa59 = new SubstNode(a59, n2, new LG(), 0.02);

        bn = new BNet();
        bn.add(sn0, sa18, sn1, sa57, sn2, sa25, sa59);
        // assign values to some variables, what we observe at the leaves of the tree
        sa18.setInstance('N');
        sa57.setInstance('D');
        sa25.setInstance('E');
        sa59.setInstance('E');

        // Option 0: 8,000 possibilities? yikes!
        // Also, we now need to nest three loops (with four variables, nest four loops, not great...)
        // Option half-way: make factor tables
        // To make a factor table we need to provide a map of instantiated variables
        Map<Variable,Object> evidence = new HashMap<>();
        for (BNode node : bn.getOrdered()) {
            Variable var = node.getVariable();
            Object val = node.getInstance(); // will be null if not instantiated
            evidence.put(var, val);
        }
        AbstractFactor f0 = sn0.makeDenseFactor(evidence);
        f0.display();
        AbstractFactor f1 = sn1.makeDenseFactor(evidence);
        f1.display();
        AbstractFactor f2 = sn2.makeDenseFactor(evidence);
        f2.display();
        AbstractFactor f18 = sa18.makeDenseFactor(evidence);
        f18.display();
        AbstractFactor f57 = sa57.makeDenseFactor(evidence);
        f57.display();
        AbstractFactor f25 = sa25.makeDenseFactor(evidence);
        f25.display();
        AbstractFactor f59 = sa59.makeDenseFactor(evidence);
        f59.display();
        // now just multiply them;
        // in theory it does not matter in what order; in practice it influences the time and space complexity hugely...
        // the factors can get very big
        AbstractFactor f0p1 = Factorize.getProduct(f0, f1);
        AbstractFactor f0p1p2 = Factorize.getProduct(f0p1, f2);
        AbstractFactor f0p1p2pa18 = Factorize.getProduct(f0p1p2, f18);
        AbstractFactor f0p1p2pa18pa57 = Factorize.getProduct(f0p1p2pa18, f57);
        AbstractFactor f0p1p2pa18pa57pa25 = Factorize.getProduct(f0p1p2pa18pa57, f25);
        AbstractFactor f0p1p2pa18pa57pa25pa59 = Factorize.getProduct(f0p1p2pa18pa57pa25, f59);
        // extract the entry with the greatest value, and normalise to make a probability of it
        double fsum = f0p1p2pa18pa57pa25pa59.getSum();
        double maxp = 0;
        int maxidx = -1;
        for (int idx : f0p1p2pa18pa57pa25pa59) {
            double p = f0p1p2pa18pa57pa25pa59.getValue(idx);
            if (p > maxp) {
                maxidx = idx;
                maxp = p;
            }
        }
        Object[] maxkey = f0p1p2pa18pa57pa25pa59.getKey(maxidx);
        EnumVariable[] evars = f0p1p2pa18pa57pa25pa59.getEnumVars();
        System.out.println("Key with greatest probability: " + evars[0] + "=" + maxkey[0] + ", " + evars[1] + "=" + maxkey[1] + ", " + evars[2] + "=" + maxkey[2] + " P=" + maxp);

        // Option 1: instantiate inference engine on BN
        ve = new VarElim();
        ve.instantiate(bn);
        q = ve.makeMPE();  // the inference engine will figure out what variables that have been instantiated to form the condition
        qr = (CGTable) ve.infer(q);
        assigned = qr.getMPE(); // extract the answer; for MPE this means all un-instantiated variable assignments
        for (Variable.Assignment a : assigned) {
            System.out.println(a.var + " = " + a.val.toString());
        }

    }


    
}
