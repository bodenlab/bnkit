/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bn.ctmc;

import bn.*;
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
    private final SubstModel model;
    private final Object[] values;
    
    private double time = 1.0;
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

    @Override
    public String getName() {
        return var.getName();
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
        throw new UnsupportedOperationException("Not supported."); 
    }

    @Override
    public Variable getVariable() {
        return var;
    }

    @Override
    public List<EnumVariable> getParents() {
        return this.parentAsList;
    }

    @Override
    public EnumTable getTable() {
        throw new UnsupportedOperationException("Not supported."); 
    }

    @Override
    public Distrib getDistrib(Object[] key) {
        throw new UnsupportedOperationException("Not supported."); 
    }

    @Override
    public Distrib getDistrib() {
        throw new UnsupportedOperationException("Not supported."); 
    }

    @Override
    public void print() {
        System.out.print("Idx ");
        System.out.print("Parent ");
        for (Object value : this.values)
            System.out.print(String.format("[%6s]", value.toString()));
        System.out.println();
        for (int i = 0; i < this.values.length; i++) {
            System.out.print(String.format("%3d ", i));
            Object Y = this.values[i];
            System.out.print(String.format(" %-6s ", Y.toString()));
            for (Object X : this.values) {
                System.out.print(String.format(" %-5.3f ", get(X, Y)));
            }
            System.out.println();
        }
    }

    @Override
    public String getType() {
        return "SubstNode";
    }

    @Override
    public String getStateAsText() {
        StringBuilder sbuf = new StringBuilder("\n");
        sbuf.append(time);
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean setState(String dump) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isRoot() {
        return false; // cannot be root
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

    @Override
    public Factor makeFactor(BNet bn) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Factor makeFactor(Map<Variable, Object> relevant) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void tieTo(BNode source) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public BNode getTieSource() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    
}
