/*
    bnkit -- software for building and using probabilistic models
    including Bayesian networks.
    Copyright (C) 2014-2016  M. Boden et al.

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
package rbm;

import bn.prob.EnumDistrib;
import dat.EnumVariable;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

/**
 *
 * @author mikael
 */
public abstract class AbstractRBM {
    
    /** Random number seed; set to change generation of weights. */
    public static int RANDOM_SEED = 1;
    protected Random rand;
    
    protected EnumVariable[] h; // hidden nodes (all Boolean)
    protected EnumVariable[] v; // visible nodes (must be enumerable)
    protected EnumDistrib[] Ph; // distributions at hidden nodes (all Boolean)
    protected EnumDistrib[] Pv; // distributions at visible nodes (must be enumerable)
    protected boolean[][] lnk;  // flag to indicate if nodes are linked
    
    public boolean isLinked(int idxVis, int idxHid) {
        return lnk[idxHid][idxVis];
    }
    
    public int getNVisible() {
        return v.length;
    }
    
    public EnumVariable[] getVisibleVars() {
        return v;
    }

    public EnumDistrib[] getVisibleDistribs() {
        return Pv;
    }
    
    public int getNHidden() {
        return h.length;
    }
    
    public EnumVariable[] getHiddenVars() {
        return h;
    }
    
    public EnumDistrib[] getHiddenDistribs() {
        return Ph;
    }

    private HashMap<EnumVariable, Integer> visvarmap = null;
    private HashMap<EnumVariable, Integer> hidvarmap = null;

    private static HashMap<EnumVariable, Integer> mapVarsToIndices(EnumVariable[] vars) {
        HashMap<EnumVariable, Integer> varmap = new HashMap<>();
        for (int i = 0; i < vars.length; i ++)
            varmap.put(vars[i], i);
        return varmap;
    }

    /**
     * Get the index of the specified input variable
     * @param var sought variable
     * @return the index of the variable, -1 if it is not found
     */
    public int getVisibleIndex(EnumVariable var) {
        if (v == null)
            return -1;
        if (visvarmap == null)
            visvarmap = mapVarsToIndices(v);
        Integer index = visvarmap.get(var);
        if (index != null)
            return index.intValue();
        return -1;
    }

    /**
     * Get the indices of the specified input variables
     * @param vars sought variables
     * @return array with the indices of the variables; value of -1 indicates missing
     */
    public int[] getVisibleIndex(EnumVariable[] vars) {
        if (visvarmap == null)
            visvarmap = mapVarsToIndices(v);
        int[] indices = new int[vars.length];
        for (int i = 0; i < indices.length; i ++)
            indices[i] = getVisibleIndex(vars[i]);
        return indices;
    }

    /**
     * Get the index of the specified hidden variable
     * @param var sought variable
     * @return the index of the variable, -1 if it is not found
     */
    public int getHiddenIndex(EnumVariable var) {
        if (h == null)
            return -1;
        if (hidvarmap == null)
            hidvarmap = mapVarsToIndices(h);
        Integer index = hidvarmap.get(var);
        if (index != null)
            return index.intValue();
        return -1;
    }

    /**
     * Get the indices of the specified hidden variables
     * @param vars sought variables
     * @return array with the indices of the variables; value of -1 indicates missing
     */
    public int[] getHiddenIndex(EnumVariable[] vars) {
        if (hidvarmap == null)
            hidvarmap = mapVarsToIndices(h);
        int[] indices = new int[vars.length];
        for (int i = 0; i < indices.length; i ++)
            indices[i] = getHiddenIndex(vars[i]);
        return indices;
    }

    public void setLinked(boolean state) {
        for (int j = 0; j < h.length; j ++) {
            for (int i = 0; i < v.length; i ++)
                this.lnk[j][i] = state;
        }
    }

    public void setLinked(int[] visidx, int[] hididx, boolean state) {
        for (int j = 0; j < hididx.length; j ++) {
            for (int i = 0; i < visidx.length; i ++) {
                this.lnk[hididx[j]][visidx[i]] = state;
            }
        }
    }

    public void setLinked(int from_visidx, int to_visidx, int hididx, boolean state) {
        for (int i = from_visidx; i < to_visidx; i ++) {
            this.lnk[hididx][i] = state;
        }
    }

    public void setLinkedWindow(int viswidth, boolean state) {
        int stepsize = (getNVisible() - viswidth) / (getNHidden() - 1);
        for (int j = 0; j < getNHidden(); j ++) {
            int start = j * stepsize;
            setLinked(start, start + viswidth, j, state);
        }
    }

    /**
     * Designate the links from the specified input variables as "in use".
     * @param inputvars
     */
    public void setLinked(EnumVariable[] inputvars) {
        int[] visidx = getVisibleIndex(inputvars);
        // assume all hidden nodes are linked
        for (int j = 0; j < h.length; j ++) {
            Arrays.fill(lnk[j], false);
            for (int i = 0; i < visidx.length; i ++) {
                if (visidx[i] != -1)
                    lnk[j][visidx[i]] = true;
            }
        }
    }

    /**
     * Designate the links from the specified input variables to the specified hidden variables as "in use".
     * @param inputvars input variables
     * @param hidvars hidden variables
     */
    public void setLinked(EnumVariable[] inputvars, EnumVariable[] hidvars) {
        int[] visidx = getVisibleIndex(inputvars);
        int[] hididx = getVisibleIndex(hidvars);
        for (int j = 0; j < hididx.length; j ++) {
            if (hididx[j] != -1) {
                for (int i = 0; i < visidx.length; i++) {
                    if (visidx[i] != -1)
                        lnk[hididx[j]][visidx[i]] = true;
                }
            }
        }
    }

    public EnumDistrib assignVisible(int index, Object value) {
        Object[] values = Pv[index].getDomain().getValues();
        double[] probs = new double[values.length];
        for (int i = 0; i < probs.length; i ++)
            probs[i] = (values[i].equals(value)) ? 1 : 0;
        Pv[index].set(probs);
        return Pv[index];
    }

    public static double logistic(double x) {
        double y = 1.0 / (1.0 + Math.exp(-x));
        return y;
    }

    public static double softmax(double[] x, int j) {
        double denom = 0;
        for (int i = 0; i < x.length; i ++)
            denom += Math.exp(x[i]);
        return Math.exp(x[j]) / denom;
    }

    public abstract Object[] encode(Object[] input, Object[] hidden);
    
    public abstract Object[] decode(Object[] hidden, Object[] visible);

    /**
     * Calculate a new hidden state vector from a given visible/input state vector;
     * note that this function updates the (internal) probabilities of the hidden state.
     * @param input input state (array)
     * @return hidden state (array); allocated internally so do what you want with it, but consider memory use...
     */
    public Object[] encode(Object[] input) {
        Object[] hinst = new Object[this.h.length];
        return encode(input, hinst);
    }

    /**
     * Calculate the visible state from the hidden;
     * note that this function updates the (internal) probabilities of the visible state.
     * @param hidden hidden state (array)
     * @return visible state (array); allocated internally so do what you want with it, but consider memory use...
     */
    public Object[] decode(Object[] hidden) {
        Object[] vinst = new Object[this.v.length];
        return decode(hidden, vinst);
    }


    public double err = -Double.MAX_VALUE;

    public abstract Double[][] getCDGradient(Object[][] minibatch, int niter);
    public abstract void setCDGradient(Double[][] delta);

    public Object[] encode_decode_clamped(Object[] clamp) {
        return encode_decode_clamped(clamp, 0);
    }

    public Object[] encode_decode_clamped(Object[] clamp, int niter) {
        Object[] hinst = new Object[getNHidden()];
        hinst = encode(clamp, hinst);
        Object[] vinst = new Object[getNVisible()];
        vinst = decode(hinst, vinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i++)
            decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(decoded, hinst);
            vinst = decode(hinst, vinst);
            for (int i = 0; i < vinst.length; i++)
                decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        }
        return decoded;
    }

    public Object[] encode_decode_restricted(Object[] input) {
        return encode_decode_restricted(input, 0);
    }

    public Object[] decode_restricted(Object[] hinst, Object[] input) {
        Object[] vinst = new Object[getNVisible()];
        vinst = decode(hinst, vinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (input[i] == null ? null : vinst[i]);
        return decoded;
    }

    /**
     * Encode and decode input states repeatedly, leaving inputs that are null un-instantiated throughout
     * @param input
     * @param niter number of repeats, 0 for performing the encode/decode once, 1 for repeating it, 2 for two repeats, etc.
     * @return
     */
    public Object[] encode_decode_restricted(Object[] input, int niter) {
        Object[] hinst = new Object[getNHidden()];
        hinst = encode(input, hinst);
        Object[] vinst = new Object[getNVisible()];
        vinst = decode(hinst, vinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (input[i] == null ? null : vinst[i]);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(decoded, hinst);
            vinst = decode(hinst, vinst);
            for (int i = 0; i < vinst.length; i++)
                decoded[i] = (input[i] == null ? null : vinst[i]);
        }
        return decoded;
    }

    public Object[] encode_decode_full(Object[] input) {
        return encode_decode_full(input, 0);
    }

    /**
     * Encode and decode input states repeatedly
     * @param input
     * @param niter number of repeats, 0 for performing the encode/decode once, 1 for repeating it, 2 for two repeats, etc.
     * @return
     */
    public Object[] encode_decode_full(Object[] input, int niter) {
        Object[] hinst = new Object[getNHidden()];
        hinst = encode(input, hinst);
        Object[] vinst = new Object[getNVisible()];
        vinst = decode(hinst, vinst);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(vinst, hinst);
            vinst = decode(hinst, vinst);
        }
        return vinst;
    }


}
