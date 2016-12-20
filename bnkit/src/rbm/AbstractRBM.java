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

    public int getVisibleIndex(EnumVariable var) {
        if (v == null)
            return -1;
        if (visvarmap == null) {
            visvarmap = new HashMap<>();
            for (int i = 0; i < v.length; i ++)
                visvarmap.put(var, i);
        }
        Integer index = visvarmap.get(var);
        if (index != null)
            return index.intValue();
        return -1;
    }

    public int getHiddenIndex(EnumVariable var) {
        if (h == null)
            return -1;
        if (hidvarmap == null) {
            hidvarmap = new HashMap<>();
            for (int i = 0; i < h.length; i ++)
                hidvarmap.put(var, i);
        }
        Integer index = hidvarmap.get(var);
        if (index != null)
            return index.intValue();
        return -1;
    }

    public void setLinked(EnumVariable[] inputvars) {
        // assume all hidden nodes are linked
        int[] visidx = new int[inputvars.length];
        for (int i = 0; i < inputvars.length; i ++) {
            int idx = getVisibleIndex(inputvars[i]);
            if (idx != -1)
                visidx[i] = idx;
        }
        for (int j = 0; j < h.length; j ++) {
            Arrays.fill(lnk[j], false);
            for (int i = 0; i < visidx.length; i ++)
                lnk[j][visidx[i]] = true;
        }
    }

    public void setLinked(EnumVariable[] inputvars, EnumVariable[] hidvars) {
        throw new RuntimeException("Not yet implemented");
    }

    public EnumDistrib assignVisible(int index, Object value) {
        Object[] values = Pv[index].getDomain().getValues();
        double[] probs = new double[values.length];
        for (int i = 0; i < probs.length; i ++)
            probs[i] = (values[i].equals(value)) ? 1 : 0;
        Pv[index].set(probs);
        return Pv[index];
    }

    public abstract Object[] encode(Object[] input);
    
    public abstract Object[] decode(Object[] hidden);

    public double err = -Double.MAX_VALUE;

    public abstract Double[][] getCDGradient(Object[][] minibatch, int niter);
    public abstract void setCDGradient(Double[][] delta);

    public Object[] encode_decode_clamped(Object[] clamp) {
        return encode_decode_clamped(clamp, 0);
    }

    public Object[] encode_decode_clamped(Object[] clamp, int niter) {
        Object[] hinst = encode(clamp);
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i++)
            decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(decoded);
            vinst = decode(hinst);
            for (int i = 0; i < vinst.length; i++)
                decoded[i] = (clamp[i] == null ? vinst[i] : clamp[i]);
        }
        return decoded;
    }

    public Object[] encode_decode_restricted(Object[] input) {
        return encode_decode_restricted(input, 0);
    }

    public Object[] decode_restricted(Object[] hinst, Object[] input) {
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (input[i] == null ? null : vinst[i]);
        return decoded;
    }

    public Object[] encode_decode_restricted(Object[] input, int niter) {
        Object[] hinst = encode(input);
        Object[] vinst = decode(hinst);
        Object[] decoded = new Object[this.v.length];
        for (int i = 0; i < vinst.length; i ++)
            decoded[i] = (input[i] == null ? null : vinst[i]);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(decoded);
            vinst = decode(hinst);
            for (int i = 0; i < vinst.length; i++)
                decoded[i] = (input[i] == null ? null : vinst[i]);
        }
        return decoded;
    }

    public Object[] encode_decode_full(Object[] input) {
        return encode_decode_full(input, 0);
    }

    public Object[] encode_decode_full(Object[] input, int niter) {
        Object[] hinst = encode(input);
        Object[] vinst = decode(hinst);
        for (int n = 0; n < niter; n ++) {
            hinst = encode(vinst);
            vinst = decode(hinst);
        }
        return vinst;
    }


}
