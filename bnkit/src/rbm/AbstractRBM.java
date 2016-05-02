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
    
    public abstract Object[] encode(Object[] input);
    
    public abstract Object[] decode(Object[] hidden);
    
}
