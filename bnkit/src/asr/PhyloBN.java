/*
 * bnkit -- software for building and using Bayesian networks
 * Copyright (C) 2014  M. Boden et al.
 * 
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
package asr;

import bn.BNet;
import bn.BNode;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.ctmc.SubstModel;
import bn.ctmc.SubstNode;
import bn.prob.EnumDistrib;
import dat.EnumVariable;
import dat.Variable;
import dat.phylo.*;

import java.util.*;

/**
 * Class for a Bayesian network that represents branch points,
 * and their relationships as extracted from a phylogenetic tree.
 * Based on deprecated class PhyloBNet.
 * @author mikael
 */
public class PhyloBN {

    private final BNet bn;
    private final SubstModel model;
    private double rate = 1.0;
    BNode[] bnodes = null;

    private PhyloBN(SubstModel model) {
        bn = new BNet();
        this.model = model;
    }

    public BNet getBN() {
        return bn;
    }

    /**
     * "Joint" reconstruction view of Bayesian network.
     */
    public class Joint implements TreeDecor<Object> {

        final private Object[] values;

        public Joint() {
             values = new Object[bnodes.length];
        }

        @Override
        public Object getDecoration(int idx) {
            return values[idx];
        }

        @Override
        public void decorate(TreeInstance ti) {
            Map<String, Integer> quick = new HashMap<>();
            // instantiate all nodes for which there are values, i.e. leaf nodes most probably
            for (int i = 0; i < bnodes.length; i ++) {
                quick.put(bnodes[i].getVariable().getName(), i);
                Object y = ti.getInstance(i);
                if (y != null)
                    bnodes[i].setInstance(y);
                else
                    bnodes[i].resetInstance();
            }
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(bn);
            Query q_mpe = ve.makeMPE();
            CGTable r1 = (CGTable)ve.infer(q_mpe);
            Variable.Assignment[] assign = r1.getMPE();
            for (Variable.Assignment assign1 : assign) {
                Integer idx = quick.get(assign1.var.getName());
                if (idx != null)
                    values[idx] = assign1.val;
            }
        }
    }

    /**
     * "Marginal" reconstruction view of Bayesian network.
     */
    public class Marginal implements TreeDecor<EnumDistrib> {

        final private EnumDistrib[] values;
        private int ancidx = 0;

        public Marginal(String ancestorID) {
            values = new EnumDistrib[bnodes.length];
            for (int i = 0; i < bnodes.length; i++) {
                if (bnodes[i].getVariable().getName().equals(ancestorID)) {
                    ancidx = i;
                    break;
                }
            }
        }

        public Marginal(int ancidx) {
            values = new EnumDistrib[bnodes.length];
            this.ancidx = ancidx;
        }

        @Override
        public EnumDistrib getDecoration(int idx) {
            return values[idx];
        }

        @Override
        public void decorate(TreeInstance ti) {
            // instantiate all nodes for which there are values, i.e. leaf nodes most probably
            for (int i = 0; i < bnodes.length; i ++) {
                Object y = ti.getInstance(i);
                if (y != null)
                    bnodes[i].setInstance(y);
                else
                    bnodes[i].resetInstance();
            }
            // set-up the inference engine
            VarElim ve = new VarElim();
            ve.instantiate(bn);
            Query q = ve.makeQuery(bnodes[ancidx].getVariable());
            CGTable r1 = (CGTable)ve.infer(q);
            values[ancidx] = (EnumDistrib)r1.getJDF().getDistrib(bnodes[ancidx].getVariable());
        }
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     * @param tree indexed phylogenetic tree
     * @param model evolutionary model
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBN create(IdxTree tree, SubstModel model) {
        return create(tree, model, 1.0);
    }

    /**
     * Construct a BN for specified phylogenetic tree using a single, supplied model.
     * @param tree phylogenetic tree
     * @param model evolutionary model
     * @param rate the evolutionary rate to be applied
     * @return the phylogenetic Bayesian network
     */
    public static PhyloBN create(IdxTree tree, SubstModel model, double rate) {
        PhyloBN pbn = new PhyloBN(model);
        pbn.rate = rate;
        EnumVariable[] bvars = new EnumVariable[tree.getSize()];
        pbn.bnodes = new BNode[tree.getSize()];
        // first, create variables...
        for (int idx : tree) {
            BranchPoint bp = tree.getBranchPoint(idx);
            EnumVariable rvar = new EnumVariable(model.getDomain(), bp.getID().toString());
            bvars[idx] = rvar;
        }
        // second, create nodes...
        for (int idx : tree) {
            BranchPoint bp = tree.getBranchPoint(idx);
            EnumVariable rvar = bvars[idx];
            int parent = tree.getParent(idx);
            if (parent < 0) { // there's no parent
                pbn.bnodes[idx] = new SubstNode(rvar, model);
            } else { // there's a parent
                EnumVariable prvar = bvars[parent];
                pbn.bnodes[idx] = new SubstNode(rvar, prvar, model, bp.getDistance());
            }
        }
        pbn.bn.add(pbn.bnodes);
        return pbn;
    }

}
