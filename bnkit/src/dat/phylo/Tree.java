/*
 * binfkit -- software for bioinformatics
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
package dat.phylo;

import asr.Parsimony;
import bn.prob.GammaDistrib;
import dat.EnumSeq;
import dat.EnumSeq.Alignment;
import dat.Enumerable;
import dat.file.Newick;

import java.io.IOException;
import java.util.*;

/**
 * Class to represent a single phylogenetic tree, refactored from old PhyloTree (now deprecated).
 * Rooted, multifurcating tree for representing phylogenetic relationships.
 * Functionality includes labeling and traversing branchpoints; reading and writing to Newick format;
 * Programmers should note that almost all functionality is implemented through recursion.
 *
 * The current design separates the tree topology (with branch points and their labels, represented by this class)
 * from instantiations (values assigned to tips and internal branch points, represented by TreeInstance,
 * of which several can be based on the same topology).
 *
 * @author mikael
 */
public class Tree extends IdxTree {

    private final BranchPoint root; // the root of the tree, all descendants linked by pointer in BranchPoint class

    /**
     * Constructor for tree from a root with all nodes connected off that.
     * Use factory methods to construct trees by assembling BranchPoints relative a root.
     */
    public Tree(BranchPoint root) {
        super(straightenTree(root));
        this.root = root;        // assume all branch points are OK and store them
    }

    /**
     * Create the array of branch points from a single branch point, which in turn is linked to its children
     * @param root the ultimate branch point
     * @return an array of all branch points
     */
    private static BranchPoint[] straightenTree(BranchPoint root) {
        List<BranchPoint> branchPoints = root.getSubtree();
        BranchPoint[] bpoints = new BranchPoint[branchPoints.size()];
        branchPoints.toArray(bpoints);
        return bpoints;
    }

    /**
     * Load a tree file with given name and of given format
     * @param filename file
     * @param format format
     * @return instance of tree
     * @throws IOException
     */
    public static Tree load(String filename, String format) throws IOException {
        if (format.equalsIgnoreCase("newick") || format.equalsIgnoreCase("nwk"))
            return Newick.load(filename);
        else
            throw new IOException("Unknown format: " + format);
    }

    /**
     * Save the current tree to a file with given name and of given format
     * @param filename file
     * @param format format
     * @throws IOException if an error happens during the write operation
     */
    public void save(String filename, String format) throws IOException {
        if (format.equalsIgnoreCase("newick") || format.equalsIgnoreCase("nwk"))
            Newick.save(this, filename, Newick.MODE_DEFAULT);
        else if (format.equalsIgnoreCase("ancestor") || format.equalsIgnoreCase("anwk"))
            Newick.save(this, filename, Newick.MODE_ANCESTOR);
        else
            throw new IOException("Unknown format: " + format);
    }

    /**
     * String representation in the Newick format.
     *
     * @return string representation of tree
     */
    @Override
    public String toString() {
        return root.toString();
    }

    /**
     * Find the node with the specified label.
     *
     * @param content label or label
     * @return matching node, or null if not found
     */
    public BranchPoint find(Object content) {
        return root.find(content);
    }

    /**
     * Label the internal branch points
     * by the convention "N" then the number incrementing from 0 for root,
     * by depth-first search. Overwrites the labels at internal branch points.
     */
    public void setInternalLabels() {
        root.setInternalLabels(0);
    }

    /**
     * Create a random tree
     * @param nLeaves number of leaves
     * @param gamma_k parameter for gamma specifying distance between branchpoints ("shape" aka alpha, or "a")
     * @param gamma_lambda parameter for gamma specifying distance between branchpoints ("scale" or "b", also referred by its inverse 1/beta, where beta is "rate")
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree Random(int nLeaves, long seed, double gamma_k, double gamma_lambda, int max_desc, int min_desc) {
        String[] leaves = new String[nLeaves];
        for (int i = 0; i < leaves.length; i ++)
            leaves[i] = "A" + (i + 1);
        return Random(leaves, seed, gamma_k, gamma_lambda, max_desc, min_desc);
    }

    /**
     * Create a random tree
     * @param leafLabels labels to be assigned to leaves
     * @param gamma_k parameter for gamma specifying distance between branchpoints ("shape" aka alpha)
     * @param gamma_lambda parameter for gamma specifying distance between branchpoints ("scale", also referred to as 1/beta, where beta is "rate")
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */

    public static Tree Random(String[] leafLabels, long seed, double gamma_k, double gamma_lambda, int max_desc, int min_desc) {
        Random rand = new Random(seed);
        GammaDistrib gd = new GammaDistrib(gamma_k, gamma_lambda);
        gd.setSeed(seed);
        List<BranchPoint> nodes = new ArrayList<>(leafLabels.length);
        // the initial, complete set of nodes that need ancestors
        for (String label : leafLabels) {
            BranchPoint bp = new BranchPoint(label);
            bp.setDistance(Math.max(gd.sample(), 0.0000001));
            nodes.add(bp);
        }
        // create ancestor nodes by picking two children, connecting them via a new ancestor,
        // then updating the list of nodes yet to be allocated an ancestor
        int M = nodes.size();
        BranchPoint bp = M > 0 ? nodes.get(0) : null;
        while (M > 1) {
            bp = new BranchPoint("");
            bp.setDistance(Math.max(gd.sample(), 0.0000001));
            // how many nodes to pick?
            int N = Math.min(rand.nextInt(max_desc - min_desc + 1) + min_desc, nodes.size());
            for (int j = 0; j < N; j ++) {
                int pick = rand.nextInt(nodes.size());
                BranchPoint node = nodes.get(pick);
                bp.addChild(node);
                node.setParent(bp);
                nodes.remove(pick);
                M -= 1;
            }
            nodes.add(bp);
            M += 1;
        }
        if (bp != null) {
            Tree t = new Tree(bp);
            t.setInternalLabels();
            return t;
        }
        throw new TreeRuntimeException("Could not create tree: invalid input");
    }

    public static void main(String[] args) {
        Tree phyloTree = Newick.parse("((A:0.6,((B:3.3,(C:1.0,D:2.5)cd:1.8)bcd:5,((E:3.9,F:4.5)ef:2.5,G:0.3)efg:7)X:3.2)Y:0.5,H:1.1)I:0.2");
        System.out.println(phyloTree.root);
        phyloTree.setInternalLabels();
        System.out.println(phyloTree.root);
        try {
            Tree edge1 = Newick.load("/Users/mikael/simhome/ASR/edge1.nwk");
            System.out.println(edge1);
            Alignment aln = new Alignment(EnumSeq.Gappy.loadClustal("/Users/mikael/simhome/ASR/gap1.aln", Enumerable.aacid));
            TreeInstance ti = edge1.getInstance(aln.getNames(), aln.getGapColumn(1));
            Parsimony tip = new Parsimony(ti.getTree());
            System.out.println(tip);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }


}

