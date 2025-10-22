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
import bn.Distrib;
import bn.prob.EnumDistrib;
import bn.prob.GammaDistrib;
import bn.prob.GaussianDistrib;
import dat.EnumSeq;
import dat.EnumSeq.Alignment;
import dat.Enumerable;
import dat.file.Newick;
import smile.stat.distribution.ExponentialFamilyMixture;
import smile.stat.distribution.GammaDistribution;
import smile.stat.distribution.Mixture;
import stats.RateModel;

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
     * Create a random tree topology, i.e. a tree without distances
     * @param nLeaves number of leaves
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree topology with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree RandomTopology(int nLeaves, long seed, int max_desc, int min_desc) {
        String[] leaves = new String[nLeaves];
        for (int i = 0; i < leaves.length; i ++)
            leaves[i] = "A" + (i + 1);
        return RandomTopology(leaves, seed, max_desc, min_desc);
    }

    /**
     * Create a random tree topology, i.e. a tree without distances
     * @param leafLabels labels to be assigned to leaves
     * @param seed random seed
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree topology with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree RandomTopology(String[] leafLabels, long seed, int max_desc, int min_desc) {
        Random rand = new Random(seed);
        List<BranchPoint> nodes = new ArrayList<>(leafLabels.length);
        // the initial, complete set of nodes that need ancestors
        for (String label : leafLabels) {
            BranchPoint bp = new BranchPoint(label);
            nodes.add(bp);
        }
        // create ancestor nodes by picking N children, connecting them via a new ancestor,
        // then updating the list of nodes yet to be allocated an ancestor
        int M = nodes.size();
        BranchPoint bp = M > 0 ? nodes.get(0) : null;
        while (M > 1) {
            bp = new BranchPoint("");
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
        throw new TreeRuntimeException("Could not create tree topology: invalid input");
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
     * @param seed random seed
     * @param gamma_k parameter for gamma specifying distance between branchpoints ("shape" aka alpha)
     * @param gamma_lambda parameter for gamma specifying distance between branchpoints ("scale", also referred to as 1/beta, where beta is "rate")
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree Random(String[] leafLabels, long seed, double gamma_k, double gamma_lambda, int max_desc, int min_desc) {
        GammaDistrib gd = new GammaDistrib(gamma_k, gamma_lambda);
        gd.setSeed(seed);
        return Random(leafLabels, gd, max_desc, min_desc, seed);
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
    public static Tree Random(int nLeaves, RateModel distrib, int max_desc, int min_desc, long seed) {
        String[] leaves = new String[nLeaves];
        for (int i = 0; i < leaves.length; i ++)
            leaves[i] = "A" + (i + 1);
        return Random(leaves, distrib, max_desc, min_desc, seed);
    }


    public static Tree Random(String[] leafLabels, RateModel distrib, int max_desc, int min_desc, long seed) {
        Random rand = new Random(seed);
        List<BranchPoint> nodes = new ArrayList<>(leafLabels.length);
        // the initial, complete set of nodes that need ancestors
        for (String label : leafLabels) {
            BranchPoint bp = new BranchPoint(label);
            bp.setDistance(Math.max(distrib.sample(), 0.0000001));
            nodes.add(bp);
        }
        // create ancestor nodes by picking N children, connecting them via a new ancestor,
        // then updating the list of nodes yet to be allocated an ancestor
        int M = nodes.size();
        BranchPoint bp = M > 0 ? nodes.get(0) : null;
        while (M > 1) {
            bp = new BranchPoint("");
            bp.setDistance(Math.max(distrib.sample(), 0.0000001));
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

    public IdxTree getIdxTree() {
        BranchPoint[] bps = Tree.straightenTree(root);
        return new IdxTree(bps);
    }

    public static void main0(String[] args) {
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

    private static void usage() {
        System.out.println("phylotree typically loads or synthesises one tree, then tries to mimic the \n(pseudo-biological) properties of that first tree when synthesising a target tree.");
        System.out.println("The first tree can be provided as a Newick file or by specifying a range of \nparameters.");
        System.out.println("The second tree is based on fitting a mixture of Gamma distributions to the \nfirst tree, calibrated to topological biases of branch lengths.");
        usage(0, null);
    }
    private static void usage(int errno, String errmsg) {
        if (errmsg != null && errno > 0)
            System.err.println("Error " + errno + ": " + errmsg);
        System.out.println("Usage: phylotree [options]");
        System.out.println("Options:");
        System.out.println("  -l, --load <String>>      Specify the file name from which a source tree is loaded");
        System.out.println("  -s, --save <String>       Specify the file name to which the target tree is saved");
        System.out.println("  -a, --alpha <Double>      Specify the alpha parameter for the source Gamma distribution");
        System.out.println("  -b, --beta <Double>       Specify the beta parameter for the source Gamma distribution");
        System.out.println("  -n, --nleaves <Integer>   Specify the number of leaves in the source tree");
        System.out.println("  -d, --meandist <Double>   Specify the mean distance to root in source tree");
        System.out.println("      --seed <Integer>      Specify the random seed");
        System.out.println("      --ncomp <Integer>     Specify the number of components in target Gamma mixture");
        System.out.println("      --iter <Integer>      Iterations to calibrate topological branch length biases");
        //System.out.println("      --res <Integer>       Specify the resolution decomposing mixture weights");
        System.out.println("  -v, --verbose             Enable verbose mode");
        System.out.println("      --matlab              Print out MATLAB commands for plotting");
        System.out.println("  -h, --help                Show this help message");
        System.exit(errno);
    }

    static boolean VERBOSE = false;
    static boolean MATLAB = false;
    static long SEED = 0;
    static Double alpha = null, beta = null;
    static Integer nLeaves = null;
    static Double MEANDIST = null;
    static String savetree = null;
    static int TREERESOLUTION = 5; // number of bins for collecting different gamma mixture weights
    static Integer NCOMP = 3;
    static Integer NITER = 100;


    public static void main(String[] args) {
        Tree tree1 = null, tree2 = null;
        for (int i = 0; i < args.length; i ++) {
            switch (args[i]) {
                case "-h":
                case "--help":
                    usage();
                    break;
                case "-a":
                case "--alpha":
                    if (i + 1 < args.length) {
                        alpha = Double.parseDouble(args[++i]);
                    } else {
                        usage(5, "Missing ALPHA[double] after " + args[i]);
                    }
                    break;
                case "-b":
                case "--beta":
                    if (i + 1 < args.length) {
                        beta = Double.parseDouble(args[++i]);
                    } else {
                        usage(6, "Missing BETA[double] after " + args[i]);
                    }
                    break;
                case "-n":
                case "--nleaves":
                    if (i + 1 < args.length) {
                        nLeaves = Integer.parseInt(args[++i]);
                    } else {
                        usage(7, "Missing NLEAVES[int] after " + args[i]);
                    }
                    break;
                case "-d":
                case "--meandist":
                    if (i + 1 < args.length) {
                        MEANDIST = Double.parseDouble(args[++i]);
                    } else {
                        usage(8, "Missing MEANDIST[double] after " + args[i]);
                    }
                    break;
                case "--ncomp":
                    if (i + 1 < args.length) {
                        NCOMP = Integer.parseInt(args[++i]);
                    } else {
                        usage(9, "Missing NCOMP[int] after " + args[i]);
                    }
                    break;
                case "--iter":
                    if (i + 1 < args.length) {
                        NITER = Integer.parseInt(args[++i]);
                    } else {
                        usage(10, "Missing ITER[int] after " + args[i]);
                    }
                    break;
                case "-l":
                case "--load":
                    if (i + 1 < args.length) {
                        String filename = args[++i];
                        try {
                            tree1 = Newick.load(filename);
                        } catch (IOException ex) {
                            usage(1, "Could not load file " + filename);
                        }
                    } else {
                        usage(2, "Missing FILE[String] after " + args[i]);
                    }
                    break;
                case "-s":
                case "--save":
                    if (i + 1 < args.length) {
                        savetree = args[++i];
                    } else {
                        usage(2, "Missing FILE[String] after " + args[i]);
                    }
                    break;
                case "--seed":
                    if (i + 1 < args.length) {
                        SEED = Long.parseLong(args[++i]);
                    } else {
                        usage(4, "Missing SEED[long] after " + args[i]);
                    }
                    break;
                case "--res":
                    if (i + 1 < args.length) {
                        TREERESOLUTION = Integer.parseInt(args[++i]);
                    } else {
                        usage(11, "Missing RESOLUTION[int] after " + args[i]);
                    }
                    break;
                case "-v":
                case "--verbose":
                    VERBOSE = true;
                    break;
                case "--matlab":
                    MATLAB = true;
                    break;
                default:
                    usage(3, "Aborting because unknown argument " + args[i]);
                    break;
            }
        }

        // First we process an already loaded tree or synthesise a new tree
        GaussianDistrib gds1 = null;
        if (tree1 != null) { // if tree was loaded already...
            gds1 = IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix());
            if (VERBOSE)
                System.out.println("Loaded tree: size = " + tree1.getSize() + "\tLeaf-to-root distance distribution: " + IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix()));
        } else if (alpha != null && nLeaves != null) {
            if (beta == null) // if beta is not set...
                beta = alpha; // mean of gamma is 1.0
            tree1 = Tree.Random(nLeaves, SEED, alpha, beta, 2, 2);
            gds1 = IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix());
            if (VERBOSE)
                System.out.println("Synthesised tree: size = " + tree1.getSize() + "\tLeaf-to-root distance distribution: " + gds1);
            if (MEANDIST != null) {
                gds1.setMean(MEANDIST);
                tree1.fitDistances(NITER, gds1, SEED);
                gds1 = IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix());
                if (VERBOSE)
                    System.out.println("Distance-adjusted tree: size = " + tree1.getSize() + "\tLeaf-to-root distance distribution: " + gds1);
            }
        }

        if (tree1 == null)
            usage(12, "No tree to analyse");

        // Now we are estimating a mixture of Gamma distributions from the first/source tree (loaded or synthesised)
        // 1. ONE mixture of a specified number of Gamma distributions
        GammaDistrib.Mixture gdm = IdxTree.getGammaMixture(tree1, NCOMP, SEED); // find mixture and weights for each component
        if (gdm == null) {
            System.out.println("Could not find a mixture of " + NCOMP + " components");
            System.exit(1);
        }
        if (MATLAB)
            System.out.println("dx = " + Arrays.toString(tree1.getValidDistances()) + "; dy = zeros(size(dx, 2));");
        int compwmax = 0;
        double[] shapes = new double[NCOMP];
        double[] scales = new double[NCOMP];
        double[] priors = new double[NCOMP];
        GammaDistrib[] gds = new GammaDistrib[NCOMP];
        for (int i = 0; i < NCOMP; i++) {
            double a = shapes[i] = gdm.distribs[i].getShape();
            double b = scales[i] = gdm.distribs[i].getScale();
            gds[i] = new GammaDistrib(shapes[i], 1 / scales[i]);
            double weight = priors[i] = gdm.priors[i];
            if (weight > gdm.priors[compwmax])
                compwmax = i;
            if (VERBOSE)
                System.out.println("Gamma " + i +"\tShape (k, a in Matlab)= " + a + "\tScale (theta, b in Matlab)= " + b + " * " + weight);
            if (MATLAB)
                System.out.println("y" + (i + 1) + " = gampdf(x, " + a + ", " + b + ") * " + weight + ";");
        }
        // 2. Generate a new tree based on the above mixture of Gamma distributions
        //   a) assume that branch lengths are uniformly distributed across topology
        tree2 = Tree.Random(tree1.getNLeaves(), gdm, 2, 2, SEED);
        if (VERBOSE) {
            System.out.print("Generated intermediate tree (seed " + SEED + "): size = " + tree2.getSize() + "\t");
            Double[][] dmat2 = tree2.getDistance2RootMatrix();
            GaussianDistrib gd2 = IdxTree.getLeafDistanceDistribution(dmat2);
            System.out.println("Leaf-to-root distance distribution: " + gd2);
        }
        //   b) adjust the placement of distances to better fit the original distribution
        tree2.fitDistances(NITER, gds1, SEED + 202);
        if (VERBOSE) {
            System.out.print("Generated final tree (seed " + SEED + "): size = " + tree2.getSize() + "\t");
            Double[][] dmat2 = tree2.getDistance2RootMatrix();
            GaussianDistrib gd2 = IdxTree.getLeafDistanceDistribution(dmat2);
            System.out.println("Leaf-to-root distance distribution: " + gd2);
        }
        if (savetree != null) {
            try {
                tree2.save(savetree, "newick");
            } catch (IOException ex) {
                System.err.println("Could not save tree to file " + savetree);
            }
        }
        System.out.println("Done.");
    }
}

