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
        // create ancestor nodes by picking N children, connecting them via a new ancestor,
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

    /**
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture
     * @return a sample from the mixture
     */
    public static double sampleMixture(GammaDistrib[] mixture, double[] priors) {
        double r = Math.random();
        double sum = 0.0;
        for (int i = 0; i < mixture.length; i++) {
            sum += priors[i];
            if (r < sum)
                return mixture[i].sample();
        }
        return mixture[mixture.length - 1].sample();
    }

    /**
     * @param mixture a mixture of Gamma distributions, a uniform prior is assumed
     * @param r a value to evaluate the mixture at
     * @return distribution of the components that generated the value r
     * */
    public static EnumDistrib getComponentDistribution(GammaDistrib[] mixture, double r) {
        double[] lhoods = new double[mixture.length];
        EnumDistrib ed = new EnumDistrib(new Enumerable(mixture.length));
        for (int i = 0; i < mixture.length; i++)
            lhoods[i] = mixture[i].get(r) + 1e-10;
        ed.set(lhoods);
        return ed;
    }
    /**
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param r a value to evaluate the mixture at
     * @return distribution of the components that generated the value r
     * */
    public static EnumDistrib getComponentDistribution(GammaDistrib[] mixture, double[] priors, double r) {
        if (priors == null)
            return getComponentDistribution(mixture, r);
        double[] lhoods = new double[mixture.length];
        EnumDistrib ed = new EnumDistrib(new Enumerable(mixture.length));
        for (int i = 0; i < mixture.length; i++)
            lhoods[i] = mixture[i].get(r) * priors[i] + 1e-10;
        ed.set(lhoods);
        return ed;
    }

    /**
     * Samples the component that has generated the given value from a mixture of Gamma distributions.
     *
     * This method evaluates the mixture of Gamma distributions at the given value and samples
     * the component that generated the value and returns its index.
     *
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param r a value to evaluate the mixture at
     * @return sample (the index of) the component that has generated the value r
     * */
    public static int sampleComponent(GammaDistrib[] mixture, double[] priors, double r) {
        return (Integer)getComponentDistribution(mixture, priors, r).sample();
    }

    /**
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param r a value to evaluate the mixture at
     * @return sample the component that has generated the value r
     * */
    public static int getMaxComponent(GammaDistrib[] mixture, double[] priors, double r) {
        return (Integer)getComponentDistribution(mixture, priors, r).getMaxIndex();
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     *
     * This method samples the components of a mixture of Gamma distributions based on the given rates, assuming uniform priors.
     * It performs a specified number of samples and returns an estimated distribution of the components.
     *
     * @param mixture a mixture of Gamma distributions
     * @param rates the values to evaluate the mixture for
     * @param nsamples number of samples from which the estimate is based
     * @return an estimate of the distribution of the components that generated the rates
     * */
    public static EnumDistrib estComponentDistribution(GammaDistrib[] mixture, double[] rates, int nsamples) {
        return estComponentDistribution(mixture, null, rates, nsamples);
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     *
     * This method samples the components of a mixture of Gamma distributions based on the given rates and priors.
     * It performs a specified number of samples and returns an estimated distribution of the components.
     *
     * @param mixture a mixture of Gamma distributions
     * @param priors prior probabilities for the components of the mixture; set to null for a uniform prior
     * @param rates the values to evaluate the mixture for
     * @param nsamples number of samples from which the estimate is based
     * @return an estimate of the distribution of the components that generated the rates
     * */
    public static EnumDistrib estComponentDistribution(GammaDistrib[] mixture, double[] priors, double[] rates, int nsamples) {
        if (rates == null || rates.length == 0)
            throw new RuntimeException("Invalid rates: cannot be empty");
        int draws = 0;
        int[] samples = new int[mixture.length];
        while (draws < nsamples) {
            for (int i = 0; i < rates.length; i++) {
                int comp = sampleComponent(mixture, priors, rates[i]);
                samples[comp] += 1;
                draws += 1;
            }
        }
        double[] p = new double[mixture.length];
        for (int k = 0; k < samples.length; k++)
            p[k] = (double)samples[k] / (double)draws;
        EnumDistrib ed = new EnumDistrib(new Enumerable(mixture.length));
        ed.set(p);
        return ed;
    }

    /**
     * Estimates the distribution of the components that generated the given rates.
     * This method uses expectation maximization to estimate the distribution of the components that generated the rates.
     *
     * @param ncomponents number of components in the mixture
     * @return a mixture of Gamma distributions fitted to the distances in the tree; null if EM fails
     */
    public Mixture getGammaMixture(int ncomponents) {
        double[] data = getValidDistances();
        Arrays.sort(data);
        double[] subset = new double[data.length/ncomponents];
        GammaDistribution[] gamma = new GammaDistribution[ncomponents];
        Mixture.Component[] components = new Mixture.Component[ncomponents];
        for (int i = 0; i < ncomponents; i++) {
            System.arraycopy(data, i*subset.length, subset, 0, subset.length);
            gamma[i] = GammaDistribution.fit(subset);
            for (int j = 0; j < subset.length; j++) {
                // System.out.print(subset[j] + " ");
            }
            components[i] = new Mixture.Component(1.0/ncomponents, gamma[i]);
        }
        try {
            Mixture mixture = ExponentialFamilyMixture.fit(data, components, 0, 200, 1E-4);
            return mixture;
        } catch (StackOverflowError e) { // EM will occasionally throw a StackOverflowError exception
            System.out.println("Stack Overflow during EM fitting. Falling back to single component.");
            return getSingleGammaMixture(data);
        }
    }

    private Mixture getSingleGammaMixture(double[] data) {
        GammaDistribution single = GammaDistribution.fit(data);
        Mixture.Component[] singleComponent = {
                new Mixture.Component(1.0, single)
        };
        return new Mixture(singleComponent);
    }

    public static GammaDistrib[] getGammaDistribs(Mixture mixture) {
        int ncomponents = mixture.components.length;
        GammaDistrib[] gds = new GammaDistrib[ncomponents];
        for (int i = 0; i < ncomponents; i++) {
            double a = ((GammaDistribution)(mixture.components[i].distribution)).k;
            double b = ((GammaDistribution)(mixture.components[i].distribution)).theta;
            gds[i] = new GammaDistrib(a, 1/b);
        }
        return gds;
    }

    public static double[] getGammaPriors(Mixture mixture) {
        int ncomponents = mixture.components.length;
        double[] priors = new double[ncomponents];
        for (int i = 0; i < ncomponents; i++) {
            priors[i] = mixture.components[i].priori;
        }
        return priors;
    }

    public int[] getMixtureLabels(Mixture mixture) {
        int[] labels = new int[getSize()];
        for (int i = 0; i < getSize(); i++)
            labels[i] = this.sampleComponent(getGammaDistribs(mixture), getGammaPriors(mixture), distance[i]);
        return labels;
    }

    private static int getBin(double[] thresholds, double value) {
        for (int i = 0; i < thresholds.length - 1; i++) {
            if (value <= thresholds[i])
                return i;
        }
        return thresholds.length - 1;
    }

    /**
     * Create a random tree
     * @param nLeaves number of leaves
     * @param gds Gamma distributions making up a mixture specifying distance between branchpoints
     * @param priors prior probabilities for the components of the mixture
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree RandomMixture(int nLeaves, long seed, GammaDistrib[] gds, double[] priors, int max_desc, int min_desc) {
        String[] leaves = new String[nLeaves];
        for (int i = 0; i < leaves.length; i ++)
            leaves[i] = "A" + (i + 1);
        return RandomMixture(leaves, seed, gds, priors, max_desc, min_desc);
    }

    /**
     * Create a random tree
     * @param nLeaves number of leaves
     * @param gds Gamma distributions making up a mixture specifying distance between branchpoints
     * @param priors prior probabilities for the components of the mixture indexed by log2 widths at branchpoint
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree RandomMixture(int nLeaves, long seed, GammaDistrib[] gds, double[][] priors, double[] thresholds, int max_desc, int min_desc) {
        String[] leaves = new String[nLeaves];
        for (int i = 0; i < leaves.length; i ++)
            leaves[i] = "A" + (i + 1);
        return RandomMixture(leaves, seed, gds, priors, thresholds, max_desc, min_desc);
    }

    /**
     * Create a random tree
     * @param leafLabels labels to be assigned to leaves
     * @param seed random seed
     * @param gds Gamma distributions making up a mixture specifying distance between branchpoints
     * @param priors prior probabilities for the components of the mixture
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree RandomMixture(String[] leafLabels, long seed, GammaDistrib[] gds, double[] priors, int max_desc, int min_desc) {
        Random rand = new Random(seed);
        for (GammaDistrib gd : gds)
            gd.setSeed(seed);
        List<BranchPoint> nodes = new ArrayList<>(leafLabels.length);
        Map<BranchPoint, Integer> widths = new HashMap<>();
        // the initial, complete set of nodes that need ancestors
        for (String label : leafLabels) {
            BranchPoint bp = new BranchPoint(label);
            bp.setDistance(Math.max(sampleMixture(gds, priors), 0.0000001));
            nodes.add(bp);
            widths.put(bp, 1); // width at leaf is 1
        }
        // create ancestor nodes by picking N children, connecting them via a new ancestor,
        // then updating the list of nodes yet to be allocated an ancestor
        int M = nodes.size();
        BranchPoint bp = M > 0 ? nodes.get(0) : null;
        while (M > 1) {
            bp = new BranchPoint("");
            // how many nodes to pick?
            int N = Math.min(rand.nextInt(max_desc - min_desc + 1) + min_desc, nodes.size());
            int mywidth = 0;
            for (int j = 0; j < N; j ++) {
                int pick = rand.nextInt(nodes.size());
                BranchPoint node = nodes.get(pick);
                Integer width = widths.get(node);
                if (width != null) {
                    mywidth += width;
                    widths.remove(node);
                }
                bp.addChild(node);
                node.setParent(bp);
                nodes.remove(pick);
                M -= 1;
            }
            widths.put(bp, mywidth);
            bp.setDistance(Math.max(sampleMixture(gds, priors), 0.0000001));
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

    /**
     * Create a random tree
     * @param leafLabels labels to be assigned to leaves
     * @param seed random seed
     * @param gds Gamma distributions making up a mixture specifying distance between branchpoints
     * @param priors prior probabilities for the components of the mixture indexed by log2 widths at branchpoint
     * @param max_desc maximum number of descendants of an ancestor
     * @param min_desc minimum number of descendants
     * @return a tree where distances are drawn from a specified Gamma density and with branching factors drawn from a specified uniform distribution
     * @throws TreeRuntimeException if something is wrong with the input labels
     */
    public static Tree RandomMixture(String[] leafLabels, long seed, GammaDistrib[] gds, double[][] priors, double[] thresholds, int max_desc, int min_desc) {
        Random rand = new Random(seed);
        for (GammaDistrib gd : gds)
            gd.setSeed(seed);
        List<BranchPoint> nodes = new ArrayList<>(leafLabels.length);
        Map<BranchPoint, Integer> widths = new HashMap<>();
        // the initial, complete set of nodes that need ancestors
        for (String label : leafLabels) {
            BranchPoint bp = new BranchPoint(label);
            bp.setDistance(Math.max(sampleMixture(gds, priors[0]), 0.0000001));
            nodes.add(bp);
            widths.put(bp, 1); // width at leaf is 1
        }
        // create ancestor nodes by picking N children, connecting them via a new ancestor,
        // then updating the list of nodes yet to be allocated an ancestor
        int M = nodes.size();
        BranchPoint bp = M > 0 ? nodes.get(0) : null;
        while (M > 1) {
            bp = new BranchPoint("");
            // how many nodes to pick?
            int N = Math.min(rand.nextInt(max_desc - min_desc + 1) + min_desc, nodes.size());
            int mywidth = 0;
            for (int j = 0; j < N; j ++) {
                int pick = rand.nextInt(nodes.size());
                BranchPoint node = nodes.get(pick);
                Integer width = widths.get(node);
                if (width != null) {
                    mywidth += width;
                    widths.remove(node);
                }
                bp.addChild(node);
                node.setParent(bp);
                nodes.remove(pick);
                M -= 1;
            }
            widths.put(bp, mywidth);
            double log2width = Math.log(mywidth)/Math.log(2);
            bp.setDistance(Math.max(sampleMixture(gds, priors[getBin(thresholds, log2width)]), 0.0000001));
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

    /**
     * Fits the root-to-leaf distances to match a target distribution by swapping their order,
     * not changing the composition.
     *
     * This method iteratively adjusts the distances in the tree to better fit the specified target distribution.
     * It performs a specified number of iterations, where in each iteration it samples pairs of leaves and swaps
     * their distances to improve the fit to the target distribution.
     *
     * @param niter the number of iterations to perform
     * @param target the target distribution to fit the distances to
     */
    public void fitDistances(int niter, Distrib target, long seed) {
        Random rand = new Random(seed);
        Double[][] dmat = getDistance2RootMatrix();
        int[][] connected = new int[this.getNLeaves()][];
        for (int i = 0; i < connected.length; i++)
            connected[i] = IdxTree.getConnected(i, 1, dmat);
        for (int i = 0; i < niter; i++) {
            // Step 1: sample a pair of leaves j1 v j2, (TODO: currently "greedy" not stochastic)
            // where j1 is "good" (root distance does not need to be changed) and
            // j2 is "bad" (root distance should be changed, as a consequence of swapping distances in j1)
            double[] leafdists = IdxTree.getLeafDistances(dmat);
            double[][] logodds = new double[leafdists.length][];
            double sumodds = 0;
            int nomj1 = 0, nomj2 = 0;
            for (int j1 = 0; j1 < leafdists.length; j1++) {
                double p1 = target.get(leafdists[j1]);
                int[] conn = connected[j1];
                logodds[j1] = new double[conn.length];
                int cnt = 0;
                for (int j2 : conn) {
                    double p2 = target.get(leafdists[j2]);
                    logodds[j1][cnt] = Math.max(Math.log(p1 / p2), 0);
                    sumodds += logodds[j1][cnt];
                    cnt ++;
                }
            }
            double toss = rand.nextDouble();
            double accum = 0;
            boolean found = false;
            for (int j1 = 0; j1 < logodds.length; j1 ++) {
                for (int cnt = 0; cnt < logodds[j1].length; cnt ++) {
                    int j2 = connected[j1][cnt];
                    double normalised = logodds[j1][cnt] / sumodds;
                    if (toss < (normalised + accum)) {
                        nomj1 = j1;
                        nomj2 = j2;
                        found = true;
                        break;
                    }
                    accum += normalised;
                }
                if (found)
                    break;
            }

            // Step 2: identify what branches b1 v b2 that should be swapped in j1,
            // so that ...
            Set<Integer> shared = new HashSet<>(); // branches shared between j1 and j2
            Set<Integer> unique = new HashSet<>(); // branches unique to j1
            for (int b = 1; b < dmat[nomj1].length; b ++) {
                if (dmat[nomj1][b] != null) {
                    if (dmat[nomj2][b] != null)
                        shared.add(b);
                    else
                        unique.add(b);
                }
            }
            // Step 3: make a choice,
            Integer favb1 = null, favb2 = null;
            // two possibilities...
            if (leafdists[nomj2] < leafdists[nomj1]) { // increase j2 relative j1:
                // pick largest unique to j1
                for (Integer b1 : unique)
                    if (favb1 == null)
                        favb1 = b1;
                    else if (dmat[nomj1][b1] > dmat[nomj1][favb1])
                        favb1 = b1;
                // pick smallest shared between j1 and j2
                for (Integer b2 : shared)
                    if (favb2 == null)
                        favb2 = b2;
                    else if (dmat[nomj1][b2] < dmat[nomj1][favb2])
                        favb2 = b2;
            } else { // decrease j2 relative j1
                // pick smallest unique to j1
                for (Integer b1 : unique)
                    if (favb1 == null)
                        favb1 = b1;
                    else if (dmat[nomj1][b1] < dmat[nomj1][favb1])
                        favb1 = b1;
                // pick largest shared between j1 and j2
                for (Integer b2 : shared)
                    if (favb2 == null)
                        favb2 = b2;
                    else if (dmat[nomj1][b2] > dmat[nomj1][favb2])
                        favb2 = b2;
            }
            // Step 4: make the swap in the distance matrix
            ;
            double x1 = dmat[nomj1][favb1];
            double x2 = dmat[nomj1][favb2];
            for (int leaf = 0; leaf < dmat.length; leaf ++) {
                if (dmat[leaf][favb1] != null)
                    dmat[leaf][favb1] = x2;
                if (dmat[leaf][favb2] != null)
                    dmat[leaf][favb2] = x1;
            }
        }
        setDistance2RootMatrix(dmat);
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


    public static Tree generateTreeFromMixture(Tree tree1,  int NCOMP, long SEED, int NITER) {
        // First we process an already loaded tree or synthesise a new tree
        GaussianDistrib gds1 = null;
        gds1 = IdxTree.getLeafDistanceDistribution(tree1.getDistance2RootMatrix());


        // Now we are estimating a mixture of Gamma distributions from the first/source tree (loaded or synthesised)
        // 1. ONE mixture of a specified number of Gamma distributions
        Mixture mixture = null;
        mixture = tree1.getGammaMixture(NCOMP); // find mixture and weights for each component

        int compwmax = 0;
        int realCompCount = mixture.components.length;
        double[] shapes = new double[realCompCount];
        double[] scales = new double[realCompCount];
        double[] priors = new double[realCompCount];
        GammaDistrib[] gds = new GammaDistrib[realCompCount];
        for (int i = 0; i < realCompCount; i++) {
            GammaDistribution dist = (GammaDistribution) (mixture.components[i].distribution);
            double a = shapes[i] = dist.k;
            double b = scales[i] = dist.theta;
            gds[i] = new GammaDistrib(a, 1.0 / b);
            double weight = priors[i] = mixture.components[i].priori;

            if (weight > mixture.components[compwmax].priori) {
                compwmax = i;
            }
        }


        // 2. Generate a new tree based on the above mixture of Gamma distributions
        //   a) assume that branch lengths are uniformly distributed across topology
        Tree tree2 = Tree.RandomMixture(tree1.getNLeaves(), SEED, gds, priors, 2, 2);


        //   b) adjust the placement of distances to better fit the original distribution
        tree2.fitDistances(NITER, gds1, SEED + 202);

        return tree2;
    }


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
        Mixture mixture = null;
        mixture = tree1.getGammaMixture(NCOMP); // find mixture and weights for each component
        if (mixture == null) {
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
            double a = shapes[i] = ((GammaDistribution) (mixture.components[i].distribution)).k;
            double b = scales[i] = ((GammaDistribution) (mixture.components[i].distribution)).theta;
            gds[i] = new GammaDistrib(shapes[i], 1 / scales[i]);
            double weight = priors[i] = mixture.components[i].priori;
            if (weight > mixture.components[compwmax].priori)
                compwmax = i;
            if (VERBOSE)
                System.out.println("Gamma " + i +"\tShape (k, a in Matlab)= " + a + "\tScale (theta, b in Matlab)= " + b + " * " + weight);
            if (MATLAB)
                System.out.println("y" + (i + 1) + " = gampdf(x, " + a + ", " + b + ") * " + weight + ";");
        }
        // 2. Generate a new tree based on the above mixture of Gamma distributions
        //   a) assume that branch lengths are uniformly distributed across topology
        tree2 = Tree.RandomMixture(tree1.getNLeaves(), SEED, gds, priors, 2, 2);
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

