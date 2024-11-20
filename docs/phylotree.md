# phylotree

phylotree is a tool to experiment with the generation of phylogenentic trees that mimic natural or synthetic ones.

It completes two major steps; first, phylotree fits a mixture of Gamma distributions on the distances with no 
regard to where they are topologically, or relative to one another. Phylogenetic trees that have been generated from
sequence alignment, in turn following a data collection process, are likely to have a bias in the branch lengths.

The second step is to account for these topological biases of the branch lengths, often observable in terms of the 
leaf-to-root distances that can deviate strongly between a tree that assumes no bias and one based on real data. 
This is done by calibrating the order of the branch lengths in the unbiased tree. 
The result is a tree that is more similar to the original tree in terms of mean leaf-to-root distances.

How is this done? The mixture of Gamma distributions is used to tentatively assign distances to all branches in 
the target tree independently of their position. Indeed, the mean distance of branches in the source tree is 
generally well-replicated by the mean in the target tree. phylotree estimates a (Gaussian) distribution of the 
leaf-to-root distances in the source tree, which serves to evaluate the topological bias of distances. 
In fact, for many trees, the mean leaf-to-root distances in natural and synthetic trees are different. 
We note that a path in the target tree, defines a leaf-to-root order of distances. Importantly, changing the 
order in one path does not influence its leaf-to-root distance, but it does change the 
leaf-to-root distances of paths that share one of these distances. 

Using the Gaussian distribution of leaf-to-root distances from the source tree, phylotree samples pairs of 
leaves in proportion to the log odds ratio of their  
leaf-to-root distances; this implies that we pick a primary leaf which has a high likelihood of being generated 
by the source distribution, and a secondary leaf which has a low likelihood. Then distances on the path from the 
primary leaf to the root are chosen such that swapping them has the greatest effect in the direction that would 
make the secondary leaf-to-root distance improve its likelihood. The two distances are updated and the process
is repeated.

### Using phylotree 

`Usage: `

phylotree typically loads or synthesises one tree, then tries to mimic the
(pseudo-biological) properties of that first tree when synthesising a target tree.

The first tree can be provided as a Newick file or by specifying a range of
parameters.

The second tree is based on fitting a mixture of Gamma distributions to the
first tree, calibrated to topological biases of branch lengths.

`Usage: dat.phylo.Tree` or phylotree `[options]`

    -l, --load <String>>      Specify the file name from which a source tree is loaded
    
    -s, --save <String>       Specify the file name to which the target tree is saved
    
    -a, --alpha <Double>      Specify the alpha parameter for the source Gamma distribution
    
    -b, --beta <Double>       Specify the beta parameter for the source Gamma distribution
    
    -n, --nleaves <Integer>   Specify the number of leaves in the source tree
    
    -d, --meandist <Double>   Specify the mean distance to root in source tree
    
    --seed <Integer>          Specify the random seed
    
    --ncomp <Integer>         Specify the number of components in target Gamma mixture
    
    --iter <Integer>          Iterations to calibrate topological branch length biases
    
    -v, --verbose             Enable verbose mode
    
    -h, --help                Show this help message

### Notes:



### Examples

```bash