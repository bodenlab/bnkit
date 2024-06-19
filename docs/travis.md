# Tracing Ancestors Via Indels and Substitutions (TrAVIS)

TrAVIS is a phylogenetic tree and ancestor sequence genertor.
It forms part of the GRASP-suite (Foley et al., 2022).

TrAVIS samples the Gamma distribution γ(κ, θ) (where κ ∈ {0.5, 1, 2} defines the shape 
and θ = 0.2 the scale) to set distances on each branch (in turn normalised to the mean δ), 
bi-furcating each branch point until the specified number of sequences have been mapped 
as leaves.
For each tree, substitutions, insertions and deletions are randomly introduced at each 
branch point. Starting at the root with an arbitrary amino acid sequence, the sequence at 
each of its children is determined recursively. 
For a sequence, each position is considered as a possible site for a mutation as a Poisson 
process with a relative rate 2r; following the recipe given by Cartwright Cartwright (2008), 
η = e−2rt is the probability that no indel has occurred, which implies a substitution 
(determined via a model like that suggested by Le and Gasquel Le and Gascuel (2008), for 
instance); the probability of an insertion is (1 − η)/2 which is equivalent to that of a 
deletion. Finally, the length of the insertion (or deletion) is given by a sampling the 
Poisson distribution f(λ) (where mean λ = 1, which will give 0 or 1 37% of the time, with 
greater widths less frequent).

### Using TrAVIS

Usage: `asr.TrAVIS`
    `[<ancestor-seq>]`

    `[-nwk <tree-file> -out <output-file-or-dir>]`

    `{-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}`

    `{-load}`

    `{-rates <a>}`

    `{-seed <random>}`

    `{-extants <5(default)>}`

    `{-dist <mean-extant-to-root>`

    `{-shape <1.1(default)>}`

    `{-scale <0.2(default)>}`

    `{-indel <1.0(default)>}`

    `{-delprop <0.5(default)>}`

    `{-indelmodel <zero-truncated-poisson(default)|poisson>	{-lambda  <1(default)>}`

    `{-gap}`

    `{-format <FASTA(default)|CLUSTAL|DOT|TREE|RATES|DIR>}`

    `{-verbose}`

    `{-help}`

    where
    tree-file is a phylogenetic tree on Newick format
    output-file-or-dir is the filename or name of directory of results
    "-gap" means that the gap-character is included in the resulting output (default for CLUSTAL format)
    "-verbose" means that details of evolutionary events are printed out on standard-output)
    "-help" prints out the parameters and instructions of how to use this tool
    

### Notes:

Evolutionary models for proteins include Jones-Taylor-Thornton (default), Dayhoff-Schwartz-Orcutt,
Le-Gasquel and Whelan-Goldman;
DNA models include Jukes-Cantor and Yang (general reversible process model).

Tree is set to have specified extants and gets distances from the Gamma distribution, with parameters:
shape (aka a and K), and
scale (aka b, where Lambda=1/b)

mean distance to root is used to scale distances in the tree (noting that greater number of extants
indirectly amplifies the time scope of the tree).

Position-specific evolutionary rates from a Gamma distribution with specified parameter "a" and mean 1;
if rate is unspecified, a uniform rate 1 is used.

Rates for insertions and deletions are scaled by -indel <multiplier> (multiplier > 1 reduces, multiplier < 1 increases chance of indels).

The proportion of deletions relative to insertions and deletions is given by -delprop <proportion> (proportion < 0.5 means that insertions will dominate.

The mean rate of occurrence in the Poisson distribution of the indel length by lambda	
