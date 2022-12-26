# Usage TrAVIS

Usage: `asr.TrAVIS`
`[<ancestor-seq>]`
`[-nwk <tree-file> -out <output-file-or-dir>]`

`{-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}`

`{-rates <a>}`

`{-seed <random>}`

`{-extants <5(default)>}`

`{-dist <mean-extant-to-root>`

`{-shape <1.1(default)>}`

`{-scale <0.2(default)>}`

`{-indel <param>}`

`{-gap}`

`{-format <FASTA(default)|CLUSTAL|DOT|TREE|RATES|DIR>}`

where `tree-file` is a phylogenetic tree on Newick format,
`output-file-or-dir` is the filename or name of directory of results, and
`-gap` means that the gap-character is included in the resulting output (default for CLUSTAL format).

### Notes:

Evolutionary models for proteins include Jones-Taylor-Thornton (default), Dayhoff-Schwartz-Orcutt,
Le-Gasquel and Whelan-Goldman;
DNA models include Jukes-Cantor and Yang (general reversible process model).

Tree is set to have specified extants and gets distances from the Gamma distribution, with parameters:
`shape` (aka `a` and `K`), `scale` (aka `b`, where `Lambda=1/b`);
mean distance to root is used to scale distances in the tree (noting that greater number of extants
indirectly amplifies the time scope of the tree).

Position-specific evolutionary rates from a Gamma distribution with specified parameter `a` and mean 1;
if rate is unspecified, a uniform rate 1 is used.
