# Introduction #

bnkit is a collection of classes that aim to help implement a Bayesian network in Java.
Currently it allows you to define a Bayesian network structure, use the network for inference,
and some learning from data.
It deals with missing values, both discrete and continuous variables.
There are reasonably efficient implementations of variable elimination, and expectation-maximisation
Iearning. bnkit does _not_ do structure learning.

The current version has been used for a number of studies now, but should still be regarded as work-in-progress.
It is a complete re-write of an earlier (in-house) version the Boden lab has been using for research since 2009 or so.

To generate documentation use javadoc. There are a few examples in the code and a very
brief tutorial in the "bn" package documentation.

As part of work funded by the Australian Research Council, the inference engine in bnkit has been bundled with code for phylogenetic analysis. In particular, the "reconstruction" package interfaces with services for ancestral sequence reconstruction (ASR). GRASP (https://github.com/bodenlab/ASR) implements a web server for ASR.

# GRASP local, terminal version #

Users who would like to automate ASR, there is now a limited-feature, command-line version of GRASP.

## GraspCmd: What can it do? ##

It accepts an alignment (FASTA or Clustal formats) and phylogenetic tree (Newick format) with perfectly concordant labels, to infer ancestor sequences by joint or marginal reconstruction by maximum likelihood. In the process, the program also infers the most parsimonious insertion and deletion events which are internally represented via partial-order graphs; it identifies the _most supported_ path of sequence inclusions at each ancestor as it currently has no function to save the (potentially more complex) partial-order graph. Instead, the program saves all sequences (in the case of joint reconstruction, again as FASTA or Clustal) or one sequence (in the case of marginal reconstrution; optionally with amino acid distributions as a TSV file). It can also re-save the tree with assigned ancestor labels.

## GraspCmd: How do I run it? ##

First, you will need Java version 8 or newer. Any operating system with Java should work, including Mac OS, MS Windows and Linux.

Then, you have a choice: you can clone/download bnkit in its entirety; you may need JUnit 5 testing to get everything working, but this is only required if you want to run software tests--but unless you are a developer you would not need this.

Alternatively, just download the pre-compiled JAR file.

# Wiki #
https://github.com/bodenlab/bnkit/wiki
