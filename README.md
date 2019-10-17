# Introduction

`bnkit` is a collection of classes that aim to help implement a Bayesian network in Java (at least version 8 is required).
Currently it allows you to define a Bayesian network structure, use the network for inference,
and some learning from data.
It deals with missing values, both discrete and continuous variables.
There are reasonably efficient implementations of variable elimination, and expectation-maximisation
Iearning. `bnkit` does _not_ do structure learning.

The current version has been used for a number of studies now, but should still be regarded as work-in-progress.
It is a complete re-write of an earlier (in-house) version the Boden lab has been using for research since 2009 or so.

To generate documentation use javadoc. There are a few examples in the code and a very
brief tutorial in the `bn` package documentation.

# GRASP-suite

As part of work funded by the Australian Research Council, the inference engine in bnkit has been bundled with code for phylogenetic analysis. In particular, the `reconstruction` package interfaces with services for ancestral sequence reconstruction (ASR). GRASP (https://github.com/bodenlab/ASR) implements a web server for ASR. The portal for everything around GRASP is [https://github.com/bodenlab/GRASP-suite])(https://github.com/bodenlab/GRASP-suite).

Users who would like to automate ASR, there is now [a limited-feature, command-line version of GRASP](https://github.com/bodenlab/GRASP-suite).
