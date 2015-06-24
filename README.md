# Introduction #

bnkit is a collection of classes that aim to help implement a Bayesian network in Java.
Currently it allows you to define a Bayesian network structure, use the network for inference,
and some learning from data.
It deals with missing values, both discrete and continuous variables.
There are basic implementations of variable elimination, and expectation-maximisation
Iearning. bnkit does _not_ do structure learning.

The current version is very much a work-in-progress.
It is a complete re-write of an earlier (in-house) version we have been using for research.
So far only basic functionality has been moved across--so some attractive features are missing.
We also expect to extend the functionality and address efficiency concerns.

To generate documentation use javadoc. There are a few examples in the code and a very
brief tutorial in the "bn" package documentation.

The GUI is at a very early stage and is not suited to serious use.
To compile the GUI, you need JGraphX (v.6 or later; https://github.com/jgraph/jgraphx).

# Details #

To be written
