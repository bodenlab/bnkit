<div >

# bnkit - Bayesian Network Toolkit

</div>

---

## Table of Contents

- [Introduction](#introduction)
- [What is bnkit used for?](#what-is-bnkit-used-for)
- [Quick Start](#quick-start)

---

## Introduction

`bnkit` is a collection of classes for implementing Bayesian networks in Java (**Java 11+ required**). It lets you:

- Define a Bayesian network structure
- Run inference over that network
- Learn from data — including data with missing values, both discrete and continuous

Under the hood it includes reasonably efficient implementations of **variable elimination** and **expectation-maximisation (EM) learning**.

> `bnkit` does *not* do structure learning, and does *not* implement dynamic Bayesian networks.

The current version has been used across a number of published studies, but is still best regarded as **work-in-progress**. It's a complete rewrite of an earlier in-house version the Boden lab has used for research since around 2009.

> 📚 **Docs:** Generate documentation with `javadoc`. A handful of code examples and a brief tutorial live in the `bn` package documentation.

---

## What is bnkit used for?

Two tools currently build on `bnkit`, for two different kinds of analysis:

| Tool | What it does | Docs |
|---|---|---|
| **[GRASP](https://github.com/bodenlab/GRASP)** | Graphical Representation of Ancestral Sequence Predictions — ancestral sequence reconstruction by maximum likelihood, scalable to very large datasets. Also implements several algorithms for inferring insertions and deletions | [`docs/graspcmd.md`](docs/graspcmd.md) |
| **TreeGazer** | Annotates ancestor (internal) and extant (leaf) nodes on a phylogenetic tree, using known properties (discrete or continuous) at a subset of nodes. Also estimates prediction uncertainty and identifies which nodes are most informative about others | [`docs/treegazer.md`](docs/treegazer.md) |

---

## Quick Start

1. **Install Java 11+** — any OS works.
2. **Get a jar** — download pre-built jars for GRASP and TreeGazer from the [releases page](https://github.com/bodenlab/bnkit/releases).

Want the latest code instead? Clone this repo and build it yourself:

- We recommend [Maven](https://maven.apache.org/) — `pom.xml` files are provided for both GRASP and TreeGazer.
- Full build steps (Maven and IntelliJ) are in each tool's docs: [GRASP](docs/graspcmd.md) · [TreeGazer](docs/treegazer.md)

---

## bnkit is part of [GRASP-suite](https://github.com/bodenlab/GRASP-suite)

Funded by the Australian Research Council, `bnkit`'s inference engine has been bundled with phylogenetic analysis code. The `asr` package in particular interfaces with services for **ancestral sequence reconstruction (ASR)** — [GRASP](https://github.com/bodenlab/GRASP) implements a web server around it.

The portal for everything GRASP-related: **[github.com/bodenlab/GRASP-suite](https://github.com/bodenlab/GRASP-suite)**

---
