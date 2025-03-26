# Introduction

`bnkit` is a collection of classes that aim to help implement a Bayesian network in Java (at least version 11 is required).
Currently, it allows you to define a Bayesian network structure, use the network for inference,
and some learning from data. It deals with missing values, both discrete and continuous variables.
There are reasonably efficient implementations of variable elimination, and expectation-maximisation
learning. `bnkit` does _not_ do structure learning. It is _not_ implementing dynamic Bayesian networks. 

The current version has been used for a number of studies now, but should still be regarded as work-in-progress.
It is a complete re-write of an earlier (in-house) version the Boden lab has been using for research since 2009 or so.

To generate documentation use javadoc. There are a few examples in the code and a very
brief tutorial in the `bn` package documentation.

# bnkit is part of [GRASP-suite](https://github.com/bodenlab/GRASP-suite)

As part of work funded by the Australian Research Council, the inference engine in bnkit has been bundled with code 
for phylogenetic analysis. In particular, the `asr` package interfaces with services for ancestral sequence 
reconstruction (ASR). GRASP (https://github.com/bodenlab/GRASP) implements a web server for ASR. 
The portal for everything around GRASP is [https://github.com/bodenlab/GRASP-suite](https://github.com/bodenlab/GRASP-suite).

## Command-line interface to GRASP

In `bnkit` there is a command-line interface of GRASP `asr.GRASP`. This Java application can prove useful if you want 
to automate tasks, run reconstructions on your own dedicated hardware, and/or access the latest features. 
This version is essentially a command-line interface to the backend features of the web-based service. 
It is worth noting that the web-based version has the advantage of a visual user interface, but that it 
may lack the latest functionality.

If you want to use `asr.GRASP`, go [here](docs/graspcmd.md).

## JSON API to GRASP

Another way to interact with GRASP is to use the JSON API of `asr.GServer`, 
documented [here](docs/json-api.md).

## Configuring GRASP to run in IntelliJ
- Open IntelliJ IDEA
- Select File > New > Project From Version Control
- For URL type http://github.com/bodenlab/bnkit
- For Directory specify where you want to save the project
- Under Run > Edit Configurations
- Use '+' to add Application
- Give it a name, e.g. 'GRASP'
- Specify 'asr.GRASP' as the main class
- Now you can run GRASP (green arrow / ctrl-R / Run > Run GRASP )
- You may need to upgrade your Java version / SDK in your settings (Version 11 or above)
- If you want to add command line arguments, you can add these under Run > Edit Configurations and add them to the 'Program arguments' textbox

## Generate an external JAR using Maven
- Click on Maven (on the far right vertical bar)
- Click on bnkit > Lifecycle > package
- A jar will be created in the <directory>/target/bnkit.jar
