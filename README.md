# Framework for testing calcium-based synaptic plasticity models
Hello there! This is a Matlab framework for testing calcium-based model of plasticity both at the level of a synapse and at the level of a network of neurons.

## Introduction
This repository contains a framework for testing models of calcium-based synaptic plasticty.
Voltage-dependent calcium ions influx into both the presynaptic and postsynaptic sides have been shown to be both necessary and sufficient for induction of synaptic plasticity.

## How to run this?
- To test a plasticity rule at the synapse level, parameterize and run simuSynapse.m
- To implement a calcium-based rule at the network level (only pheno rule supported at this moment), parameterize and run simuNetwork.m

## Available models
Two models are included:
- naive model, with a single variable for synaptic plasticity. This model is a simplified version of the model introduced by Grauper & Brunel in [1].
- 3-vars phenomenological model, with a phosphorylation variable and a plasticity variable, dependent on the former by a noise-dependent transform

## Available networks
Right now the type of network implemented is that studied by N. Brunel in [2]. Additional network architectures will be added in the futture

## Available tests
The framework allows for the following tests, all launched form Models/simu.m:
- run a single instance of plasticity, starting from chosen initial conditions and provided an history of calcium spikes
- compute STDP curves for the set of parameters defined in the "environment" section
- compute STDP dependency on frequency
- compute STDP dependency on number of pairings
- plotting the dependency on frequency of number of pairings against experimental data

## Defining a calcium history
### How is calcium let into the synapse?
Calcium is let into the synapse both on the presynaptic and on the postsynaptic side.
Following an action potential's arrival at the synapse on the presynaptic side, a burst of calcium is released through voltage-dependent ion channels, then the calcium is eliminated.
A similar phenomenon of calcium influx occurs when the postsynaptic neuron fires and its action potential backpropagates to the synapse (bAP).
This causes NMDA channels to become actionnable by presynaptic stimulation.
As we considered simplified models, we assume that both the presynaptic and postsynaptic calcium history play a role in plasticity.
Here they act linearly through their sum, but generalized models can be created easily.

To account for the sudden influx of calcium, we use a burst model where the calcium concentration exponentially decays back to zero.

### How do I provide a calcium-history to test my model?
The different methods take as argument a presynaptic and a postsynaptic calcium history, in the form of a set of burst times.
Note that the model does not yet account for variability in spiking (for instance through Poisson processes).

## Coming soon
- Modules for easy analysis at network level and pluing of new plasticity rules and network architectures

## References
[1] Michael Graupner & Nicolas Brunel, Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location, March 6th 2012, PNAS Vol. 109, N. 10, 3991 3996

[2] Nicolas Brunel, Dynamics of Sparsely Connected Networks of Excitatory and Inhibitory Spiking Neurons, June 9th 1999, https://link.springer.com/article/10.1023/A:1008925309027
