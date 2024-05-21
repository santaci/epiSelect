# SOS
## SimOutbreakSelection
A simulation-based framework for exploring and detecting epidemic-driven selection

Dependencies
-------------
`SimOutbreakSelection` (SOS) relies on SLiM (), python (3.x.x), and R (v).
An appropriate environment can me installed using conda or micromamba.

Installation
------------
Preparing Inputs
==============
Input I. Simulating the demographic history of the host
--------------
This can be a well-defined model such as the Gravel model of human evoluton or can be a custom designed one based on curated information from relevant literature or in-house demographic inferencing.
For example purposes, the demographic history of African Cape buffalo and the Medieval Swedish human population have been provided to simulate in SLiM. 


Input II. Simulating the epidemic of interest with selection using SLiM
--------------
This is the most bespoke part of the framework---designing the trajectory and implementing the key components into the simulation:
- Bottleneck induced by the epidemic
- Viability of those who are homozygous for the advantageous allele (i.e. selected variant)
- The starting frequency of the variant under selection
The selection-driven simulations used in our manuscript can be found in the `inputs/` folder of each respective epidemic. Note that each mode of inheritance (mode of selection) is a separate SLiM script.


Input III. Save generations of interest
--------------
Within the SLiM model of epidemic-driven selection, make sure to save the generations of interest for the output tree sequence. The command to do this in SLiM would be:
`sim.treeSeqRememberIndividuals(p2.individuals)`
In this example, `p2` is the population of interest. Doing this in the early part ensures you are keeping those that have died.
