# `SOS`: SimOutbreakSelection


## !!! ``SOS`` is still a WIP.  Please be mindful that there are paths that require the user to manually change in some scripts. !!!


Dependencies
-------------
`SimOutbreakSelection` or `SOS` relies on [`SLiM`](https://github.com/MesserLab/SLiM/releases/download/v3.3.1/SLiM.zip) (3.3.1) , python (3.7.12), and R (4.2.2). 

Crucially `SOS` uses the following modules within python to handle `SLiM` tree sequence outputs: 
* `tskit` (0.2.3) 
* `msprime` (0.7.4)
* `pyslim` (0.403)

An appropriate environment can be installed using conda or micromamba. Micromamba is highly recommended.

Installation
------------
```
# Create a conda environment
micromamba env create -f environment.yml
```

Preparing Inputs
==============
Input I. Simulating the demographic history of the host
--------------
This can be a well-defined model such as the Gravel model of human evoluton or can be a custom designed one based on curated information from relevant literature or in-house demographic inferencing.

For example purposes, the demographic history of African Cape buffalo (`rinderpest/demographies/`) and the Medieval Swedish human population (`plague/demographies/`) have been provided to simulate in SLiM. 


Input II. Simulating the epidemic of interest with selection using SLiM
--------------
This is the most bespoke part of the framework—designing the trajectory and implementing the key components into the simulation:
- Bottleneck induced by the epidemic
- Viability of those who are homozygous for the advantageous allele (i.e. selected variant)
- The starting frequency of the variant under selection
The selection-driven simulations used in our manuscript can be found in the `inputs/` folder of each respective epidemic. Note that each mode of inheritance (mode of selection) is a separate SLiM script.


Input III. Save generations of interest
--------------
Within the SLiM model of epidemic-driven selection, make sure to keep track of Pedigree IDs of the individuals that survive in the generation* of interest. 

This would be done within the `intialize() {}` callback with the following command: 
`initializeSLiMOptions(keepPedigrees=T);`

Survivor pedigree IDs will be used to subset apart the dead individuals from the same `generation`* when we write those IDs out into a separate file. Those IDs will then be used by `SOS` downstream.

In order to keep the dead from the generations of interest within the output tree sequence, we save the individuals in the early part of the generation*

The command to do this in SLiM would be:
`sim.treeSeqRememberIndividuals(p2.individuals);`

In this example, `p2` is the population of interest.

****This is now `cycle` or `tick` for the newer SLiM versions.***

Running `SOS` using Snakemake workflows
==============
Snakemake workflows can be found within the `rules/` folders of each example epidemic in this repository. There are rules to run the demographic history (if you haven't done so already) as well as the epidemic part of the simulation.

The epidemic simulation rules require a configuration file (.YAML) with the necessary parameters such as starting allele frequency (`f`), viability of the homozygous advantageous allele(`v`), and corresponding inheritance `model`* (i.e. `add` or `rec`). Other information such as which generations will be compared (`gen1` and `gen2`), the generation where there are dead to compare to survivors (`dead`), number of samples per grouping(`n`), number of simulations to perform (`sim`), and number of subsamplings per simulation (`resamp`) can be found within the configuration file. All parameters can be adjusted by the user. 
```
## Example configuration YAML

model: "add"
recovery: "grad"
v: 1.0
f: 0.1
sim: 1
n: 100
resamp: 100
deme: 1
gen1: 40001
gen2: 40003
midchrom: 33759515
chrom: [21, 22]
OUTMAIN: "output"
dead:
```


More examples can be found in the `configs/` folder of each respective epidemic.

****NOTE: At the moment only SLiM code for additive and recessive selection have been provided as examples. An equation to implement a dominant model can be found in the Materials & Methods section of our manuscript.***

Quick usage
-------------
Activate the SOS environment
```
micromamba activate SOS
```

An example usage using the demo data (`rinderpest/`) in this repository would be to first run the demographic history workflow:
```
cd rinderpest/
snakemake all -p -s rules/african_buffalo_deme.snakemake --cores 10
```

Once that has completed and we have our tree sequence outputs of the host's demographic history, we can move on to the our epidemic-driven selection workflow.

In this example below, the Snakemake file contains rules to use the demopgraphic history tree sequences as input for our nonWF SLiM simulations. The number of simulations and subsamplings can be found in the configuration file. The `--configfile` path lets us know we are running an additive model with n=100 per group where Dead and Alive are being compared starting at an allele frequency of 0.1 (`DA01`). 
```
snakemake all -p -s rules/buffalo_runs.snakemake --configfile configs/add/100/n100_grad_DA01.yaml --cores 20 --debug-dag
```

Further information can be adjused within the configuration file and the Snakemake file.

Once your epidemic-selection workflow is done, gather the results for those simulations into one mastertable using the command below. You will multiple mastertables if you explore different parameter combinations.
```
## Usage: get_results.sh <model> <generation1> <generation2> <sample_size> [demeGeneration] [midchrom] [nsims] [nresamp] [workdir]

bash get_results.sh add_grad_v1.0_f0.1_sim 40002 40002 100
```

Then calculate power for all the performed simulations of the same viability (`v`).
This script points to the directory where mastertables for the simulations of VAA=1 are stored (`outputs/plots/add/VAA1/`). `3` signifies the number of top candidates to use when considering outliers for all power calculations (this is the default). 

```
## Usage: Rscript calculate_power.R <path to candidate tables> [# of top candidates. Default = 3]

Rscript calculate_power.R outputs/plots/add/VAA1/ 3
```
