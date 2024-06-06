# Activate necessary modules
import subprocess
import argparse
import os
import os.path
import time
from datetime import datetime
import sys
import numpy as np
import random
from art import *

print("Loading menu...")

def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--tree",
        type=str,
        dest='tree',
        help="Path to the tree file without mutations.",
    )

    parser.add_argument(
        "--muTree",
        type=str,
        dest='muTree',
        help="Path to the mutated and recapitated tree file."
    )

    parser.add_argument(
        "--bottle",
        type=int,
        default=58002,
        dest='bottle',
        help="SLiM generation(s) of the bottleneck."
    )

    parser.add_argument(
        "--select",
        type=int,
        dest='select',
        help="Site of selected variant",
        default=13811692
    )

    parser.add_argument(
        "--generations",
        type=int,
        nargs=2,
        default=[58002, 58012],
        dest='generations',
        help="Generations to be compared. (e.g. 58002 58012)."
    )
    parser.add_argument(
        "--chromosomes",
        type=int,
        nargs='+',
        default=[21, 22],
        dest='chromosomes',
        help="Chromosomes simulated in SLiM. (e.g. 21 22)."
    )
    parser.add_argument(
        "--dead",
        type=int,
        help="Provide the generation for the dead of interest. (e.g. 58002)"
    )

    parser.add_argument(
        "--mixCem",
        type=float,
        default=None,
        dest='mixCem',
        help="Proportion of diseased dead in cemetery mixed with unknown survivors (e.g. 0.5)."
    )

    parser.add_argument(
        "--putparents",
        type=str,
        help="Provide the putparent CSV with both dead and alive individuals."
    )

    parser.add_argument(
        "--sampling",
        type=int,
        nargs='+',
        default=100,
        help="Provide the number of samples to compare for each generation/grouping. If comparing with the dead, input that sample size first."
    )

    parser.add_argument(
        "--plots",
        action='store_true',
        help="Plot out results."
    )
    
    parser.add_argument(
        "--winfst",
        action='store_true',
        help="Windowed FST 100kb."
    )

    parser.add_argument(
        "--vcf",
        action='store_true',
        help="Creates files to build VCF."
    )

    parser.add_argument(
        "--nosel",
        action='store_true',
        help="Does not run selscan."
    )

    parser.add_argument(
        "--pmap",
        action='store_true',
        help="Does not create genetic map for iHS. Uses physical positions."
    )

    parser.add_argument(
        "--recmaps",
        type=str,
        help="Provide the directory where recombination maps per chromosome can be found. Files must be named as 'slim_chr#_recode.map'."
    )

    parser.add_argument(
        "--realhet",
        action='store_true',
        help="Counts true heterozygous sites for generation (After)."
    )

    parser.add_argument(
        "--midchrom",
        type=int,
        default=33759515,
        help="Position to divide the chromosome made in SLiM for plots and VCF."
    )

    parser.add_argument(
        "--ihs",
        action='store_true',
        help="Produces inputs for selscan."
    )

    parser.add_argument(
        "--gwas",
        action='store_true',
        help="Produces inputs for GWAS."
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Provide the number of threads for selscan."
    )

    parser.add_argument(
        "--demegen",
        type=int,
        default=58000,
        dest='demegen',
        help="Last SLiM generation of demographic history."
    )

    parser.add_argument(
        "--mu",
        type=float,
        default=2.36e-8,
        dest='mu',
        help="Mutation rate for neutral mutations."
    )

    parser.add_argument(
        "--rec",
        type=float,
        default=1e-8,
        dest='rec',
        help="Uniform recombination rate for recapitation."
    )

    parser.add_argument(
        "--ne",
        type=int,
        default=7310,
        dest='ne',
        help="Ancestral effective population size."
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=int(time.time()),
        dest='seed',
        help="Seed for simulation."
    )

    parser.add_argument(
        "--dir",
        type=str,
        default=None,
        dest='dir',
        help="Working directory if different."
    )
    tprint("SOS",font='defleppard')
    print("SimOutbreakSelection v.1.0")
    print("Cindy Santander and Ida Moltke.")
    args = parser.parse_args()
    print("Using " + str(args.threads) + " thread(s).\n")

    # Create log-file of arguments
    full = vars(parser.parse_args())
    deaf = vars(parser.parse_args([]))
    global bfsamp
    global afsamp
    if type(args.sampling)==list and len(args.sampling) > 1:
        bfsamp = int(args.sampling[0])
        afsamp = int(args.sampling[1])
    else:
        bfsamp = int(args.sampling[0])
        afsamp = int(args.sampling[0])


    with open("nonWF_"+str(args.generations[0])+"_"+str(args.generations[1])+ "_n" + str(bfsamp) +"_n" + str(afsamp) + "_r" + str(args.seed) + ".args", "w") as f:
            f.write("SimOutbreakSelection v.1.0\n")
            f.write("Time: " + datetime.now().strftime("%d/%m/%Y %H:%M:%S") + "\n")
            f.write("Directory: " + str(os.getcwd()) + "\n")
            f.write("Options:\n")
            for key in full:
                    if full[key] != deaf[key]:
                            if type(full[key]) is bool:
                                    f.write("\t--" + str(key) + "\n")
                            else:
                                    f.write("\t--" + str(key) + " " + str(full[key]) + "\n")
    del full, deaf
    return parser.parse_args(args=None if argv else ['--help'])


# Read in sampling function
def gen_samples(gen, num_samples):
    time_inds = set(mutated.individuals_alive_at(gen))
    # Determine all samples that are in a population
    samples_of_pop = dict()
    for u in time_inds:
        pop = mutated.individual(u).population
        samples_of_pop.setdefault(pop, []).append(u)

    # Create new sample IDs
    sample_names = []
    sampled_nodes = []
    for pop in sorted(samples_of_pop):
        samples = sorted(np.random.choice(samples_of_pop[pop], num_samples, replace=False))
        for sam in samples:
            sampled_nodes.append(sam)
    actual_samples = np.asarray([mutated.individual(pep).nodes for pep in sampled_nodes]).flatten().tolist()
    return actual_samples, sampled_nodes


def get_gen_inds(target, num_samples):
    gen = slim_gens - target
    sample_inds, gen_inds = gen_samples(gen, num_samples)
    return sample_inds, gen_inds


def bottle_samples(target, status):
    gen = slim_gens - target
    # Take the samples of interest:
    # Extract individuals who survived the epidemic
    np.warnings.filterwarnings('ignore')
    file = open(putparents)
    f = file.readlines()
    bot_line = int((target-deme_gen)*2)
    f = f[bot_line]
    f = f.split('\t')
    for i in range(1, len(f)):
        f[i] = int(f[i])
    survivors = set(f[1:])
    bottle = mutated.individual_ages_at(gen)
    bottle_alive = bottle[bottle >= 0]
    print("There are these many ages in the bottleneck: ")
    print(len(bottle))
    print("There are these many that existed in the bottleneck:")
    print(len(bottle_alive))

    print("Start loop now.")
    num_alive = len(mutated.individuals_alive_at(gen)[bottle_alive == 0])
    print("There are these many alive: ")
    print(num_alive)
    allBottle = set(mutated.individuals_alive_at(gen)[bottle_alive == 0])
    print("Getting pedIds...")

    print("Making dictionary...")
    bottle_inds = {mutated.individual(i).metadata.pedigree_id:i for i in allBottle}
    print("Grabbing survivor IDs...")
    survivor_ids = [bottle_inds.get(k) for k in survivors]

    print("Grabbing survivor nodes...")
    survivor_nodes = np.asarray([mutated.individual(j).nodes for j in survivor_ids]).flatten().tolist()
    print("Grabbing victims using sets...")
    victims = set(bottle_inds.keys())-set(bottle_inds.keys()).intersection(survivors)
    print("Grabbing victim ids...")
    dead_ids = [bottle_inds.get(k) for k in victims]

    print("Grabbing victim nodes...")
    dead_nodes = np.asarray([mutated.individual(j).nodes for j in dead_ids]).flatten().tolist()
    print("Loop finished.")
    bottle_samples = dead_ids
    bottle_nodes = dead_nodes
    other_samples = survivor_ids
    other_nodes = survivor_nodes
    if status == 'dead':
        print("DEAD")
        print(len(dead_ids),len(dead_nodes))
    else:
        print("ALIVE")
    print(len(survivor_ids), len(survivor_nodes))
    return bottle_samples, bottle_nodes, other_samples, other_nodes

def mixed_cemetery(target, num_samples, status, ratio):
    bottle_sam, bottle_nodes, other_sam, other_nodes = bottle_samples(target, status)
    dead_mix = random.sample(bottle_sam, k=int(ratio*num_samples))
    survivor_mix = random.sample(other_sam, k=int(num_samples-(ratio*num_samples)))
    bottle_ids = dead_mix + survivor_mix
    print(bottle_ids)
    bottle_gen_inds = np.asarray([mutated.individual(p).nodes for p in bottle_ids]).flatten().tolist()
    return bottle_gen_inds, bottle_ids, dead_mix, survivor_mix

def get_bottle_inds(target, num_samples, status, alt_samples):
    bottle_sam, bottle_nodes, other_sam, other_nodes = bottle_samples(target, status)
    bottle_ids = random.sample(bottle_sam, k=num_samples)
    other_ids = random.sample(other_sam, k=alt_samples)
    bottle_gen_inds = np.asarray([mutated.individual(p).nodes for p in bottle_ids]).flatten().tolist()
    other_gen_inds = np.asarray([mutated.individual(p).nodes for p in other_ids]).flatten().tolist()
    return bottle_gen_inds, bottle_ids, other_gen_inds, other_ids

def merge_gen_gt(before, before_inds, after, after_inds):
    # Obtain simplified tree for samples (across generations) that will be compared
    merged = mutated.simplify(samples = before + after)
    merged_pos = [int(x) for x in merged.tables.sites.asdict()['position']]
    merged_inds = before_inds + after_inds

    # Obtains number of variants for length of array
    num_var = merged.num_sites
    count_of_allele_1 = np.zeros(num_var, dtype=int)
    count_of_allele_2 = np.zeros(num_var, dtype=int)
    #Fills array with genotypes
    i = 0
    for variant in merged.variants():
        tmp = variant.genotypes
        num_gtbf = bfsamp * 2
        count_of_allele_1[i] = int(np.sum(tmp[:num_gtbf]))
        num_gtaf = afsamp * 2
        count_of_allele_2[i] = int(np.sum(tmp[num_gtbf:]))
        i += 1
    return count_of_allele_1, count_of_allele_2, merged_pos, merged_inds, num_gtbf, num_gtaf


def matmul_r2(genotypes):
    norm = (genotypes - genotypes.mean(1, keepdims=True)) / genotypes.std(1, ddof=1, keepdims=True)
    # num_snps, num_indivs =
    r2 = (np.dot(norm, norm.T) / (norm.shape[1] - 1)) ** 2
    return (r2)


def get_LD_at_site(genotypes, site_index):
    # normalize the genotype matrix
    norm = np.zeros(genotypes.shape)
    idx = np.sum(genotypes, axis=1) != 0
    norm[idx] = (genotypes[idx] - genotypes[idx].mean(1, keepdims=True)) / genotypes[idx].std(1, ddof=1, keepdims=True)
    r2 = ((np.dot(norm[site_index], norm.T) / (norm.shape[1] - 1)) ** 2)
    return (r2)

def jonas(oldG):
    m, n = oldG.shape
    newG = np.empty((m, n//2, 2), dtype=np.int8)
    for j in range(m):
        c = 0
        for i in range(n//2):
            for d in range(2):
                newG[j,i,d] = oldG[j,c]
                c += 1
    return newG

def main(argv):
    # Parsing command-line arguments
    args = parse_args(argv)
    # Work directory
    if args.dir != None:
        os.chdir(args.dir)

    # Import heavy libraries
    import gzip
    import h5py
    import scipy as sp
    from scipy import interpolate
    import warnings
    import pandas as pd
    import msprime
    import pyslim
    import allel
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style('white')
    sns.set_style('ticks')
    sns.set_context('notebook')
    import tskit
    import collections
    import glob

    # Saving command-line arguments
    np.random.seed(args.seed)
    file_name = args.tree
    file_mut = args.muTree
    bottle_gen = int(args.bottle)
    selected = args.select
    before_gen = args.generations[0]
    after_gen = args.generations[1]
    chromo = args.chromosomes
    midchrom = args.midchrom
    dead = args.dead
    mix_cem = args.mixCem
    sampling = args.sampling
    vcf = args.vcf
    nosel = args.nosel
    realhet = args.realhet
    plots = args.plots
    winfst = args.winfst
    ihs = args.ihs
    gwas = args.gwas
    global putparents
    putparents = args.putparents
    global mutated
    global threads
    threads = args.threads
    mu = args.mu
    rec = args.rec
    eff_pop = args.ne
    global deme_gen
    deme_gen = args.demegen
    pmap = args.pmap
    recmaps = args.recmaps

    start = datetime.now()
    if type(args.sampling)==list and len(args.sampling) > 1:
        bfsamp = int(args.sampling[0])
        afsamp = int(args.sampling[1])
    else:
        bfsamp = int(args.sampling[0])
        afsamp = int(args.sampling[0])


    if os.path.exists(file_mut):
        print("\nReading recapitated tree...")
        mutated = pyslim.load(file_mut)
    else:
        print("\nRecapitating tree...")
        ts = pyslim.load(file_name)
        # Recapitate
        recap = ts.recapitate(recombination_rate=rec, Ne=eff_pop)

        # Mutate
        mutated = pyslim.SlimTreeSequence(msprime.mutate(recap, rate=mu, random_seed=62, keep=True))
        mutated.dump(file_mut)
        mutated = pyslim.load(file_mut)
        mu_pos = [int(x) for x in mutated.tables.sites.asdict()['position']]
        while mu_pos.count(selected) > 1:
            print(mu_pos.count(selected))
            mutated = pyslim.SlimTreeSequence(msprime.mutate(recap, rate=rec, keep=True))
            mutated.dump(file_mut)
            mutated = pyslim.load(file_mut)
            break
    end = datetime.now()
    time_taken = end - start
    print('\nTime: ',time_taken)

    # SLiM generations and tskit generations
    global slim_gens
    slim_gens = mutated.slim_generation
    sim_gens = slim_gens - deme_gen
    print("SlimGen", "tskitGen", "GenSize" )
    for gen, j in zip(range(sim_gens, -1, -1), range(0, sim_gens + 1)):
        mutated.individuals_alive_at(gen)
        print(j + deme_gen, gen, len(mutated.individuals_alive_at(gen)))
    if realhet==True:
        all=len(mutated.individuals_alive_at(slim_gens-after_gen))
        hetsamps,hetnodes=gen_samples(slim_gens-after_gen,all)
        realhetgt=mutated.simplify(hetsamps)
        simpsamps=realhetgt.individuals_alive_at(slim_gens-after_gen)
        lalala=jonas(realhetgt.genotype_matrix())
        realhetgt = allel.GenotypeArray(lalala)
        del lalala
        realhet_count=[realhetgt[:,i].count_het() for i in range(0,len(simpsamps),1)]
        np.savetxt('nonWF_%s_het.out' % (after_gen), np.c_[realhet_count], delimiter='\t', comments='',fmt=['%i'])
        del realhetgt



    status_inds=[]
    print(before_gen, dead, bottle_gen, after_gen)

    if before_gen == bottle_gen and after_gen == bottle_gen:
        print("Extracting the dead and survivors from bottleneck...")

        before, before_inds, after, after_inds = get_bottle_inds(bottle_gen, bfsamp, 'dead', afsamp)
        status = "Dead"
        for s in range(1, bfsamp + 1): status_inds.append(status)
        print(afsamp)
        status = "Survivors"
        for s in range(1, afsamp + 1): status_inds.append(status)

        print("You are comparing Dead and Survivors of the same generation.")

    else:
        # Before individuals
        if before_gen == bottle_gen:
            print("before = bottle")
            if dead == before_gen and mix_cem==None:
                print("Getting dead from first input generation...")
                before, before_inds, _, _ = get_bottle_inds(bottle_gen, bfsamp, 'dead', afsamp)
                status = "Dead"
                for s in range(1,bfsamp+1): status_inds.append(status)
                # Requested mixed cemetery of diseased and survivors
            elif dead == before_gen and mix_cem!=None:
                print("Getting dead from first input generation...")
                print("Requested mixed cemetery...")
                before, before_inds, dead_mix, survivor_mix = mixed_cemetery(bottle_gen, bfsamp, 'dead', mix_cem)
                status = "Dead"
                for s in range(1,bfsamp+1): status_inds.append(status)
            else:
                print("dead != before")
                _, _, before, before_inds = get_bottle_inds(bottle_gen, bfsamp, 'alive', afsamp)
                status = "Survivors"
                for s in range(1, bfsamp + 1): status_inds.append(status)
        else:
            print("Extracting individuals before bottleneck...")
            before, before_inds = get_gen_inds(before_gen, bfsamp)
            status = "Alive"
            for s in range(1, bfsamp + 1): status_inds.append(status)

        # After individuals
        if after_gen == bottle_gen:
            if dead == after_gen and mix_cem==None:
                print("Getting dead from second input generation...")
                after, after_inds, _, _ = get_bottle_inds(bottle_gen, afsamp, 'dead', bfsamp)
                status="Dead"
                for s in range(1, afsamp + 1): status_inds.append(status)
                # Requested mixed cemetery of diseased and survivors
            elif dead == after_gen and mix_cem!=None:
                print("Getting dead from second input generation...")
                print("Requested mixed cemetery...")
                after, after_inds, dead_mix, survivor_mix = mixed_cemetery(bottle_gen, afsamp, 'dead', mix_cem)
                status="Dead"
                for s in range(1, afsamp + 1): status_inds.append(status)

            else:
                other, other_inds, after, after_inds = get_bottle_inds(bottle_gen, afsamp, 'alive', bfsamp)
                del other
                del other_inds
                status = "Survivors"
                for s in range(1, afsamp + 1): status_inds.append(status)
        else:
            print("Extracting individuals after bottleneck...")
            after, after_inds = get_gen_inds(after_gen, afsamp)
            status = "Alive"
            for s in range(1, afsamp + 1): status_inds.append(status)

    # Allele counts, positions, individual IDs for both generations
    print("Merging variants for allele count...")
    print(len(before), len(before_inds),len(after),len(after_inds))
    start = datetime.now()
    count_of_allele_1, count_of_allele_2, merged_pos, merged_inds, num_gtbf, num_gtaf = merge_gen_gt(before,before_inds,after,after_inds)
    status_inds = np.array(status_inds)
    end = datetime.now()
    time_taken = end - start
    print('\nTime: ',time_taken)

    # LD Decay Estimation
    print("Estimating LD decay...")
    nac1 = count_of_allele_1 / num_gtbf
    nac2 = count_of_allele_2 / num_gtaf
    print("Beginning to simplify sub-sampled tree...")
    merged = mutated.simplify(samples = before + after)
    merged_pos = [int(x) for x in merged.tables.sites.asdict()['position']]
    positions = merged_pos
    if selected in merged_pos:
        site = merged_pos.index(selected)
    else:
        print("WARNING: Your selected variant is not in the list. Taking closest variant now...")
        selected = min(merged_pos, key=lambda x:abs(x-selected))
        site = merged_pos.index(selected)

    # BEFORE
    print("LD Before Processing...")
    np.seterr(divide='ignore', invalid='ignore')
    biG = merged.genotype_matrix()
    genotypes = biG[:, :num_gtbf]
    del biG
    r2_site = get_LD_at_site(genotypes=genotypes, site_index=site)
    gen1_r2 = r2_site

    # AFTER
    print("LD After Processing...")
    biG = merged.genotype_matrix()
    genotypes = biG[:, num_gtbf:]
    del biG
    r2_site = get_LD_at_site(genotypes=genotypes, site_index=site)
    gen2_r2 = r2_site

    # OUTPUT Fst WITH LD VALUES for usage in R
    ac1 = np.stack([num_gtbf - count_of_allele_1, count_of_allele_1], axis=1)
    ac2 = np.stack([num_gtaf - count_of_allele_2, count_of_allele_2], axis=1)
    pos = {}
    y = {}
    pos = merged_pos
    if winfst==True:
        blen = 100000
    ##ONLY FOR WINFST 
        y, windows, _ = allel.windowed_hudson_fst(pos, ac1, ac2, size=blen, step=int(blen/2))
        x = windows[:,0]
        gen2_r2=[0]*len(y)
    else:
        blen = 1
        windows = allel.moving_statistic(pos, statistic=lambda v: [v[0], v[-1]], size=blen)
        x = np.asarray(windows).mean(axis=1)
        y, _, _ = allel.windowed_hudson_fst(pos, ac1, ac2, windows=windows)
    
    x2 = [int(l)-midchrom if int(l) >= midchrom else int(l) for l in x]
    chr_fst = [chromo[1] if int(l) >= midchrom else chromo[0] for l in x]
    snp_fst = ["chr%s_%s" % (chromo[1],int(l)-midchrom) if int(l) >= midchrom else "chr%s_%s" % (chromo[0],int(l))for l in x]
    
    if winfst==True:
        np.savetxt('nonWF_%s_%s_%s_%s_FST.out.gz' % (bfsamp,afsamp,before_gen,after_gen), np.c_[chr_fst,x2,y,gen2_r2,snp_fst], header="CHR\tPOS\tFST\tR2\tSNP", delimiter='\t', comments='',fmt=['%s','%s','%3s','%s','%s'])
    else:
        np.savetxt('nonWF_%s_%s_%s_%s_FST.out.gz' % (bfsamp,afsamp,before_gen,after_gen), np.c_[chr_fst,x2,y,count_of_allele_1,count_of_allele_2,gen2_r2,snp_fst], header="CHR\tPOS\tFST\tAC1\tAC2\tR2\tSNP", delimiter='\t', comments='',fmt=['%s','%s','%3s','%s','%s','%3s','%s'])
    
    del x2
    del chr_fst
    del snp_fst


    #PLOTTING CODE with LD colouring
    if plots==True:
        mpl.use('Agg')
        # Windowed FST - Full - LD color
        ac1 = np.stack([num_gtbf - count_of_allele_1, count_of_allele_1], axis=1)
        ac2 = np.stack([num_gtaf - count_of_allele_2, count_of_allele_2], axis=1)
        chrom = '%s & %s' % (chromo[0],chromo[1])
        fig = plt.figure(figsize=(20, 10))
        i = 1
        pos = {}
        y = {}
        pos = merged_pos
        blen = 1
        windows = allel.moving_statistic(pos, statistic=lambda v: [v[0], v[-1]], size=blen)
        x = np.asarray(windows).mean(axis=1)
        y, _, _ = allel.windowed_hudson_fst(pos, ac1, ac2, windows=windows)

        site = np.where(x == selected)[0]
        site = site.tolist()[0]

        # plot
        #dingdong = gen2_r2 > 0.1
        ax1 = plt.scatter(x, y, c=gen2_r2, cmap=plt.cm.jet, alpha=0.5);
        plt.scatter(x[site], y[site], c='deeppink', marker="*", s=100);
        plt.vlines([midchrom], 0, 1.0, colors='r', linestyles='dashed')
        plt.colorbar(ax1, label=r'$r^{2}$')
        plt.ylim(0, 1.0);
        plt.xlim(0, max(pos))
        plt.ylabel("$F_{ST}$", fontsize=20)
        plt.xlabel('Chromosome %s position\n (%s bp bins)' % (chrom, blen), fontsize=20);
        plt.xticks(fontsize=18);
        plt.yticks(fontsize=18);
        # LD coloured FST
        plt.savefig('FST_LD_%s_%s_%s_%s.png' % (before_gen, after_gen, bfsamp, afsamp));

        # LD coloured JSFS
        fig = plt.figure(figsize=(20, 10))
        site = merged_pos.index(selected)
        dingdong = gen2_r2 > 0.1
        plt.scatter((nac1)[~dingdong], (nac2)[~dingdong], color="darkblue", zorder=1, alpha=0.5);
        plt.scatter((nac1)[dingdong], (nac2)[dingdong], c=gen2_r2[dingdong],
                    cmap=plt.cm.get_cmap('jet'), zorder=2, alpha=0.6);

        plt.clim(0, 1);
        plt.colorbar(label=r'$r^{2}$');

        selected_dingdong = np.abs(np.array(merged_pos) - selected) == 0
        plt.scatter((nac1)[selected_dingdong], (nac2)[selected_dingdong], c='deeppink',
                    marker="*",
                    s=120, zorder=3);

        plt.ylabel('Pop. of Generation After: %s (%s)' % (after_gen,status_inds[-1]), fontsize=20);
        plt.xlabel('Pop. of Generation Before: %s (%s)' % (before_gen, status_inds[0]), fontsize=20);
        plt.xticks(fontsize=18);
        plt.yticks(fontsize=18);

        plt.savefig('SFS_LD_%s_%s_%s_%s.png' % (before_gen, after_gen, bfsamp, afsamp))

        # Regular SFS
        before_sub = count_of_allele_1[count_of_allele_1>0]
        after_sub = count_of_allele_2[count_of_allele_2>0]
        before_sub = before_sub[before_sub < num_gtbf]
        after_sub = after_sub[after_sub < num_gtaf]

        fig, axs = plt.subplots(1, 2, sharey=False, tight_layout=False, figsize=(20, 10))

        if num_gtbf <= 100:
            axs[0].bar(np.arange(1,num_gtbf),np.bincount(before_sub, minlength=num_gtbf)[1:]/len(before_sub[before_sub > 0]), label='%s (%s)' % (before_gen, status_inds[0]), width=0.5)
            axs[0].bar(np.arange(1,num_gtaf)+0.5,np.bincount(after_sub, minlength=num_gtaf)[1:]/len(after_sub[after_sub > 0]), label='%s (%s)' % (after_gen,status_inds[-1]), width=0.5)
        else:
            axs[0].bar(np.arange(1,101),np.bincount(before_sub, minlength=num_gtbf)[1:101]/len(before_sub[before_sub > 0]), label='%s (%s)' % (before_gen, status_inds[0]), width=0.5)
            axs[0].bar(np.arange(1,101)+0.5,np.bincount(after_sub, minlength=num_gtaf)[1:101]/len(after_sub[after_sub > 0]), label='%s (%s)' % (after_gen,status_inds[-1]), width=0.5)
        axs[0].set_ylabel('site frequency / total # of sites');
        axs[0].set_xlabel('allele frequency class - no 0s');
        axs[0].legend();

        axs[1].plot(np.arange(1,num_gtbf),np.bincount(before_sub, minlength=num_gtbf)[1:]/len(before_sub[before_sub>0]), marker='o',linestyle='',label='%s (%s)' % (before_gen, status_inds[0]))
        axs[1].plot(np.arange(1,num_gtaf),np.bincount(after_sub, minlength=num_gtaf)[1:]/len(after_sub[after_sub>0]), marker='*',linestyle='',label='%s (%s)' % (after_gen,status_inds[-1]))

        axs[1].set_ylabel('site frequency / total # of sites');
        axs[1].set_xlabel('allele frequency class - no 0s - full spectra');
        axs[1].legend();
        plt.savefig('SFS_%s_%s_%s_%s.png' % (before_gen, after_gen, bfsamp, afsamp))

    # OUTPUT VCF for GWAS or iHS or nSL
    dirname=os.path.dirname(os.path.abspath(__file__)) + "/scripts/"
    gt_vcf_gwas=dirname + "gt_vcf_gwas.sh"
    if type(sampling)==list and len(sampling) > 1:
        sampling=str('%s_%s') % (bfsamp, afsamp)
    else:
        sampling=int(sampling[0])
    if vcf==True:
        if gwas==True and dead!=None and ihs==True:
            # For GWAS VCF
            np.savetxt('nonWF_%s_gwas.inds' % sampling, np.c_[merged_inds, status_inds], delimiter='\t', fmt='%s')
            gwas_name = str('nonWF_%s_gwas' % sampling)
            chrom = int(chromo[0])
            simplified_inds = []
            for q in merged.individuals_alive_at(slim_gens-before_gen):
                ind = merged.individual(q)
                if merged.node(ind.nodes[0]).is_sample():
                    simplified_inds.append(q)
                    assert merged.node(ind.nodes[1]).is_sample()
            print(len(simplified_inds))
            if before_gen!=after_gen or mix_cem!=None:
                for q in merged.individuals_alive_at(slim_gens-after_gen):
                    ind = merged.individual(q)
                    if merged.node(ind.nodes[0]).is_sample():
                        simplified_inds.append(q)
                        assert merged.node(ind.nodes[1]).is_sample()

            merged_inds.sort()
            merged_inds = [str(x) for x in merged_inds]

            print(len(simplified_inds), len(merged_inds))
            print(simplified_inds)
            print(merged_inds)
            with gzip.open('%s.vcf.gz' % (gwas_name), "wt") as vcf_file:
                merged.write_vcf(vcf_file,contig_id="chr"+str(chrom),individuals=simplified_inds,individual_names=merged_inds)
            process = subprocess.Popen(['bash', gt_vcf_gwas, gwas_name, str(chrom), str(midchrom), str(before_gen), str(after_gen)])
            process.wait()

            # For IHS VCFs
            merged = mutated.simplify(samples=before)
            merged_pos = [int(x) for x in merged.tables.sites.asdict()['position']]
            merged_inds = before_inds

            np.savetxt('nonWF_%s_%s_ihs.inds' % (bfsamp, before_gen), np.transpose(merged_inds), delimiter='\t', fmt='%s')

            print("Creating iHS VCFs now...")
            ihs_name = str('nonWF_%s_%s_ihs' % (bfsamp, before_gen))
            chrom = int(chromo[0])
            simplified_inds = []
            for q in merged.individuals_alive_at(slim_gens-before_gen):
                ind = merged.individual(q)
                if merged.node(ind.nodes[0]).is_sample():
                    simplified_inds.append(q)
                    assert merged.node(ind.nodes[1]).is_sample()
            merged_inds = [str(x) for x in merged_inds]
            with gzip.open('%s.vcf.gz' % (ihs_name), "wt") as vcf_file:
                merged.write_vcf(vcf_file,contig_id="chr"+str(chrom),individuals=simplified_inds,individual_names=merged_inds)
            process = subprocess.Popen(['bash', gt_vcf_gwas, ihs_name, str(chrom), str(midchrom), str(before_gen), str(after_gen)])
            process.wait()

            merged = mutated.simplify(samples=after, reduce_to_site_topology=True)
            merged_pos = [int(x) for x in merged.tables.sites.asdict()['position']]
            merged_inds = after_inds

            np.savetxt('nonWF_%s_%s_ihs.inds' % (afsamp, after_gen), np.transpose(merged_inds), delimiter='\t', fmt='%s')
            ihs_name = str('nonWF_%s_%s_ihs' % (afsamp, after_gen))
            chrom = int(chromo[0])
            simplified_inds = []
            for q in merged.individuals_alive_at(slim_gens-after_gen):
                ind = merged.individual(q)
                if merged.node(ind.nodes[0]).is_sample():
                    simplified_inds.append(q)
                    assert merged.node(ind.nodes[1]).is_sample()
            merged_inds = [str(x) for x in merged_inds]
            print(len(simplified_inds), len(merged_inds))
            with gzip.open('%s.vcf.gz' % (ihs_name), "wt") as vcf_file:
                merged.write_vcf(vcf_file,contig_id="chr"+str(chrom),individuals=simplified_inds,individual_names=merged_inds)
            process = subprocess.Popen(['bash', gt_vcf_gwas, ihs_name, str(chrom), str(midchrom), str(before_gen), str(after_gen)])
            process.wait()

            all_gens = [before_gen,after_gen]
            if status_inds[0]=="Dead" and status_inds[-1]=="Survivors":
                stat=["d","a"]
                for g,h in zip(all_gens,stat):
                    if g==before_gen:
                        sampling=bfsamp
                    else:
                        sampling=afsamp
                    ihs_name = str('nonWF_%s_%s_ihs' % (sampling, g))
                    chrom = int(chromo[0])
            else:
                h=None
                for g in all_gens:
                    if g==before_gen:
                        sampling=bfsamp
                    else:
                        sampling=afsamp

                    ihs_name = str('nonWF_%s_%s_ihs' % (sampling, g))
                    chrom = int(chromo[0])

                if pmap!=True:
                    for n in chromo:
                        chrom = int(n)
                        if h!=None:
                            callset = allel.read_vcf('%s_chr%s_%s.vcf.gz' % (ihs_name, chrom, h), region='chr%s' % (chrom))
                        else:
                            callset = allel.read_vcf('%s_chr%s.vcf.gz' % (ihs_name, chrom), region='chr%s' % (chrom))
                        vcf_pos = callset['variants/POS'].tolist()
                        #print(len(vcf_pos))
                        map_rate = [x.split('\t')[2] for x in open(recmaps + 'slim_chr%s_recode.map' % (chrom)).readlines()]
                        map_rate = [float(i) for i in map_rate]
                        map_cm = [x.split('\t')[3] for x in open(recmaps + 'slim_chr%s_recode.map' % (chrom)).readlines()]
                        map_cm = [float(i) for i in map_cm]

                        map_pos = [x.split('\t')[1] for x in open(recmaps + 'slim_chr%s_recode.map' % (chrom)).readlines()]
                        map_pos = [int(i) for i in map_pos]

                        extrap_cm = interpolate.interp1d(map_pos, map_cm)
                        fin_chr = ['chr%s' % chrom] * len(vcf_pos)
                        fin_snp = []
                        vcf_map = []

                        for val in vcf_pos:
                            if val < map_pos[1]:
                                # CALCULATE USING RATE
                                new_cm = float(map_cm[1] - (map_pos[1] - val) * (map_rate[1] / 1000000))
                                vcf_map.append(format(new_cm, 'f'))
                            elif val > map_pos[-1]:
                                new_cm = float(map_cm[-1] + (val - map_pos[-1]) * (map_rate[-1] / 1000000))
                                vcf_map.append(format(new_cm, 'f'))
                            else:
                                new_cm = extrap_cm(val)
                                vcf_map.append(format(new_cm, 'f'))
                            fin_snp.append('chr%s_%s' % (chrom,val))
                        if h!=None:
                            np.savetxt(str('nonWF_%s_%s_chr%s_%s.map' % (sampling, g, chrom, h)), np.c_[fin_chr, fin_snp, vcf_map, vcf_pos], delimiter='\t', fmt=['%s','%s','%.8s','%s'])
                        else:
                            np.savetxt(str('nonWF_%s_%s_chr%s.map' % (sampling, g, chrom)), np.c_[fin_chr, fin_snp, vcf_map, vcf_pos], delimiter='\t', fmt=['%s','%s','%.8s','%s'])

        elif gwas==True and dead==None:
            print("You don't have dead people.")

        elif gwas==True and dead!=None:
            # For GWAS VCF
            np.savetxt('nonWF_%s_gwas.inds' % sampling, np.c_[merged_inds, status_inds], delimiter='\t', fmt='%s')
            gwas_name = str('nonWF_%s_gwas' % sampling)
            chrom = int(chromo[0])
            simplified_inds = []
            for q in merged.individuals_alive_at(slim_gens-before_gen):
                ind = merged.individual(q)
                if merged.node(ind.nodes[0]).is_sample():
                    simplified_inds.append(q)
                    assert merged.node(ind.nodes[1]).is_sample()
            #print(len(simplified_inds))
            if before_gen!=after_gen or mix_cem!=None:
                for q in merged.individuals_alive_at(slim_gens-after_gen):
                    ind = merged.individual(q)
                    if merged.node(ind.nodes[0]).is_sample():
                        simplified_inds.append(q)
                        assert merged.node(ind.nodes[1]).is_sample()
            merged_inds.sort()
            merged_inds = [str(x) for x in merged_inds]
            print(len(simplified_inds),len(merged_inds))
            with gzip.open('%s.vcf.gz' % (gwas_name), "wt") as vcf_file:
                merged.write_vcf(vcf_file,contig_id="chr"+str(chrom),individuals=simplified_inds,individual_names=merged_inds)
            process = subprocess.Popen(['bash', gt_vcf_gwas, gwas_name, str(chrom), str(midchrom), str(before_gen), str(after_gen)])
            process.wait()


        elif ihs==True:
            # For IHS VCFs
            merged = mutated.simplify(samples=before)
            merged_pos = [int(x) for x in merged.tables.sites.asdict()['position']]
            merged_inds = before_inds

            np.savetxt('nonWF_%s_%s_ihs.inds' % (bfsamp, before_gen), np.transpose(merged_inds), delimiter='\t', fmt='%s')

            print("Creating iHS VCFs now...")
            ihs_name = str('nonWF_%s_%s_ihs' % (bfsamp, before_gen))
            chrom = int(chromo[0])
            simplified_inds = []
            for q in merged.individuals_alive_at(slim_gens-before_gen):
                ind = merged.individual(q)
                if merged.node(ind.nodes[0]).is_sample():
                    simplified_inds.append(q)
                    assert merged.node(ind.nodes[1]).is_sample()
            merged_inds = [str(x) for x in merged_inds]
            with gzip.open('%s.vcf.gz' % (ihs_name), "wt") as vcf_file:
                merged.write_vcf(vcf_file,contig_id="chr"+str(chrom),individuals=simplified_inds,individual_names=merged_inds)
            process = subprocess.Popen(['bash', gt_vcf_gwas, ihs_name, str(chrom), str(midchrom), str(before_gen), str(after_gen)])
            process.wait()

            merged = mutated.simplify(samples=after, reduce_to_site_topology=True)
            merged_pos = [int(x) for x in merged.tables.sites.asdict()['position']]
            merged_inds = after_inds

            np.savetxt('nonWF_%s_%s_ihs.inds' % (afsamp, after_gen), np.transpose(merged_inds), delimiter='\t', fmt='%s')
            ihs_name = str('nonWF_%s_%s_ihs' % (afsamp, after_gen))
            chrom = int(chromo[0])
            simplified_inds = []
            for q in merged.individuals_alive_at(slim_gens-after_gen):
                ind = merged.individual(q)
                if merged.node(ind.nodes[0]).is_sample():
                    simplified_inds.append(q)
                    assert merged.node(ind.nodes[1]).is_sample()
            merged_inds = [str(x) for x in merged_inds]
            print(len(simplified_inds), len(merged_inds))
            with gzip.open('%s.vcf.gz' % (ihs_name), "wt") as vcf_file:
                merged.write_vcf(vcf_file,contig_id="chr"+str(chrom),individuals=simplified_inds,individual_names=merged_inds)
            process = subprocess.Popen(['bash', gt_vcf_gwas, ihs_name, str(chrom), str(midchrom), str(before_gen), str(after_gen)])
            process.wait()

            
            all_gens = [before_gen,after_gen]
            for g in all_gens:
                if g==before_gen:
                    sampling=bfsamp
                else:
                    sampling=afsamp

                ihs_name = str('nonWF_%s_%s_ihs' % (sampling, g))
                chrom = int(chromo[0])

                if pmap!=True:
                    for n in chromo:
                        chrom = int(n)
                        callset = allel.read_vcf('%s_chr%s.vcf.gz' % (ihs_name, chrom), region='chr%s' % (chrom))
                        vcf_pos = callset['variants/POS'].tolist()
                        #print(len(vcf_pos))
                        map_rate = [x.split('\t')[2] for x in open(recmaps + 'slim_chr%s_recode.map' % (chrom)).readlines()]
                        map_rate = [float(i) for i in map_rate]
                        map_cm = [x.split('\t')[3] for x in open(recmaps + 'slim_chr%s_recode.map' % (chrom)).readlines()]
                        map_cm = [float(i) for i in map_cm]

                        map_pos = [x.split('\t')[1] for x in open(recmaps + 'slim_chr%s_recode.map' % (chrom)).readlines()]
                        map_pos = [int(i) for i in map_pos]

                        extrap_cm = interpolate.interp1d(map_pos, map_cm)
                        fin_chr = ['chr%s' % chrom] * len(vcf_pos)
                        fin_snp = []
                        vcf_map = []

                        for val in vcf_pos:
                            if val < map_pos[1]:
                                # CALCULATE USING RATE
                                new_cm = float(map_cm[1] - (map_pos[1] - val) * (map_rate[1] / 1000000))
                                vcf_map.append(format(new_cm, 'f'))
                            elif val > map_pos[-1]:
                                new_cm = float(map_cm[-1] + (val - map_pos[-1]) * (map_rate[-1] / 1000000))
                                vcf_map.append(format(new_cm, 'f'))
                            else:
                                new_cm = extrap_cm(val)
                                vcf_map.append(format(new_cm, 'f'))
                            fin_snp.append('chr%s_%s' % (chrom,val))
                        np.savetxt(str('nonWF_%s_%s_chr%s.map' % (sampling, g, chrom)), np.c_[fin_chr, fin_snp, vcf_map, vcf_pos], delimiter='\t', fmt=['%s','%s','%.8s','%s'])


    else:
        return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
