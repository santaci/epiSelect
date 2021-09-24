import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-in1")
parser.add_argument("-in2")
parser.add_argument("-out")
args = parser.parse_args()

file1 = pd.read_csv(args.in1, dtype=int, header=None, usecols=[1], sep="\t").values.squeeze()
file2 = pd.read_csv(args.in2, dtype=int, header=None, names=["gen", "id", "geno"], sep="\t")
c = np.zeros(file2.shape[0], dtype=int)

vals2 = file2["id"].values
n = vals2.shape[0]

for i in range(n):
	c[i] = np.sum(file1 == vals2[i])

file2['counts'] = c
file2.to_csv(args.out, sep="\t", header=False, index=False)