import sys
import csv
import math
import random as rd
import numpy as np

lower = int(sys.argv[1])
upper = int(sys.argv[2]) + 1
T_skip = int(sys.argv[3])

max_time = 1000000
loop_time = int(max_time / T_skip)
tree_name = "TD_"
space_name = "SD_"
X = 1
Y = 40000
N_ratio = 0.1
head_length = 5


# Functions
def get_dir(file_name, sub_dir):
    """
    This function will return the correct path to the relevant data files for the LTEE data.
    This file is designed to work with the processed barcode counts from LTEE fitness assays. For other users,
    please format the data as seen in the regular documentation.

    Parameters:
        file_name(str) - name of the data file including file type
        sub_dir(str) - name of the folder in which the desired file resides

    Returns:
        data_path(str) - the exact path to where the named file resides
    """
    start = os.path.abspath(__file__)
    start_dir = os.path.dirname(start)
    data_path = os.path.join(start_dir, sub_dir, file_name)

    return data_path


def Generate_DNA(tree, sample):
    sites = 0
    for i in sample:
        if i > sites:
            sites = i

    chromosomes = {}
    for s in sample:
        gene = np.zeros(sites + 1)
        branch = tree[s]
        if branch[2] > 0:
            if branch[2] == 1:
                gene[s] = 1
            else:
                snps = [s]
                for b in range(0, branch[2] - 1):
                    parent = tree[snps[b]]
                    snps.append(parent[1])
                for b in range(0, len(snps)):
                    gene[snps[b]] = 1
        chromosomes[s] = gene

    return (chromosomes, sites)


def dNdS(tree, sample, chromosomes):
    n_total = 0
    s_total = 0
    for s in sample:
        gene = chromosomes[s]
        pop = sample[s]
        for m in range(0, len(gene)):
            if gene[m] == 1:
                branch = tree[m]
                if branch[3] == 1:
                    n_total += pop
                else:
                    s_total += pop
    return (n_total, s_total)


def SFS(sample, chromosomes, sites, n):
    site_count = np.zeros(sites + 1)
    for s in sample:
        gene = chromosomes[s]
        pop = sample[s]
        for m in range(1, sites + 1):
            if gene[m] == 1:
                site_count[m] += pop
    print("sites done")

    sfs = np.zeros(n - 1)
    fixed_list = []
    for seg in range(0, sites + 1):
        if site_count[seg] > 0:
            if int(site_count[seg]) == n:
                fixed_list.append(seg)
            else:
                ind = int(site_count[seg] - 1)
                sfs[ind] += 1

    return (sfs, fixed_list)


def Theta_Pi(sfs, n):
    denom = float(n * (n - 1)) / 2.0
    numer = 0.0
    for i in range(1, n):
        numer += i * (n - i) * float(sfs[i - 1])
    theta = numer / denom

    return theta


def TajimasD(sfs, n):
    a1 = 0.0
    a2 = 0.0
    for i in range(1, n):
        a1 += 1.0 / float(i)
        a2 += + 1.0 / float(i * i)
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)

    S = math.fsum(sfs)

    if S > 1:
        denom = math.sqrt(e1 * S + e2 * S * (S - 1))
        D = (Theta_Pi(sfs, n) - S / a1) / denom
    else:
        D = 100.0

    return D


def Stats(tree, sample, n, local_name, local_num, time):
    DNA = Generate_DNA(tree, sample)
    chromosomes = DNA[0]
    sites = DNA[1]
    """
    dna_file = open(str(local_num) + '_' + local_name + '_dna.csv', "w")
    for s in sample:
        line_holder = str(s) + ',' + str(sample[s])
        gene = chromosomes[s]
        for m in range(0, len(gene)):
            line_holder += ',' + str(gene[m])
        dna_file.write(line_holder)
        dna_file.write("\n")
    dna_file.close()
    """
    print(str(local_num) + " " + local_name + " dna done")

    spectrum = SFS(sample, chromosomes, sites, n)
    sfs = spectrum[0]
    fixed = spectrum[1]
    sfs_file = open(str(local_num) + '_' + local_name + str(time) +'_sfs.csv', "w")
    sfs_file.write(str(n))
    sfs_file.write("\n")
    for ind in range(0, len(sfs)):
        line_holder = str(ind + 1) + ',' + str(sfs[ind])
        sfs_file.write(line_holder)
        sfs_file.write("\n")
    sfs_file.close()
    if local_name == "omni":
        fixed_file = open("FD_" + str(local_num) + str(time) + '.csv', "w")
        fixed_file.write(str(len(fixed)))
        fixed_file.write("\n")
        for ind in range(0, len(fixed)):
            fixed_file.write(str(fixed[ind]))
            fixed_file.write("\n")
        fixed_file.close()
    print(str(local_num) + " " + local_name + " sfs done")

    counts = dNdS(tree, sample, chromosomes)
    ratio = (1 - N_ratio) / N_ratio
    if counts[1] == 0:
        ns = -1
    else:
        ns = ratio * counts[0] / counts[1]
    D = TajimasD(sfs, n)

    return (D, ns, len(fixed))


# Code
for r in range(lower, upper):
    for gen in range(1, loop_time + 1):
        T = gen * T_skip
        T_D = [[" ", "Tajimas D", "Tajimas D", "dN/dS", "dN/dS", "Mutations Fixed", "Mutations Fixed"], ["Simulation #", "Whole Sampling", "Random Sampling", "Whole Sampling", "Random Sampling", "Whole Sampling", "Random Sampling"]]
        tree_file = tree_name + str(r) + ".csv"
        space_file = space_name + str(r) + "_" + str(T) + ".csv"
        with open(tree_file, 'rt') as f:
            csv_reader = csv.reader(f)
            header = next(csv_reader)
            clone_num = int(header[head_length + 1])
            tree_data = {}
            for c in range(0, clone_num):
                clones = next(csv_reader)
                tkey = int(clones[0])
                tval = [int(clones[0]), int(clones[1]), int(clones[2]), int(clones[3])]
                tree_data[tkey] = tval
        f.close()

        space_data = np.zeros((Y, X))
        with open(space_file, 'rt') as g:
            csv_reader = csv.reader(g)
            for y in range(0, Y):
                row = next(csv_reader)
                for x in range(0, X):
                    space_data[y][x] = int(row[x])
        g.close()

        pool = []
        ss = 100
        for y in range(0, Y):
            for x in range(0, X):
                if space_data[y][x] >= 0:
                    pool.append(int(space_data[y][x]))

        omni_sample = {}
        # regn_sample = {}
        for p in range(0, len(pool)):
            if pool[p] in tree_data:
                if pool[p] not in omni_sample:
                    omni_sample[pool[p]] = 1
                else:
                    omni_sample[pool[p]] = omni_sample[pool[p]] + 1

        rd.shuffle(pool)
        rand_sample = {}
        for p in range(0, ss):
            if pool[p] in tree_data:
                if pool[p] not in rand_sample:
                    rand_sample[pool[p]] = 1
                else:
                    rand_sample[pool[p]] = rand_sample[pool[p]] + 1

        omni_stats = Stats(tree_data, omni_sample, len(pool), "omni", r, T)
        print(str(r) + " omni stats done")
        rand_stats = Stats(tree_data, rand_sample, ss, "rand", r, T)
        print(str(r) + " rand stats done")
        temp = open("Stat_" + str(r) + "_" + str(T) + ".csv", "w")
        res = [r, omni_stats[0], rand_stats[0], omni_stats[1], rand_stats[1], omni_stats[2], rand_stats[2]]
        line_holder = str(res[0])
        for k in range(1, len(res)):
            line_holder += ',' + str(res[k])
        temp.write(line_holder)
        temp.write("\n")
        temp.close()
        T_D.append([r, omni_stats[0], rand_stats[0], omni_stats[1], rand_stats[1], omni_stats[2], rand_stats[2]])
