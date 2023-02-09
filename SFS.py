import csv
import numpy as np
import matplotlib.pyplot as plt
import os

dim = "0DN"
low = 1
high = 50
sample = "rand"
t_skip = 10000
genome = 200


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


def plot_sfs(sfs, freq, dimension, filename):
    plt.plot(freq, sfs, "k.")
    neutral = np.zeros(len(sfs))
    for p in range(0, len(sfs)):
        neutral[p] = 0.15 / (p + 1)
    plt.plot(freq, neutral, "b--")
    plt.xscale("log")
    plt.ylim([0, 0.4])
    plt.title(dimension + " Site Frequency Spectrum")
    plt.xlabel('frequency of individuals')
    plt.ylabel('proportion of the mutations')
    out_file = get_dir(filename, dimension + "_sfs")
    plt.savefig(out_file)
    plt.close()


def get_sfs(low_end, high_end, sample_type, size, gen, dimension, plot=True):
    sfs = np.zeros(size)
    name = "_" + sample_type + str(gen) + "_sfs.csv"
    for r in range(low_end, high_end + 1):
        data_file = str(r) + name
        fold = dimension + "_"
        if r < 26:
            fold += "1-25"
        else:
            fold += "26-50"
        list_file = get_dir(data_file, fold)
        with open(list_file, 'rt') as g:
            csv_reader = csv.reader(g)
            header = next(csv_reader)
            fixed_num = int(header[0])
            for n in range(1, fixed_num):
                counts = next(csv_reader)
                xval = int(counts[0])
                yval = float(counts[1])
                sfs[xval - 1] += int(yval)
        g.close()

    if plot:
        norm_sfs = sfs / sfs.sum()
        freq = np.zeros(size)
        for i in range(0, size):
            freq[i] = (i + 1) / float(size)
        plot_sfs(norm_sfs, freq, dimension, dimension + "_" + sample_type + str(gen) + "_sfs.png")

    out_name = "SFS_" + dimension + sample_type + str(gen) + "_" + str(low_end) + "-" + str(high_end) + ".csv"
    out_file = get_dir(out_name, dimension + "_sfs")
    fix_file = open(out_file, "w")
    fix_file.write(str(len(sfs)))
    fix_file.write("\n")
    for i in range(0, len(sfs)):
        line_holder = str(i + 1) + "," + str(sfs[i])
        fix_file.write(line_holder)
        fix_file.write("\n")
    fix_file.close()


loop = int(1000000 / t_skip)
for time in range(1, loop + 1):
    point = time * t_skip
    print(point)
    get_sfs(low, high, sample, genome, point, dim)
