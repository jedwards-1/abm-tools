import sys
import csv
import numpy as np

lower_bound = int(sys.argv[1])
upper_bound = int(sys.argv[2])
stat_name = "Stat_"
T_D = [[" ", "Tajimas D", "Tajimas D", "dN/dS", "dN/dS", "Mutations Fixed", "Mutations Fixed"], ["Simulation #", "Whole Sampling", "Random Sampling", "Whole Sampling", "Random Sampling", "Whole Sampling", "Random Sampling"]]
for r in range(lower_bound, upper_bound + 1):
    stat_file = stat_name + str(r) + ".csv"
    with open(stat_file, 'rt') as f:
        csv_reader = csv.reader(f)
        read1 = next(csv_reader)
        data1 = []
        for i in range(0, len(read1)):
            data1.append(float(read1[i]))
        T_D.append(data1)
    f.close()

averages = ["Average:"]
deviations = ["Standard Deviation:"]
for q in range(1, len(T_D[2])):
    var = np.zeros(len(T_D) - 2)
    for a in range(2, len(T_D)):
        var[a - 2] = T_D[a][q]
    averages.append(np.mean(var))
    deviations.append(np.std(var))
T_D.append(averages)
T_D.append(deviations)

stats = open('Stats_Summary_' + str(lower_bound) + '-' + str(upper_bound) + '.csv', "w")
for u in range(0, len(T_D)):
    line_holder = str(T_D[u][0])
    for l in range(1, len(T_D[u])):
        line_holder += ',' + str(T_D[u][l])
    stats.write(line_holder)
    stats.write("\n")
stats.close()
