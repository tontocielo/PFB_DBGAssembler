#!/usr/bin/env python3

import matplotlib.pyplot as plt

file = open('kmer_contig.txt', 'r')

x_axis = []
y_axis = []
for line in file:
    line = line.rstrip()
    line = line.split('\t')
    x_axis.append(int(line[0]))
    y_axis.append(int(line[1]))

print(x_axis,y_axis)

plt.scatter(x_axis, y_axis)
plt.show()
