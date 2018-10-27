#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

filename = sys.argv[1]
kmer_length = int(sys.argv[2])

def kmer_count(filename, kmer_length):
   FILE = open(filename, 'r')
   if kmer_length % 2 == 0: raise ValueError('Kmer length cannot be an even number')
   kmers = {}
   for line in FILE:
    header = line
    seq = next(FILE)
    seq = seq.rstrip()
    plus = next(FILE)
    quality = next(FILE)
    #sliding window
    for i in range(0, len(seq) - kmer_length + 1):
       kmer = seq[i:i + kmer_length]
       if kmer in kmers:  kmers[kmer] += 1
       else:  kmers[kmer] = 1
   return(kmers)

kmers = kmer_count(filename, kmer_length)

kmers_sorted = sorted(kmers.items(), key = lambda x: x[1], reverse=True)
print(kmers_sorted)
plt.hist(kmers.values())
plt.plot()

plt.savefig('graph1.png')


#call as
#print(kmer_count(filename, kmer_length))















