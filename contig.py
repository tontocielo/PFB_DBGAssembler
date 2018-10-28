#!/usr/bin/env python3

import sys
import random

filename = sys.argv[1]
kmer_length = int(sys.argv[2])

def rev_comp(x):
  output =  x.lower().replace('a', 'T').replace('t', 'A').replace('g','C').replace('c','G')[::-1]
  return output

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
       if kmer in kmers:  
          kmers[kmer] += 1
          kmers[rev_comp(kmer)] += 1
          
       else: 
          kmers[kmer] = 1
          kmers[rev_comp(kmer)] = 1
   return(kmers)

kmers_dict = kmer_count(filename, kmer_length)
print(len(kmers_dict), 'kmers gathered in the dictionary')

# remove low frequency kmers
kmers_dict = {key:val for key, val in kmers_dict.items() if val >= 6}
   
 
print(len(kmers_dict), 'kmers left in the dictionary after cut-off')

kmers_used = set()
alphabet = ['A', 'G', 'C', 'T']
contigs = []


for kmer in kmers_dict:
     kmers_used.add(kmer)
     kmers_used.add(rev_comp(kmer))

     direction = 'right'
     contig = kmer

     # extend to the right           
     while direction == 'right':           
         prefix = contig[-kmer_length + 1:]
         extensions = {}
         counter = 0    
          
         for i in alphabet:
             new_kmer = prefix + i
             extensions[new_kmer] = 0
             if new_kmer in kmers_dict and new_kmer not in kmers_used:
                counter += 1
                extensions[new_kmer] += 1

         if counter == 1:
            for record in extensions:
                if extensions[record] == 1:
                   contig = contig + record[-1:]
                   kmers_used.add(record)
                   kmers_used.add(rev_comp(record))
                
         else:
             direction = 'left'

     # proceed by extending in the opposite direction
     while direction == 'left':
         prefix = contig[0:kmer_length - 1]
         extensions = {}
         counter = 0

         for i in alphabet:
             new_kmer = i + prefix
             extensions[new_kmer] = 0
             if new_kmer in kmers_dict and new_kmer not in kmers_used:
                counter += 1
                extensions[new_kmer] += 1

         if counter == 1:
            for record in extensions:
                if extensions[record] == 1:
                   contig = record[0] + contig
                   kmers_used.add(record)
                   kmers_used.add(rev_comp(record))
         else:
            direction = 'stop'


     contigs.append(contig)
         




contigs = sorted(contigs, key = len, reverse = True)
for i in range(10):
  print(len(contigs[i]), '\t', contigs[i])























 
  



