#!/usr/bin/env python3

# Usage:
# contig.py <ARGV1> <ARGV2> <ARGV3> <ARGV4>
# Where:
#    ARGV1 = input 'fastq' or 'fq' file containing sequencing reads
#    ARGV2 = kmer length. Default is 31
#    ARGV3 = kmer frequency cut-off. Default is 6
#    ARGV4 = output file  

import sys

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
kmers_dict = {key:val for key, val in kmers_dict.items() if val >= int(sys.argv[3])}
   
print(len(kmers_dict), 'kmers left in the dictionary after removing low frequency kmers')

kmers_used = set()
alphabet = ['A', 'G', 'C', 'T']
contigs = []

for kmer in kmers_dict:
     if kmer in kmers_used:
         continue
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
 #            for key in extensions: kmers_used.add(key)


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
#            for key in extensions: kmers_used.add(key)


     contigs.append(contig)
         
contigs = sorted(contigs, key = len, reverse = True)

# Calculate quick statistics for contig assebmly
#assembly_size
genome_size = 50000
counter = 0
n50 = 0
l50 = 0
for i in range(len(contigs)):
  counter += len(contigs[i])
  if counter >= genome_size / 2:
     n50 = len(contigs[i])
     l50 = i + 1
     break
#print('N50', n50)
#print('L50', l50)

assembly_size = 0
for i in range(len(contigs)):
  print('>contig' + (str(i)))
  print(contigs[i])
  assembly_size += len(contigs[i]) 

print(assembly_size)

output = open('output.txt', 'a')

output.write(str(kmer_length) + '\t' + str(len(contigs)) + '\n')






















 
  



