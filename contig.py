#!/usr/bin/env python3

import sys
import random

filename = sys.argv[1]
kmer_length = int(sys.argv[2])
OUTPUT = open('output.txt', 'w')
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
kmers_used = set()
contigs = []

def extend_right(kmers_dict, contig = ''):
   if len(contig) == 0:
      key = random.choice(list(kmers_dict.keys()))
      prefix = key[1:]
      temp1 = key
      kmers_used.add(key)
      kmers_used.add(rev_comp(key))
      kmers_dict.pop(key)
      kmers_dict.pop(rev_comp(key))
   else:
      prefix = contig[-kmer_length + 1:]
      temp1 = contig

   alphabet = ['A', 'G', 'C', 'T']
   count = 0
   four = {}
      
   for i in alphabet:
         new_kmer = prefix + i
         four[new_kmer] = 0
         if new_kmer in kmers_dict:
            count += 1
            four[new_kmer] += 1
   if count == 1:
         for record in four:
             if four[record] == 1:
                temp1 = temp1 + record[-1:]
                kmers_used.add(record)
                kmers_used.add(rev_comp(record))
                kmers_dict.pop(record)
                kmers_dict.pop(rev_comp(record))
                return extend_right(kmers_dict, contig = temp1)
                
   else: 
     return extend_left(kmers_dict, temp1)

def extend_left(kmers_dict, contig = ''):
 temp1 = contig
# print('Left function starts here')
# print(len(kmers_dict.keys()))
 OUTPUT.write(str(len(kmers_dict.keys())) + '\n')
 if len(kmers_dict.keys()) != 0:
   prefix = contig[0:kmer_length - 1]
   alphabet = ['A', 'G', 'C', 'T']
   count = 0
   four = {}
      
   for i in alphabet:
         new_kmer = i + prefix
         four[new_kmer] = 0
         if new_kmer in kmers_dict:
            count += 1
            four[new_kmer] += 1

   if count == 1:
         for record in four:
             if four[record] == 1:
                temp1 = record[0] + temp1
                kmers_used.add(record)
                kmers_used.add(rev_comp(record))
                kmers_dict.pop(record)
                kmers_dict.pop(rev_comp(record))
                return extend_left(kmers_dict, contig = temp1)
                
   else: 
     return temp1
 else:
   return temp1
#test = extend_right(kmers_dict)
#print(test)


while len(kmers_dict.keys()) != 0:
  contig = extend_right(kmers_dict)
  contigs.append(contig)


contigs = sorted(contigs, key = len, reverse = True)
for contig in contigs:
#  print(str(len(contig)) + '\n' + contig)
  OUTPUT.write(str(len(contig)) + '\n' + contig + '\n')

OUTPUT.close()
  



























 
  



