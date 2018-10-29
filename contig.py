#!/usr/bin/env python3
# Usage:
# contig.py <-i> <-k> <-c> 
# Options:
#    -i = input 'fastq' or 'fq' file containing sequencing reads
#    -k = kmer length. Default is 31
#    -c = kmer frequency cut-off. If no value given, the cutoff will be selected by the assembler
#    -o = output contigs in a fasta format
#    -r = write report and statistics in a file
    
import sys
import getopt
import re

#--------------- Get options from the command line------------------------
inputfile = ''
kmer_length = 31
cutoff = 0
report = ''
outfile = ''
opts, argv = getopt.getopt(sys.argv[1:], 'hi:k:c:o:r:', ['help','ifile','kmer','cutoff', 'output', 'report'])
for opt, arg in opts:
    if opt in ('-h', '--help'):
       # Print help 
       print('Program usage:')
       print('assembler.py -i <inputfile> -k <kmer_length> -c <cutoff>')
       print('\t' + '-i/--ifile' + '\t' + 'input file containing sequencing reads in the fastq format')
       print('\t' + '-k/--kmer' + '\t' + 'kmer length. Default is 31')
       print('\t' + '-c/--cutoff' + '\t' + 'kmer frequency cut-off. If no value provided, the cutoff will be selected by the assembler')
       print('\t' + '-o/--output' + '\t' + 'output file for assembled contigs')
       print('\t' + '-r/--report' + '\t' + 'write report to a specified file')
       sys.exit()
    elif opt in ('-i','--ifile'): inputfile = arg
    elif opt in ('-k','--kmer'): kmer_length = int(arg)
    elif opt in ('-c', '--cutoff'): cutoff = int(arg)
    elif opt in ('-o', '--output'): outfile = arg
    elif opt in ('-r', '--report'): report = arg  

#-----------------Convert to reverse complement---------------------------
def rev_comp(x):
  output =  x.lower().replace('a', 'T').replace('t', 'A').replace('g','C').replace('c','G')[::-1]
  return output

#-----------------Count unique kmers---------------------------------------
def kmer_count(inputfile, kmer_length):
   FILE = open(inputfile, 'r')
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

######################## End of function definitions #################################
if report != '':  REPORT =  open(report, 'w')
if outfile != '': ASSEMBLY = open(outfile, 'w')

kmers_dict = kmer_count(inputfile, kmer_length)
print(len(kmers_dict), 'kmers gathered in the dictionary')
if report != '': REPORT.write(str(len(kmers_dict)) + '\t' + 'kmers gathered in the dictionary' + '\n')

# remove low frequency kmers
kmers_dict = {key:val for key, val in kmers_dict.items() if val >= cutoff}
print(len(kmers_dict), 'kmers left in the dictionary after removing low frequency kmers')
if report != '': REPORT.write(str(len(kmers_dict)) + '\t' + 'kmers left in the dictionary after removing low frequency kmers' + '\n')

kmers_used = set()
alphabet = ['A', 'G', 'C', 'T']
contigs = []

for kmer in kmers_dict:
     if kmer in kmers_used: continue  # Don't reuse kmers

     kmers_used.add(kmer)             # Put kmers in the 'used' category so they would not be counted or used for contig extension later
     kmers_used.add(rev_comp(kmer))

     direction = 'right'
     # Begin with the random kmer as a seed for extension
     contig = kmer

     # Extend contig in the right direction           
     while direction == 'right':           
         prefix = contig[-kmer_length + 1:]
         extensions = {}
         counter = 0    
          
         for i in alphabet:            # Construct 4 possible kmers to extend the contig
             new_kmer = prefix + i
             extensions[new_kmer] = 0
             if new_kmer in kmers_dict and new_kmer not in kmers_used:
                counter += 1
                extensions[new_kmer] += 1

         if counter == 1:              # If there is only 1 way to extend the contig, do it
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
print(len(contigs))
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
print('N50', n50)
print('L50', l50)
if report != '' : REPORT.write('N50:' + '\t' + str(n50) + '\n') # Write stats to the report
if report != '' : REPORT.write('L50:' + '\t' + str(l50) + '\n')

for i in range(len(contigs)):
  print('>contig' + str(i) + '\t' + 'length: ' + str(len(contigs[i])))
  record =  re.findall(r'\w{1,100}', contigs[i])                         #
  record = '\n'.join(record)                                    #
  if outfile != '' : ASSEMBLY.write('>contig' + str(i) + '\n')  # Write contigs to the file, wrap line if they are longer than 100 nt
  if outfile != '' : ASSEMBLY.write(record + '\n')          #























 
  



