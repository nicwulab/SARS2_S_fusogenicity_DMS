#!/usr/bin/python
import os
import sys
import operator
from Bio import SeqIO
from collections import Counter

def ProcessMultilib(Rfile):
  print ("Reading %s" % Rfile)
  records = SeqIO.parse(Rfile,"fastq")
  variants = [] 
  record_count = 0
  for record in records:
    record_count += 1
    Rseq  = record.seq
    Rroi = Rseq
    if ((Rroi[0:24] == "ACATCTGCCCTGCTGGCCGGCACA") and (Rroi[-24:] == "CAGAGCAAGAGAGTGGACTTTTGC") and (len(Rroi[24:-24]) == 462)): # Only include those that have the correct forward primer sequence, correct reverse primer sequence, and the correct number of nucleotides between the primers
      Rroi = Rroi[24:-24] # Trim forward and reverse primers
      variants.append(Rroi)
    #if record_count == 1000: break
  return Counter(variants)

def Output(input_dict, bin_0_rep1_dict, bin_1_rep1_dict, bin_2_rep1_dict, bin_3_rep1_dict, mNG2_neg_rep1_dict, mNG2_pos_rep1_dict, bin_0_rep2_dict, bin_1_rep2_dict, bin_2_rep2_dict, bin_3_rep2_dict, mNG2_neg_rep2_dict, mNG2_pos_rep2_dict, bin_0_rep3_dict, bin_1_rep3_dict, bin_2_rep3_dict, bin_3_rep3_dict, outfile):
  print ("Compiling results into %s" % outfile)
  outfile = open(outfile,'w')
  muts = list(set(list(input_dict.keys())+
                  list(bin_0_rep1_dict.keys())+
                  list(bin_1_rep1_dict.keys())+
                  list(bin_2_rep1_dict.keys())+
                  list(bin_3_rep1_dict.keys())+
                  list(mNG2_neg_rep1_dict.keys())+
                  list(mNG2_pos_rep1_dict.keys())+
                  list(bin_0_rep2_dict.keys())+
                  list(bin_1_rep2_dict.keys())+
                  list(bin_2_rep2_dict.keys())+
                  list(bin_3_rep2_dict.keys())+
                  list(mNG2_neg_rep2_dict.keys())+
                  list(mNG2_pos_rep2_dict.keys())+
                  list(bin_0_rep3_dict.keys())+
                  list(bin_1_rep3_dict.keys())+
                  list(bin_2_rep3_dict.keys())+
                  list(bin_3_rep3_dict.keys())))
  outfile.write("\t".join(['mut','input','bin0_rep1','bin1_rep1', 'bin2_rep1', 'bin3_rep1', 'fus_neg_rep1', 'fus_pos_rep1', 'bin0_rep2', 'bin1_rep2', 'bin2_rep2', 'bin3_rep2', 'fus_neg_rep2', 'fus_pos_rep2', 'bin0_rep3', 'bin1_rep3', 'bin2_rep3', 'bin3_rep3'])+"\n")
  for mut in muts:
    ipt_count   = input_dict[mut]
    bin_0_rep1_count = bin_0_rep1_dict[mut]
    bin_1_rep1_count = bin_1_rep1_dict[mut]
    bin_2_rep1_count = bin_2_rep1_dict[mut]
    bin_3_rep1_count = bin_3_rep1_dict[mut]
    mNG2_neg_rep1_count = mNG2_neg_rep1_dict[mut]
    mNG2_pos_rep1_count = mNG2_pos_rep1_dict[mut]
    bin_0_rep2_count = bin_0_rep2_dict[mut]
    bin_1_rep2_count = bin_1_rep2_dict[mut]
    bin_2_rep2_count = bin_2_rep2_dict[mut]
    bin_3_rep2_count = bin_3_rep2_dict[mut]
    mNG2_neg_rep2_count = mNG2_neg_rep2_dict[mut]
    mNG2_pos_rep2_count = mNG2_pos_rep2_dict[mut]
    bin_0_rep3_count = bin_0_rep3_dict[mut]
    bin_1_rep3_count = bin_1_rep3_dict[mut]
    bin_2_rep3_count = bin_2_rep3_dict[mut]
    bin_3_rep3_count = bin_3_rep3_dict[mut]
    outfile.write("\t".join(map(str,[mut,ipt_count,bin_0_rep1_count,bin_1_rep1_count, bin_2_rep1_count, bin_3_rep1_count, mNG2_neg_rep1_count, mNG2_pos_rep1_count, bin_0_rep2_count, bin_1_rep2_count, bin_2_rep2_count, bin_3_rep2_count, mNG2_neg_rep2_count, mNG2_pos_rep2_count, bin_0_rep3_count, bin_1_rep3_count, bin_2_rep3_count, bin_3_rep3_count]))+"\n")
  outfile.close()

def main():
  outfile = 'result/S2_HR1_DMS_count_nuc.tsv'
  input_dict  = ProcessMultilib('fastq_merged/Input.assembled.fastq')
  bin_0_rep1_dict  = ProcessMultilib('fastq_merged/Bin_0_rep1.assembled.fastq')
  bin_1_rep1_dict  = ProcessMultilib('fastq_merged/Bin_1_rep1.assembled.fastq')
  bin_2_rep1_dict  = ProcessMultilib('fastq_merged/Bin_2_rep1.assembled.fastq')
  bin_3_rep1_dict  = ProcessMultilib('fastq_merged/Bin_3_rep1.assembled.fastq') 
  mNG2_neg_rep1_dict  = ProcessMultilib('fastq_merged/mNG2_neg_rep1.assembled.fastq')
  mNG2_pos_rep1_dict  = ProcessMultilib('fastq_merged/mNG2_pos_rep1.assembled.fastq')
  bin_0_rep2_dict  = ProcessMultilib('fastq_merged/Bin_0_rep2.assembled.fastq')
  bin_1_rep2_dict  = ProcessMultilib('fastq_merged/Bin_1_rep2.assembled.fastq')
  bin_2_rep2_dict  = ProcessMultilib('fastq_merged/Bin_2_rep2.assembled.fastq')
  bin_3_rep2_dict  = ProcessMultilib('fastq_merged/Bin_3_rep2.assembled.fastq') 
  mNG2_neg_rep2_dict  = ProcessMultilib('fastq_merged/mNG2_neg_rep2.assembled.fastq')
  mNG2_pos_rep2_dict  = ProcessMultilib('fastq_merged/mNG2_pos_rep2.assembled.fastq')
  bin_0_rep3_dict  = ProcessMultilib('fastq_merged/Bin_0_rep3.assembled.fastq')
  bin_1_rep3_dict  = ProcessMultilib('fastq_merged/Bin_1_rep3.assembled.fastq')
  bin_2_rep3_dict  = ProcessMultilib('fastq_merged/Bin_2_rep3.assembled.fastq')
  bin_3_rep3_dict  = ProcessMultilib('fastq_merged/Bin_3_rep3.assembled.fastq') 
  Output(input_dict, bin_0_rep1_dict, bin_1_rep1_dict, bin_2_rep1_dict, bin_3_rep1_dict, mNG2_neg_rep1_dict, mNG2_pos_rep1_dict, bin_0_rep2_dict, bin_1_rep2_dict, bin_2_rep2_dict, bin_3_rep2_dict, mNG2_neg_rep2_dict, mNG2_pos_rep2_dict, bin_0_rep3_dict, bin_1_rep3_dict, bin_2_rep3_dict, bin_3_rep3_dict, outfile)

if __name__ == "__main__":
  main()
