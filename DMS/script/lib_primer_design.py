#!/usr/bin/python
import sys
import operator
import itertools
from Bio import SeqIO

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def nucleotide_to_codon_seq(nucleotide_seq):
  codon_seq = [nucleotide_seq[n*3:n*3+3] for n in range(0,int(len(nucleotide_seq)/3))]
  return codon_seq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "NNK":"X"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def reverse_translation(nucleotide_seq):
  protein_seq = translation(nucleotide_seq)
  AA_codon = {
     'C': ['TGT', 'TGC'],
     'D': ['GAT', 'GAC'],
     'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
     'Q': ['CAA', 'CAG'],
     'M': ['ATG'],
     'N': ['AAC', 'AAT'],
     'P': ['CCT', 'CCG', 'CCA', 'CCC'],
     'K': ['AAG', 'AAA'],
     '_': ['TAG', 'TGA', 'TAA'],
     'T': ['ACC', 'ACA', 'ACG', 'ACT'],
     'F': ['TTT', 'TTC'],
     'A': ['GCA', 'GCC', 'GCG', 'GCT'],
     'G': ['GGT', 'GGG', 'GGA', 'GGC'],
     'I': ['ATC', 'ATA', 'ATT'],
     'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
     'H': ['CAT', 'CAC'],
     'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
     'W': ['TGG'],
     'V': ['GTA', 'GTC', 'GTG', 'GTT'],
     'E': ['GAG', 'GAA'],
     'Y': ['TAT', 'TAC'] }
  codon_seq = nucleotide_to_codon_seq(nucleotide_seq)
  codon_list_of_list = []
  for aa, codon in zip(protein_seq, codon_seq):
    assert(aa==translation(codon))
    possible_codons = [new_codon for new_codon in AA_codon[aa] if hamming(new_codon,codon) <= 1 and new_codon[0:2]==codon[0:2]]
    codon_list_of_list.append(possible_codons)
  seqs = []
  for seq in itertools.product(*codon_list_of_list):
    seq = ''.join(list(seq))
    assert(translation(seq)==protein_seq)
    seqs.append(seq)
  return (seqs)

def checking_bc_hm_dist(codon_pos, seq_rand_codon, cassettes_barcoded, cassettes_selected):
  for pos in cassettes_selected.keys():
    cas_rand_codon = cassettes_barcoded[cassettes_selected[pos]]
    diff = list(set(seq_rand_codon).symmetric_difference(set(cas_rand_codon)))
    if len(diff)-diff.count(codon_pos)-diff.count(pos) < 2: return "not okay"
  return "okay"

def select_barcode(cassettes_barcoded, aa_per_casette):
  cassettes_selected = {}
  for n in range(0,aa_per_casette):
    codon_pos = n + 1
    count = 0
    for seq in sorted(cassettes_barcoded.keys(), key=lambda x:max(cassettes_barcoded[x])):
      seq_rand_codon = cassettes_barcoded[seq]
      if codon_pos not in seq_rand_codon:
        count += 1
        dist_check = checking_bc_hm_dist(codon_pos, seq_rand_codon, cassettes_barcoded, cassettes_selected)
        if dist_check == 'okay':
          cassettes_selected[codon_pos] = seq
          break
  assert(len(cassettes_selected.keys())==aa_per_casette)
  return cassettes_selected

def barcoding(casette_seq):
  protein_seq = translation(casette_seq)
  possible_seqs = reverse_translation(casette_seq)
  barcoded_casette_seqs = [seq for seq in possible_seqs if hamming(seq, casette_seq) >= 2]
  casette_seq_codons = nucleotide_to_codon_seq(casette_seq)
  barcoded_casette_seq_dict = {}
  for barcoded_casette_seq in barcoded_casette_seqs:
    barcoded_casette_seq_codons = nucleotide_to_codon_seq(barcoded_casette_seq)
    codon_pos = 0
    codon_rand = []
    for casette_codon, barcoded_codon in zip(casette_seq_codons, barcoded_casette_seq_codons):
      codon_pos += 1
      if casette_codon != barcoded_codon:
        codon_rand.append(codon_pos)
    if len(codon_rand) >= 2:
      barcoded_casette_seq_dict[barcoded_casette_seq] = codon_rand
  return barcoded_casette_seq_dict

def design_cassettes_Fprimer(output_file,template,cassettes):
    print ("writing: %s" % output_file)
    outfile = open(output_file,'w')
    aa_per_casette = 8
    for n in cassettes:
        seq = template[n*24:n*24+60]
        cassette_original = seq[21:45]
        cassettes_barcoded = barcoding(cassette_original) #generate WT sequences with internal barcode (synonymous mutation)
        cassettes_selected = select_barcode(cassettes_barcoded, aa_per_casette) #selection 8 internal barcoded WT sequences
        id = '>Cassette'+str(n+1)
        for aa in range(0,aa_per_casette,1):
            primerID = id+ '_'+ str(aa+1)
            constant_5 = seq[0:21]
            constant_3 = seq[45::]
            codon_pos  = aa + 1
            cassette_seq = cassettes_selected[codon_pos]
            mut_region = cassette_seq[0:aa*3]+"NNK"+cassette_seq[aa*3+3::] 
            outfile.write(primerID+"\n"+constant_5+mut_region+constant_3+"\n")
    outfile.close()

def design_cassettes_Rprimer(output_file,template, cassettes):
    print ("writing: %s" % output_file)
    outfile = open(output_file,'w')
    for n in cassettes:
        seq = template[n*24:n*24+21]
        complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
        rcseq = seq.translate(complements)[::-1]
        id = '>Cassette'+str(n+1)
        outfile.write(id+ '_Rprimer'+  "\n"+ rcseq + "\n")
    outfile.close()

def read_ref(reffile):
  records = SeqIO.parse(reffile,"fasta")
  ref_dict = {}
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    ref_dict[ID] = seq
  return ref_dict

def main():
    #Reference sequence should contain 21 nt upstream of the mutation of interest
    ref_dict = read_ref('Fasta/sars-cov-2-spike-hr1-linker.fa')
    template = ref_dict['SARS2-HR1-linker']
    cassettes = range(0,19,1) #cassettes = range(0,N,1), where N is the number of cassattes needed
    F_output_file = 'Fasta/SARS-2-Spike-HR1_lib_Fprimer_bc.fa'
    R_output_file = 'Fasta/SARS-2-Spike-HR1_lib_Rprimer_bc.fa'
    design_cassettes_Fprimer(F_output_file,template,cassettes)
    design_cassettes_Rprimer(R_output_file,template,cassettes)

if __name__ == "__main__":
  main()
