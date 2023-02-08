import os
from Bio import SeqIO


def concatenate (fileName, OutPutfileName): #combine all the fasta sequence into one fasta file
    with open(OutPutfileName, 'w') as w_file:
        for filen in fileName:
         with open(filen, 'rU') as o_file:
                seq_records = SeqIO.parse(o_file, 'fasta')
                SeqIO.write(seq_records, w_file, 'fasta')

my_file = open("Genbank_Id_List.txt", "r")
data = my_file.read()
id_list = data.replace('\n', ',').split(",")
print(id_list)
my_file.close()
filename = [s + ".fasta" for s in id_list]

concatenate (filename, "sarbecovirus_sequences.fasta")