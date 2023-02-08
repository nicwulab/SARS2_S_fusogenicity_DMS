from email import header
from Bio import SeqIO
import csv   
import pandas as pd

for record in SeqIO.parse('sarbecovirus_S_S2_alignment.txt', 'fasta'):
    print(record.id, record.seq)
def split(word):
    return list(word)
align_seq = []
align_seq_id= []
for record in SeqIO.parse('sarbecovirus_S_S2_alignment.txt', 'fasta'):
    align_seq_sub = list(record.seq)
    align_seq_id.append(record.id)
    align_seq.append(align_seq_sub)

query_seq = align_seq[0]
freq_result = []

blank_index = 0
for res in range(len(query_seq)):
    if query_seq[res] != '-':
        freq_result_sub = []
        actual_res_num = res - blank_index + 883
        freq_result_sub.append(actual_res_num)
        freq_result_sub.append(query_seq[res])
        count = 0
        for subject in align_seq:
            if query_seq[res]==subject[res]:
                count = count + 1
        freq_result_sub.append(count)
        freq_result_sub.append(count/27)
        freq_result.append(freq_result_sub)
    else:
        blank_index = blank_index + 1
print(freq_result)

data_frame = pd.DataFrame(freq_result)

data_frame.to_csv("sarbecovirus_S_S2_sequence_conservation.csv", index = False, header = ["pos","resi","conservation_count","align_freq"])