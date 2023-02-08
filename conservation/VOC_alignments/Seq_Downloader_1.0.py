import os
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "wenhaoo2@illinois.edu"  # Always tell NCBI who you are
my_file = open("VOC_Genbank_ID.txt", "r")
data = my_file.read()
id_list = data.replace('\n', ',').split(",")
print(id_list)
my_file.close()
filename = [s + ".fasta" for s in id_list] #idListNew comes from idList with ".fasta"
index = 0
for x in filename:
    if not os.path.isfile(x):
     # Downloading...
        net_handle = Entrez.efetch( db="nucleotide", id= id_list[index], rettype="fasta", retmode="text" )
        out_handle = open(x, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        index = index + 1
        print("Saved")
    print("Parsing...")
    record = SeqIO.read(x, "fasta")
    print(record)
