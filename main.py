from seq_input import read_fasta
from read_orf import readingFrames
from read_orf import ORFData
from read_orf import printORFs
from read_orf import sequenceORFs
from seq_output import maxLenIndex
from seq_output import outputORF

print("This program will find ORFs for sequences.\n When prompted, please enter a FASTA file with 1 or more sequences.")

file_path = input("Enter a FASTA file: ")

min_bp = int(input("Enter minimum length in bp for ORFs: "))

seq = read_fasta(file_path)
# Test
#import sys
#sys.path.append('D:/share/divergene/script/final')
#seq = read_fasta("D:/share/divergene/script/final/sequence.fasta")
#min_bp = 300

rf = {}
ORF_dict = {}
ORF_pair_dict = {}
ORF_seq_dict = {}
min_aa = round(min_bp/3)


for i in seq.keys():
    rf[i] = readingFrames(seq[i])
    ORF_dict[i] = ORFData(rf[i])
    ORF_pair_dict[i] = printORFs(ORF_dict[i][0],ORF_dict[i][1],min_aa)
    ORF_seq_dict[i] = sequenceORFs(ORF_pair_dict[i],rf[i])
    outputORF(ORF_pair_dict[i],ORF_seq_dict[i],max_codon=15,seq_name=i,max_only=True)

