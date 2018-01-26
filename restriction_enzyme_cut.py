"""
Restriction Enzyme Fragment Lengths
by Tovah Markowitz
runs with python 2.7
"""

##################################################################################
# Modules

from Bio import SeqIO
import optparse
import re
import matplotlib.pyplot as plt

##################################################################################
# Functions

def re_cut(fasta_filename,re_seq,out_root):
    # step1: load in genome sequence
    sequence = list(SeqIO.parse(fasta_filename, 'fasta'))

    lengths = [ ]
    # step2: for each chromosome:
    # 1) find recognition sequence
    # 2) calculate lengths of fragments
    # 3) output list of lengths
    for i in range(len(sequence)):
        chr_seq = str(sequence[i].seq)
        chr_len = len(chr_seq)
        positions = [m.start() for m in re.finditer(re_seq, chr_seq)]
        positions.insert(0,1)
        positions.append(chr_len)
        L = [positions[n] - positions[n-1] for n in range(len(positions))]
        lengths.extend(L)


    # step3: create histogram output
    plt.hist(lengths,bins=20)
    plt.title("Genome Cut With Sequence " + re_seq)
    plt.xlabel("Length of cut fragment (bp)")
    plt.ylabel("Count")
    plt.savefig(out_root + "_cut_" + re_seq + '_hist.png')
    plt.close()

    # step4: calculate number of fragments greater than 1000
    # and produce output document
    large = sum(j > 1000 for j in lengths)
    maximum = max(lengths)
    f = open((out_root + "_CutOver1000.txt"),'a')
    f.write(re_seq + "\t" + str(large) + "\t" + str(maximum) + "\n")
    f.close()
    

def all_4mers():
    '''Function to create a list of all possible 4mers'''
    letters = ["A","C","G","T"]
    fourmers = [0]*256
    count=0
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    fourmers[count] = (letters[i] + letters[j]
                                           + letters[k] +letters[l])
                    count += 1
    return fourmers

##################################################################################
# Main

desc="""The input is a fasta file and a restriction enzyme sequence. 
The output is a histogram of output fragment lengths. Also included is a
text file of the number of fragments over 1kb in length and the maximum
fragment length.
"""

# parse object for managing input options
parser = optparse.OptionParser(description=desc)

parser.add_option('-f', dest='fastafile', default='', help='Name and path \
of the fasta file of the reference genome.')
parser.add_option('-o', dest='outroot', default='', help='Identifying \
information to add to every output file.')

(options,args) = parser.parse_args()
fasta_filename = options.fastafile
out_root = options.outroot

fourmers = all_4mers()
for k in range(len(fourmers)):
    re_cut(fasta_filename, fourmers[k], out_root)
