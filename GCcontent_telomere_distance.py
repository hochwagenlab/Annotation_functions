# GCcontent_telomere_distance.py
# made by Tovah Markowitz
# uses python2.7
# 1/3/18
########################
# MODULES

from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import optparse

########################
# FUNCTIONS

def divide_equal_steps(chr_lens, tel_dist, num_max_groups):
    """Divide genome into equal-sized chunks based upon distance from telomeres.
    """
    dividers = { }
    for chrom in chr_lens:
        num_groups = chrom[1] / (tel_list*2)
        if num_groups > num_max_groups:
            num_groups = num_max_groups
        mid_group = chr_len - (2 * num_groups * tel_dist) 
        tmp = tel_dist*num_groups + mid_group + 1
        dividers[chrom[0]] = (range(tel_dist, tel_dist*num_groups + 1, tel_dist)
                        + range(tmp, (tmp + (tel_dist*num_groups)), tel_dist))
    return(dividers)

def divide_unequal_groups(chr_lens, cut_points):
    """Divide genome into unequal-sized chunks based upon distance from 
    telomeres.
    """
    dividers = { }
    for chrom in chr_lens:
        tmp = [chrom[1] - cut_pt for cut_pt in cut_points]
        tmp.sort()
        dividers[chrom[0]] = cut_points + tmp
    return(dividers)

def import_genes(genesGFF):
    """Read in genes from GFF file.
    """
    # read in GFF file with genic positions listed
    f = open(genesGFF, 'r')
    gff_in = f.readlines()
    f.close
    genes = [row.strip().split('\t') for row in gff_in]
    # make start and end positions into integers
    for i in range(len(genes)):
        genes[i][3] = int(genes[i][3])
        genes[i][4] = int(genes[i][4])
    genes.sort()
    return(genes)

def check_gene_overlaps(genes):
    """Ensure that none of the genes overlap.
    """
    for i in range (len(genes) - 1):
        if (genes[i][0] == genes[i+1][0]):
                if (genes[i][4] >= genes[i+1][3]):
                        print "Error: some genes in gff overlap."

def determine_bin(dividers, chromosome, position):
    """Determine which bin a specific position belongs to.
    """
    bin_num = 0
    for div in dividers[chromosome]:
           if div <= position:
                bin_num += 1
    return bin_num

def calcGC_genic_intergenic( genicGFF, dividers, fasta_dict):
    """This function splits the regions of the genome by both bins and
    intergenic vs. genic sequence. For each short fragment of the genome,
    it then calculates number of bases and the number that are G/C.
    """
    genic = import_genes(genicGFF)
    check_gene_overlaps(genes)
    # preassign keys and vectors to add GC counts and base counts to each group
    # groups determined from number of dividers per chromosome plus 1
    chrom = ''
    genGC = { }
    interGC = { }
    for key in dividers.keys():
        genGC[key] = [ [[],[]] for _ in range(len(dividers[key]) + 1) ]
        interGC[key] = [ [[],[]] for _ in range(len(dividers[key]) + 1) ]
    for row in genic:
        
        # if new chromosome
        if row[0] != chrom:
            if chrom is not '':
                if binB == len(dividers[chrom]):
                    seqB = str(fasta_dict[chrom].seq[end:])
                elif (binB + 1) == len(dividers[chrom]):
                    div = dividers[chrom][binB]
                    seqB = str(fasta_dict[chrom].seq[end:div])
                    seqA = str(fasta_dict[chrom].seq[div:])
                    interGC[chrom][binB + 1][0].append(seqA.count("G")
                                                           + seqA.count("C"))
                    interGC[chrom][binB + 1][1].append(len(seqA))
                interGC[chrom][binB][0].append(seqB.count("G")
                                                   + seqB.count("C"))
                interGC[chrom][binB][1].append(len(seqB))
            end = 1
            chrom = row[0]
            binB = 0
            
        # start with intergenic region
        binA = determine_bin(dividers, chrom, row[3])
        if binA == binB:
            seqA = str(fasta_dict[chrom].seq[end:row[3]])
        elif (binA - 1) == binB:
            div = dividers[chrom][binA-1]
            seqB = str(fasta_dict[chrom].seq[end:div])
            seqA = str(fasta_dict[chrom].seq[div:row[3]])
            interGC[chrom][binB][0].append(seqB.count("G") + seqB.count("C"))
            interGC[chrom][binB][1].append(len(seqB))
        interGC[chrom][binA][0].append(seqA.count("G") + seqA.count("C"))
        interGC[chrom][binA][1].append(len(seqA))
        
        # then analyze genic region
        end = row[4]
        binB = determine_bin(dividers, chrom, row[4])
        # if both ends of genic region are in the same bin
        if binA == binB:
            seqA = str(fasta_dict[chrom].seq[row[3]:row[4]])            
        elif (binA + 1) == binB:
            div = dividers[chrom][binA]
            seqA = str(fasta_dict[chrom].seq[row[3]:div])
            seqB = str(fasta_dict[chrom].seq[div:row[4]])
            genGC[chrom][binB][0].append(seB.count("G") + seqB.count("C"))
            genGC[chrom][binB][1].append(len(seqB))
        genGC[chrom][binA][0].append(seqA.count("G") + seqA.count("C"))
        genGC[chrom][binA][1].append( len(seqA) )
    return(interGC, genGC)

def sumGC(GCstructure, dividers):
    """Calculate the total number of bases and the number of G/C bases
    within each large region. Input should be output of 
    calc_GC_genic_intergenic. 
    """
    GCout = { }
    for key in dividers.keys():
        GCout[key] = [ [ ] for _ in range(len(dividers[key]) + 1)]
    for key in GCstructure.keys():
        for i in range(len(dividers[key]) + 1):
            GC = sum(GCstructure[key][i][0])
            seq_len = sum(GCstructure[key][i][1])
            GCout[key][i] = [GC, seq_len]
    return GCout

def write_sumGC_results(genicGC_sum, intergenicGC_sum, file_name):
    f = open(file_name, 'w')
    f.write('chromosome' + '\t' + 'group_num' + '\t' + 'genic_GC_count' + '\t'
                + 'genic_num_bases' + '\t' + 'intergenic_GC_count' + '\t'
                + 'intergenic_num_bases'+ '\n' )
    for key in genicGC_sum.keys():
        for group in range(len(genicGC_sum[key])):
            tmp1 = '\t'.join([str(i) for i in genicGC_sum[key][group]])
            tmp2 = '\t'.join([str(i) for i in intergenicGC_sum[key][group]])
            f.write(key + '\t' + str(group+1) + '\t' + tmp1 + '\t' + tmp2
                        + '\n')
    f.close()

def group_summarize(input):
    """Calculate G/C ratio over each large region. Regions are also
    merged based upon distance from telomere on each chromosome.
    """
    num_groups = max([len(input[key]) for key in input.keys()]) / 2 + 1
    groupGC = [ [ ] for _ in range(num_groups) ]
    group_len = [ [ ] for _ in range(num_groups) ]
    group_ratio = [ [ ] for _ in range(num_groups) ]
    for key in input.keys():
        for i in range(num_groups):
            # start by grabbing the paired bins
            # all of which will be equally sized per group
            if len(input[key])/2 > i:
                tmpGC = input[key][i][0] + input[key][:i+1:-1][0][0]
                tmp_len = input[key][i][1] + input[key][:i+1:-1][0][1]
                groupGC[i].append(tmpGC)
                group_len[i].append(tmp_len)
                group_ratio[i].append(float(tmpGC) / float(tmp_len))
            # if central bin, just grab bin and add into groupGC and group_len
            # for group_ratio, separate slightly from other values in group
            elif len(input[key])/2 == i:
                tmpGC = input[key][i][0]
                tmp_len = input[key][i][1]
                groupGC[i].append(tmpGC)
                group_len[i].append(tmp_len)
                group_ratio[i].append([float(tmpGC) / float(tmp_len)])
    outRatio = [ [ ] for _ in range(num_groups) ]
    for i in range(num_groups):
        out_ratio[i] = float(sum(groupGC[i])) / float(sum(group_len[i]))
    return(out_ratio, group_ratio)

def write_ratio_results(genic_ratio, intergenic_ratio, file_name):
    f = open(file_name, 'w')
    f.write('groupNum' + '\t' + 'genicVal' + '\t' + 'intergenicVal' + '\n')
    for row in range(len(genicRatio)):
        if genic_ratio[row][0] == 0:
            genic_ratio[row].pop(0)
            intergenic_ratio[row].pop(0)
        for pos in range(len(genic_ratio[row])):
            if isinstance(genic_ratio[row][pos], float):
                genic_val = str(genic_ratio[row][pos])
                ingenic_val = str(intergenic_ratio[row][pos])
            else:
                genic_val = str( genic_ratio[row][pos][0])
                intergenic_val = str(intergenic_ratio[row][pos][0])
            f.write(str(row+1) + '\t' + genic_val + '\t' + ingenic_val + '\n')
    f.close()
    
########################
# MAIN

desc="""
This script is designed to split the genome into chunks and calculate
the proportion that is G/C. Note: if the reference genome contains N, the
A/T content cannot be calculated directly from the output of this function.
First, the genome is split into large bins based upon their distance from
telomeres. Second, each bin is split into smaller windows based upon whether
they are genic or intergenic regions. (Caution: this only works when genes do 
not overlap.) Third, for each independent window, the number of G/C bases and 
window length are calculated. Fourth, total genic G/C content and intergenic 
content G/C content is calculated separately per bin. For the summarized output
files, bins are combined based upon distance from telomere end (i.e., if a 
chromosome has three bins, the outermost two bins will be combined). Ratios are
calculated last.
Example 1: 
`python GCcontent_telomere_distance.py -f "Yue.SK1.genome.fa" 
-g "SK1.gene.gff" -c "20000,110000" -r "YueSK1_GC_20000_100000_raw.txt"`
Example 2: 
`python GCcontent_telomere_distance.py -f "Yue.SK1.genome.fa" 
-g "SK1.gene.gff" -c 30000 -m 9 -s "YueSK1_GC_30000_9cuts_summary.txt"`
"""

# parse object for managing input options
parser = optparse.OptionParser(description=desc)

parser.add_option('-c', dest='cutpoints', default= '30000', help='Indication \
of fragment size. If cutting the genome into equal size fragments, give one \
value. If cutting the genome into unequal size fragments, give a quoted \
comma-delimited list. Default: 30000')
parser.add_option('-f', dest='fastafile', default='', help='Name and path \
of the fasta file of the reference genome.')
parser.add_option('-g', dest='gfffile', default='', help='Name and path of \
the list of genes for the reference genome in GFF format.')
parser.add_option('-m', dest='maxcuts', default='3', help='Maximum number of \
cuts per chromosome when using equal size fragments. Default: 3')
parser.add_option('-r', dest='rawout', default='', help='Output file name. \
This file, if called, will contain the G/C content and total sequence length \
of each fragment analyzed.')
parser.add_option('-s', dest='sumout', default='', help='Output file name. \
This file, if called, will contain the summarized ratio values from the \
analysis.')

(options,args) = parser.parse_args()
cut_points = options.cutpoints
fasta_file = options.fastafile
gff_file = options.gfffile
max_cuts = options.maxcuts
out_file_raw = options.rawout
out_file_summary = options.sumout

FSA = SeqIO.to_dict(SeqIO.parse( fasta_file, 'fasta'))
chr_lens = [[key,len(FSA[key].seq)] for key in FSA.keys()]

cut_points = map(int,cut_points.split(','))
cut_points.sort()

if len(cut_points2) != 1:
    dividers = divide_unequal_groups(chr_lens, cut_points)
else:
    dividers = divide_equal_steps(chr_lens, cut_points, int(max_cuts))

intergenicGC, genicGC = calcGC_genic_intergenic(gff_file, dividers, FSA)
intergenicGCsum = sumGC(intergenicGC, dividers)
genicGCsum = sumGC(genicGC, dividers)

if out_file_raw:
    write_sumGC(genicGCsum, intergenicsumGC, out_file_raw)

if out_file_summary:
    intergenic_ratio1, intergenic_ratio2 = group_summarize(intergenicGCsum)
    genic_ratio1, genic_ratio2 = group_summarize(genicGCsum)
    write_ratio_results(genic_ratio2, intergenic_ratio2, out_file_summary)

