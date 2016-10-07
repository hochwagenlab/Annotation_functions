"""
Sequence Extraction
by Tovah
(Based on Sunny's sequence extraction scripts)

Inputs: genome fasta file, file with positions, query (optional), range (optional)
Outputs: fasta file of specific regions
Updated: 09/21/16
Change: updated function to work with SK1 and S288C chromosome naming systems
Updated: 10/7/16 to improve error messages and documentation
"""

##################################################################################
# Modules

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import optparse
import re
import sys

##################################################################################
# Functions

def main(fastaFile, inputFile, output, r_up, r_down, query):
    # read in fasta sequence
    genome = SeqIO.to_dict( SeqIO.parse( fastaFile, 'fasta' ) )
    # run readBedGff
    input = readBedGff(inputFile)

    r_up = int(r_up)
    r_down = int(r_down)
    
    # note chromosome number starts with 1 and python indexing starts with 0
    # takes into account range
    # if query, find query within input file
    if query is not "None":
        queryInfo = [ i for i in input if query in i[3] ]
        if not queryInfo:
            print( "Query not found in " + inputFile + ". Try systematic name." )
            sys.exit()
    else:
        queryInfo = input

    # extract sequences and write output FASTA file
    output_fasta = [0] * len(queryInfo)
    for i in range( len(queryInfo) ):
        output_fasta[i] = seq_extract( queryInfo[i], genome, r_up, r_down )
    # write output file
    f = open( output, 'w' )
    f.write( '\n'.join( i for i in output_fasta ) )
    f.close()

def seq_extract( row, fasta, r_up, r_down ):
    # take "row" from readBedGff file and extract sequence using biopython
    # put into FASTA format
    if row[4] == '+':
        start = row[1] - r_up - 1
        end = row[2] + r_down
        sequence = str( fasta[ row[0] ].seq[start:end] )
    else:
        start = row[2] + r_up
        end = row[1] - r_down - 1
        sequence = fasta[ row[0] ].seq[end:start]
        sequence = str( sequence.reverse_complement() )
    if "=" in row[3]:
        name = ( re.search(r'ID=(.*);', row[3] ).group(1)
            + ',' + re.search('Name=(.*)', row[3] ).group(1) )
    else:
        name = row[i][3]
    out_fasta = ( '>' + name + '\t' + fasta[row[0]].id
        + ':' + str(start) + '-' + str(end) + '('
        + row[4] + ')' + '\n' + sequence )
    return out_fasta
             
def readBedGff(inputName):
    # function to read in Gff or Bed files and convert into a universal format
    # read in input file
    f = open(inputName,'r')
    input = f.readlines()
    f.close()
    
    # split into a table
    input = [ row.strip().split('\t') for row in input if not row.startswith("#") ]
        
    # determine type of input file and make four column table
    # column 1: chr#, column 2: start, column 3: end, column 4: name, column 5: orientation
    if len(input[0]) == 9:          # gff file
        print( "Input file is a GFF" )
        for i in range( len(input) ):
            chr = input[i][0]
            input[i] = [ chr, int(input[i][3]), int(input[i][4]), input[i][8], input[i][6] ]
    elif len(input[0]) == 3:          # simplest bed file
        print( "Input file is a 3-column BED" )
        for i in range( len(input) ):
            chr = input[i][0]
            input[i] = [ chr, int(input[i][1]), int(input[i][2]), 'id'+i, '+']
    elif len(input[0]) == 5:          # bed file made for summits
        print( "Input file is a 5-column BED" )
        for i in range( len(input) ):
            chr =input[i][0]
            input[i] = [ chr, int(input[i][1]), int(input[i][2]), input[i][3], '+']
    else:
        print (inputFile + " cannot be read. Make sure you are using a GFF file, a \
        3 column BED file, or a 5 column BED file.")
        sys.exit()
    
    return (input)

################################################################################
# Main

# parse object for managing input options

desc="""
Purpose: to extract sequence(s) from a FASTA file based upon sequence coordinates
given from a separate GFF or BED file.
If no specific query is given, the output file will contain sequences for every
row in the input file. This function will produce an output as long as the
chromosome naming system is the same in both FASTA and input file. R_up and r_down
expands the coordinates of the sequence to be extracted beyond those in the input
file. Example:
python FastaSeqExtract_biopython.py -f sk1_MvO_V1.fasta -i SK1_annotation_modified.gff
-q "YIL072W" -o "YIL072W.fa"
"""

parser = optparse.OptionParser(description=desc)

# essential data, defines commandline options
parser.add_option ('-f', dest = 'fastaFile', default = '', help = "This input \
is the name of the FASTA file of the entire genome.")
parser.add_option ('-i', dest = 'inputFile', default = '', help = "This input \
is the name of the file that contains positions of sequences to be extracted.")
parser.add_option ('-u', dest = 'r_up', default = '0', help = "This input \
is the length of sequence to be extracted upstream of the designated feature. \
Default = 0bp")
parser.add_option ('-d', dest = 'r_down', default = '0', help = "This input \
is the length of sequence to be extracted downstream of the designated feature. \
Default = 0bp")
parser.add_option ('-q', dest = 'query', default = 'None', help = "This input \
is the name of a specific query to search for. This includes: GeneID, Gene Name, \
and summit name. This is case-sensitive. Default: entire input file")
parser.add_option ('-o', dest = 'output', default = '', help = "This input \
is the name of the newly created output FASTA file.")

# load the inputs
(options,args) = parser.parse_args()

# reads the inputs from commandline
fastaFile = options.fastaFile
inputFile = options.inputFile
r_up = options.r_up
r_down = options.r_down
query = options.query
output = options.output

a = main(fastaFile, inputFile, output, r_up, r_down, query)


#fastaFile="sk1_MvO_V1.fasta"
#inputFile="RPLS.gff"
