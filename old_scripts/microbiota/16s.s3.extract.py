#!/usr/bin/env /ldfssz1/ST_META/share/flow/anaconda3/bin/python3 

# ==============================================================================
# 16s s2 extract:
#               1. select manually alignment result by AliView software on windows
#                   a. 515F - 806R
#                   b. 515F - 926R
#               2. delete unaliagnment fasta in the start and end
#               3. delete the base gap in each fasta and calculate each fasta length
#               4. simulate 16s rDNA PCR in silico by FastPCR software
#                   a. generate forward and reverse primer on 515F,806R,926R
#                   b. do manully simulation in windows platform
#                   c. summary simulational result
#
#   Notice: Program to run according to parameters
#
# Authors: ZouHua
#
# Please type "./16s.s3.extract.py -h" for usage help
#
# ==============================================================================


__author__ = ('ZouHua (zouhua@genomics.cn)')
__version__ = '0.1'
__date__ = '20190110'

import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


try:
    import argparse as ap
except ImportError:
    sys.exit("Unable to find argparse module")


def parse_arguments(args):
    """
    1. -a: filter alignment file in AliView
    2. -t: trim outcome on gap base
    3. -n: calculate each fasta result
    4. -c: degenerate codons file
    5. -p: produce forward and reverse primer based on degenerate codons
    6. -f: FastPCR result
    7. -s: summary PCR simulational result
    """

    parser = ap.ArgumentParser(
        description="DESCRIPTION\n"
        "16s.s3.extract version "+__version__+" ("+__date__+"): \n"
        "16s extract \n"
        "AUTHORS: "+__author__+"\n\n"
        "",
        formatter_class=ap.RawTextHelpFormatter,
        prog="16s.rDNA.extract.py")
    parser.add_argument(
        '-a', '--alig', metavar='<file>', type=str,
        help="trim alignment fasta file\n")
    parser.add_argument(
        '-t', '--trim', metavar='<file>', type=str,
        help="trim gap base fasta\n")
    parser.add_argument(
        '-n', '--numb', metavar='<file>', type=str,
        help="number of base in each fasta\n")
    parser.add_argument(
        '-c', '--code', metavar='<file>', type=str,
        help="degenerate codons file\n")
    parser.add_argument(
        '-p', '--prim', metavar='<file>', type=str,
        help="primer outcome\n")
    parser.add_argument(
        '-f', '--fpcr', metavar='<file>', type=str,
        help="outcome of FastPCR\n")
    parser.add_argument(
        '-s', '--sum', metavar='<file>', type=str,
        help="summary FastPCR outcome\n")
    parser.add_argument(
        '-d', '--dir', metavar='<directory>', type=str,
        help="output directory of extract\n",
        required=True)
    parser.add_argument(
        '-v', '--version', action='version',
        version="16s.s3.extract version {} ({})".format(__version__, __date__),
        help="Prints the current 16s.s3.extract version and exit")
    return parser.parse_args()


# forward primer 515F
forward_primer_515F = "GTGYCAGCMGCCGCGGTAA"
# reverse primer 806R & 926R
reverse_primer_806R = "GGACTACNVGGGTWTCTAAT"
reverse_primer_926R = "CCGYCAATTYMTTTRAGTTT"


def trim_fasta(alignment_file, trim_file, calculate_file):
    """
    1. trim alignment outcomes
    2. calculate trim result
    """
    with open(trim_file, "w") as trim_handle, open(calculate_file, "w") as calculate_handle:
        records = SeqIO.parse(alignment_file, "fasta")
        for record in records:
            sequence = record.seq
            if not sequence.startswith("-") and not sequence.endswith('-'):
                copy_sequence = (str(sequence) + '.')[:-1]
                new_sequence = copy_sequence.replace('-', '')
                new_sequence_length = len(new_sequence)
                if new_sequence_length > 0 :
                    new_record = SeqRecord(
                                        Seq(new_sequence, generic_dna),
                                        id=record.id, name=record.name,
                                        description=record.description)
                    SeqIO.write(new_record, trim_handle, "fasta")
                    calculate_handle.write(record.id + "\t" + str(new_sequence_length) + "\n")


def codon_into_dictionary(codon_file):
    """
    storage codons into dict
    """
    codon_dictionary = {}
    with open(codon_file, 'r') as codon_handle:
        for codon in codon_handle:
            codon = codon.strip()
            samples = re.split(r'\s+', codon)
            # key -> degenerate codons; value -> base
            codon_dictionary[samples[0]] = samples[1]
    return(codon_dictionary)


def primer_produce(codon_file, primer_file):
    """
    A degenerate nucleotide sequence and translate it into its multiple possible oligos
    """
    codon = codon_into_dictionary(codon_file)
    with open(primer_file, "w") as primer_handle:
        f_515F_nest = recursive_nucleotides(forward_primer_515F, codon)
        f_515_list = primer_type_output(f_515F_nest, "515   f")
        r_806R_nest = recursive_nucleotides(reverse_primer_806R, codon)
        r_806R_list = primer_type_output(r_806R_nest, "806   r")
        r_926R_nest = recursive_nucleotides(reverse_primer_926R, codon)
        r_926R_list = primer_type_output(r_926R_nest, "926   r")

        for seq_515 in f_515_list:
            primer_handle.write(seq_515 + "\n")
        for seq_806 in r_806R_list:
            primer_handle.write(seq_806 + "\n")
        for seq_926 in r_926R_list:
            primer_handle.write(seq_926 + "\n")


def recursive_nucleotides(sequence, codon_dictionary):
    """
    Recursive function that replaces degenerate nucleotides with all combinations
    """
    sequence_set = []
    if re.search(r'(.*)([RYSWMKBDHVN])(.*)', sequence, re.M):
        search = re.search(r'(.*)([RYSWMKBDHVN])(.*)', sequence, re.M)
        former = search.group(1)
        codon = search.group(2)
        latter = search.group(3)
        degenerate_nucletides = codon_dictionary[codon]
        bases = re.split(r'/', degenerate_nucletides)
        for base in bases:
            new_sequence = "".join([former, base, latter])
            sequence_set.append(recursive_nucleotides(new_sequence, codon_dictionary))
        return(sequence_set)
    else:
        return(sequence)


def flatten_nest_list(nest_list):
    """
    release nested list by recursive_nucleotides
    """
    for element in nest_list:
        if not isinstance(element, list):
            yield element
        else:
            for elem in flatten_nest_list(element):
                yield elem


def primer_type_output(nest_list, types):
    """
    return list -> types \t sequence
    """
    sequence_type_primer = []
    for seq in flatten_nest_list(nest_list):
        sequence_type = types + "   " + seq
        sequence_type_primer.append(sequence_type)
    return(sequence_type_primer)


def summary_fastpcr(trim_file, fastpcr_file, fastpcr_summary):
    """
    1. total sequence number
    2. amplified sequences
    3. only forward primer sequences
    4. only reverse primer sequences
    5. no matching both-side sequences
    """
    count_total = 0
    count_both = 0
    count_forward = 0
    count_reverse = 0
    # total sequence
    records = SeqIO.parse(trim_file, "fasta")
    for record in records:
        count_total = count_total + 1

    # split sign in fastpcr file
    fastpcr_file_handle = open(fastpcr_file, "r")
    lines = fastpcr_file_handle.read()
    lines = re.split(r'In silico Primer', lines)
    for line in lines:
        if re.findall(r'f  5', line) and re.findall(r'r  5', line):
            count_both = count_both + 1
        elif re.findall(r'f  5', line):
            count_forward = count_forward + 1
        elif re.findall(r'r  5', line):
            count_reverse = count_reverse + 1
    fastpcr_file_handle.close()

    # summary output
    count_non = count_total - (count_both + count_forward + count_reverse)
    fastpcr_summary_handle = open(fastpcr_summary, "w")
    fastpcr_summary_handle.write(
        "#sequences: " + str(count_total) + "\n" +
        "#sequences could be amplified: " + str(count_both) + "\n" +
        "#sequences only match forward primer: " + str(count_forward) + "\n" +
        "#sequences only match reverse primer: " + str(count_reverse) + "\n" +
        "#sequences don't match both-side primers: " + str(count_non) + "\n")
    fastpcr_summary_handle.close()


def makedir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def main():
    args = parse_arguments(sys.argv)

    # all the result directory
    cwd = os.getcwd()
    result_folder = "/".join([cwd, args.dir])
    makedir(result_folder)
    # trim sequence
    if args.alig is not None:
        trim_output_file = "/".join([result_folder, args.trim])
        length_sequence = "/".join([result_folder, args.numb])
        trim_fasta(args.alig, trim_output_file, length_sequence)

    # primer
    if args.code is not None:
        primer_name = os.path.basename(args.prim)
        primer_codon = "/".join([result_folder, primer_name])
        primer_produce(args.code, primer_codon)

    # summary fastPCR
    if args.fpcr is not None:
        trim_name = "/".join([result_folder, args.trim])
        fastPCR_summary = "/".join([result_folder, args.sum])
        summary_fastpcr(trim_name, args.fpcr, fastPCR_summary)


main()
