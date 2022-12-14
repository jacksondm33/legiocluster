#!/usr/bin/env python


"""Parse vcf output."""


import csv
import logging
import matplotlib.pyplot as plt
import numpy as np
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def translate_cigar(cigar):
    """
    Translates a CIGAR score from alphanumeric format to all-alphabetical for
      easy counting of all 'X', 'D', and 'I'
    helper function to read_VCF_file()
    param: str CIGAR = a CIGAR score in ther form '#M#X#I#D', where # is the
           count and Match, eXchanged (mismatch), Insertion, Deletion; e.g.:
           '1X2M2X1I3M'
    return: a cigar string with only alphabetical characters, e.g. 'XMMXXIMMM'
    """

    number = ''
    letter = ''
    tl_cigar = ''          # translated CIGAR string: 'XXX' instead of '3X'
    for c in cigar:        # character by character
        if c.isnumeric():  # if number
            number += c    # reassemble the number
        elif c.isalpha():  # if letter
            letter = c
            subcigar = (int(number) * letter)  # replaces '3X' with 'XXX'
            tl_cigar += subcigar       # reassemble complete CIGAR string
            number = ''    # reset
            letter = ''    # rest
    return tl_cigar


def parse_fb_vcf_line(vcf_line):
    """
    Parses a line from a FreeBayes-style VCF-file.
    helper function to read_VCF_file()
    NOTE: CHROM is the header of each contig or genome and has been edited for
          each reference genome to match the SPAdes output for contigs, which
          includes the length of the contig, coverage is set to 1.0 if unknown
    param: str vcf_line = a line of data from a VCF-file
    return: (contig name, contig position, variant POSition, REFerence
            sequence, ALTernative sequence, QUALity score, TYPE of mutation
            (snp, mnp, ins, del, or complex), and the CIGAR string (Match,
            eXchange, Deletion, Insertion))
    example line (split for readability):
        CHROM:  S-paucimobilis-NBRC-13935_NZ-BBJS00071.2_length_10828_cov_1.000
        POS:    32577
        ID:     .
        REF:    A
        ALT:    T
        QUAL:   1729.33
        FILTER: .
        INFO:   AB=0;ABP=0;AC=1;AF=1;AN=1;AO=52;CIGAR=1X;DP=52;DPB=52;DPRA=0;
                EPP=4.51363;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;
                NUMALT=1;ODDS=398.192;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;
                PRO=0;QA=1954;QR=0;RO=0;RPL=26;RPP=3.0103;RPPR=0;RPR=26;RUN=1;
                SAF=23;SAP=4.51363;SAR=29;SRF=0;SRP=0;SRR=0;TYPE=snp
        FORMAT: GT:DP:AD:RO:QR:AO:QA:GL
        OTHER:  1:52:0,52:0:0:52:1954:-176.103,0
    """

    # the entries in a line of a vcf-file
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, OTHER \
    = vcf_line.split()
    # extract the contig length from the header
    contig_name = CHROM.split('_')[1]
    contig_len = CHROM.split('_')[3]
    lo_infos = INFO.split(';')
    for info in lo_infos:
        if info.startswith('CIGAR'):
            CIGAR = info.split('=')[1]
        elif info.startswith('TYPE'):
            TYPE = info.split('=')[1]
    return (contig_name, int(contig_len), int(POS), REF, ALT, float(QUAL),
            TYPE, CIGAR)


def read_vcf_file(vcf_file):
    """
    Reads a VCF-file generated by FreeBayes, 'freebayes.vcf' and returns a
      list of the position of mutations.
    return: lo_variant_posns = list of positions of SNPs and indels, used to
            plot their distributions
    return: int V1 = mutation events (SNPs + indel events)
    return: int V2 = number of indel events, regardless of length
    return: int V3 = number of bases in indels
    return: int V4 = SNP count
    """

    # lists of numbers, where each number is a position in the reference genome
    lo_variant_posns = []  # list containing all SNPs and bases in indels
                           #  used for plotting distribution of SNPs and indels
    lo_mutation_posns = [] # list containing all SNPs and indel events
    cum_length  = 0        # cumulative position (contig lengths summed up)
    prev_contig = 'none'   # name of previous contig
    prev_length = 0        # length of previous contig
    SNP_count   = 0        # counts all individual SNPs ('X' in CIGAR)
    event_count = 0        # counts all indel events (where D = DDD = I = III)

    # extract data from the vcf-file
    with open(vcf_file, 'r') as in_file:
        for line in in_file:
            line = line.rstrip('\\n')
            if line.startswith('#'): # ignore header or comment rows
                continue
            else:
                contig_name, contig_len, POS, REF, ALT, QUAL, TYPE, CIGAR \
                = parse_fb_vcf_line(line)

                # if new contig, add the length from the previous contig
                # to the cumulative length, then update name and length
                if contig_name != prev_contig:
                    cum_length += prev_length
                    prev_length = contig_len
                    prev_contig = contig_name

                # get the position of the variant within the translated
                # CIGAR string, calculate the position of the variant within
                #  the genome, and add it to the list
                tl_cigar = translate_cigar(CIGAR)

                last_string = ''

                for index, string in enumerate(tl_cigar):

                    # variant count
                    # each SNP, ins, or del counted as one individual mutation
                    if string in ['X','I','D']:
                        lo_variant_posns.append(cum_length + POS + index)

                    # SNP count, individually
                    if string == 'X':
                        lo_mutation_posns.append(cum_length + POS + index)
                        SNP_count += 1
                    # indel event count
                    # only the first inserted or deleted base of an indel is
                    # counted if more of the same follows; e.g.: I and IIIIII
                    # are counted both only as one mutation event; two events
                    # in the same CIGAR, e.g. DDDDIIII, are counted as two
                    # for consistency with compare_SNP_files.py
                    elif string == 'I' and last_string != 'I':
                        lo_mutation_posns.append(cum_length + POS + index)
                        event_count += 1
                        last_string = 'I'
                    elif string == 'D' and last_string != 'D':
                        lo_mutation_posns.append(cum_length + POS + index)
                        event_count += 1
                        last_string = 'D'

    # these should be the same values as obtained by compare_SNP_files.py
    V1 = len(lo_mutation_posns)             # SNPs + indel events
    V2 = event_count                        # indel events
    V3 = len(lo_variant_posns) - SNP_count  # bases in indels
    V4 = SNP_count                          # SNPs

    return lo_variant_posns, V1, V2, V3, V4


def write_to_file(report_file, reference_file, isolate, to_mutations, ref_seq_len,
                  SNP_THRESHOLD):
    """
    Writes the summary data from the FreeBayes VCF file to report.txt.
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: tup to_mutations = (V1, V2, V3, V4)
    param: int ref_seq_len = length of the reference genome
    param: int SNP_THRESHOLD = number of mutation events (formerly SNPs only)
           above which an isolote will be added to the list of candidate
           reference genomes
    output: summary statistics written to the report file
    """

    # turn the int into str
    V1, V2, V3, V4 = [str(V) for V in to_mutations]
    data_str = V1 + ' (' + V2 + ', ' +  V3 + ', ' + V4 + ')'

    # write results to the report
    with open(report_file, 'a') as report:
        print('\\n\\nSNPs and INDEL events between ' + isolate\
              + ' and reference ' + str(reference_file)[:-3] + ' (FreeBayes):',
              file=report)

        print('\\nFound', data_str, 'SNPs and INDEL events compared to a ' \
              + 'reference genome of', ref_seq_len, 'bp.', file=report)
        print('(Note that the indel event count might be slightly lower in'\
              ' the SNP-matrix.)\\n', file=report)

        # if too many SNPs/INDELs, make the query it's own reference
        if to_mutations[0] >= SNP_THRESHOLD:
            print('\\nNOTE:\\nNumber of SNPs and INDEL events above threshold. '\
                  + 'Adding the isolate to the list of reference genomes.',
                  file=report)

        print('Figure: SNP/INDEL distribution', file=report)


def plot_mutation_dist(mutation_dist_file, isolate, lo_variant_posns, ref_seq_len):
    """
    Produces a histogram showing the distribution of SNPs in the genome.
      FreeBayes combines SNPs and indels into MNPs (multi-nucleotide
      polymorphisms) or complex events (composite insertion and substitution
      events), which have been unravelled into single mutations.
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: int ref_seq_len = length of the reference genome
    param: list lo_variant_posns = list of positions of SNPs and indels
    output: a histogram, 'mutation_dist.png', that shows the number of SNPs
            and number of bases in indels in intervals (bins) of 5000 bases
    """

    fig, ax = plt.subplots()
    max_x_val = int(np.ceil(ref_seq_len/500000)) + 1  # highest value on x-axis
    # print a histogram to file, were each bin covers about 5000 bp of the
    #   genome over the length of the genome
    plt.hist(lo_variant_posns, bins=int(ref_seq_len/5000), \
             range=(0, ref_seq_len))
    plt.title('SNP/indel distribution for ' + isolate + ' in 5 kb intervals')
    plt.xlabel('Position [million bases]')
    plt.ylabel('Number of mutations per 5kb')
    # change the tick labels: for a genome of 3.5 Mb, need a tick label every
    #   0.5 Mb or 500000 bp => make 8 tick labels (0, 0.5, 1, ...)
    # choose which x locations to have ticks: every 500 kb
    ax.set_xticks([500000 * x for x in range(0, max_x_val)])
    # set the labels to display at those ticks: 500 kb intervals
    ax.set_xticklabels([0.5 * x for x in range(0, max_x_val)])
    plt.savefig(mutation_dist_file)
    plt.close()


def seq_len(reference_file):
    """
    Returns the length of the reference genome in bp.
    return: int length of the sequence
    """

    seq = ''
    with open(reference_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            if not line.startswith('>'):
                seq = seq + line
    return len(seq)


def parse_vcf_output(vcf_file, reference_file, mutation_dist_file, output_file, report_file, isolate, SNP_THRESHOLD):
    """Parses the vcf file and writes the results to the report."""

    # reads the freebayes VCF file
    lo_variant_posns, V1, V2, V3, V4 = read_vcf_file(vcf_file)

    # return a tuple of values for SNPs and indel counts
    to_mutations = (V1, V2, V3, V4)

    # returns the length of the refernce genome in bp
    ref_seq_len = seq_len(reference_file)

    # makes a plot of the SNP distribution
    plot_mutation_dist(mutation_dist_file, isolate, lo_variant_posns, ref_seq_len)

    # adds text to the report file
    write_to_file(report_file, reference_file, isolate, to_mutations, ref_seq_len,
                  SNP_THRESHOLD)

    with open(output_file, 'a', newline='') as output:
        output_writer = csv.writer(output)
        output_writer.writerow(to_mutations)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(parse_vcf_output("$vcf", "$fasta", "$mutation_dist", "$output", "$report", "$meta.id", int("$snp_threshold")))
