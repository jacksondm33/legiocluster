#!/usr/bin/env python


"""Make SNP cons fa."""


import csv
import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def read_ref_file(reference_file):
    """
    Reads a fasta file with the sequence of the reference genome (consisting
      of one or more contigs) and returns a list of headers and sequences.
    return: list lo_contigs = list of headers and sequences, e.g.:
        [[NODE_1_length_6526_cov_26.4, 'ACTTGTACTAATTGGCTGATTGTTGACATAA...'],
         [NODE_2_length_5226_cov_30.2, 'GTACTAATTGGCTGATTGTCTTCCAACATAA...'],
         ...]
    """

    lo_contigs = []
    contig = ''
    seq = ''

    with open(reference_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            # extracts the contig name and adds a new list to lo_contigs that
            # includes the new contig name and '' (default) for the sequence
            if line.startswith('>'):
                contig = line.split()[0][1:]
                lo_contigs.append([contig, seq])
            # takes a line representing a sequence and adds it to the last
            # sequence in the list of lists
            else:
                lo_contigs[-1][1] += line

    return lo_contigs


def convert_to_base_list(lo_contigs):
    """
    Takes a list of [[header_1, sequence_1], ...] and converts it to a list
      of [[contig, position, base], ...]
    param: list lo_contigs = [[header_1, sequence_1], ...]
    return: list lo_bases = list of [contig-name, position, base] for all
            contig sequences
    """

    lo_bases = []
    for contig in lo_contigs:
        for i, base in enumerate(contig[1]):
            lo_bases.append([contig[0], str(i+1), base])

    return lo_bases


def write_ref_seq(snp_cons_file, reference, lo_bases):
    """
    Takes a list of [[contig-name, position, base], ...] and writes the base
      to file such that each base occupies it's own line.
    param: list lo_bases = list of [contig-name, position, base] for the entire
           reference sequence
    output: a '_SNP_cons.txt' file for the reference strain
    """

    with open(snp_cons_file, 'w') as outfile:
        print('# sequence for reference ' + reference, file=outfile)
        for base in lo_bases:
            print(base[2], file=outfile)


def make_snp_cons_fa(reference_file, snp_cons_file, reference):
    """
    main function
    param: str isolate = isolate name, e.g.: 'IDR001234'
    output: a '_SNP_cons.txt' file added to the /VCF_folder
    """

    # returns a list of contigs [header, sequence] for the reference
    lo_contigs = read_ref_file(reference_file)

    # returns list of [contig, posn, base] for the reference
    lo_bases = convert_to_base_list(lo_contigs)

    # generates a <ref>_SNP_cons.txt file
    write_ref_seq(snp_cons_file, reference, lo_bases)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(make_snp_cons_fa("$fasta", "$snp_cons", "$meta.ref"))

