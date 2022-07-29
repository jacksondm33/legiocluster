#!/usr/bin/env python


"""Compare SNPs."""


import csv
import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def get_seq_data(snp_cons_file):
    """
    Returns a list where each genome position is a list item; compared to the
      reference genome, an item can be a single base, multiples bases (INS),
      '-' (DEL), 'N' (unmapped), or 'n' (ambiguous)
    param: str path_file = path and full name of the '_SNP_cons.txt' file
    return: lo_data = genome sequence as a list, e.g.:
            ['A','A','T','G','CAGGA','C','-','-','-','T','N','n','G','A', ...]
    """

    lo_bases = []
    with open(snp_cons_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('''\n''')
            if not line.startswith('#'):
                lo_bases.append(line)
    return lo_bases


def get_base_count(str1, str2):
    """
    Counts the number of mutations if two strains have insertions relative to
      the reference genome and one of the insertions is larger than the other.
      Returns a close approximation of the number of bases that are different;
      getting an exact number would require something like Smith-Waterman,
      which is too complex for this application. Alternative: counting bases.
      Start at str[1:] because str[0] is the base before the actual INS.
    param: str str1 = a string of bases, where len(str1) >= len (str2)
    param: str str2 = a string of bases
    return: number of bases that are different between the strings
    helper function to compare_two_genomes()
    """

    count = abs(str1[1:].count('A') - str2[1:].count('A')) \
          + abs(str1[1:].count('C') - str2[1:].count('C')) \
          + abs(str1[1:].count('G') - str2[1:].count('G')) \
          + abs(str1[1:].count('T') - str2[1:].count('T'))
    # at most, there can be only len(str1)-1 changes: 'AGG' versus 'AT' is
    # one insertion and one mismatch (relative to the reference, they are
    # two insertions of 2 and 1 bases, respectively)
    if count >= (len(str1) - 1):
        count = len(str1) - 1
    return count


def get_pairwise_count(str1, str2):
    """
    Counts the number of mutations if two strains that have insertions relative
      to the reference genome and both of the insertions are of equal length.
      Returns an approximation of the number of bases that are different after
      a side-by-side comparison.
      helper function to compare_two_genomes()
    param: str str1 = a string of bases
    param: str str2 = a string of bases, where one str is larger than the other
    return: int count of number of bases that are different
    """

    count = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            count += 1
    return count


def compare_two_genomes(lo_fst, lo_snd):
    """
    Compares two lists with sequence data, where each item is either a base,
      an insertion (MI+), a deletion ('-'), ambiguous ('n'), or unmapped ('N').
      The indexing is the same as that of the reference genome, an INS is
      listed as MI+ (e.g. 'ATT'), hence no disruption of the index.
    param: list lo_fst = list of bases for the first genome
    param: list lo_snd = list of bases for the second genome
    return: variant_count = the number of SNPs, INS, and DEL, while ignoring
            'n' or 'N'
    """

    SNP_count        = 0
    INDEL_base_count = 0

    # one list element at a time
    for i in range(len(lo_fst)):
        # extract the bases, INS, or DEL after removing residual white spaces
        fst = lo_fst[i].split()[0]
        snd = lo_snd[i].split()[0]
        # no difference (same base or same INS, DEL, 'n', 'N')
        if fst == snd:
            continue
        # ignore ambiguous or unmapped bases
        elif (fst in ['n','N']) or (snd in ['n','N']):
            continue
        # the two items are different: must be a SNP, INS or DEL
        elif fst != snd:
            # a simple, single SNP
            if fst in ['A','C','G','T'] and snd in ['A','C','G','T']:
                SNP_count += 1
            else:
                # a DEL in one, but not the other, sequence
                if fst == '-' or snd == '-':
                    INDEL_base_count += 1
                # both have insertions and the first INS is larger
                elif len(fst) > len(snd):
                    INDEL_base_count += get_base_count(fst, snd)
                # both have insertions and the second INS is larger
                elif len(fst) < len(snd):
                    INDEL_base_count += get_base_count(snd, fst)
                # both have insertions of the same length, but are
                # different from each other
                elif len(fst) == len(snd):
                    INDEL_base_count += get_pairwise_count(fst, snd)

    return SNP_count, INDEL_base_count


def get_indels(lo_bases):
    """
    Converts a list of bases into two lists: one for indels, one for
      ambiguous bases.
    param: list lo_bases = list of single bases, multiple bases (insertions),
           '-' (deletions), 'n' or 'N' (ambiguous); Note that the position of
           item corresponds to it's position in the reference genome
    return: list lo_indels = list of [position, base] for indels, where
            each item is either one or more '-' or two or more bases. e.g.:
            [[3018, 'CGATTT'], [29548, '-'], [117852, '-------------'], ... ]
    return: list lo_nNs = list of positions that are either 'n' or 'N', e.g.:
            [182, 248, 611, 761, 963, 1122, 1579, ...]
    helper function to indel_comp_manager()
    """

    lo_indels = []  # (position, bases) of insertions (MI+) or deletions (D+)
    lo_nNs    = []  # position of 'n' or 'N'

    posn = 0         # base count = position relative to the reference genome
    last_posn = -1   # keeps track if deletions are consecutive

    for line in lo_bases:

        posn += 1  # update position

        # ambiguous bases or gaps: add posn to list
        if line in ['n','N']:
            lo_nNs.append(posn)

        # insertions, which have the form MI+, where M is the same base
        # found in the reference and I+ represents >= 1 inserted bases
        elif len(line) > 1:
            lo_indels.append([posn, line])

        # deletion (single character, '-', combine to multi-deletion,
        # '------', as applicable):
        elif line == '-':
            # extend existing deletion if the posn is next to a
            # previous deletion
            if posn == last_posn + 1:
                lo_indels[-1][1] += line
                last_posn = posn
            # a new deletion
            else:
                lo_indels.append([posn, line])
                last_posn = posn

    return lo_indels, lo_nNs


def compare_events(isolate1, isolate2, lo_indels1, lo_indels2, lo_nNs2):
    """
    Compares indel sequences of two isolates. Need to run this function twice:
      once for A versus B, then B versus A.
      Note 1: deletions are one or more per position (e.g.: '1234 -') while
      insertions are two or more bases per position (e.g.: '5678 ACGT', where
      'A' isthe match to the corresponding base in the reference sequence
      Note 2: insertions in isolate1 and isolate2 have to be identical, else,
      they will be counted as separate events, e.g.: ('1234 atttttttt') and
      ('1234 atttgtttt') will be considered as two events.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate1 = name of the first isolate (= isolate_A or isolate_B)
    param: str isolate2 = name of the second isolate (= isolate_B or isolate_A)
    param: list lo_indels1 = [position, insertion or deletion] for isolate1
           e.g.: [[4, 'CAGA'], [6, '---'], [10, '-'], [11, 'GA']]
    param: list lo_indels2 = [position, insertion or deletion] for isolate2
    param: list lo_nNs2 = ambiguous bases in isolate2
    param: bool write_ic_file = if True, write results to file (for QC and
           development)
    helper function to indel_comp_manager()
    """

    # an "event" is a single insertion or deletion of one or more bases
    lo_events = []

    # event_no is the count of indels
    for event_no, indel1 in enumerate(lo_indels1):
        # the indel is unique to that isolate
        if indel1 not in lo_indels2:
            # split up the (posn, base) tuple
            posn1, base1 = indel1
            # if it's a deletion: e.g.: (1234, '-')
            if '-' in base1:
                # checks if any deleted base in isolate1 is ambiguous in
                # isolate2 by checking if any position of the deletion in iso1
                # is also on the list of ns or Ns from iso2
                for i in range(posn1, posn1 + len(base1)):
                    # check if the base at that position is ambiguous
                    if i not in lo_nNs2:
                        lo_events.append(event_no)

            # an insertion: e.g.: '567 ACG' where 'CG' is inserted, but not 'A'
            else:
                # check if the base at that position is ambiguous
                if posn1 not in lo_nNs2:
                    lo_events.append(event_no)

    # removing duplicates and reporting the final number
    no_events = len(set(lo_events))

    return no_events


def indel_comp_manager(isolate_A, isolate_B, lo_bases_A, lo_bases_B):
    """
    Mananges the comparison of two SNP_cons.txt files with each other to
      determine the number of indel events.
    Note: indels versus 'n' or 'N' are not counted; indels that are similar,
      but of different length, are counted as separate events
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate_A = name of one isolate
    param: str isolate_B = name of another isolate
    param: list lo_bases_A = list of 'ACGT', 'ACGT'+ (insertion), '-', 'n',
           'N' for isolate_A
    param: list lo_bases_B = list of 'ACGT', 'ACGT'+ (insertion), '-', 'n',
           'N' for isolate_B
    return: int = the sum of the number of indel events for A:B and B:A
    """

    # convert list of bases into lo_indels and lo_nNs
    lo_indels_A, lo_nNs_A = get_indels(lo_bases_A)
    lo_indels_B, lo_nNs_B = get_indels(lo_bases_B)

    # compare A versus B, then B versus A
    no_events_A = compare_events(isolate_A, isolate_B, lo_indels_A,
                                 lo_indels_B, lo_nNs_B)
    no_events_B = compare_events(isolate_B, isolate_A, lo_indels_B,
                                 lo_indels_A, lo_nNs_A)

    # return sum of indel events
    return no_events_A + no_events_B


def compare_snps(snp_cons_file, lo_cluster_snp_cons, pairwise_diffs_file):
    """
    Organizes the pairwise comparison of '_SNP_cons.txt' files, one per
      isolate in a cluster. Isolate_A, the query, will be compared to all other
      files in lo_files.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str VCF_path_ref = cluster-specific folder with mutation data
    param: str isolate_A = name of the query isolate, user supplied
    param: list lo_files = list of the names of '_SNP_cons.txt' files for
           a cluster of isolates
    return: list lo_pairwise_diffs = [(G1, G2, V1, V2, V3, V4), ...], where
            G1 and G2 are the names of the two genomes,
            V1 = indel events + SNPs = mutation events,
            V2 = number of indel events,
            V3 = count of bases in indels,
            V4 = SNPs
            e.g.:  [('iso1', 'iso2', 19, 13, 41, 6),
                    ('iso2', 'iso1', 19, 13, 41, 6)]
    """

    lo_pairwise_diffs = []

    isolate_A = snp_cons_file.split('_SNP_cons.txt')[0]

    # extract the genome sequence
    lo_bases_A = get_seq_data(snp_cons_file)

    for cluster_snp_cons_file in lo_cluster_snp_cons:

        # name of the comparator genome
        isolate_B = cluster_snp_cons_file.split('_SNP_cons.txt')[0]

        # don't compare the isolate to itself
        if isolate_B != isolate_A:
            # extract the genome sequence
            lo_bases_B = get_seq_data(cluster_snp_cons_file)
            # check that both lists have same number of loci
            if len(lo_bases_A) != len(lo_bases_B):
                break
            # count of SNPs and bases in INDELs
            cSNP, cINDEL = compare_two_genomes(lo_bases_A, lo_bases_B)
            # get the number of indel events
            no_events = indel_comp_manager(isolate_A, isolate_B,
                                           lo_bases_A, lo_bases_B)
            # add (G1, G2, V1, V2, V3, V4) to the list
            lo_pairwise_diffs.append((isolate_A, isolate_B, no_events + cSNP,
                                     no_events, cINDEL, cSNP))

    with open(pairwise_diffs_file, 'w', newline='') as outfile:  # write to csv
        csv_writer = csv.writer(outfile)
        for row in lo_pairwise_diffs:
            csv_writer.writerow(row)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(compare_snps("$snp_cons", "$cluster_snp_cons".split(), "$pairwise_diffs"))
