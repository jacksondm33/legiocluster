#!/usr/bin/env python


"""Parse kraken output."""


import logging
import numpy as np
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def sort_kraken_res(kraken_file, GENUS):
    """
    Sort the headers of contigs into two lists, depending on if they match
      the GENUS or not.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str GENUS = taxonimic genus to which the isolate belongs
    return: list lo_headers = headers of contigs matching the genus
    return: list lo_wrong_genera = full kraken taxonomic information for
            contigs not matching the genus
    """

    lo_headers = []
    lo_wrong_genera = []
    with open(kraken_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            if GENUS in line:
                header = line.split()[0]
                lo_headers.append(header)
            else:
                lo_wrong_genera.append(line)

    return lo_headers, lo_wrong_genera


def write_good_contigs_to_fasta(contigs_file, good_contigs_file, isolate, lo_headers,
                                MIN_COV=3.0, MIN_LEN=250):
    """
    Copies contigs from a 'SPAdes_contigs.fa' file if they meet criteria
      for Genus, minimum length, and coverage. Returns count of all contigs
      processed, list of lengths, and list of contigs that are of good quality
      but wrong genus.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: list lo_headers = headers of contigs matching the genus
    param: int MIN_COV = minimum read coverage, default 3.0 fold
    param: int MIN_LEN = minimum contig length, default 250 bp
    return: int count = count of all contigs
    return: list lo_lens = list of the lengths of contigs that are above
            MIN_COV and above MIN_LEN
    return: list lo_bad_contigs = list of headers from contigs that passed
            QC requirements, but that are of the wrong genus
    output: isolate + '_cc.fasta' file
    """

    write_to_file = False
    count = 0
    lo_lens = []
    lo_bad_contigs = []

    with open(contigs_file, 'r') as infile:
        with open(good_contigs_file, 'a') as outfile:

            # if header is on the list, write it and all sequence lines
            # to the new fasta file, else ignore
            for line in infile:
                line = line.rstrip('\\n')

                # check the header
                if line.startswith('>'):
                    count += 1
                    # check that contig meets minimal criteria
                    n, contig, l, length, c, coverage = line.split('_')
                    meets_QC = False
                    if int(length) >= MIN_LEN and float(coverage) >= MIN_COV:
                        meets_QC = True
                        lo_lens.append(int(length))

                    # write to file if contig is of right Genus and meets QC
                    # each line starts with '>', the headers in the list do not
                    if line[1:] in lo_headers and meets_QC:
                        write_to_file = True
                        print(line, file=outfile)
                    # keep track of contigs that are of wrong species but that
                    # are large and of sufficient coverage
                    # 1000bp and 7.5x is used by the Lpn-pipeline to filter
                    # out low quality contigs
                    elif line[1:] not in lo_headers and int(length) >= MIN_LEN\
                    and float(coverage) >= MIN_COV:
                        lo_bad_contigs.append(line)
                        write_to_file = False
                    else:
                        write_to_file = False

                # write sequence to file only if header checks out
                else:
                    if write_to_file:
                        print(line, file=outfile)

    return count, lo_lens, lo_bad_contigs


def write_bad_contigs_to_fasta(contigs_file, bad_contigs_file, lo_bad_contigs):
    """
    Copies contigs from a 'SPAdes_contigs.fa' to a new fasta file, but
      only those contigs that were of the wrong genus: allows easy search
      with BLASTN (e.g. to look for plasmids).
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: list lo_bad_contigs = list of headers from contigs that passed
           QC requirements, but that are of the wrong genus
    output: 'wrong_genus_contigs.fasta' file
    """

    write_to_file = False

    with open(contigs_file, 'r') as infile:
        with open(bad_contigs_file, 'a') as outfile:

            # if header that's on the list, write it and all sequence lines
            # to the new fasta file, else ignore
            for line in infile:
                line = line.rstrip('\\n')

                # check the header
                if line.startswith('>'):
                    # write to file if contig is of wrong Genus
                    if line in lo_bad_contigs:
                        write_to_file = True
                        print(line, file=outfile)
                    else:
                        write_to_file = False

                # write sequence to file only if header checks out
                else:
                    if write_to_file:
                        print(line, file=outfile)


def combine_failed_contigs(lo_wrong_genera, lo_bad_contigs):
    """
    Combines two lists to include the data in the report.
    param: list lo_wrong_genera = Kraken data for all contigs of wrong genus
    param: list lo_bad_contigs = headers from SPAdes for contigs with wrong
           genus but good data quality
    return: list lo_failed_contigs = list of suspicious contigs
    """

    lo_failed_contigs = []

    for item in lo_wrong_genera:
        header = '>' + item.split()[0]
        if header in lo_bad_contigs:
            lo_failed_contigs.append(item)

    return lo_failed_contigs


def calc_N50(lo_lens):
    """
    Takes a list of contig lengths, sorted in descending order, and
      returns the contig length corresponding to the N50.
    helper function to write_report()
    param: lo_lens = list of the lengths of contigs that are above
            MIN_COV and above MIN_LEN
    return: int or float = N50
    """

    if lo_lens != []:
        slo_length = sorted(lo_lens, reverse=True)
        total_len = sum(slo_length)
        cum_len = 0
        i = 0
        while cum_len <= total_len/2:
            cum_len += slo_length[i]
            i += 1
        # the last contig to be cum_len <= total_len/2 is the N50
        return slo_length[i-1]
    else:
        return -1


def write_report(report_file, isolate, lo_lens, count, lo_failed_contigs):
    """
    Writes some summary statistics to the report.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: list lo_lens = list of the lengths of contigs that are above
           MIN_COV and above MIN_LEN
    param: int count = count of all contigs assembled by SPAdes
    param: list lo_failed_contigs = list of suspicious contigs
    output: text added to report.txt
    """

    with open(report_file, 'a') as outfile:
        print('\\n\\nContigs overview:\\n(after Kraken species ID and removal'\
              + ' of questionable contigs)', file=outfile)
        print('- Folder:                  {:21s}'.format(work_dir),
              file=outfile)
        print('- Isolate:                 {:21s}'.format(isolate),
              file=outfile)
        print('- Input contigs [n]:       {:12,d}'.format(count),
              file=outfile)
        print('- Acceptable contigs [n]:  {:12,d}'.format(len(lo_lens)),
              file=outfile)
        print('- Total length [bp]:       {:12,d}'.format(sum(lo_lens)),
              file=outfile)
        print('- N50 [bp]:                {:12,d}'.format(calc_N50(lo_lens)),
              file=outfile)
        print('- Largest contig [bp]:     {:12,d}'.format(max(lo_lens)),
              file=outfile)
        print('- Mean contig [bp]:        {:14,.1f}'.format(np.mean(lo_lens)),
              file=outfile)
        print('- Median contig [bp]:      {:14,.1f}'.format(np.median(lo_lens)),
              file=outfile)
        print('- Smallest contig [bp]:    {:12,d}'.format(min(lo_lens)),
              file=outfile)
        if lo_failed_contigs != []:
            print('- Contigs of concern (wrong genus, but good quality data):',
                  file=outfile)
            for failed_contig in lo_failed_contigs:
                if ';' in failed_contig:
                    print('\t', failed_contig.split('\t')[0], '\t',
                          failed_contig.split(';')[-1], file=outfile)
                else:
                    print('\t', failed_contig, file=outfile)
        print('\\n', file=outfile)


def parse_kraken_output(kraken_file, contigs_file, good_contigs_file, bad_contigs_file, report_file, isolate, genus):
    """
    Main function: run Kraken on SPAdes output files
    param: list lo_phylo_tree_data = list of:
           str sp_abbr = three letter species abbreviation, e.g.: 'Lpn'
           str isolate = isolate name, e.g.: 'IDR001234'
           str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
           str ref = name of a reference strain (not used here)
    output: text added to report.txt
    output: Kraken data files
    output: several FASTA files
    """

    lo_headers, lo_wrong_genera = sort_kraken_res(kraken_file, genus)

    logger.info('\\n', len(lo_headers), 'contigs matched the Genus')
    logger.info(len(lo_wrong_genera), 'contigs did not match the Genus')

    # correct genus, good quality contigs: write to fasta file
    count, lo_lens, lo_bad_contigs = \
    write_good_contigs_to_fasta(contigs_file, good_contigs_file, isolate, lo_headers)

    # wrong genus: write to fasta file for easy blast search
    write_bad_contigs_to_fasta(contigs_file, bad_contigs_file, lo_bad_contigs)

    # combine info from lo_wrong_genera and lo_bad_contigs into one:
    # contigs of poor quality or wrong genus
    lo_failed_contigs = combine_failed_contigs(lo_wrong_genera, lo_bad_contigs)

    # writes some statistics to file
    write_report(report_file, isolate, lo_lens, count, lo_failed_contigs)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(parse_kraken_output("$kraken", "$contigs", "$good_contigs", "$bad_contigs", "$report", "$meta.id", "$genus"))
