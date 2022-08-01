#!/usr/bin/env python


"""Calculate coverage."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def calculate_perc_ge_q30(reads_file, fastqc_results, report_file):
    """
    Calculates the percent of bases with a quality score greater than Q30
    param: str reads_file = name of file with forward or reverse reads
           processed by Trimmomatic
    param: str fastqc_results = name of directory with fastqc results
    param: str report_file = output report file
    """

    n_ge_Q30 = 0  # number of bases with quality score >= Q30
    n_all = 0
    consider = False

    # extracts the data from the FastQC results file, 'fastqc_data.txt'
    with open(Path(fastqc_results) / 'fastqc_data.txt', mode='r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            if line.startswith('>>Per sequence quality scores'):
                consider = True
            elif line.startswith('#Quality'):
                continue
            elif consider and len(line.split()) == 2:
                score, n = line.split()
                # sums up the number of all bases
                n_all += float(n)
                # sums up the number of bases with a quality score >= 30
                if int(score) >= 30:
                    n_ge_Q30 += float(n)
            elif line.startswith('>>END_MODULE'):
                consider = False

    # calculates the percentage of base with quality score >= Q30
    if n_ge_Q30 > 0 and n_all > 0:
        perc_ge_Q30 = round(((n_ge_Q30 * 100) / n_all), 3)
    else:
        perc_ge_Q30 = -1

    # write to report
    with open(report_file, 'a') as report:
        print('\\nPercentage of bases with quality score >= Q30 ('\
              + str(n_ge_Q30) + ' * 100) / ' + str(n_all) + ' = '\
              + str(perc_ge_Q30) + '\\n', file=report)


def calculate_coverage(reads_file, fastqc_results, report_file, med_genome_len):
    """
    Calculate coverage.
    param: str reads_file = name of file with forward or reverse reads
           processed by Trimmomatic
    param: str fastqc_results = name of directory with fastqc results
    param: int med_genome_len = median genome length for that species
    param: str report_file = output report file
    """

    # extracts the data from the FastQC results file, 'fastqc_data.txt'
    with open(Path(fastqc_results) / 'fastqc_data.txt', mode='r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            if line.startswith('Total Sequences'):
                total_seqs = int(line.split()[-1])
            elif line.startswith('Sequence length'):
                seq_range = line.split()[-1]
                if '-' in seq_range:
                    max_seq_len = int(seq_range.split('-')[1])
                else:
                    max_seq_len = int(seq_range)

    # calculation based on PulseNet SOP
    coverage = round((total_seqs * max_seq_len * 2) / med_genome_len, 3)

    with open(report_file, 'a') as report:
        print('\\nCoverage: (' + str(total_seqs) + ' * ' + str(max_seq_len)\
              + ' * 2) / ' + str(med_genome_len) + ' = ' + str(coverage),
              file=report)

    calculate_perc_ge_q30(reads_file, fastqc_results, report_file)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(calculate_coverage("$reads", "$fastqc_results", "$report", int("$med_genome_len")))
