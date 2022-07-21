#!/usr/bin/env python


"""Calculate coverage."""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def calculate_coverage(reads_file, fastqc_results, med_genome_len, report_file):
    """
    Calculate coverage.
    param: str reads_file = name of file with forward or reverse reads
           processed by Trimmomatic
    param: str fastqc_results = name of directory with fastqc results
    param: int med_genome_len = median genome length for that species
    param: str report_file = output report file
    """

    # extracts the data from the FastQC results file, 'fastqc_data.txt'
    with open(fastqc_results / 'fastqc_data.txt', mode='r') as infile:
        for line in infile:
            line = line.rstrip('\n')
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
        print('\nCoverage: (' + str(total_seqs) + ' * ' + str(max_seq_len)\
              + ' * 2) / ' + str(med_genome_len) + ' = ' + str(coverage),
              file=report)


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
    with open(fastqc_results / 'fastqc_data.txt', mode='r') as infile:
        for line in infile:
            line = line.rstrip('\n')
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
        print('\nPercentage of bases with quality score >= Q30 ('\
              + str(n_ge_Q30) + ' * 100) / ' + str(n_all) + ' = '\
              + str(perc_ge_Q30) + '\n', file=report)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fastqc-results",
        metavar="FASTQC_RESULTS",
        type=Path,
        help="Fastqc results directory",
        required=True,
    )
    parser.add_argument(
        "--med-genome-len",
        metavar="MED_GENOME_LEN",
        type=int,
        help="Median genome length",
        required=True,
    )
    parser.add_argument(
        "--reads-file",
        metavar="READS_FILE",
        type=Path,
        help="Reads file",
        required=True,
    )
    parser.add_argument(
        "--report-file",
        metavar="REPORT_FILE",
        type=Path,
        help="Output report file",
        required=True,
    )
    parser.add_argument(
        "--log-level",
        metavar="LOG_LEVEL",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        help="The desired log level",
        default="INFO",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.reads_file.is_file():
        logger.error(f"The given input file {args.reads_file} was not found!")
        sys.exit(2)
    if not args.fastqc_results.is_dir():
        logger.error(f"The given input directory {args.fastqc_results} was not found!")
        sys.exit(2)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    calculate_coverage(args.reads_file, args.fastqc_results, args.med_genome_len, args.report_file)
    calculate_perc_ge_q30(args.reads_file, args.fastqc_results, args.report_file)


if __name__ == "__main__":
    sys.exit(main())
