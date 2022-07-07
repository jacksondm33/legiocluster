#!/usr/bin/env python


"""Extract fastqc results."""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def extract_fastqc_results(reads_file, fastqc_results, report_file):
    """
    Extracts basic statistics and summary results from the 'fastqc_data.txt'
      and 'summary.txt' files and writes them to the report.
    param: str reads_file = (path and) name of file with the raw forward or
           reverse reads, e.g.: "IDR200001234.fastq.gz"
    param: str fastqc_results = name of directory with fastqc results
    param: str report_file = output report file
    output: writes basic statistics to the report file
    """

    # extract selected data from the 'fastqc_data.txt' file
    with open(report_file, 'a') as report:
        print('\nRead quality control (FastQC results):', file=report)
        # Original name of the file with the raw reads
        print('Results for processed reads from:', reads_file, file=report)
        with open(fastqc_results / 'fastqc_data.txt', mode='r') as infile_1:
            for line in infile_1:
                line = line.rstrip('\n')
                if line.startswith('Filename') \
                or line.startswith('Total Sequences') \
                or line.startswith('Sequences flagged') \
                or line.startswith('Sequence length') \
                or line.startswith('%GC'):
                    print(line, file=report)

    # extract all data from the 'summary.txt' file
    lo_qc_results = []
    with open(report_file, 'a') as report:
        with open(fastqc_results / 'summary.txt', mode='r') as infile_2:
            for line in infile_2:
                line = line.rstrip('\n')
                # remove read_file name in each line
                qc_result, what = line.split('	')[:2]
                lo_qc_results.append(qc_result)
                print(qc_result + '   ' + what, file=report)


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
        help="The desired log level (default WARNING).",
        default="WARNING",
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
    extract_fastqc_results(args.reads_file, args.fastqc_results, args.report_file)


if __name__ == "__main__":
    sys.exit(main())
