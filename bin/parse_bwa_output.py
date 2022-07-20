#!/usr/bin/env python


"""Parse BWA output."""


import argparse
import csv
import logging
import sys
from pathlib import Path
from statistics import fmean, median, pstdev


logger = logging.getLogger()


def calculate_frag_len(sam_file, report_file):
    """
    Extracts the fragment lengths from a SAM file and writes the mean and
    other stats to the report.
    param: str sam_file = input sam file
    output: text added to report
    """

    lo_frag_lens = []

    with open(sam_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            if line.startswith('@'):
                continue
            else:
                frag_len = line.split()[8]
                if frag_len.startswith('-'):
                    continue
                else:
                    lo_frag_lens.append(int(frag_len))

    with open(report_file, 'a') as report:
        print('\n\nGenomic fragments:', file=report)
        print('Smallest fragment:\t', min(lo_frag_lens), file=report)
        print('Mean length:\t\t', round(fmean(lo_frag_lens), 2), file=report)
        print('S.D.:\t\t\t', round(pstdev(lo_frag_lens), 2), file=report)
        print('median:\t\t', median(lo_frag_lens), file=report)
        print('Largest fragment:\t', max(lo_frag_lens), file=report)


def parse_bwa_output(flagstat_file, idxstats_file, output_file, reference_file, sam_file, report_file, MAPPED_THRESHOLD):
    """Parse flagstat file."""

    with open(report_file, 'a') as report:
        print('\n\nMapping the query against strain ' + str(reference_file)\
              + ' (BWA MEM):', file=report)

        with open(flagstat_file, 'r') as flagstat:
            flagstat_data = ''.join(flagstat.readlines())
            percent_mapped = float(flagstat_data.split('mapped (')[1].split('%')[0])

        print('\nAlignment QC (Samtools flagstat):', file=report)
        print(flagstat_data, file=report)
        print('\n\nPercentage of mapped reads:', percent_mapped, file=report)

        if percent_mapped <= MAPPED_THRESHOLD:
            print('\nNOTE:\nPercentage of mapped reads below threshold.\n'\
                + 'Adding the isolate to the list of candidate reference '\
                + 'genomes.',
                  file=report)

        with open(idxstats_file, 'r') as idxstats:
            idxstats_data = ''.join(idxstats.readlines())

        print('\n\nAlignment QC (Samtools idxstats):', file=report)
        print('ref_fa_file\tlen\tmapped\tunmapped', file=report)
        print(idxstats_data, file=report)

    calculate_frag_len(sam_file, report_file)

    with open(output_file, 'a', newline='') as output:
        output_writer = csv.writer(output)
        output_writer.writerow([percent_mapped])


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--flagstat-file",
        metavar="FLAGSTAT_FILE",
        type=Path,
        help="Input flagstat file",
        required=True,
    )
    parser.add_argument(
        "--idxstats-file",
        metavar="IDXSTATS_FILE",
        type=Path,
        help="Input idxstats file",
        required=True,
    )
    parser.add_argument(
        "--output-file",
        metavar="OUTPUT_FILE",
        type=Path,
        help="Output file",
        required=True,
    )
    parser.add_argument(
        "--reference-file",
        metavar="REFERENCE_FILE",
        type=Path,
        help="Input reference file",
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
        "--sam-file",
        metavar="SAM_FILE",
        type=Path,
        help="Input sam file",
        required=True,
    )
    parser.add_argument(
        "--log-level",
        metavar="LOG_LEVEL",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        help="The desired log level (default WARNING).",
        default="WARNING",
    )
    parser.add_argument(
        "--mapped-threshold",
        metavar="MAPPED_THRESHOLD",
        type=float,
        help="Mapped threshold",
        default=0.0,
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.flagstat_file.is_file():
        logger.error(f"The given input file {args.flagstat_file} was not found!")
        sys.exit(2)
    if not args.idxstats_file.is_file():
        logger.error(f"The given input file {args.idxstats_file} was not found!")
        sys.exit(2)
    if not args.reference_file.is_file():
        logger.error(f"The given input file {args.reference_file} was not found!")
        sys.exit(2)
    if not args.sam_file.is_file():
        logger.error(f"The given input file {args.sam_file} was not found!")
        sys.exit(2)
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    parse_bwa_output(args.flagstat_file, args.idxstats_file, args.output_file, args.reference_file, args.sam_file, args.report_file, args.mapped_threshold)


if __name__ == "__main__":
    sys.exit(main())
