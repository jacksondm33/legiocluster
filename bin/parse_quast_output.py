#!/usr/bin/env python


"""Parse quast output."""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def parse_quast_output(contigs_file, quast_report_file, report_file, CONTIG_THRESHOLD):
    """
    Parses the Quast output and writes the results to the report.
    param: str contigs_file = name of a sequence file to be QC'd
    output: info from the Quast file will be added to 'report.txt',
            issues a Warning if too many contigs
    return int contigs = the number of all contigs
    """

    contigs = 0

    # adds a header to '_report.txt'
    with open(report_file, 'a') as report:
        print('\n\nAssembly quality check (Quast results) for', \
              str(contigs_file) + ':', file=report)

        # opens the report from Quast and writes data to report file
        with open(quast_report_file, mode='r') as in_file:
            for line in in_file:
                line = line.rstrip('\n')
                print(line, file=report)
                # extracts the number of contigs and warns if too many
                if line.startswith('# contigs (>= 0 bp)'):
                    contigs = int(line.split()[-1])

        if contigs > CONTIG_THRESHOLD:
            print('\nWARNING:', file=report)
            print(contigs, 'contigs', file=report)
            print('THE SAMPLE MIGHT BE CONTAMINATED\n\n', file=report)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contigs-file",
        metavar="CONTIGS_FILE",
        type=Path,
        help="Input contigs file",
        required=True,
    )
    parser.add_argument(
        "--quast-report-file",
        metavar="QUAST_REPORT_FILE",
        type=Path,
        help="Input quast report file",
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
    parser.add_argument(
        "--contig-threshold",
        metavar="CONTIG_THRESHOLD",
        type=int,
        help="Contig threshold",
        default=300,
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.contigs_file.is_file():
        logger.error(f"The given input file {args.contigs_file} was not found!")
        sys.exit(2)
    if not args.quast_report_file.is_file():
        logger.error(f"The given input file {args.quast_report_file} was not found!")
        sys.exit(2)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    parse_quast_output(args.contigs_file, args.quast_report_file, args.report_file, args.contig_threshold)


if __name__ == "__main__":
    sys.exit(main())
