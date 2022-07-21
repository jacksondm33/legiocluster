#!/usr/bin/env python


"""Parse qualimap output."""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def parse_qualimap_output(qualimap_file, report_file):
    """
    Extracts results from the Qualimap output and writes it to report
    output: text added to report
    """

    # opens the report from Qualimap and writes selected lines to report
    with open(qualimap_file, 'r') as in_file:

        # opens the report file
        with open(report_file, 'a') as report:

            # adds a header to report
            print('\n\nMapping quality check (Qualimap results):', \
                  file=report)

            # writing specific data to the report
            for line in in_file:
                line = line.rstrip('\n')
                if line.startswith('     number of bases')\
                or line.startswith('     number of contigs')\
                or line.startswith('     number of reads')\
                or line.startswith('     number of mapped reads')\
                or line.startswith('     number of mapped bases')\
                or line.startswith('     mean mapping quality'):
                    print(line.replace('     ',''), file=report)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--qualimap-file",
        metavar="QUALIMAP_FILE",
        type=Path,
        help="Input qualimap file",
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
    if not args.qualimap_file.is_file():
        logger.error(f"The given input file {args.qualimap_file} was not found!")
        sys.exit(2)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    parse_qualimap_output(args.qualimap_file, args.report_file)


if __name__ == "__main__":
    sys.exit(main())
