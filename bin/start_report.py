#!/usr/bin/env python


"""Starts the report."""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def start_report(sp_abbr, isolate, reads, report_file):
    """
    Creates a report file and writes the header.
    output: a report file with header information
    """

    lo_data = [('Species:\t\t', sp_abbr),
               ('Isolate name:\t\t', isolate),
               ('Forward reads:\t\t', reads[0]),
               ('Reverse reads:\t\t', reads[1])]

    with open(report_file, 'a') as report:
        print('REPORT\n', file=report)
        for datum in lo_data:
            print(datum[0], datum[1], file=report)
        print('\n', file=report)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--isolate",
        metavar="ISOLATE",
        help="Isolate",
    )
    parser.add_argument(
        "--reads",
        metavar="READS",
        nargs=2,
        type=Path,
        help="Reads files",
    )
    parser.add_argument(
        "--report-file",
        metavar="REPORT_FILE",
        type=Path,
        help="Output report file",
    )
    parser.add_argument(
        "--sp-abbr",
        metavar="SP_ABBR",
        help="Species abbreviation",
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
    if not args.reads[0].is_file():
        logger.error(f"The given input file {args.reads[0]} was not found!")
        sys.exit(2)
    if not args.reads[1].is_file():
        logger.error(f"The given input file {args.reads[1]} was not found!")
        sys.exit(2)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    start_report(args.sp_abbr, args.isolate, args.reads, args.report_file)


if __name__ == "__main__":
    sys.exit(main())
