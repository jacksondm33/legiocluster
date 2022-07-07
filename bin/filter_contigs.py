#!/usr/bin/env python


"""Filter contigs."""


import argparse
import logging
import random
import sys
from pathlib import Path


logger = logging.getLogger()


def filter_contigs(contigs_in, contigs_out, MIN_CONTIG_LEN, MIN_CONTIG_COV):
    """
    Goes through the spades contigs file and writes the contigs above
      MIN_CONTIG_LEN and MIN_CONTIG_COV to a new file.
    param: str contigs_in = input contigs file
    param: str contigs_out = output contigs file
    param: int MIN_CONTIG_LEN = min lengths of contig in bases
    param: float MIN_CONTIG_COV = minimal contig coverage per base
    output: fasta file named after the isolate with contigs that meet the
            selection criteria
    """

    with open(contigs_in, 'r') as in_file:
        with open(contigs_out, 'a') as out_file:
            for line in in_file:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    data = line.split('_')
                    length = int(data[3])
                    cov = float(data[5])
                    if length > MIN_CONTIG_LEN\
                    and cov > MIN_CONTIG_COV:
                        do_copy = True
                    else:
                        do_copy = False
                if do_copy:
                    print(line, file=out_file)

    logger.info('\nfilter_contigs:')
    logger.info('minimal contig length:', MIN_CONTIG_LEN)
    logger.info('minimal contig coverage:', MIN_CONTIG_COV)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contigs-in",
        metavar="CONTIGS_IN",
        type=Path,
        help="Input contigs file",
        required=True,
    )
    parser.add_argument(
        "--contigs-out",
        metavar="CONTIGS_OUT",
        type=Path,
        help="Output contigs file",
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
        "--min-contig-cov",
        metavar="MIN_CONTIG_COV",
        type=float,
        help="Minimum contig coverage",
        default=1.0,
    )
    parser.add_argument(
        "--min-contig-len",
        metavar="MIN_CONTIG_LEN",
        type=int,
        help="Minimum contig length",
        default=1,
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.contigs_in.is_file():
        logger.error(f"The given input file {args.contigs_in} was not found!")
        sys.exit(2)
    args.contigs_out.parent.mkdir(parents=True, exist_ok=True)
    filter_contigs(args.contigs_in, args.contigs_out, args.min_contig_len, args.min_contig_cov)


if __name__ == "__main__":
    sys.exit(main())
