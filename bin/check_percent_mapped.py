#!/usr/bin/env python


"""Check percent mapped."""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def check_percent_mapped(percent_mapped, MIN_PERCENT_MAPPED):
    """Check percent mapped."""
    # QC check in case there are too many unmapped reads
    if percent_mapped < MIN_PERCENT_MAPPED:
        logger.error('There were only ' + str(percent_mapped) \
            + '% mapped reads, which is far too few.')
        sys.exit(2)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--percent-mapped",
        metavar="PERCENT_MAPPED",
        type=float,
        help="Percent mapped",
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
        "--min-percent-mapped",
        metavar="MIN_PERCENT_MAPPED",
        type=float,
        help="Minimum percent mapped",
        default=0.0,
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    check_percent_mapped(args.percent_mapped, args.min_percent_mapped)


if __name__ == "__main__":
    sys.exit(main())
