#!/usr/bin/env python


"""
Opens both read files simultaneuously, checks if either read sequence
  contains 25 or more Gs in a row (xG=25), and writes the reads to new
  files if neither one does.
"""


import argparse
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


def get_header_symbol(file):
    """
    Returns the first character of the first line, which identifies the
      header of a read, ususally a "@".
    helper function to remove_poly_Gs()
    param: str file = name of the read file
    return: str CHAR = first character of the header
    """
    with open(file, 'r') as infile:
        for line in infile:
            CHAR = line[0]
            return CHAR


def read_line(path_file):
    """
    The yield converts the read file into an iterable item (like a list),
      returning one line at a time.
    helper function to remove_poly_Gs()
    param: str path_file = path and filename of read file
    yield: one line at a time from the read file
    """
    with open(path_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            yield line


def write_line(path_file, read):
    """
    Write one read at a time, spread over four lines:
        header, sequence, header, quality score.
    helper function to remove_poly_Gs()
    param: str path_file = path and filename of the new read file
    param: list read = [header, sequence, header, quality score] = one read
    output: a new read file
    """
    with open(path_file, 'a') as outfile:
        for data in read:
            print(data, file=outfile)


def remove_poly_Gs(reads_in_1, reads_in_2, reads_out_1, reads_out_2, log_file, xG):
    """
    Opens both read files simultaneuously, checks if either read sequence
      contains 25 or more Gs in a row (xG=25), and writes the reads to new
      files if neither one does.
    param: str reads_in_1 = input reads file 1
    param: str reads_in_2 = input reads file 2
    param: str reads_out_1 = output reads file 1
    param: str reads_out_2 = output reads file 2
    param: str log_file = log file
    param: int xG = remove reads with this many 'G's in a row
    output: two new read files and a log file
    """
    bad_read_count = 0

    # get the character that identifies the start of a read header, usually a
    #  '@', but not always
    CHAR = get_header_symbol(reads_in_1)

    # zip() reads both files, one line at a time, returns a tuple
    for new_lines in zip(read_line(reads_in_1),
                         read_line(reads_in_2)):

        line_1, line_2 = new_lines

        # one list of data per read, start with [header]
        if line_1.startswith(CHAR):
            F_read = [line_1]
            R_read = [line_2]
        # add [seq, +, qual]
        else:
            F_read.append(line_1)
            R_read.append(line_2)

        # if [header, seq, +, qual], then check for poly-Gs and write to file,
        #  but only if no poly-Gs are found
        if len(F_read) == 4 and len(R_read) == 4:
            if not ('G' * xG in F_read[1])\
            and not ('G' * xG in R_read[1])\
            and not ('C' * xG in F_read[1])\
            and not ('C' * xG in R_read[1]):
                write_line(reads_out_1, F_read)
                write_line(reads_out_2, R_read)
            else:
                bad_read_count += 1

    print('\nDiscarded', bad_read_count, 'read pairs that contained >= '\
          + str(xG) + ' Gs.\n')

    with open(log_file, 'a') as log:
        print('\nDiscarded', bad_read_count, 'read pairs that contained >= '\
              + str(xG) + ' Gs.\n', file=log)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Remove poly Gs.",
        epilog="Example: python remove_poly_gs.py reads_1.fq reads_2.fq reads_nog_1.fq reads_nog_2.fq sample.log",
    )
    parser.add_argument(
        "reads_in_1",
        metavar="READS_IN_1",
        type=Path,
        help="Input reads file 1",
    )
    parser.add_argument(
        "reads_in_2",
        metavar="READS_IN_2",
        type=Path,
        help="Input reads file 2",
    )
    parser.add_argument(
        "reads_out_1",
        metavar="READS_OUT_1",
        type=Path,
        help="Output reads file 1",
    )
    parser.add_argument(
        "reads_out_2",
        metavar="READS_OUT_2",
        type=Path,
        help="Output reads file 2",
    )
    parser.add_argument(
        "log_file",
        metavar="LOG_FILE",
        type=Path,
        help="Log file",
    )
    parser.add_argument(
        "--xg",
        metavar="XG",
        type=int,
        help="Remove reads with this many 'G's in a row",
        default=25,
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.reads_in_1.is_file():
        logger.error(f"The given input file {args.reads_in_1} was not found!")
        sys.exit(2)
    if not args.reads_in_2.is_file():
        logger.error(f"The given input file {args.reads_in_2} was not found!")
        sys.exit(2)
    args.reads_out_1.parent.mkdir(parents=True, exist_ok=True)
    args.reads_out_2.parent.mkdir(parents=True, exist_ok=True)
    args.log_file.parent.mkdir(parents=True, exist_ok=True)
    remove_poly_Gs(args.reads_in_1, args.reads_in_2, args.reads_out_1, args.reads_out_2, args.log_file, args.xg)


if __name__ == "__main__":
    sys.exit(main())
