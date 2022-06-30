#!/usr/bin/env python


"""
This module reduces the number of reads (at random) after trimming with
Trimmomatic if there are too many reads as specified in the pipeline.
"""


import argparse
import logging
import random
import sys
from pathlib import Path


logger = logging.getLogger()


def fq_reader(file):
    """
    Takes the path and file name of a fastq file and returns a list of
      (header1, sequence, header2, quality score) tuples, one per read.
    param: str file = name of the forward or reverse read file
    return: list lo_reads = list of (header1, seq, header2, qual)
    """

    lo_reads = []
    count = 0
    with open(file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n')
            count += 1
            if count == 1:
                header1 = line
            elif count == 2:
                seq = line
            elif count == 3:
                header2 = line
            elif count == 4:
                qual = line
                lo_reads.append((header1, seq, header2, qual))
                count = 0
    return lo_reads


def make_lo_random_indices(N, k):
    """
    Selects k numbers drawn from a population of N.
    param int N = population size
    param int k = unique numbers to draw out of N
    return list lo_indices = k numbers drawn from N, without replacement
    """
    lo_indices = random.sample([i for i in range(N)], k)
    return sorted(lo_indices)


def fq_writer(outfile, lo_reads, lo_indices):
    """
    Extracts those reads specified by the list of indices from the
    list of reads and writes them to a new file.
    param: str outfile = output file, either 'paired_reads_1.fq' or
           'paired_reads_2.fq'
    lo_reads = list of (header1, seq, header2, qual) tuples
    lo_indices = list of k numbers drawn from a population of size N
    """
    with open(outfile, 'a') as write_file:
        for i in lo_indices:
            for j in  lo_reads[i]:
                print(j, file=write_file)


def read_reducer(reads_in_1, reads_in_2, reads_out_1, reads_out_2, random, k, START, STOP):
    """
    param: str reads_in_1 = input reads file 1
    param: str reads_in_2 = input reads file 2
    param: str reads_out_1 = output reads file 1
    param: str reads_out_2 = output reads file 2
    param: bool random: if True, selects k reads at random (w/o replacement);
           if False, use all reads between START and STOP
    param: int k = number of reads to select in random mode; if random=False,
           k!=0, then k reads will be chosen from the end of the file; to
           select k reads from the start of the file, set START=0, STOP=k
    param: int START = lower limit, index of first read to be included
    param: int STOP = upper limit, index of first read to be excluded
    output: a new file with fewer reads as the input file
    """
    # extracting reads from the forward read file and get total number of reads
    lo_F_reads = fq_reader(reads_in_1)
    N = len(lo_F_reads)
    text_F_file = 'There are ' + str(N) + ' reads in the F-read file.'

    # limit k and STOP to the size of N if either one is larger than N
    if k > N:
        k = N
    if STOP > N:
        STOP = N
    text_input = 'User input\nk = '+str(k) + '\nSTART = '+str(START) \
    + '\nSTOP = '+str(STOP)

    # lo_indices will be used to select reads from lo_F_reads and lo_R_reads
    # random choice of k reads
    if random:
        lo_indices = make_lo_random_indices(N, k)
        text_indices = 'generated ' + str(len(lo_indices)) + ' indices at random'
    # non-random, requires two vaklues out of k, START, STOP
    else:
        # use START - STOP as range, which is 0 to N (= all) by default
        if k == 0:
            lo_indices = [i for i in range(START, STOP)]
            text_indices = 'generated ' + str(len(lo_indices))\
            + ' indices from ' + str(START) + ' to ' + str(STOP)
        # select k reads from the end
        else:
            lo_indices = [i for i in range(N - k, N)]
            text_indices = 'generated ' + str(len(lo_indices))\
            + ' indices from ' + str(N-k) + ' to ' + str(N)

    # write selected reads to file
    fq_writer(reads_out_1, lo_F_reads, lo_indices)

    # extracting and writing the reverse reads using the same indices
    lo_R_reads = fq_reader(reads_in_2)
    text_R_file = 'There are ' + str(len(lo_R_reads))\
    + ' reads in the R-read file.'

    # There should be the same number of reads in the F- and R-read file
    if len(lo_R_reads) == N:
        fq_writer(reads_out_2, lo_R_reads, lo_indices)
        text_final = 'Writing new read files compete.'
    else:
        text_final = 'Could not complete writing files.'

    # logging text
    for text in ['\n\nRead reduction:', text_input, text_F_file, text_R_file,
                 text_indices, text_final]:
        logger.info(text)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
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
        "--random",
        metavar="RANDOM",
        type=bool,
        help="Random",
        default=False,
    )
    parser.add_argument(
        "--k",
        metavar="K",
        type=int,
        help="k",
        default=0,
    )
    parser.add_argument(
        "--start",
        metavar="START",
        type=int,
        help="Start",
        default=0,
    )
    parser.add_argument(
        "--stop",
        metavar="STOP",
        type=int,
        help="Stop",
        default=999999999,
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
    read_reducer(args.reads_in_1, args.reads_in_2, args.reads_out_1, args.reads_out_2,
                 args.random, args.k, args.start, args.stop)


if __name__ == "__main__":
    sys.exit(main())
