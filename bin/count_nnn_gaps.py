#!/usr/bin/env python


"""Count nnn gaps."""


import argparse
import csv
import logging
import matplotlib.pyplot as plt
import re
import sys
from numpy import mean, std, ceil
from pathlib import Path


logger = logging.getLogger()


def parse_file(depth_file):
    """
    Extracts the read depths values from a 'samtools_depth.txt' file.
    return: list lo_depths = read depth for each base in the genome, sorted by
            position; e.g.: [0,0,1,1,4,5,7,19,45, ...]
    """

    lo_depths = []
    with open(depth_file, 'r') as in_file:
        for line in in_file:
            line = line.rstrip('\n')
            depth = line.split()[2]
            lo_depths.append(int(depth))
    return lo_depths


def translate_low_coverage(lo_depths, MIN_DEPTH):
    """
    Translation of a list of numbers into a string of 'a' and 'B', where an
      'a' indicates that the corresponding number was Above MIN_DEPTH or
      'B' Below it
    param: list lo_depths = read depths for each base in the genome
    param: int MIN_DEPTH = minimal value to be sufficiently mapped by reads
    return str ab_string = string of 'a' and 'B'
           e.g.: [1,1,5,6,4,1], MIN_DEPTH=3 => 'BBaaaB'
    """

    ab_string = ''
    for depth in lo_depths:
        if depth >= MIN_DEPTH:
            ab_string += 'a'
        else:
            ab_string += 'B'
    return ab_string


def depth_by_group(lo_depths, MIN_DEPTH):
    """
    Splits lo_depths into three lists (depth == 0, 0 < depth < MIN_DEPTH,
      depth >= MIN_DEPTH), and calculates the mean, standard deviation and
      length of each list.
    param: list lo_depths = read depths for each base in the genome
    param: int MIN_DEPTH = minimal value to be sufficiently mapped by reads
    return: tuple of tuples with mean, std and length for each list
    """

    lo_above = []
    lo_below = []
    lo_zeros = []

    for depth in lo_depths:
        if depth >= MIN_DEPTH:
            lo_above.append(depth)
        elif 0 < depth < MIN_DEPTH:
            lo_below.append(depth)
        elif depth == 0:
            lo_zeros.append(depth)

    if len(lo_above) > 0:
        above = (mean(lo_above), std(lo_above), len(lo_above))
    else:
        above = (0,0,0)
    if len(lo_below) > 0:
        below = (mean(lo_below), std(lo_below), len(lo_below))
    else:
        below = (0,0,0)
    if len(lo_zeros) > 0:
        zero =  (mean(lo_zeros), std(lo_zeros), len(lo_zeros))
    else:
        zero = (0,0,0)

    return above, below, zero


def count_gaps(ab_string, GAP_LENGTH):
    """
    Counts the number of all gaps >= GAP_LENGTH.
    param: str ab_string = string of 'a' and 'B', e.g. 'BBaaaB...'
    param: int GAP_LENGTH = minimal gap length to be counted
    return: int no_gaps = number of gaps >= GAP_LENGTH
    """

    query = 'B' * GAP_LENGTH + '+'
    # re.findall('BBBB+', str) finds all substrings of 4 or more 'B' in str
    lo_gaps = re.findall(query, ab_string)
    no_gaps = len(lo_gaps)
    lo_gap_lens = [len(x) for x in lo_gaps]
    return lo_gap_lens, no_gaps


def write_report(report_file, lo_depth_stats, MIN_DEPTH, GAP_LENGTH):
    """
    Writes summary data to report.txt.
    param: list lo_depth_stats = depth statistics to be writted to report.txt
    param: int MIN_DEPTH = minimal value to be sufficiently mapped by reads
    param: int GAP_LENGTH = minimal gap length
    """

    depth_mean, depth_sd, count_all, count_below, count_above, \
        depth_above, depth_below, depth_zero, lo_gap_lens, no_gaps = lo_depth_stats

    percent_below = round((count_below * 100 / count_all), 2)
    percent_above = 100 - percent_below

    with open(report_file, 'a') as report:

        print('\nAlignment QC (Samtools depth):', file=report)

        print('Total number of bases: ' + str(count_all), file=report)
        print('Number (percent) of bases with read depth < '\
              + str(MIN_DEPTH) + ': \t' + str(count_below)\
              + ' (' + str(percent_below) + '%)', file=report)
        print('Number (percent) of bases with read depth >= '\
              + str(MIN_DEPTH) + ':\t' + str(count_above)\
              + ' (' + str(percent_above) + '%)', file=report)

        print('Average read depth (S.D.): ' + str(depth_mean) + ' ('\
              + str(depth_sd) + ')', file=report)

        print('Average read depth (S.D., count) for bases with read depth'\
              + ' >= ' + str(MIN_DEPTH) + ':\t\t'\
              + str(round(depth_above[0], 2)) + '\t('\
              + str(round(depth_above[1], 2)) + ',\t'\
              + str(round(depth_above[2], 2)) + ')', file=report)
        print('Average read depth (S.D., count) for bases with read depth'\
              + ' > 0 and < ' + str(MIN_DEPTH) + ':\t'\
              + str(round(depth_below[0], 2)) + '\t('\
              + str(round(depth_below[1], 2)) + ',\t'\
              + str(round(depth_below[2], 2)) + ')', file=report)
        print('Average read depth (S.D., count) for bases with read depth'\
              + ' == 0:\t\t'\
              + str(round(depth_zero[0], 2)) + '\t('\
              + str(round(depth_zero[1], 2)) + ',\t'\
              + str(round(depth_zero[2], 2)) + ')', file=report)

        print('Number of gaps >= ' + str(GAP_LENGTH) + ' bases:', no_gaps,\
              file=report)
        print('List of gaps >= ' + str(GAP_LENGTH) + ' bases:',\
              sorted(lo_gap_lens), file=report)

        print('Total number of bases in gaps >= ' + str(GAP_LENGTH)\
              + ' bases:', sum(lo_gap_lens), file=report)

        # placeholders for images in html file
        suffix = '' # TODO: replace this
        print('Figure: Read depth per base' + suffix + ' (plot)', file=report)
        print('Figure: Read depth per base' + suffix + ' (histogram)',\
              file=report)


def write_log(MIN_DEPTH, GAP_LENGTH, INTERVAL):
    """
    Writes settinigs to log.txt.
    param: int MIN_DEPTH = minimal value to be sufficiently mapped by reads
    param: int GAP_LENGTH = minimal gap length
    param: int INTERVAL = size of subsections of the genome, e.g. 5000 bp
    """

    logger.info('\ncount_nnn_gaps.py settings:')
    logger.info('minimal read depth:', MIN_DEPTH)
    logger.info('minimal gap length:', GAP_LENGTH)
    logger.info('counting interval:', INTERVAL)


def calc_n_per_interval(ab_string, INTERVAL):
    """
    Counts the number of low depth bases within a region of length INTERVAL.
    param: str ab_string = string of 'a' and 'B', e.g. 'BBaaaB...'
    param: int INTERVAL = size of subsections of the genome, e.g. 5000 bp
    return: list lo_depth_per_interval = list of counts of low coverage bases
            per INTERVAL, e.g.: [0,0,167,5000,5000,321,0,0,0,...]
    """

    lo_depth_per_interval = [ab_string[i:i+INTERVAL].count('B') for i in\
                             range(0, len(ab_string), INTERVAL)]
    return lo_depth_per_interval


def histo_read_depth_distr(histo_depths_file, lo_depths):
    """
    Plots a histogram of the read depth per base distribution.
    param: list lo_depths = read depths for each base in the genome
    output: histogram of the read depth per base distribution, including the
            mean and mean +/- 3 StDev
    """

    average = mean(lo_depths)
    SD = std(lo_depths)

    plt.hist(lo_depths, bins=20, color='brown')
    # Draw a default (v-)line at x that spans the y-range
    plt.axvline(x= average, color='blue', linewidth=3)
    plt.axvline(x= (average + 3 * SD), color='blue')
    plt.axvline(x= (average - 3 * SD), color='blue')
    plt.title('Read depth per base (average +/- 3 SD)')
    plt.xlabel('Read depth per base')
    plt.ylabel('Number of bases')
    plt.savefig(histo_depths_file)
    plt.close()


def plot_read_depth_distr(plot_depths_file, lo_depths):
    """
    Plots the read depth per base distribution.
    param: list lo_depths = read depths for each base in the genome
    output: histogram of the read depth per base distribution, including
            lines indicating the mean +/- 3 * StDev
    """

    average = mean(lo_depths)
    SD = std(lo_depths)
    max_x_val = int(ceil(len(lo_depths)/500000)) + 1  # highest value on x-axis
    fig, ax = plt.subplots()

    plt.plot(lo_depths, color='blue')
    # Draw a line at y that spans the x-range
    plt.axhline(y= average, color='red', linewidth=3)
    plt.axhline(y= (average + 3 * SD), color='orange')
    plt.axhline(y= (average - 3 * SD), color='orange')
    plt.title('Read depth per base (average +/- 3 SD)')
    plt.xlabel('Read depth per base [Mb]')
    plt.ylabel('Number of reads')

    # set the tick labels to 500000 bp intervals
    ax.set_xticks([500000 * x for x in range(0, max_x_val)])
    # change the tick labels from absolute values to intervals of 0.5 Mb
    ax.set_xticklabels([0.5 * x for x in range(0, max_x_val)])

    plt.savefig(plot_depths_file)
    plt.close()


def count_nnn_gaps(depth_file, histo_depths_file, plot_depths_file,
                   output_file, report_file, percent_mapped, MIN_DEPTH,
                   GAP_LENGTH, INTERVAL, MAX_NO_NS, MAX_NO_GAPS, MAPPED_THRESHOLD):
    """
    Main function: parses the 'samtools_depth.txt' file, writes statistics
      to the report, and plots the distribution of poorly covered regions.
    param: int MIN_DEPTH = minimal value to be sufficiently mapped by reads
    param: int GAP_LENGTH = minimal gap length
    param: int INTERVAL = size of subsections of the genome, e.g. 5000 bp
    return: float depth_mean = average read depth per base
    return: float depth_sd = standard deviation read depth per base
    """

    # list of depths values, where the first element is base number 1
    lo_depths = parse_file(depth_file)

    # string of 'a' and 'B' for above/below MIN_DEPTH
    ab_string = translate_low_coverage(lo_depths, MIN_DEPTH)

    # counts the number of gaps larger than GAP_LENGTH
    lo_gap_lens, no_gaps = count_gaps(ab_string, GAP_LENGTH)

    # mean and standard deviation for all bases
    depth_mean = round(mean(lo_depths), 2)
    depth_sd = round(std(lo_depths), 3)
    # count number of bases above and below MIN_DEPTH
    count_all = len(ab_string)
    count_below = ab_string.count('B')
    count_above = count_all - count_below
    # average read depth, standard deviation and number for bases with
    # depth == 0, 0 < depth < MIN_DEPTH, and depth >= MIN_DEPTH
    depth_above, depth_below, depth_zero = depth_by_group(lo_depths, MIN_DEPTH)
    # combine data to one list for passing to write_report()
    lo_depth_stats = [depth_mean, depth_sd, count_all, count_below, \
                      count_above, depth_above, depth_below, depth_zero, \
                      lo_gap_lens, no_gaps]

    # add data to the report and log file
    write_report(report_file, lo_depth_stats, MIN_DEPTH, GAP_LENGTH)
    write_log(MIN_DEPTH, GAP_LENGTH, INTERVAL)

    # plot the distribution of read depths per base
    histo_read_depth_distr(histo_depths_file, lo_depths)
    plot_read_depth_distr(plot_depths_file, lo_depths)

    # if there are too many unmapped bases, abort unless it might be a
    # candidate reference
    if (count_below > MAX_NO_NS) and not (percent_mapped < MAPPED_THRESHOLD):
        logger.error('There are ' + str(n_count) \
                     + ' unmapped bases, which is far too many.')
        sys.exit(2)

    # if there are too many gaps, abort unless it might be a candidate reference
    if (no_gaps > MAX_NO_GAPS) and not (percent_mapped < MAPPED_THRESHOLD):
        logger.error('There are ' + str(gaps) \
                     + ' gaps compared to the reference genome,'\
                     + ' which is far too many.')
        sys.exit(2)

    with open(output_file, 'a', newline='') as output:
        output_writer = csv.writer(output)
        output_writer.writerow([depth_mean])
        output_writer.writerow([depth_sd])


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--depth-file",
        metavar="DEPTH_FILE",
        type=Path,
        help="Input depth file",
        required=True,
    )
    parser.add_argument(
        "--histo-depths-file",
        metavar="HISTO_DEPTHS_FILE",
        type=Path,
        help="Output histo depths file",
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
        "--percent-mapped",
        metavar="PERCENT_MAPPED",
        type=float,
        help="Percent mapped",
        required=True,
    )
    parser.add_argument(
        "--plot-depths-file",
        metavar="PLOT_DEPTHS_FILE",
        type=Path,
        help="Output plot depths file",
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
        "--gap-length",
        metavar="GAP_LENGTH",
        type=int,
        help="Gap length",
        default=100,
    )
    parser.add_argument(
        "--interval",
        metavar="INTERVAL",
        type=int,
        help="Interval",
        default=5000,
    )
    parser.add_argument(
        "--log-level",
        metavar="LOG_LEVEL",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        help="The desired log level",
        default="INFO",
    )
    parser.add_argument(
        "--mapped-threshold",
        metavar="MAPPED_THRESHOLD",
        type=float,
        help="Mapped threshold",
        default=0.0,
    )
    parser.add_argument(
        "--max-no-gaps",
        metavar="MAX_NO_GAPS",
        type=int,
        help="Maximum number of gaps",
        default=999999999,
    )
    parser.add_argument(
        "--max-no-ns",
        metavar="MAX_NO_NS",
        type=int,
        help="Maximum number of unmapped bases",
        default=999999999,
    )
    parser.add_argument(
        "--min-depth",
        metavar="MIN_DEPTH",
        type=int,
        help="Minimum depth",
        default=1,
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.depth_file.is_file():
        logger.error(f"The given input file {args.depth_file} was not found!")
        sys.exit(2)
    args.histo_depths_file.parent.mkdir(parents=True, exist_ok=True)
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    args.plot_depths_file.parent.mkdir(parents=True, exist_ok=True)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    count_nnn_gaps(args.depth_file, args.histo_depths_file, args.plot_depths_file,
                   args.output_file, args.report_file, args.percent_mapped, args.min_depth,
                   args.gap_length, args.interval, args.max_no_ns, args.max_no_gaps, args.mapped_threshold)


if __name__ == "__main__":
    sys.exit(main())
