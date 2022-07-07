#!/usr/bin/env python


"""Parse spades output."""


import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

logger = logging.getLogger()


def parse_spades_output(contigs_file, contig_len_dist_file,
                        contig_cov_dist_file, contig_len_x_cov_dist_file,
                        contig_ind_len_file, contig_ind_cov_file,
                        min_contig_len, min_contig_cov, max_no_contigs,
                        report_file):
    """
    Writes information from the SPAdes output file to report.txt and returns
    info about each contig.
    param: str contigs_file = contigs file
    param: str report_file = report file
    return: list lo_contig_data = [(contig number, length, coverage), ...];
            e.g.: [(1, 238256, 41.824755), (2, 208256, 8.247), ...]
    """

    lo_contig_data = []

    with open(report_file, 'a') as report:
        print('\n\nDe novo assembly (SPAdes):', file=report)
        print('contig\tlength (bp)\tcoverage', file=report)

        with open(contigs_file, 'r') as contigs:
            for line in contigs:
                if line.startswith('>'):
                    line = line.rstrip('\n')
                    data = line.split('_')
                    # contig number, length, coverage; e.g.:
                    # >NODE_1_length_238256_cov_41.824755
                    print(data[1], '\t', data[3], '\t', data[5], file=report)
                    # collect contig info
                    lo_contig_data.append((int(data[1]), int(data[3]),
                                           float(data[5])))

        # Write placeholders for figures to the report
        print('\nFigure: contigs vs length', file=report)
        print('\nFigure: contigs vs coverage', file=report)
        print('\nFigure: contig length distribution', file=report)
        print('\nFigure: contig coverage distribution', file=report)
        print('\nFigure: contig length * coverage distribution', file=report)

    plot_sum_dist(contig_len_dist_file, contig_cov_dist_file, lo_contig_data, min_contig_len, min_contig_cov)

    contig_stats = plot_len_x_cov_dist(contig_len_x_cov_dist_file, lo_contig_data, min_contig_len, min_contig_cov)

    write_to_file(report_file, contig_stats, min_contig_len, min_contig_cov)

    plot_ind_dist(contig_ind_len_file, contig_ind_cov_file, lo_contig_data)

    n_contigs = len(lo_contig_data)

    if n_contigs > max_no_contigs:
        logger.error(f"There were {n_contigs} contigs, which is far too many.")
        sys.exit(2)


def write_to_file(report_file, contig_stats, MIN_CONTIG_LEN, MIN_CONTIG_COV):
    """
    Writes results from the contig analysis to the report.
    param: str report_file = report file
    param: list contig_stats = list of four lists with contigs that
          1) passed none of the thresholds for contig length or coverage
          2) passed one threshold
          3) passed both thresholds
          4) had a very high coverage
    param: int MIN_CONTIG_LEN = minimum contig length (1000)
    param: float MIN_CONTIG_COV = minimum contig coverage (7.5)
    output: writes results to report file
    """

    a,b,c,d    = contig_stats
    no_contigs = sum(contig_stats)

    with open(report_file, 'a') as report:
        print('\nContig analysis:', file=report)
        print('(min length: ' + str(MIN_CONTIG_LEN) + ' bp, min coverage: '
              + str(MIN_CONTIG_COV) + 'x)', file=report)
        print('contigs that fail both thresholds: ',
              round(a*100/no_contigs, 2), '%', file=report)
        if (a*100/no_contigs) > 20:
            print('WARNING: The sample seems to be contaminated!', file=report)
        elif (a*100/no_contigs) > 5:
            print('NOTE: The sample might be contaminated!', file=report)
        print('contigs that are too short or have a low coverage: ',
              round(b*100/no_contigs, 2), '%', file=report)
        print('contigs that meet both thresholds: ',
              round(c*100/no_contigs, 2), '%', file=report)
        print('contigs with a high coverage (> 250x): ',
              round(d*100/no_contigs, 2), '%', file=report)
        if (d*100/no_contigs) > 0.5:
            print('NOTE: There might be a plasmid!', file=report)


def plot_it_1(output_file, data, MIN_CONTIG_VALUE, MIN, MAX, COLOR1, COLOR2,
              TITLE, XLABEL, YLABEL):
    """
    Plots a histogram of distributions and saves it to file.
    helper function to plot_length_dist() and plot_coverage_dist()
    param: str output_file = output file
    param: data = lo_len_data or lo_cov_data
    param: MIN_CONTIG_VALUE = either min contig length threshold (e.g.: 1000)
           or min contig coverage threshold (e.g.: 7.5)
    param: int or float MIN = minimal data range to plot
    param: int MAX = maximal data range to plot
    param: str COLOR1 = bar color: 'blue' or 'brown'
    param: str COLOR2 = line color: 'red' or 'blue'
    param: str TITLE = title of the chart
    param: str XLABEL = label for the x-axis
    param: str YLABEL = label for the y-axis
    param: str FILE_NAME = name of the file
    """

    fig, ax = plt.subplots()
    # print a histogram to file, were each bin covers about 5000 bp of the
    #   genome over the length of the genome
    plt.hist(data, color=COLOR1,
             bins=np.logspace(np.log10(MIN),np.log10(MAX), 30), alpha=0.5)
    plt.title(TITLE)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.gca().set_xscale('log')  # set the scale on the x-axis to log
    plt.axvline(MIN_CONTIG_VALUE, color=COLOR2)  # solid line, min length cutoff
    plt.savefig(output_file)
    plt.close()


def plot_sum_dist(contig_len_dist_file, contig_cov_dist_file, lo_contig_data, MIN_CONTIG_LEN, MIN_CONTIG_COV):
    """
    Plots histograms of contig-lengths and contig-coverage distributions and
      saves them to disk.
    param: str contig_len_dist_file = contig len dist output file
    param: str contig_cov_dist_file = contig cov dist output file
    param: list lo_contig_data = [(contig number, length, coverage), ...]
    param: int MIN_CONTIG_LEN = min contig length threshold (e.g.: 1000)
    param: float MIN_CONTIG_COV = min contig coverage threshold (e.g.: 7.5)
    output: histogram with contig lengths distributions
    output: histogram with contig coverages distributions
    """

    N = str(len(lo_contig_data))
    lo_len_data = [datum[1] for datum in lo_contig_data]
    lo_cov_data = [datum[2] for datum in lo_contig_data]

    # plot the contig-lengths distribution
    plot_it_1(contig_len_dist_file, lo_len_data, MIN_CONTIG_LEN, 100, 1000000,
              'blue', 'red', 'Contig length distribution (n=' + N + ')',
              'Length [log10]', 'Number of contigs')

    # plot the contig-coverage distribution
    plot_it_1(contig_cov_dist_file, lo_cov_data, MIN_CONTIG_COV, 0.1, 1000,
              'brown', 'blue', 'Contig coverage distribution (n=' + N + ')',
              'Coverage [log10]', 'Number of contigs')


def plot_len_x_cov_dist(contig_len_x_cov_dist_file, lo_contig_data, MIN_CONTIG_LEN,
                        MIN_CONTIG_COV):
    """
    Plots a histogram of contig-lengths * contig-coverage distributions and
      saves it to file. Bars that meet the thresholds for min coverage and min
      length are shown in green, those that meet one of the two are shown in
      orange, and those that meet none are shown in red. A large red area
      suggests problematic data, such as contaminations.
    param: str contig_len_x_cov_dist_file = contig len x cov dist output file
    param: list lo_contig_data = [(contig number, length, coverage), ...]
    param: int MIN_CONTIG_LEN = min contig length threshold (e.g.: 1000)
    param: float MIN_CONTIG_COV = min contig coverage threshold (e.g.: 7.5)
    return tuple = the number of contigs that meet none, one, or both minimal
           criteria or that have a very high coverage
    output: histogram
    """

    lo_none     = []   # meets no thresholds
    lo_one      = []   # meets one of two thresholds
    lo_both     = []   # meets both thresholds
    lo_high_cov = []   # very high coverage, maybe a plasmid

    for datum in lo_contig_data:
        # very high coverage indicates plasmid
        if datum[1] > MIN_CONTIG_LEN and datum[2] > 250:
            lo_high_cov.append(datum[1] * datum[2])
        # meets both thresholds
        elif datum[1] > MIN_CONTIG_LEN and datum[2] > MIN_CONTIG_COV:
            lo_both.append(datum[1] * datum[2])
        # meets no thresholds
        elif datum[1] <= MIN_CONTIG_LEN and datum[2] <= MIN_CONTIG_COV:
            lo_none.append(datum[1] * datum[2])
        # meets one of two thresholds
        else:
            lo_one.append(datum[1] * datum[2])

    no_contigs = len(lo_none) + len(lo_one) + len(lo_both) + len(lo_high_cov)

    fig, ax = plt.subplots()
    BINS = np.logspace(np.log10(10),np.log10(100000000), 30)

    # print a histogram to file:
    # bins=np.logspace(np.log10(10),np.log10(100000000), 30) = generate 30
    #  bins, ranging from 10 to 100000000
    # red: contigs that are too short and too low coverage
    plt.hist(lo_none, color='red', bins=BINS, alpha=0.5)
    # orange: contigs that are too short or too low coverage
    plt.hist(lo_one, color='orange', bins=BINS, alpha=0.5)
    # green: contigs with sufficient length and coverage
    plt.hist(lo_both, color='green', bins=BINS, alpha=0.5)
    # black: contigs with very high coverage
    plt.hist(lo_high_cov, color='black', bins=BINS, alpha=0.75)

    plt.title('Contig Length * Coverage distribution (n='\
              + str(no_contigs) + ')')
    plt.xlabel('Length * Coverage [log10]')
    plt.ylabel('Number of contigs')
    plt.gca().set_xscale('log')  # change x-axis to log10 scale
    plt.savefig(contig_len_x_cov_dist_file)
    plt.close()

    return len(lo_none), len(lo_one), len(lo_both), len(lo_high_cov)


def plot_it_2(output_file, lo_covs, Color, contig_1k, Title, IS_COV,
              med_cov, x_label):
    """
    Plots contig coverage or contig lengths distributions
      helper function to plot_ind_dist()
    param: str output_file = output file
    param: list lo_covs = list of data to be plotted, either contig coverage
           or contig lengths
    param: str Color = 'maroon' or 'seagreen'
    param: int contig_1k = the number of the smallest contig that is >= 1000 bp
    param: str Title = title of the chart
    param: bool IS_COV = if True, the lo_covs is contig coverage data,
           else: contig lengths
    param: float med_cov = median coverage
    param: str x_label = label for the x-axis
    output: one of four possible charts saved to file
    """

    # plot a bar chart that is as wide as the list of data, as high as the
    # data, in maroon, align bars to the edge, make bars 1 pixel wide, and
    # tone down the color intensity
    plt.bar(range(len(lo_covs)), lo_covs, color=Color, align='edge',
            width=1, alpha=0.7)

    # Draw horizonal lines at y that spans the x-range
    if IS_COV:
        plt.axhline(y=1, color='orange', linewidth=2)
        plt.axhline(y=10, color='orange', linewidth=2)
        plt.axhline(y=100, color='orange', linewidth=2)
        plt.axhline(y=7.5, color='blue', linewidth=2)
        plt.axhline(y= med_cov, color='red', linewidth=3)
    else:
        plt.axhline(y=100, color='orange', linewidth=2)
        plt.axhline(y=1000, color='blue', linewidth=2)
        plt.axhline(y=10000, color='orange', linewidth=2)
        plt.axhline(y=100000, color='orange', linewidth=2)

    # Draw vertical lines at x that span the y-range
    plt.axvline(x=contig_1k, color='blue', linewidth=2)
    plt.title(Title)
    plt.xlabel('Contig')
    plt.ylabel(x_label)
    plt.yscale('log')  # apply log10 scale on y-axis
    plt.savefig(output_file)
    plt.close()


def plot_ind_dist(contig_ind_len_file, contig_ind_cov_file, lo_contig_data):
    """
    Organizes the plotting of four graphs using plot_it_2()
    param: str contig_ind_len_file = contig ind len output file
    param: str contig_ind_cov_file = contig ind cov output file
    param: list lo_contig_data = = [(contig number, length, coverage), ...]
    output: manages plotting four charts, see also plot_it_2()
    """

    # Number of the smallest contig that is >= 1000 bp
    contig_1k = len([datum[2] for datum in lo_contig_data if datum[1] >= 1000])

    # Length of all contigs
    lo_len = [(datum[1]) for datum in lo_contig_data]
    plot_it_2(contig_ind_len_file, lo_len, 'maroon', contig_1k,
              'Contig length distribution (median=' \
              + str(round(np.median(lo_len), 2)) \
              + ')\nsmallest contig >=1000 bp = ' + str(contig_1k),
              False, 0, 'Length [log10]')

    # Coverage for all contigs
    lo_cov = [datum[2] for datum in lo_contig_data]
    med_cov = np.median(lo_cov)  # median coverage
    plot_it_2(contig_ind_cov_file, lo_cov, 'seagreen', contig_1k,
              'Contig coverage distribution\n(median='  \
              + str(round(med_cov, 2)) + ')',
              True, med_cov, 'Coverage [log10]')


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contig-cov-dist-file",
        metavar="CONTIG_COV_DIST_FILE",
        type=Path,
        help="Contig cov dist output file",
        required=True,
    )
    parser.add_argument(
        "--contigs-file",
        metavar="CONTIGS_FILE",
        type=Path,
        help="Input contigs file",
        required=True,
    )
    parser.add_argument(
        "--contig-ind-cov-file",
        metavar="CONTIG_IND_COV_FILE",
        type=Path,
        help="Contig ind cov output file",
        required=True,
    )
    parser.add_argument(
        "--contig-ind-len-file",
        metavar="CONTIG_IND_LEN_FILE",
        type=Path,
        help="Contig ind len output file",
        required=True,
    )
    parser.add_argument(
        "--contig-len-dist-file",
        metavar="CONTIG_LEN_DIST_FILE",
        type=Path,
        help="Contig len dist output file",
        required=True,
    )
    parser.add_argument(
        "--contig-len-x-cov-dist-file",
        metavar="CONTIG_LEN_X_COV_DIST_FILE",
        type=Path,
        help="Contig len x cov dist output file",
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
    parser.add_argument(
        "--max-no-contigs",
        metavar="MAX_NO_CONTIGS",
        type=int,
        help="Maximum number of contigs",
        default=999999999,
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
    if not args.contigs_file.is_file():
        logger.error(f"The given input file {args.contigs_file} was not found!")
        sys.exit(2)
    args.contig_cov_dist_file.parent.mkdir(parents=True, exist_ok=True)
    args.contig_ind_cov_file.parent.mkdir(parents=True, exist_ok=True)
    args.contig_ind_len_file.parent.mkdir(parents=True, exist_ok=True)
    args.contig_len_dist_file.parent.mkdir(parents=True, exist_ok=True)
    args.contig_len_x_cov_dist_file.parent.mkdir(parents=True, exist_ok=True)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    parse_spades_output(args.contigs_file, args.contig_len_dist_file,
                        args.contig_cov_dist_file, args.contig_len_x_cov_dist_file,
                        args.contig_ind_len_file, args.contig_ind_cov_file,
                        args.min_contig_len, args.min_contig_cov,
                        args.max_no_contigs, args.report_file)


if __name__ == "__main__":
    sys.exit(main())
