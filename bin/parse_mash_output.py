#!/usr/bin/env python


"""Parse mash output."""


import argparse
import csv
import logging
import sys
from pathlib import Path


logger = logging.getLogger()


class Mash_result(object):
    """
    The output from Mash dist is represented as an object for easy comparison
    and ranking.
    """

    # initializes a Mash_result object
    # removes paths and file extentions from reference and query and converts
    #  text to numbers where appropriate
    def __init__(self, reference='', query='', distance=1, pvalue=1,
                 hashes=''):
        self._reference = reference.split('/')[-1].split('.')[0]
        self._query     = query.split('/')[-1].split('.')[0]
        self._distance  = float(distance)
        self._pvalue    = float(pvalue)
        self._hashes    = hashes

    # getter functions to return object attributes
    def get_all(self):
        return (self._reference, self._distance, self._pvalue, self._hashes)
    def get_reference(self):
        return self._reference
    def get_query(self):
        return self._query
    def get_distance(self):
        return self._distance
    def get_pvalue(self):
        return self._pvalue
    def get_hashes(self):
        return self._hashes

    # functions needed for comparisons
    def __lt__(self, other):
        return self._distance < other._distance
    def __le__(self, other):
        return self._distance <= other._distance
    def __eq__(self, other):
        return self._distance == other._distance
    def __ne__(self, other):
        return self._distance != other._distance
    def __gt__(self, other):
        return self._distance > other._distance
    def __ge__(self, other):
        return self._distance >= other._distance


def read_species_file(species_file):
    """
    Extracts species dict from the input species file.
      helper function to sort_mash_output()
    param: str species_file = input species file
    param: dict do_species = dict of species abbreviations mapped to names
    """

    do_species = dict()
    with open(species_file, 'r', newline='') as species:
        species_reader = csv.reader(species)
        for row in species_reader:
            do_species[row[0]] = row[1]
    return do_species


def read_distance_file(dist_file):
    """
    Extracts data from the distances file created by 'Mash dist'.
      helper function to sort_mash_output()
    param: str dist_file = name of TAB file produced by mash dist
    return: list lo_refs = Mash_result objects with reference name, query name,
            Mash distance, p-value, and number of matching hashes, sorted by
            distance
    """

    # retrieve the data and write them to a sorted list as a Mash_result object
    lo_refs = []
    with open(dist_file, 'r') as in_file:
        for line in in_file:
            line = line.rstrip('\n')
            if line != '':
                # reference, query, Mash distance, p-value, matching hashes:
                REF, QRY, DIS, PVL, HSH = line.split()
                # creating a Mash_result object
                mr = Mash_result(REF, QRY, DIS, PVL, HSH)
            else:  # stand-in if Mash failed for some reason
                mr = Mash_result('-fail-', '-fail-', 1, 1, '-fail-')
            lo_refs.append(mr)
    return sorted(lo_refs)


def write_references_file(references_file, lo_min_dist_refs):
    with open(references_file, 'a', newline='') as references:
        references_writer = csv.writer(references)
        for ref in lo_min_dist_refs:
            references_writer.writerow([ref])


def write_to_file(report_file, lo_refs, lo_sm_dist, header_text, do_species):
    """
    Write Mash results to report file.
      helper function to parse_mash_output()
    param: str report_file = output report file
    param: list lo_refs = Mash_result objects sorted by distance
    param: list lo_sm_dist = list of references with the smallest Mast distances
    param: str header_text = header text for the report.txt
    param: dict do_species = dict of species abbreviations mapped to names
    output: text written to file
    """

    # write results to the report file
    with open(report_file, 'a') as report:
        # header
        print(header_text, file=report)
        # if more than one reference
        if len(lo_sm_dist) > 1:

            # print reference(s) with smallest distance, then the runner-up
            for i in range(len(lo_sm_dist)):
                if i < len(lo_sm_dist) - 1:
                    print('\nReference with the shortest distance', file=report)
                else:
                    print('\nRunner up', file=report)
                # extract data from the Mash_result object
                print('Strain name:\t\t', lo_sm_dist[i].get_reference(),
                      file=report)
                print('Mash distance:\t\t', lo_sm_dist[i].get_distance(),
                      file=report)
                print('P-value:\t\t', lo_sm_dist[i].get_pvalue(),
                      file=report)
                print('Matching hashes:\t', lo_sm_dist[i].get_hashes(),
                      file=report)

                # extract the species name in case of the contamination check
                if do_species is not None:
                    tent_species = do_species.get(lo_sm_dist[i].get_reference(),
                                                  'UNKNOWN SPECIES')
                    print('These reads seem to have come from:', tent_species, \
                          'or a related species.', file=report)
        # in case there is only one reference
        else:
            print('\nReference with the shortest distance', file=report)
            print('Strain name:\t\t', lo_sm_dist[0].get_reference(),
                  file=report)
            print('Mash distance:\t\t', lo_sm_dist[0].get_distance(),
                  file=report)
            print('P-value:\t\t', lo_sm_dist[0].get_pvalue(),
                  file=report)
            print('Matching hashes:\t', lo_sm_dist[0].get_hashes(),
                  file=report)
            print('\nRunner up: none', file=report)


def check_quality(ref_object, report_file):
    """
    Checks if distance and p-value for the species-reference with the smallest
      distance are below thresholds. The sample might be contaminated if not.
    param: str ref_object = reference object with smallest distance
    param: str report_file = output report file
    return: bool passed_qc = True if distance and p-value are below thresholds
    output: writes reference_ID, Mash_distance, P_value, Matching_hashes
            to file
    """

    MAX_MASH_PVAL = 0.0000001  # maximum acceptable P value
    MAX_MASH_DIST = 0.1        # maximum acceptable value for Mash distance

    # unpack data
    REF, DIS, PVL, HSH = ref_object.get_all()

    # checks that the Mash distance and p-value are below acceptable values
    if PVL <= MAX_MASH_PVAL and DIS <= MAX_MASH_DIST:
        qc_text = 'PASSED QC'
        passed_qc = True
    else:
        qc_text = 'WARNING, THIS STRAIN MIGHT BE LESS THAN IDEAL'
        passed_qc = False

    # writes results to the report file
    with open(report_file, 'a') as report:
        print('\nMash QC results:', qc_text, file=report)

    return passed_qc


def parse_mash_output(dist_file, sp_abbr, references_file, report_file, species_file):
    """
    Parses a distances file and returns a list of Mash_result objects of
      those references that have the shortest distance (one or more), and the
      runner-up. All (except for the runner-up) will be tested by BWA. The one
      with the highest percentage of mapped reads will then become the
      reference. This was done since two isolates with the same Mash distance
      can have very different percentage of mapped reads.
      The while-loop will add references to the list until a reference is
      found that is outside the reasonable range of the smallest distance;
      this reference will also be printed to file for comparison.
    param: str dist_file = name of TAB file produced by mash dist
    param: str sp_abbr = species abbreviation
    param: str report_file = output report file
    param: str species_file = input species file
    output: writes reference_ID, Mash_distance, P_value, Matching_hashes to
            file, e.g.: [('IDR00123', 0.0293323, 0.0, '182/400')
                         ('IDR00456', 0.0293323, 0.0, '182/400')
                         ('IDR00789', 0.034976,  0.0, '160/400')]
    """

    lo_sm_dist = []

    # Extracts data from the distances file created by 'mash dist',
    #  where each reference is a Mash_result object
    lo_refs = read_distance_file(dist_file)
    logger.info('Mash list of references:', [mr.get_all() for mr in lo_refs])

    # if there is more than 1 reference
    if len(lo_refs) > 1:
        # adds the first two references with the smallest distances to the list
        lo_sm_dist.extend(lo_refs[0:2])

        n = 2   # first two items on list sorted by distance

        # Make a list of references with similar, short Mash distance to be
        #  checked later with BWA:
        # While second to last and last item in the list have equal distance
        # or the last item in the list has a distance similar to the first
        # item, add a new reference to the list.
        # The function x+0.001+x*0.05 was determined emprically based on
        # actual Mash distances: it runs parallel to the Hashes vs Mash
        # distance curve. At small distances, the "+0.001" will keep the two
        # curves parallel, at large distances, the "x*0.05" becomes more
        # important), e.g.: [0.00173, 0.00173, 0.00183, 0.00348]
        # 1. and 2. are equal, 3. is similar to 1., and 4. is the runner-up
        while lo_sm_dist[-2].get_distance() == lo_sm_dist[-1].get_distance() \
        or lo_sm_dist[-1].get_distance() < (lo_sm_dist[0].get_distance()+0.001 \
                     +lo_sm_dist[0].get_distance()*0.05):

            # add references that are within range to lo_sm_dist
            if n < len(lo_refs):
                lo_sm_dist.append(lo_refs[n])
                n += 1
            # if end of lo_refs is reached, add last list element again
            # (the duplicate will be removed when lo_min_dist_refs is created)
            else:
                lo_sm_dist.append(lo_refs[n-1])
                break

        # makes a list of those reference fasta files that have the smallest
        #  Mash distances
        lo_min_dist_refs = [ref.get_reference() + '.fa' for ref in \
                            lo_sm_dist[:-1]]
        logger.info('Mash list of min distance references:', lo_min_dist_refs)

    # only one reference
    else:
        lo_sm_dist.extend(lo_refs)

        # makes a list of those reference fasta files that have the smallest
        #  Mash distances
        lo_min_dist_refs = [ref.get_reference() + '.fa' for ref in \
                            lo_sm_dist]

    # writes results to the log and report files
    if species_file is not None:
        do_species = read_species_file(species_file)
        write_to_file(report_file, lo_refs, lo_sm_dist,
                      '\nContamination check (Mash):', do_species)
    else:
        write_to_file(report_file, lo_refs, lo_sm_dist,
                      '\nFinding a reference strain (Mash):', None)

    # checks that the distance and p-value are below threshold; if not, the
    #  sample might be contaminated or the reference a poor choice
    passed_qc = check_quality(lo_sm_dist[0], report_file)
    logger.info('Mash passed QC:', passed_qc)

    if species_file is not None:
        # the species based on the Mash data
        mash_species = lo_min_dist_refs[0][:-3]

        # over-riding the species check
        # sometimes the lab only provides the info that an isolate belongs to a
        #  "complex" without giving the exact species, the code below allows for
        #  this kind of uncertainty when doing the species check with Mash
        if sp_abbr == 'Spy' and mash_species == 'Sdy':
            mash_species = 'Spy'
            passed_qc = True
        elif sp_abbr == 'Cro' and mash_species in ['Cco','Cdu','Cma','Cmu',
                                                'Csa','Ctu','Cun']:
            mash_species = 'Cro'
            passed_qc = True
        elif sp_abbr == 'Ecl' and mash_species in ['Eas','Eca','Ecl','Ehh',
                                                'Eho','Eko','Elu','Ero']:
            mash_species = 'Ecl'
            passed_qc = True
        elif sp_abbr == 'Cbo':
            mash_species = 'Cbo'
            passed_qc = True

        if not passed_qc:
            logger.error("The reads did not pass the Mash QC check.")
            sys.exit(2)

        if mash_species != sp_abbr:
            logger.error("The reads did not pass the Mash species check.")
            sys.exit(2)
    else:
        if not passed_qc:
            logger.error("The reference did not pass the Mash QC check.")
            sys.exit(2)

    write_references_file(references_file, lo_min_dist_refs)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dist-file",
        metavar="DIST_FILE",
        type=Path,
        help="Distance file",
        required=True,
    )
    parser.add_argument(
        "--references-file",
        metavar="REFERENCES_FILE",
        type=Path,
        help="Output references file",
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
        "--sp-abbr",
        metavar="SP_ABBR",
        help="Species abbreviation",
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
        "--species-file",
        metavar="SPECIES_FILE",
        type=Path,
        help="Species file",
        default=None,
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.dist_file.is_file():
        logger.error(f"The given input file {args.dist_file} was not found!")
        sys.exit(2)
    if args.species_file is not None and not args.species_file.is_file():
        logger.error(f"The given input file {args.species_file} was not found!")
        sys.exit(2)
    args.references_file.parent.mkdir(parents=True, exist_ok=True)
    args.report_file.parent.mkdir(parents=True, exist_ok=True)
    parse_mash_output(args.dist_file, args.sp_abbr, args.references_file, args.report_file, args.species_file)


if __name__ == "__main__":
    sys.exit(main())
