#!/usr/bin/env python


"""Filter contigs."""


import logging
import platform
import random
import sys
import yaml
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
                line = line.rstrip('\\n')
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

    logger.info('\\nfilter_contigs:')
    logger.info('minimal contig length:', MIN_CONTIG_LEN)
    logger.info('minimal contig coverage:', MIN_CONTIG_COV)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(filter_contigs("$contigs", "$filtered_contigs", int("$min_contig_len"), float("$min_contig_cov")))
