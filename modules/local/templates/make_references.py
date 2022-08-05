#!/usr/bin/env python


"""Make references."""


import csv
import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def make_references(references_file, lo_reference, lo_cluster, lo_fasta, lo_snp_cons, lo_bwa, lo_fai, lo_mutations_matrix):
    """Make references."""

    header = ["reference", "cluster", "fasta", "snp_cons", "bwa", "fai", "mutations_matrix"]

    with open(references_file, mode="w", newline="") as references:
        writer = csv.DictWriter(references, header, delimiter=",")
        writer.writeheader()
        for row in zip(lo_reference, lo_cluster, lo_fasta, lo_snp_cons, lo_bwa, lo_fai, lo_mutations_matrix):
            writer.writerow(row)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(make_references("$output", "$references".split(), "$clusters".split(), "$fastas".split(),
                             "$snp_cons".split(), "$bwas".split(), "$fais".split(), "$mutations_matrices".split()))
