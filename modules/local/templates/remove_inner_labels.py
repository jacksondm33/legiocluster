#!/usr/bin/env python


"""Remove inner labels."""


import logging
import platform
import re
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def remove_inner_labels(svg_file, output_file):
    """
    Inner labels are numbers in the tree that are annoying because they
      sometimes overlap with isolate names. This function finds and removes
      them from the 'parsnp_tree.svg' file.
    param: str seed = combination of sp_abbr and reference, e.g.: 'Lpn/F4468/'
    output: 'parsnp_tree.svg' without inner labels
    """

    # regular expression to find the inner labels
    label = "<text class='inner-label' x='\\d+.\\d+' y='\\d+.\\d+'>\\d+.\\d+</text>"
    lo_inner_labels = []

    # find all inner labels
    with open(svg_file, 'r') as infile:
        for line in infile:
            if line.startswith('<g'):
                lo_inner_labels = re.findall(label, line)

    # captures the text in a svg file
    lo_lines = []
    with open(svg_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            lo_lines.append(line)

    # replacing the labels with ''
    for label in lo_inner_labels:
        lo_lines[6] = lo_lines[6].replace(label, '')
    # over-riding the 'parsnp_tree.svg' with the label-free text
    with open(output_file, 'w') as outfile:
        for line in lo_lines:
            print(line, file=outfile)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(remove_inner_labels("$svg", "$output"))
