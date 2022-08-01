#!/usr/bin/env python


"""Parse qualimap output."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def parse_qualimap_output(qualimap_file, report_file):
    """
    Extracts results from the Qualimap output and writes it to report
    output: text added to report
    """

    # opens the report from Qualimap and writes selected lines to report
    with open(qualimap_file, 'r') as in_file:

        # opens the report file
        with open(report_file, 'a') as report:

            # adds a header to report
            print('\\n\\nMapping quality check (Qualimap results):', \
                  file=report)

            # writing specific data to the report
            for line in in_file:
                line = line.rstrip('\\n')
                if line.startswith('     number of bases')\
                or line.startswith('     number of contigs')\
                or line.startswith('     number of reads')\
                or line.startswith('     number of mapped reads')\
                or line.startswith('     number of mapped bases')\
                or line.startswith('     mean mapping quality'):
                    print(line.replace('     ',''), file=report)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(parse_qualimap_output("$qualimap", "$report"))
