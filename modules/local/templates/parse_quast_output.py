#!/usr/bin/env python


"""Parse quast output."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def parse_quast_output(contigs_file, quast_report_file, report_file, CONTIG_THRESHOLD):
    """
    Parses the Quast output and writes the results to the report.
    param: str contigs_file = name of a sequence file to be QC'd
    output: info from the Quast file will be added to 'report.txt',
            issues a Warning if too many contigs
    return int contigs = the number of all contigs
    """

    contigs = 0

    # adds a header to '_report.txt'
    with open(report_file, 'a') as report:
        print('\\n\\nAssembly quality check (Quast results) for', \
              str(contigs_file) + ':', file=report)

        # opens the report from Quast and writes data to report file
        with open(quast_report_file, mode='r') as in_file:
            for line in in_file:
                line = line.rstrip('\\n')
                print(line, file=report)
                # extracts the number of contigs and warns if too many
                if line.startswith('# contigs (>= 0 bp)'):
                    contigs = int(line.split()[-1])

        if contigs > CONTIG_THRESHOLD:
            print('\\nWARNING:', file=report)
            print(contigs, 'contigs', file=report)
            print('THE SAMPLE MIGHT BE CONTAMINATED\\n\\n', file=report)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(parse_quast_output("$contigs", "$quast", "$report", int("$contig_threshold")))
