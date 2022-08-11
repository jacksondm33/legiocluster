#!/usr/bin/env python


"""Parse BWA output."""


import csv
import logging
import numpy as np
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def calculate_frag_len(sam_file, report_file):
    """
    Extracts the fragment lengths from a SAM file and writes the mean and
    other stats to the report.
    param: str sam_file = input sam file
    output: text added to report
    """

    lo_frag_lens = []

    with open(sam_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            if line.startswith('@'):
                continue
            else:
                frag_len = line.split()[8]
                if frag_len.startswith('-'):
                    continue
                else:
                    lo_frag_lens.append(int(frag_len))

    with open(report_file, 'a') as report:
        print('\\n\\nGenomic fragments:', file=report)
        print('Smallest fragment:\t', min(lo_frag_lens), file=report)
        print('Mean length:\t', round(np.mean(lo_frag_lens), 2), file=report)
        print('S.D.:\t', round(np.std(lo_frag_lens), 2), file=report)
        print('median:\t', np.median(lo_frag_lens), file=report)
        print('Largest fragment:\t', max(lo_frag_lens), file=report)


def parse_bwa_output(reference_file, sam_file, flagstat_file, idxstats_file, output_file, report_file, MAPPED_THRESHOLD):
    """Parse flagstat file."""

    with open(report_file, 'a') as report:
        print('\\n\\nMapping the query against strain ' + str(reference_file) + ' (BWA MEM):', file=report)

        with open(flagstat_file, 'r') as flagstat:
            flagstat_data = ''.join(flagstat.readlines())
            percent_mapped = float(flagstat_data.split('mapped (')[1].split('%')[0])

        print('\\nAlignment QC (Samtools flagstat):', file=report)
        print(flagstat_data, file=report)
        print('\\n\\nPercentage of mapped reads:', percent_mapped, file=report)

        if percent_mapped <= MAPPED_THRESHOLD:
            print('\\nNOTE:\\nPercentage of mapped reads below threshold.\\n'\
                + 'Adding the isolate to the list of candidate reference '\
                + 'genomes.',
                  file=report)

        with open(idxstats_file, 'r') as idxstats:
            idxstats_data = ''.join(idxstats.readlines())

        print('\\n\\nAlignment QC (Samtools idxstats):', file=report)
        print('ref_fa_file\tlen\tmapped\tunmapped', file=report)
        print(idxstats_data, file=report)

    calculate_frag_len(sam_file, report_file)

    with open(output_file, 'a', newline='') as output:
        output_writer = csv.writer(output)
        output_writer.writerow([percent_mapped])


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(parse_bwa_output("$fasta", "$sam", "$flagstat", "$idxstats", "$output", "$report", float("$mapped_threshold")))
