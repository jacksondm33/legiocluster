#!/usr/bin/env python


"""Extract fastqc results."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def extract_fastqc_results(lo_reads_files, lo_fastqc_results, report_file):
    """
    Extracts basic statistics and summary results from the 'fastqc_data.txt'
      and 'summary.txt' files and writes them to the report.
    param: str reads_file = (path and) name of file with the raw forward or
           reverse reads, e.g.: "IDR200001234.fastq.gz"
    param: str fastqc_results = name of directory with fastqc results
    param: str report_file = output report file
    output: writes basic statistics to the report file
    """

    for reads_file, fastqc_results in zip(lo_reads_files, lo_fastqc_results):

        # extract selected data from the 'fastqc_data.txt' file
        with open(report_file, 'a') as report:
            print('\\nRead quality control (FastQC results):', file=report)
            # Original name of the file with the raw reads
            print('Results for processed reads from:', reads_file, file=report)
            with open(Path(fastqc_results) / 'fastqc_data.txt', mode='r') as infile_1:
                for line in infile_1:
                    line = line.rstrip('\\n')
                    if line.startswith('Filename') \
                    or line.startswith('Total Sequences') \
                    or line.startswith('Sequences flagged') \
                    or line.startswith('Sequence length') \
                    or line.startswith('%GC'):
                        print(line, file=report)

        # extract all data from the 'summary.txt' file
        lo_qc_results = []
        with open(report_file, 'a') as report:
            with open(Path(fastqc_results) / 'summary.txt', mode='r') as infile_2:
                for line in infile_2:
                    line = line.rstrip('\\n')
                    # remove read_file name in each line
                    qc_result, what = line.split('	')[:2]
                    lo_qc_results.append(qc_result)
                    print(qc_result + '   ' + what, file=report)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(extract_fastqc_results("$reads".split(), "$fastqc_results".split(), "$report"))
