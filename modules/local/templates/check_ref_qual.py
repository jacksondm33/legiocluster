#!/usr/bin/env python


"""Check ref qual."""


import csv
import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def check_ref_qual(contigs_file, output_file, report_file, MED_GENOME_LEN):
    """
    Checks if the genome in a SPAdes_contigs.fa file is of sufficient
        quality to serve as new reference genome.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: int MED_GENOME_LEN = median genome length for that species
    return: bool passed_qc = if True, then the genome assembled by SPAdes is
            of sufficient quality to serve as candidate reference genome
    output: reasons why the genome failed the QC check
    """

    MIN_LEN_FIRST_CONTIG = 150000  # minimum length of the first contig
    MIN_CONTIG_COV       =     15  # minimum read coverage per contig
    MIN_CONTIG_LEN       =   1000  # minimum contig length
    MAX_NO_TOTAL_CONTIGS =    350  # maximum number of all contigs
    MAX_NO_CONTIGS_1KB   =    200  # maximum number of contigs > 1 kb

    lo_contig_data = []
    lo_cov = []
    lo_len = []
    lo_len_roc = []
    passed_qc = True
    msg = ''

    # extracting the data
    with open(contigs_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                line = line.rstrip('\\n')
                # >NODE_1_length_238256_cov_41.824755
                data = line.split('_')
                # contig number, length, coverage; e.g.:
                lo_contig_data.append((int(data[1]), int(data[3]),
                                       float(data[5])))

    # make list of contig lengths and coverages
    for contig in lo_contig_data:
        # all contigs >= 1kb in length, regardless of coverage
        if contig[1] >= MIN_CONTIG_LEN:
            lo_len_roc.append(contig[1])
        # all contigs >= 1kb in length and >= 7.5x coverage
        if contig[1] >= MIN_CONTIG_LEN and contig[2] >= MIN_CONTIG_COV:
            lo_len.append(contig[1])
            lo_cov.append(contig[2])

    # want first contig to be at least MIN_LEN_FIRST_CONTIG long
    if lo_contig_data[0][1] <= MIN_LEN_FIRST_CONTIG:
        passed_qc = False
        msg += '\\nFAIL! The first contig is too short (<'\
        + str(MIN_LEN_FIRST_CONTIG) + ').'
    # want first contig of sufficient coverage
    if lo_contig_data[0][2] <= MIN_CONTIG_COV:
        passed_qc = False
        msg += '\\nFAIL! The coverage of the first contig is insufficient (<'\
        + str(MIN_CONTIG_COV) + '-fold).'
    # if the sum of all contigs is too small compared to the median genome
    #  length for that species
    if sum(lo_len) < (MED_GENOME_LEN * 0.9):
        passed_qc = False
        msg += '\\nFAIL! The sum of all contig lengths above min coverage is '\
        + str(sum(lo_len))\
        + ', which is less than 90% of the median size for that species. '
    # too many contigs in total
    if len(lo_contig_data) > MAX_NO_TOTAL_CONTIGS:
        passed_qc = False
        msg += '\\nFAIL! There are too many total contigs ('\
        + str(len(lo_contig_data)) + '). '
    # too many contigs in above 1kb
    if len(lo_len) > MAX_NO_CONTIGS_1KB:
        passed_qc = False
        msg += 'FAIL! There are too many contigs >1kb (' + str(len(lo_len))\
        + '). '
    # WARNING ONLY
    # check that all contigs of sufficient length have sufficient coverage
    if len(lo_len) != len(lo_len_roc):
        ## seen low cov contigs that were correct species, so don't fail
        #passed_qc = False
        msg += '\\nWARNING: Some contigs >= ' + str(MIN_CONTIG_LEN)\
        + ' bp have ' + 'less than desired coverage (< '\
        + str(MIN_CONTIG_COV) +'x). '

    # prepare text for the report and log file
    if passed_qc:
        text = '\\nNote:\\nThe isolate passed the QC check for new references. '\
        + msg + '\\nPlease run a Blast search on the isolate.fa file to make '\
        + 'sure the sample is not contaminated.'\
        + '\\nhttps://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch'
    else:
        text = '\\nThe isolate FAILED the QC check for new references,'\
        + ' here is why:\\n' + msg + '\\nNOTE: This sequence needs to be added'\
        + ' manually to a folder with similar genomes.'

    # write to report
    with open(report_file, 'a') as report:
        print(text, file=report)

    with open(output_file, 'a', newline='') as output:
        output_writer = csv.writer(output)
        output_writer.writerow([passed_qc])


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(check_ref_qual("$contigs", "$output", "$report", int("$med_genome_len")))
