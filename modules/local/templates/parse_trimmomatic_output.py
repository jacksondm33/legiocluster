#!/usr/bin/env python


"""Parse trimmomatic output."""


import csv
import logging
import numpy as np
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def parse_trimmomatic_output(trimlog_file, output_file, report_file, min_reads):
    """
    Extracts data from the trimmomatic log file and writes them to the report.
    param: str trimlog_file = trimmomatic log file
    param: str report_file = output report file
    output: data added to report file
    """

    prev_read = ('', 0)    # name and trimmed length of the previous read

    input_read_pairs = 0   # all read pairs
    both_surviving = 0     # f- and r-reads surviving
    f_only_surviving  = 0  # only the forward read survived
    r_only_surviving  = 0  # only the reverse read survived
    dropped = 0            # none of the reads surviving

    lo_f_length_distr = [] # lengths distributions for trimmed forward reads
    lo_f_trim_5 = []       # lengths bases trimmed at 5' of forward reads
    lo_f_trim_3 = []       # lengths bases trimmed at 3' of forward reads

    lo_r_length_distr = [] # lengths distributions for trimmed reverse reads
    lo_r_trim_5 = []       # lengths bases trimmed at 5' of reverse reads
    lo_r_trim_3 = []       # lengths bases trimmed at 3' of reverse reads

    # extracts data from the trimmomatic log file
    with open(trimlog_file, 'r') as log:
        for line in log:
            line = line.rstrip('\\n')
            line_content = line.split(' ')

            # example for reads created by ART:
            #   pLPP_var-55400/1 249 0 249 1
            #   pLPP_var-55400/2 250 0 250 0
            # returns as name:
            #   pLPP_var-55400/1
            #   pLPP_var-55400/2
            if len(line_content) == 5:
                name, trim_length, lost_5, loc, lost_3 = line_content

            # example for actual reads:
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711 1:N:0:7 0 0 0 0
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711 2:N:0:7 0 0 0 0
            # returns as name:
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711/1
            #   M01698:26:000000000-BD7TF:1:1101:18858:1711/2
            elif len(line_content) == 6:
                name1, name2, trim_length, lost_5, loc, lost_3 = line_content
                name = name1 + '/' + name2[0]

            # example for reads downloaded from NCBI:
            # SRR6902774.1 1 length=251 251 0 251 0
            # SRR6902774.1 1 length=251 129 0 129 122
            # SRR6902774.2 2 length=250 250 0 250 0
            elif len(line_content) == 7:
                base_name, count, in_length, trim_length, lost_5, loc, lost_3 \
                = line_content
                if prev_read[0] == base_name:
                    name = base_name + '/2'
                else:
                    name = base_name + '/1'

            base_name = name[:-2]
            trim_length = int(trim_length)
            lost_5 = int(lost_5)
            lost_3 = int(lost_3)

            # total number of read pairs going in
            input_read_pairs += 0.5  # count half a pair

            # compares two reads: if they have the same name base (without the
            # '/1' and '/2'), they are counted if their length is >0 after the
            # trimming
            if prev_read[0] == base_name:
                if (prev_read[1] > 0) and (trim_length > 0):
                    both_surviving += 1
                elif (prev_read[1] > 0) and (trim_length == 0):
                    f_only_surviving += 1
                elif (prev_read[1] == 0) and (trim_length > 0):
                    r_only_surviving += 1
                elif (prev_read[1] == 0) and (trim_length == 0):
                    dropped += 1

            # collect list with read-lengths, number of bases trimmed at 5'
            # and at 3' for forward and reverse reads
            if name.endswith('/1'):
                lo_f_length_distr.append(trim_length)
                if trim_length > 0:
                    lo_f_trim_5.append(lost_5)
                    lo_f_trim_3.append(lost_3)
            if name.endswith('/2'):
                lo_r_length_distr.append(trim_length)
                if trim_length > 0:
                    lo_r_trim_5.append(lost_5)
                    lo_r_trim_3.append(lost_3)

            # updating the read for ther next comparison
            prev_read = base_name, trim_length

    # write data to report file
    with open(report_file, 'a') as report:
        print('Read pre-processing (Trimmomatic):', file=report)
        print('Adapters removed, low quality (< Q20) regions removed, short reads (<100) removed, poly-G (>25) removed', file=report)
        print('Input read pairs:\t', int(input_read_pairs),\
              sep='', file=report)
        print('Both surviving:\t', both_surviving,\
              ' (', round(both_surviving*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Forward only surviving:\t', f_only_surviving,\
              ' (', round(f_only_surviving*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Reverse only surviving:\t', r_only_surviving, \
              ' (', round(r_only_surviving*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Dropped read pairs:\t', dropped,\
              ' (', round(dropped*100/input_read_pairs, 2), '%)',\
              sep='', file=report)
        print('Mean (SD) lengths of trimmed F reads:\t',\
              round(np.mean(lo_f_length_distr), 2), \
              ' (', round(np.std(lo_f_length_distr), 3), ')',\
              sep='', file=report)
        print('Mean (SD) lengths of trimmed R reads:\t',\
              round(np.mean(lo_r_length_distr), 2),\
              ' (', round(np.std(lo_r_length_distr), 3), ')',\
              sep='', file=report)
        print("Mean (SD) no. of bases trimmed from 5' of F reads(*):\t",\
              round(np.mean(lo_f_trim_5), 2),\
              ' (', round(np.std(lo_f_trim_5), 3), ')', sep='', file=report)
        print("Mean (SD) no. of bases trimmed from 5' of R reads(*):\t",\
              round(np.mean(lo_r_trim_5), 2),\
              ' (', round(np.std(lo_r_trim_5), 3), ')', sep='', file=report)
        print("Mean (SD) no. of bases trimmed from 3' of F reads(*):\t",\
              round(np.mean(lo_f_trim_3), 2),\
              ' (', round(np.std(lo_f_trim_3), 3), ')', sep='', file=report)
        print("Mean (SD) no. of bases trimmed from 3' of R reads(*):\t",\
              round(np.mean(lo_r_trim_3), 2),\
              ' (', round(np.std(lo_r_trim_3), 3), ')', sep='', file=report)
        print('(*) if trimmed read length > 0', file=report)

    max_read_len = max(lo_f_length_distr)

    if both_surviving < min_reads:
        logger.error("Not enough reads surviving after Trimmomatic.")
        sys.exit(2)

    with open(output_file, 'a', newline='') as output:
        output_writer = csv.writer(output)
        output_writer.writerow([both_surviving])
        output_writer.writerow([max_read_len])


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(parse_trimmomatic_output("$trimlog", "$output", "$report", int("$min_reads")))
