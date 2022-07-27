#!/usr/bin/env python


"""Cleanup freebayes."""


import csv
import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def read_mpileup(mpileup_file):
    """
    Returns the content of a mpileup.vcf file as a list.
      The mpileup.vcf file generally uses one row per base unless there is an
      insertion, in which case the first row represents the single base of the
      reference and the second row represents the insertion. There can be
      unmapped regions ("N"), which will result in missing rows/bases.
    return: dict of contig + position : [ref_base, query_base] pairs, e.g.:
            {'NODE_1_1' : ['A', '<*>'],
             'NODE_1_2' : ['T', 'C,<*>'],
             'NODE_1_3' : ['G', '<*>'],
             ...}
    """

    do_mpileup_data = {}
    with open(mpileup_file, 'r') as input_file:
        for row in input_file:
            # ignore the header
            if row[0].startswith('#'):
                  continue
            # extract the data
            else:
                data = row.split('\t')
                contig = data[0]
                posn = data[1]
                ref_base = data[3]
                query_base = data[4]
                # contig name and bp-posn make a unique key for each base
                key = contig + '_' + posn
                # in case on an insertion, there will be two entries (the 1.
                # will have the wt base, the 2. the insertion): we want the
                # first entry; the if-clause prevents the 2. from overriding
                # the first
                if key not in do_mpileup_data.keys():
                    do_mpileup_data[key] = [ref_base, query_base]

    return do_mpileup_data


def translate_cigar(CIGAR):
    """
    Translates a CIGAR string from alphanumeric format to all-alphabetical for
      easy counting of all 'X', 'D', and 'I'.
    param: str CIGAR = a CIGAR string in the form '#M#X#I#D', where # is the
           count and MXID represent Match, eXchanged (mismatch), Insertion,
           and Deletion; e.g.: '1X2M2X1I3M'
    return: a cigar string with only alphabetical characters, e.g. 'XMMXXIMMM'
    """

    number = ''
    letter = ''
    tl_cigar = ''          # translated CIGAR string: 'XXX' instead of '3X'
    for c in CIGAR:        # character by character
        if c.isnumeric():  # if number
            number += c    # replace '' with the number
        elif c.isalpha():  # if letter
            letter = c
            subcigar = (int(number) * letter)  # replaces '3X' with 'XXX'
            tl_cigar += subcigar       # reassemble complete CIGAR string
            number = ''    # reset
            letter = ''    # reset
    return tl_cigar


def get_i_str(i, s, query_seq, tl_CIGAR, DEL_count):
    """
    FreeBayes output for complex mutations can be a mix of various mutations
      This function takes a translated CIGAR string, position, and mutation,
      and returns a string representing the closest matching 5' base and the
      one or more inserted bases that follow.
    param: int i = index of the inserted base in the translated CIGAR string
    param: str s = character from the CIGAR string ('I','M','D','X') at index i
    param: str query_seq = query sequence at the site of a complex mutation
    param: str tl_CIGAR = translated CIGAR string
    param: int DEL_count = number of deletions in the CIGAR string prior to i
    return: a string of bases that includes the first match and any inserted
            bases at a given index in the query, general format: MI+
    helper function to clean_up_fb_data()
    """

    I_str = query_seq[i-1-DEL_count]
    while s == 'I':
        I_str += query_seq[i-DEL_count]
        i += 1
        s = tl_CIGAR[i]
    return I_str


def clean_up_fb_data(lo_file_content):
    """
    Takes a list of data from FreeBayes and converts a "complex" mutation
      consisting of multiple mutations into individual ones; uses a translated
      CIGAR string (e.g. 'XMMX' instead of '1X2M1X')
    param: list lo_file_content = list of (contig, posn, ref_base, query_base,
           tl_CIGAR)
    return: dict of contig_posn : (posn, ref_base, query_base, tl_CIGAR)
            each mutation will be assigned the correct position;
            deletion = '-'
            an insertion of one or more bases will be preceeded by the last
            matching base, in the form: MI+
    helper function to read_freebayes_snps()
    """

    lo_added_rows = []    # new, deconvoluted mutations, to be added
    lo_deleted_rows = []  # the old, complex mutations, to be deleted

    # goes throught the list of rows one at a time, where one row is one
    # FreeBayes mutation, including complex ones, such as: 1M3I1M1X2M
    for row in lo_file_content:

        # content of a row
        contig, posn, ref_base, query_base, tl_CIGAR = row

        # need to count the number of DEL and INS to adjust the indices
        DEL_count = 0
        INS_count = 0

        # simple SNPs are 'X' => len == 1: no action needed, everything else
        # has a larger tl_CIGAR string
        if len(tl_CIGAR) > 1:

            # goes through the translated CIGAR string one charater at a time,
            # where M = match (no action needed), X = eXchange, I = Insertion,
            # D = Deletion
            for i,s in enumerate(tl_CIGAR):

                # this process is very error prone; the try/except will keep
                # the program running and return those items that need
                # correction, usually due to indexing issues
                try:
                    # if SNP, return ref_base and query_base
                    # e.g.: 'T' 'A'
                    if s == 'X':
                        new_row = (contig,
                                   str(int(posn)+i-INS_count),
                                   ref_base[i-INS_count],
                                   query_base[i-DEL_count],
                                   tl_CIGAR)
                        lo_added_rows.append(new_row)
                    # if DEL, return ref_base and '-' for the query
                    # e.g.: 'T' '-'
                    elif s == 'D':
                        new_row = (contig,
                                   str(int(posn)+i-INS_count),
                                   ref_base[i-INS_count],
                                   '-',
                                   tl_CIGAR)
                        lo_added_rows.append(new_row)
                        DEL_count += 1
                    # if INS, return last matching base at posn[i-1] for the
                    # ref_base and the last matching base followed by the
                    # inserted bases added to it for the query_base
                    # e.g.: 'T' 'TAAA'
                    elif s == 'I':
                        # skip if more than one 'I' next to each other, since
                        # that would lead to repeat entries, such as 'AAA',
                        # 'AA', 'A', in  the rows that follow
                        if tl_CIGAR[i-1] != 'I':
                            # helper function to extract the query string to
                            # insert
                            I_str = get_i_str(i, s, query_base, tl_CIGAR,
                                              DEL_count)
                            # 'i-1' because this will be an entry at the last
                            # match posn, not at the posn of the 'I'
                            new_row = (contig,
                                       str(int(posn)+i-1-INS_count),
                                       ref_base[i-INS_count-1],
                                       I_str,
                                       tl_CIGAR)
                            lo_added_rows.append(new_row)
                        INS_count += 1
                # prints data for trouble shooting in case of a failure
                except Exception as e:
                    logger.error('An Exception has occurred in cleanup_FB.clean_up_'\
                          + '''fb_data():\n''', str(e))
                    logger.error(row)
                    logger.error(new_row)
                    logger.error(i,s, INS_count, DEL_count)
                    logger.error()

            # mark the org row for deletion to prevent duplications
            lo_deleted_rows.append(row)

    # adding the new, deconvoluted rows
    lo_freebayes_data = lo_file_content + lo_added_rows

    # deleting the original rows with the complex mutations
    for r in lo_deleted_rows:
        lo_freebayes_data.remove(r)

    # sorts the rows by position: The value of the key parameter should be a
    # function that takes a single argument and returns a key to use for
    # sorting purposes. Here, that function is lambda, which converts the
    # posn to int, then sorts by that number
    lo_freebayes_data = sorted(lo_freebayes_data, key=lambda row:int(row[1]))

    # conversion to dict with contig_posn as key and ref_base, query_base,
    # tl_CIGAR as value:
    #  {'NZ_JHGY01000001.1_6': ('6', 'T', 'C', 'X'),
    #   'NZ_JHGY01000001.1_16': ('16', 'T', 'C', 'X'), ...}
    # speeds up the comparison with the mpileup data
    do_freebayes_data = {}
    for data in lo_freebayes_data:
        do_freebayes_data[data[0] + '_' + data[1]] = data[1:5]

    return do_freebayes_data


def read_freebayes_snps(vcf_file):
    """
    Returns the content of freebayes.vcf as a dictionary.
      A line in a vcf looks like this (one line per mutation):
      NZ_JHGY1.1  171  .  ACGA  GCGT  3397.62  .  AB=0;ABP=0;AC=1;AF=1;AN=1; \
      AO=106;CIGAR=1X2M1X;DP=111;
    return: a dict of
      contig-name + bp-position : (bp-position, ref-base, query-base, CIGAR)
      e.g.:  {'NZ_JHGY05.1_2202': ('2202', 'T', 'TTGTCAC', 'MIIIIIIMXMMXMIMM'),
              ...}
    """

    file_content = []
    with open(vcf_file, 'r') as input_file:
        for row in input_file:
            # ignore header rows
            if row[0].startswith('#'):
                  continue
            else:
                data = row.split('\t')
                contig = data[0]
                posn = data[1]
                ref_base = data[3]
                query_base = data[4]
                # extract the CIGAR string from a bunch of other stuff
                CIGAR = data[7].split(';')[6].split('=')[1]
                # translate from '1X3M1X' into 'XMMMX' format
                tl_CIGAR = translate_cigar(CIGAR)
                # add selected data to list
                file_content.append((contig, posn, ref_base, query_base,
                                         tl_CIGAR))

    # use helper function to deconvolute complex mutations and transform into
    # a dict with contig_name_+_posn as keys
    do_freebayes_data = clean_up_fb_data(file_content)

    return do_freebayes_data


def combine_csv(csv_file, lo_bases, do_mpileup_data, do_freebayes_data,
                isolate, reference):
    """
    Combines the data from 'mpileup.vcf' and 'freebayes.vcf' with the data
      from the reference sequence. The latter is a list that is used to make
      dictionary keys to retrieve the mpileup and freebayes data, if available.
      The combined data are written to a csv file.
    param: list lo_bases = [[contig, posn, base], ...] for the reference genome
    param: dict do_mpileup_data = contig plus position : [ref_base, query_base]
    param: dict do_freebayes_data =  contig plus position :
                                     (posn, ref_base, query_base, tl_CIGAR)
    param: str isolate = isolate name, e.g.: 'IDR001234'
    output: a CSV file combining all input data
    """

    with open(csv_file, 'w') as output:
        # generates a tab-separated csv file
        row_writer = csv.writer(output, dialect='excel-tab')
        # write the header rows
        row_writer.writerow(['# ' + reference + ' versus ' + isolate])
        row_writer.writerow(['# contig', 'posn', 'ref_base', 'mp-ref',
                             'mp-query', 'fb-posn', 'fb-ref', 'fb-query',
                             'fb-CIGAR'])

        # one base in the ref seq at a time
        for base in lo_bases:

            # each row is a list, to which data are added
            combined_rows = base

            # the key is the contig plus posn
            key = base[0] + '_' + base[1]
            # retrieve REF and ALT data from mpileup file, or ['',''] if
            # no data available, then add to the row
            mp_data = do_mpileup_data.get(key,['',''])
            combined_rows.extend(mp_data)
            # same for freebayes data
            fb_data = do_freebayes_data.get(key,['','',''])
            combined_rows.extend(fb_data)

            # write to cvsv file
            row_writer.writerow(combined_rows)


def make_consensus(snp_cons_file, csv_file, isolate, reference, diagnostic_mode):
    """
    Takes a "combined_files.csv" file made by combine_csv() and determines a
    consensus sequence for the query from the mpileup and the FreeBayes data:
    - if both or only FreeBayes calls it a mutation, it's a mutation
    - if only mpileup calls it a mutation, it's ambiguous ('n') (the mutation
      might be below FreeBayes threshold values)
    - if no mpileup data are available for that positionn, it's unmapped ('N')
    - note that some reference genomes can include the letter 'N'
    - if deletion, insert '-'
    - if insertion, add the inserted bases behind the last matching base, e.g.:
      insertion of 'CT' after 'A': 'ACT'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    output: a '_SNP_cons.txt' file that contains the combined mutation data
    """

    with open(snp_cons_file, 'w') as outfile:
        # write a header
        print('# SNPs and INDELs after comparing strain', reference,
              'with', isolate,
              '(Based on bcftools mpileup and FreeBayes data.)',
              file=outfile)

        with open(csv_file, 'r') as infile:
            for line in infile:
                line = line.rstrip('''\n''')
                data = line.split()
                # ignore header
                if data[0] == '#':
                    continue
                # adding data
                else:
                    # 3 = contig, posn, base of ref, but the position was not
                    # mapped by mpileup or FreeBayes
                    if len(data) == 3:
                        contig, posn, ref = data
                        if diagnostic_mode:
                            print(contig, posn, ref, 'N', 'unmapped',
                                  file=outfile)
                        else:
                            print('N', file=outfile)

                    # mpileup, but no FreeBayes data, are available
                    elif len(data) == 5:
                        contig, posn, ref, mp_ref, mp_alt = data
                        # query base is the same as the reference base
                        if (ref == mp_ref) and (mp_alt == '<*>'):
                            if diagnostic_mode:
                                print(contig, posn, ref, ref, file=outfile)
                            else:
                                print(ref, file=outfile)
                        # query different from reference, but no FreeBayes data
                        # to support it => possible mutation is below threshold
                        else:
                            if diagnostic_mode:
                                print(contig, posn, ref, 'n', 'unsupported',
                                      file=outfile)
                            else:
                                print('n', file=outfile)

                    # FreeBayes, but no mpileup data, are available
                    elif len(data) == 7:
                        contig, posn, ref, fb_posn, fb_ref, fb_alt, cigar \
                        = data
                        # FreeBayes would only list something, fb_alt, if it
                        # was different from the reference
                        if diagnostic_mode:
                            print(contig, posn, ref, fb_posn, fb_ref, fb_alt,
                                  cigar, file=outfile)
                        else:
                            print(fb_alt, file=outfile)

                    # we got mpileup and FreeBayes data for that posn
                    elif len(data) == 9:
                        contig, posn, ref, mp_ref, mp_alt, \
                        fb_posn, fb_ref, fb_alt, cigar = data
                        # write down FreeBayes output for that posn
                        if diagnostic_mode:
                            print(contig, posn, ref, fb_alt, cigar,
                                  file=outfile)
                        else:
                            print(fb_alt, file=outfile)

                    # catching all else
                    else:
                        if diagnostic_mode:
                            print('unexpected number of data',
                                  file=outfile)
                        else:
                            print('N', file=outfile)


def cleanup_freebayes(mpileup_file, vcf_file, vcf_path, bases_file, csv_file, isolate, reference, diagnostic_mode):
    """
    main function
    param: str isolate = isolate name, e.g.: 'IDR001234'
    output: a '_SNP_cons.txt' file added to the /VCF_folder
    """

    snp_cons_file = vcf_path / (isolate + '_SNP_cons.txt')

    lo_bases = []
    with open(bases_file, 'r', newline='') as bases:
        bases_reader = csv.reader(bases)
        for entry in bases_reader:
            lo_bases.append(entry)

    # returns dict of mpileup data
    do_mpileup_data = read_mpileup(mpileup_file)

    # returns dict of freebayes data
    do_freebayes_data = read_freebayes_snps(vcf_file)

    # makes csv file from ref seq, mpileup, and freebayes data
    combine_csv(csv_file, lo_bases, do_mpileup_data, do_freebayes_data, isolate, reference)

    # converts csv file into consensus <isolate>_SNP_cons.txt file
    make_consensus(snp_cons_file, csv_file, isolate, reference, diagnostic_mode)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(cleanup_freebayes("$mpileup", "$vcf", Path("$vcfs"), "$bases", "$csv", "$meta.id", "$meta.ref", bool("$diagnostic_mode")))

