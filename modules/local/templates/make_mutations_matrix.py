#!/usr/bin/env python


"""Make mutations matrix."""


import csv
import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def write_csv(outfile, lo_rows):
    """
    Writes a list of rows to a CSV file, where each item in a row will be
      a cell in the table.
    param: str outfile = path to and name of the CSV file
    param: list lo_rows = list of data
    """

    with open(outfile, 'w', newline='') as outfile:  # write to csv
        csv_writer = csv.writer(outfile)
        for row in lo_rows:
            csv_writer.writerow(row)


def write_mutations_matrix(mutations_matrix_file, lo_pairwise_diffs):
    """
    Takes a list of pairwise differences, and writes a matrix to file:
        V1 (V2, V3, V4)
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str VCF_path_ref = path the reference's VCF-folder
    param: list lo_pairwise_diffs = list of tuples:
           [(G1, G2, V1, V2, V3, V4), (G2, G1, V1, V2, V3, V4), ...]
    return: the mutations_matrix
    output: the mutations_matrix written to a .csv file,
            if print_to_file == True
    """

    lo_isolates = [isolate[0] for isolate in lo_pairwise_diffs] # extract isolates
    so_isolates = set(lo_isolates)  # remove duplicates
    slo_isolates = sorted(list(so_isolates))  # sort list

    # make a dictionary of isolate:position pairs
    do_posns = {}
    for n, isolate in enumerate(slo_isolates):
        do_posns[isolate] = n+1

    x = len(slo_isolates)  # number of isolates+1 = dimension of matrix

    # make matrix of dimension x+1 * x+2 and fill with 'nd' (= not determined)
    mutations_matrix = [['nd' for i in range(x+1)] for j in range(x+1)]

    # set upper left corner, mutations_matrix[0][0], to empty
    mutations_matrix[0][0] = ''

    # fill the diagonals with '0'
    for i in range(1,x+1):
        mutations_matrix[i][i] = '0 (0, 0, 0)'

    # add row headers
    for p, isolate in enumerate(slo_isolates):
        mutations_matrix[p+1][0] = isolate

    # add column headers
    for q, isolate in enumerate(slo_isolates):
        mutations_matrix[0][q+1] = isolate

    # get positions from do_posns, then fill matrix
    for pair in lo_pairwise_diffs:
        G1, G2, V1, V2, V3, V4 = pair
        # format V1 (V2, V3, V4)
        data = str(V1)+' (' + str(V2)+', ' + str(V3)+', ' + str(V4)+')'

        fst_posn = do_posns.get(G1, None)
        snd_posn = do_posns.get(G2, None)
        mutations_matrix[fst_posn][snd_posn] = data

    # write the mutations matrix to file
    write_csv(mutations_matrix_file, mutations_matrix)
    logger.info('Made mutations_matrix.csv')

    return mutations_matrix


def get_identicals(lo_pairwise_diffs, USE_SNPs):
    """
    Takes the data from a SNP-matrix and returns a list of the names of
      isolate pairs that have zero differences (mutation events [default] or
      SNPs) between each other.
    param: list lo_pairwise_diffs = list of tuples (G1, G2, V1, V2, V3, V4)
    return: list of isolate pairs that have different names and that have a
            mutation-event count (default, else SNP count) of zero
    """

    lo_identicals = []

    for pair in lo_pairwise_diffs:
        G1, G2, V1, V2, V3, V4 = pair

        # use V1 = mutation events by default, use V4 = SNPs as needed
        V_metric = V1
        if USE_SNPs:
            V_metric = V4

        # the two isolates are the same or have more than 0 differences
        if G1 == G2 or V_metric > 0:
            continue
        # isolates G1 and G2 are not the same, but have 0 differences
        else:
            lo_identicals.append([G1, G2])

    return sorted(lo_identicals)


def combine_identicals(lo_identicals):
    """
    Takes a sorted list of isolate pairs and combines all those that have zero
      mutation events or SNPs in common.
    param: list lo_identicals = list of lists, where each sublist includes a
           pair of isolate names that share zero mutation events or SNPs
    return: list of lists, where all isolates that share zero mutation events
            or SNPs are combined into one list
    """

    # list of lists, where each sublist contains two or more isolates with
    # zero variants
    lo_comb_ident = []

    # list.pop() removes one list [pair of isolates] at a time
    while lo_identicals != []:
        G1, G2 = lo_identicals.pop()

        # at the beginning, lo_comb_ident is empty, so add the first pair
        if lo_comb_ident == []:
            lo_comb_ident.append([G1, G2])

        # enumerate(lo_comb_ident) returns the count and a sublist from
        # lo_comb_ident
        for count, lo_results in enumerate(lo_comb_ident):
            # if pair already present in the results' sublist, move on
            if G1 in lo_results and G2 in lo_results:
                break
            # if one of the two is present, add the missing one
            elif G1 in lo_results and G2 not in lo_results:
                lo_results.append(G2)
            elif G2 in lo_results and G1 not in lo_results:
                lo_results.append(G1)
            # the two items of the pair are not present in any sublist and
            # lo_comb_ident has been exhausted: add isolate pair as new list
            # to lo_comb_ident
            elif count+1 == len(lo_comb_ident):
                lo_comb_ident.append([G1, G2])

    return sorted(lo_comb_ident)


def concat_identicals(lo_pairwise_diffs, lo_comb_ident, USE_SNPs):
    """
    Combines a list, lo_pairwise_diffs, of isolate pairs with their mutation
      events counts (default, else SNPs), and a list, lo_comb_ident, of
      isolates that share zero differences:
    - converts each list in lo_comb_ident into a concatenated string of names
    - replaces each isolate name in lo_pairwise_diffs with the concatenated
       name, if applicable
    - returns only pairs of isolates that have more than zero differences,
      where groups of isolates with zero differences are represented by their
      concatenated name
    param: lo_pairwise_diffs = isolate pairs with their indel and SNP counts
    param: lo_comb_ident = list of isolates that share zero differences
    return: list of isolates that have more than zero differences
    """

    # converts the list of isolate names into a string of concatenated names,
    #  separated by a newline
    lo_str_ident = [','.join(sorted(ident)) for ident in lo_comb_ident]

    lo_comb = []

    for pair in lo_pairwise_diffs:
        G1, G2, V1, V2, V3, V4 = pair

        V_metric = V1
        if USE_SNPs:
            V_metric = V4

        # replace G1 with the concatenated name, if applicable
        for ident1 in lo_str_ident:
            if G1 in ident1:
                G1 = ident1

        # replace G2 with the concatenated name, if applicable
        for ident2 in lo_str_ident:
            if G2 in ident2:
                G2 = ident2

        # add to returned list if isolate names are not identical and the entry
        #  is not already present (prevents duplicates)
        if G1 != G2 and [G1, G2, V_metric] not in lo_comb:
            lo_comb.append([G1, G2, V_metric])
        if G1 != G2 and [G2, G1, V_metric] not in lo_comb:
            lo_comb.append([G2, G1, V_metric])

    return lo_comb


def process_data(lo_rows, MODE):
    """
    Takes the list of rows and replaces the text in each cell with a modified
      text: removes suffixes from isolate names, and restricts the data to
      either MEs or SNPs
          isolate-name_suffix -> isolate-name
          ME (IDE, BID, SNP)  -> ME
          ME (IDE, BID, SNP)  -> SNP
    param: list lo_rows = list of original data from the CSV file, were each
           row is a list and each data cell is a string
    param: str MODE = determines the output data, either 'SNP' or 'ME'
    return: list lo_mod_rows = list of rows with modified data cells
    """

    lo_mod_rows = []
    for row in lo_rows:
        lo_mod_cells = []
        for cell in row:
            # extracts ME or SNP from the data cells
            if '(' in cell and ')' in cell:
                # e.g.: 23 (5, 9, 18) -> 18
                if MODE == 'SNP':
                    cell = cell.split()[3][:-1]
                # e.g.: 23 (5, 9, 18) -> 23
                elif MODE == 'ME':
                    cell = cell.split()[0]
            # removes suffixes from the isolate names
            # e.g.: IDR2000166282-01-00_S78 -> IDR2000166282-01-00
            elif cell.startswith('IDR') and '_' in cell:
                cell = cell.split('_')[0]
            lo_mod_cells.append(cell)
        lo_mod_rows.append(lo_mod_cells)
    return lo_mod_rows


def make_mutations_matrix(lo_cluster_pairwise_diffs, mutations_matrix_file,
                          snp_matrix_file, me_matrix_file,
                          concat_pairwise_snps_file, concat_pairwise_mes_file):
    """
    Main function: Compares 'SNP_cons.txt' files in a folder and returns for
      each pair of isolates a list of two tuples: [(G1, G2, V1, V2, V3, V4),
      (G2, G1, V1, V2, V3, V4)], where:
      - G1 and G2 are the names of the two isolates,
      - V1 = mutation events (= indel events + SNPs),
      - V2 = indel events (number of INS or DEL, regardless of length),
      - V3 = count of bases in indels,
      - V4 = SNPs
      Note: the two tuples are the same except for the order of G1 and G2,
        which is needed later to make the MST
      Note: the genome length might be different for each isolate due to
        INDELs; but the number of list entries for that genome (loci) should
        always be the same as the reference genome
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str isolate = isolate name, e.g.: 'IDR001234'
    param: str ref_fa_file = name of a reference strain's FASTA file
    param: str SS_dir = species-specific directory, e.g.: 'Lpn/'
    return: list of lo_concat_pairwise_MEs = isolate names and number of
            mutation events
    return: list of lo_concat_pairwise_SNPs = isolate names and number of SNPs
    output: new or updated 'mutations_matrix.csv' file
    """

    so_pairwise_diffs = set()

    for pairwise_diffs_file in lo_cluster_pairwise_diffs:
        with open(pairwise_diffs_file, newline='') as infile:
            reader = csv.reader(infile)
            for G1, G2, V1, V2, V3, V4 in reader:
                so_pairwise_diffs.add((G1, G2, int(V1), int(V2), int(V3), int(V4)))
                so_pairwise_diffs.add((G2, G1, int(V1), int(V2), int(V3), int(V4)))

    lo_pairwise_diffs = list(so_pairwise_diffs)

    # writes the mutations matrix, [V1 (V2, V3, V4)], to the
    # 'mutations_matrix.csv' file
    mutations_matrix = write_mutations_matrix(mutations_matrix_file,
                                              lo_pairwise_diffs)
    logger.info('Added or updated the mutations matrix.')

    # generate ME- and SNP-matrices from mutations_matrix
    lo_SNP_rows = process_data(mutations_matrix, 'SNP')
    write_csv(snp_matrix_file, lo_SNP_rows)

    lo_ME_rows = process_data(mutations_matrix, 'ME')
    write_csv(me_matrix_file, lo_ME_rows)

    # formatting the data for the MST: the next three functions combine
    # isolate pairs with zero indels events + SNPs to de-clutter the MST
    # returns list of isolate pairs that have zero events + SNPs between them

    ##### 1. run once to get data for mutation events #########################
    lo_ident_isol_MEs = get_identicals(lo_pairwise_diffs, False)
    logger.info(lo_ident_isol_MEs)

    # combines all identical isolates
    lo_comb_ident_ME = combine_identicals(lo_ident_isol_MEs)
    logger.info(lo_comb_ident_ME)

    # fuses list entries for identical isolates (ME = mutation event)
    lo_concat_pairwise_MEs = concat_identicals(lo_pairwise_diffs,
                                               lo_comb_ident_ME, False)
    logger.info(lo_concat_pairwise_MEs)

    ##### 2. run again to get data for SNPs only ##############################
    lo_ident_isol_SNPs = get_identicals(lo_pairwise_diffs, True)
    logger.info(lo_ident_isol_SNPs)

    # combines all identical isolates
    lo_comb_ident_SNP = combine_identicals(lo_ident_isol_SNPs)
    logger.info(lo_comb_ident_SNP)

    # fuses list entries for identical isolates (ME = mutation event)
    lo_concat_pairwise_SNPs = concat_identicals(lo_pairwise_diffs,
                                                lo_comb_ident_SNP, True)
    logger.info(lo_concat_pairwise_SNPs)

    # write the uncluttered list of isolate pairs, where isolates with zero
    # indels + SNPs have been concatenated with newlines
    write_csv(concat_pairwise_snps_file, lo_concat_pairwise_SNPs)
    write_csv(concat_pairwise_mes_file, lo_concat_pairwise_MEs)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(make_mutations_matrix("$cluster_pairwise_diffs".split(), "$mutations_matrix",
                                   "$snp_matrix", "$me_matrix", "$concat_pairwise_snps",
                                   "$concat_pairwise_mes"))
