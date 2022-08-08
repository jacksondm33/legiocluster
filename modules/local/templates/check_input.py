#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import csv
import logging
import platform
import sys
import yaml
from collections import Counter
from pathlib import Path


logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FASTQ_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    VALID_FASTA_FORMATS = (
        ".fa",
        ".fasta",
    )

    def __init__(self, id_col, fastq_cols, fasta_cols):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").

        """
        self._id_col = id_col
        self._fastq_cols = fastq_cols
        self._fasta_cols = fasta_cols
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_id(row)
        self._validate_fastq(row)
        self._validate_fasta(row)
        self._seen.add(row[self._id_col])
        self.modified.append(row)

    def _validate_id(self, row):
        """Assert that the identifier exists and convert spaces to underscores."""
        assert len(row[self._id_col]) > 0, "Identifier is required."
        # Sanitize identifier slightly.
        row[self._id_col] = row[self._id_col].replace(" ", "_")

    def _validate_fastq(self, row):
        """Assert that the FASTQ entries are non-empty and have the right format."""
        for col in self._fastq_cols:
            assert len(row[col]) > 0, "FASTQ entries must be non-empty"
            self._validate_fastq_format(row[col])

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FASTQ_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\\n"
            f"It should be one of: {', '.join(self.VALID_FASTQ_FORMATS)}"
        )

    def _validate_fasta(self, row):
        """Assert that the FASTA entries are non-empty and have the right format."""
        for col in self._fasta_cols:
            assert len(row[col]) > 0, "FASTA entries must be non-empty"
            self._validate_fasta_format(row[col])

    def _validate_fasta_format(self, filename):
        """Assert that a given filename has one of the expected FASTA extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FASTA_FORMATS), (
            f"The FASTA file has an unrecognized extension: {filename}\\n"
            f"It should be one of: {', '.join(self.VALID_FASTA_FORMATS)}"
        )

    def validate_unique_ids(self):
        """
        Assert that the identifier is unique.
        """
        assert len(self._seen) == len(self.modified), "The identifier must be unique."


def check_input(file_in, file_out, references):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,fastq_1,fastq_2
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
            SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,

    .. _viral recon samplesheet:
        https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

    """
    if references:
        required_columns = {"sample", "reference", "fasta", "snp_cons"}
    else:
        required_columns = {"sample", "fastq_1", "fastq_2", "set_ref", "make_ref"}

    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with open(file_in, newline="") as in_handle:
        reader = csv.DictReader(in_handle)
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            logger.critical(f"The csv **must** contain the column headers: {', '.join(required_columns)}.")
            sys.exit(1)
        # Validate each row.
        if references:
            checker = RowChecker("sample", [], ["fasta"])
        else:
            checker = RowChecker("sample", ["fastq_1", "fastq_2"], [])
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_ids()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with open(file_out, mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
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

    sys.exit(check_input("$input", "$input_valid", "$references" == "true"))
