#!/usr/bin/env python


"""Check percent mapped."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def check_percent_mapped(percent_mapped, min_percent_mapped):
    """Check percent mapped."""
    # QC check in case there are too many unmapped reads
    if percent_mapped < min_percent_mapped:
        logger.error('There were only ' + str(percent_mapped) \
            + '% mapped reads, which is far too few.')
        sys.exit(2)


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(check_percent_mapped(float("$percent_mapped"), float("$min_percent_mapped")))
