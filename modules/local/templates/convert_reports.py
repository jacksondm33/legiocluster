#!/usr/bin/env python


"""Convert reports."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def convert_report(report_file, html_file):
    """
    Adds a line of text from the report.txt file to the html file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str line = line of text from the report.txt file
    """

    with open(report_file, 'r') as infile:

        with open(html_file, 'a') as html:

            use_desc = False
            use_table = False

            for line in infile:

                line = line.rstrip('\\n')

                if line.count('\t') == 1:
                    if not use_desc:
                        if use_table:
                            print('</table>', file=html)
                            use_table = False
                        print('<dl class="dl-horizontal">', file=html)
                        use_desc = True
                    name, desc = line.split('\t')
                    print('<dt>' + name + '</dt><dd>' + desc + '</dd>', file=html)
                elif line.count('\t') > 1:
                    if not use_table:
                        if use_desc:
                            print('</dl>', file=html)
                            use_desc = False
                        print('<table>', file=html)
                        print('<tr>', file=html)
                        for element in line.split('\t'):
                            print('<th>' + element + '</th>', file=html)
                        print('</tr>', file=html)
                        use_table = True
                    else:
                        print('<tr>', file=html)
                        for element in line.split('\t'):
                            print('<td>' + element + '</td>', file=html)
                        print('</tr>', file=html)
                else:
                    if use_table:
                        print('</table>', file=html)
                        use_table = False
                    elif use_desc:
                        print('</dl>', file=html)
                        use_desc = False
                    print('<p>' + line + '</p>', file=html)


def convert_reports(lo_reports):
    """Convert reports"""

    for report_file in lo_reports:
        convert_report(report_file, report_file[:-4] + '.html')


if __name__ == "__main__":
    logging.basicConfig(filename="$log_file", level="$log_level", format="[%(levelname)s] %(message)s")

    versions = {}
    versions["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }
    with open("versions.yml", "w") as f:
        yaml.dump(versions, f, default_flow_style=False)

    sys.exit(convert_reports("$reports".split()))
