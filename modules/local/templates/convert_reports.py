#!/usr/bin/env python


"""Convert reports."""


import logging
import platform
import sys
import yaml
from pathlib import Path


logger = logging.getLogger()


def make_report_html(report_file):

    use_desc = False
    use_table = False
    section_name = ''
    html = []

    with open(report_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\\n')
            if line == '':
                continue
            if section_name == '':
                section_name = line[:-1]
                continue
            if line.count('\t') == 1:
                if not use_desc:
                    if use_table:
                        html.append('</table>')
                        use_table = False
                    html.append('<dl class="dl-horizontal" style="width:100%">')
                    use_desc = True
                name, desc = line.split('\t')
                if name.endswith(':'):
                    name = name[:-1]
                html.append('<dt>' + name + '</dt><dd>' + desc + '</dd>')
                continue
            if line.count('\t') > 1:
                if not use_table:
                    if use_desc:
                        html.append('</dl>')
                        use_desc = False
                    html.append('<table class="table" style="width:100%">')
                    html.append('<tr>')
                    for element in line.split('\t'):
                        html.append('<th>' + element + '</th>')
                    html.append('</tr>')
                    use_table = True
                else:
                    html.append('<tr>')
                    for element in line.split('\t'):
                        html.append('<td>' + element + '</td>')
                    html.append('</tr>')
                continue
            if use_table:
                html.append('</table>')
                use_table = False
            elif use_desc:
                html.append('</dl>')
                use_desc = False
            if line.endswith(':'):
                html.append('<h3>' + line[:-1] + '</h3>')
            else:
                html.append('<p>' + line + '</p>')

    return section_name, "\\n".join(html)


def convert_report(report_file, mqc_file, report_name):
    """
    Adds a line of text from the report.txt file to the html file.
    param: str work_dir = isolate-specific folder, e.g.: 'WH200812_001259/'
    param: str line = line of text from the report.txt file
    """

    section_name, data = make_report_html(report_file)

    report_mqc = {
        "id": report_name,
        "section_name": section_name,
        "plot_type": "html",
        "data": data,
    }

    with open(mqc_file, 'w') as outfile:
        yaml.dump(report_mqc, outfile, default_flow_style=False)


def convert_reports(lo_reports):
    """Convert reports"""

    for report_file in lo_reports:
        report_name = report_file[:-4]
        convert_report(report_file, report_name + '_mqc.yml', report_name)


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
