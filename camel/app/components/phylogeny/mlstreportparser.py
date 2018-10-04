import json
import logging
from typing import Dict, Tuple, List

import os

from camel.app.components.phylogeny.mlsttabularparser import MlstTabularParser
from camel.app.tools.pipelines.sequence_typing.htmlreportertyping import HtmlReporterTyping


class MlstReportParser(object):
    """
    Extracts the sequence typing information from the HTML report output.
    """

    @staticmethod
    def get_scheme_dir(html_file: str, html_name: str) -> str:
        """
        Returns the directory with the output of the MLST scheme.
        :param html_file: HTML file
        :param html_name: HTML Galaxy name
        :return: Scheme directory.
        """
        files_dir = f"{html_file.split('.')[0]}_files"
        st_dir = os.path.join(files_dir, 'sequence_typing')
        if not os.path.isdir(st_dir):
            raise ValueError("Input file '{}' is not a valid sequence typing output".format(html_name))
        if len(os.listdir(st_dir)) != 1:
            raise ValueError("Multiple sequence typing outputs found")
        return os.path.join(st_dir, os.listdir(st_dir)[0])

    @staticmethod
    def get_analysis_info(scheme_dir: str) -> Dict:
        """
        Returns the analysis metadata from the given scheme directory.
        :param scheme_dir: Folder containing the typing output
        :return: Analysis info dictionary
        """
        info_path = os.path.join(scheme_dir, HtmlReporterTyping.INFO_FILENAME)
        if not os.path.isfile(info_path):
            raise ValueError("No analysis info file found (you might have to rerun the sequence typing tool)")
        with open(info_path) as handle:
            metadata = json.load(handle)
            logging.info('Scheme: {}'.format(metadata['scheme']))
            logging.info('Sample: {}'.format(metadata['sample']))
        return metadata

    @staticmethod
    def parse_typing_output(scheme_dir: str) -> List[Tuple[str, str]]:
        """
        Parses the sequence typing output file.
        :param scheme_dir: Folder containing the typing output
        :return: Parsed allele ids
        """
        try:
            tabular_file = [os.path.join(scheme_dir, x) for x in os.listdir(scheme_dir) if x.endswith('.tsv')][0]
        except IndexError:
            raise FileNotFoundError("No tabular output file found")
        return MlstTabularParser.parse_tabular_input(tabular_file)

    @staticmethod
    def parse_html_all(html_input: List[Tuple[str, str]]) -> Dict[str, List[Tuple[str, str]]]:
        """
        Parses all HTML input files that were provided trough the command line arguments.
        :param html_input: List of HTML input tuples
        :return: Parsed allele ids
        """
        allele_ids = {}
        for html_file, html_name in html_input:
            logging.info('Parsing file: {}'.format(html_name))
            scheme_dir = MlstReportParser.get_scheme_dir(html_file, html_name)
            metadata = MlstReportParser.get_analysis_info(scheme_dir)
            if metadata['sample'] in allele_ids:
                logging.warning("Duplicate sample! {}".format(metadata['sample']))
                continue
            allele_ids[metadata['sample']] = MlstReportParser.parse_typing_output(scheme_dir)
        return allele_ids
