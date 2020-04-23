#!/usr/bin/env python
import argparse
from pathlib import Path

import shutil

from typing import Optional, Sequence

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.genedetection.dbhelper import DBHelper
from camel.app.components.genedetection.genedetectionutils import GeneDetectionUtils
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources import CSS_STYLE


class MainMakeGeneDetectionDB(object):
    """
    This class is used to create databases for the gene detection tool.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes this tool.
        :param args: (Optional) arguments
        """
        self._args = MainMakeGeneDetectionDB.parse_arguments(args)
        fasta_name = self._args.fasta_name if self._args.fasta_name is not None else Path(self._args.fasta).name
        self._db_name = FileSystemHelper.make_valid(Path(fasta_name).stem)
        self._helper = DBHelper(self._db_name, self._args.working_dir)
        self._clusters = None
        self._new_name_by_header = None

    @staticmethod
    def parse_arguments(args: Optional[Sequence[str]] = None) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :return: Parsed arguments
        """
        argument_parser = argparse.ArgumentParser()
        argument_parser.add_argument('--fasta', required=True, help="Input FASTA file.")
        argument_parser.add_argument('--fasta-name', help="Name of the input FASTA file (for Galaxy input).")
        argument_parser.add_argument('--identity-cutoff', default=80, type=int)
        argument_parser.add_argument('--output-dir', required=True)
        argument_parser.add_argument('--output-html', required=True)
        argument_parser.add_argument('--working-dir', default=str(Path('.').absolute()))
        return argument_parser.parse_args(args)

    def run(self) -> None:
        """
        Runs this tool.
        :return: None
        """
        output_dir = Path(self._args.output_dir)
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        input_fasta = self._helper.standardize_fasta_headers(Path(self._args.fasta))
        self.__export_blast_db(input_fasta, output_dir)
        self.__export_srst2_db(input_fasta, output_dir)
        self._helper.export_metadata(self._db_name, output_dir)
        self.__export_report()

    def __export_blast_db(self, input_fasta: Path, output_dir: Path) -> None:
        """
        Creates and exports a gene detection BLAST database from the given FASTA file.
        :param input_fasta: Input FASTA file
        :param output_dir: Output directory
        :return: None
        """
        # Create file
        dir_indexing = self._helper.get_working_subdir('index_blast')
        new_path = dir_indexing / input_fasta.name
        shutil.copyfile(str(input_fasta), str(new_path))

        # Cluster
        self._helper.index_samtools_faidx(new_path, dir_indexing)
        self._helper.index_blast(new_path, dir_indexing)

        # Export files
        for f in dir_indexing.iterdir():
            shutil.copyfile(str(f), str(output_dir / f.name))

    def __export_srst2_db(self, input_fasta: Path, output_dir: Path) -> None:
        """
        Exports a database for SRST2.
        :param input_fasta: Input FASTA file
        :param output_dir: Output directory
        :return: None
        """
        # Cluster FASTA
        dir_clustering = self._helper.get_working_subdir('clustering')
        fasta_seq_headers = dir_clustering / 'seq_headers.fasta'
        self._new_name_by_header = self._helper.convert_fasta_headers_to_seq(input_fasta, fasta_seq_headers)
        self._clusters = self._helper.get_clusters_form_fasta(fasta_seq_headers, self._args.identity_cutoff)

        # Create SRST2 FASTA
        dir_indexing = self._helper.get_working_subdir('index_srst2')
        fasta_srst2 = dir_indexing / f'{self._db_name}-clustered_{self._args.identity_cutoff}.fasta'
        self._helper.create_srst2_fasta(fasta_seq_headers, fasta_srst2, self._clusters)

        # Index
        self._helper.index_samtools_faidx(fasta_srst2, dir_indexing)
        self._helper.index_bowtie2(fasta_srst2, dir_indexing)
        self._helper.index_blast(fasta_srst2, dir_indexing)

        # Export files
        self._helper.export_mapping(self._new_name_by_header, output_dir)
        for f in dir_indexing.iterdir():
            shutil.copyfile(str(f), str(output_dir / f.name))

    def __export_report(self) -> None:
        """
        Creates a report with some info on the database.
        :return: None
        """
        self._report = HtmlReport(self._args.output_html)
        self._report.initialize('Gene detection database', CSS_STYLE)
        self._report.add_html_object(self.__create_db_info_section())
        self._report.add_html_object(self.__create_clusters_section())
        self._report.add_html_object(SnakePipelineUtils.create_commands_section(
            self._helper.informs, self._args.working_dir))
        self._report.save()

    def __create_db_info_section(self) -> HtmlReportSection:
        """
        Creates the report section with the database info.
        :return: HTML report section
        """
        section_db_info = HtmlReportSection('Database info')
        section_db_info.add_table([
            ['Name:', self._db_name],
            ['Size:', sum(len(c.seq_ids) for c in self._clusters)],
            ['Nb. clusters:', len(self._clusters)],
            ['Clustering cutoff: ', f'{self._args.identity_cutoff}%']
        ], table_attributes=[('class', 'information')])
        return section_db_info

    def __create_clusters_section(self) -> HtmlReportSection:
        """
        Creates the report section with the cluster information.
        :return: HTML report section
        """
        section_clusters = HtmlReportSection('Clusters')
        allele_by_seq_id = {s: GeneDetectionUtils.parse_header(h)[1]['allele'] for s, h in
                            self._new_name_by_header.items()}
        table_data = [[
            cluster.name,
            len(cluster.seq_ids),
            ', '.join([allele_by_seq_id[s] for s in cluster.seq_ids])
        ] for cluster in self._clusters]
        section_clusters.add_table(table_data, ['Cluster', 'Size', 'Sequence ids'], [('class', 'data')])
        return section_clusters


if __name__ == '__main__':
    Camel.get_instance()
    main = MainMakeGeneDetectionDB()
    main.run()
