import logging

import os

from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.tools.tool import Tool


class CfsanSnpPipeline(Tool):
    """
    Runs the SNP pipeline developed by CFSAN.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super().__init__('CFSAN SNP Pipeline', '2.0.2', camel)

    def _check_input(self) -> None:
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No reference FASTA input found.")
        if len(self._tool_inputs['FASTA']) != 1:
            raise InvalidInputSpecificationError("Only one reference FASTA file is supported.")
        if 'FASTQ' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No FASTQ input found.")
        if 'VAL_name' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No sample name input found ('VAL_name')")
        if len(self._tool_inputs['FASTQ']) % 2 != 0:
            raise InvalidInputSpecificationError("Only paired end input is supported (2 entries per sample)")
        if len(self._tool_inputs['VAL_name']) != len(self._tool_inputs['FASTQ']) / 2:
            raise InvalidInputSpecificationError(
                "Number of sample names does not correspond to the number of read files.")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        reads_folder = self.__create_reads_input()
        self.__build_command(reads_folder)
        self._execute_command()
        self.__set_output()
        self.__analyze_metrics_file()

    def __get_reference(self) -> str:
        """
        Returns the reference FASTA file, a symbolic link is created so the pipeline can index it.
        :return: Path to reference FASTA
        """
        link_path = os.path.join(self._folder, self._tool_inputs['FASTA'][0].basename)
        logging.info("Creating symlink for reference FASTA file: {}".format(link_path))
        if not os.path.exists(link_path):
            os.symlink(self._tool_inputs['FASTA'][0].path, link_path)
        return link_path

    def __create_reads_input(self) -> str:
        """
        Creates the input for the pipeline.
        :return: Reads input folder
        """
        reads_folder = os.path.join(self._folder, 'reads')
        if not os.path.isdir(reads_folder):
            os.mkdir(reads_folder)
        for i in range(0, len(self._tool_inputs['FASTQ']), 2):
            forward_reads = self._tool_inputs['FASTQ'][i]
            reverse_reads = self._tool_inputs['FASTQ'][i + 1]
            sample_name = FileSystemHelper.make_valid(self._tool_inputs['VAL_name'][i // 2].value)
            logging.info("Adding sample '{}' as input".format(sample_name))
            sample_folder = os.path.join(reads_folder, sample_name)
            logging.info("Creating sample folder '{}'".format(sample_folder))
            os.mkdir(sample_folder)
            os.symlink(forward_reads.path, os.path.join(sample_folder, '{}_1.fastq'.format(sample_name)))
            os.symlink(reverse_reads.path, os.path.join(sample_folder, '{}_2.fastq'.format(sample_name)))
        return reads_folder

    def __build_command(self, reads_folder: str) -> None:
        """
        Builds the command.
        :param reads_folder: Reads folder
        :return: None
        """
        command_parts = [
            'export CLASSPATH=$GATK_JAR:$PICARD_JAR:$CLASSPATH;',
            self._tool_command,
            '-s {}'.format(reads_folder),
            '--conf {}'.format(os.path.join(os.path.dirname(__file__), 'snppipeline.conf')),
            self.__get_reference(),
            ' '.join(self._build_options())
        ]
        self._command.command = ' '.join(command_parts)

    def __set_output(self) -> None:
        """
        Sets the output.
        :return: None
        """
        self._tool_outputs['FASTA'] = [ToolIOFile(os.path.join(self._folder, 'snpma.fasta'))]
        self._tool_outputs['FASTA_Preserved'] = [ToolIOFile(os.path.join(self._folder, 'snpma_preserved.fasta'))]

        for key in ['VCF_Cons', 'VCF_Cons_preserved', 'VCF', 'VCF_preserved', 'BAM', 'Sample_names']:
            self._tool_outputs[key] = []
        with open(os.path.join(self._folder, 'sampleDirectories.txt')) as handle:
            for line in handle.readlines():
                sample_directory = line.strip()
                self._tool_outputs['VCF_Cons'].append(ToolIOFile(os.path.join(sample_directory, 'consensus.vcf')))
                self._tool_outputs['VCF_Cons_preserved'].append(
                    ToolIOFile(os.path.join(sample_directory, 'consensus_preserved.vcf')))
                self._tool_outputs['VCF'].append(ToolIOFile(os.path.join(sample_directory, 'var.flt.vcf')))
                self._tool_outputs['VCF_preserved'].append(
                    ToolIOFile(os.path.join(sample_directory, 'var.flt_preserved.vcf')))
                self._tool_outputs['BAM'].append(ToolIOFile(os.path.join(sample_directory, 'reads.sorted.bam')))
                self._tool_outputs['Sample_names'].append(ToolIOValue(os.path.basename(line.strip())))

    def __analyze_metrics_file(self) -> None:
        """
        Analyzes the metrics file and adds the information to the informs.
        :return: None
        """
        metrics_file = os.path.join(self._folder, 'metrics.tsv')
        if not os.path.isfile(metrics_file):
            raise IOError("No metrics file generated.")
        with open(metrics_file) as handle:
            lines = handle.readlines()
        column_names = lines[0].split('\t')
        for line in lines[1:]:
            line_parts = line.split('\t')
            sample_name = line_parts[0].replace('"', '')
            stats = {column_names[i]: line_parts[i] for i in range(0, len(line_parts))}
            self._informs[sample_name] = stats

    def _check_command_output(self) -> None:
        """
        Checks the command output.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Error executing {self.name}: {self.stderr}")
