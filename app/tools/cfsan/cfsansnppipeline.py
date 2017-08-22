import logging

import os

import re

from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.tools.tool import Tool


class CfsanSnpPipeline(Tool):
    """
    Runs the SNP pipeline developed by CFSAN.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: CAMEL instance
        """
        super(CfsanSnpPipeline, self).__init__('CFSAN SNP Pipeline', '0.7.0', camel)

    def _check_input(self):
        """
        Checks the tool input.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise ValueError("No reference FASTA input found.")
        if len(self._tool_inputs['FASTA']) != 1:
            raise ValueError("Only one reference FASTA file is supported.")
        if 'FASTQ' not in self._tool_inputs:
            raise ValueError("No FASTQ input found.")
        if len(self._tool_inputs['FASTQ']) % 2 != 0:
            raise ValueError("Only paired end input is supported (2 entries per sample)")
        if len(self._tool_inputs['VAL_Name']) != len(self._tool_inputs['FASTQ']) / 2:
            raise ValueError("Number of sample names does not correspond to the number of read files.")
        super(CfsanSnpPipeline, self)._check_input()

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        reads_folder = self.__create_reads_input()
        self.__build_command(reads_folder)
        self._execute_command()
        self.__set_output()
        self.__analyze_metrics_file()

    def __get_reference(self):
        """
        Returns the reference FASTA file, a symbolic link is created so the pipeline can index it.
        :return: Path to reference FASTA
        """
        link_name = os.path.join(self._folder, self._tool_inputs['FASTA'][0].basename)
        logging.info("Creating symlink for reference FASTA file: {}".format(link_name))
        os.symlink(self._tool_inputs['FASTA'][0].path, link_name)
        return link_name

    def __create_reads_input(self):
        """
        Creates the input for the pipeline.
        :return: Reads input folder
        """
        reads_folder = os.path.join(self._folder, 'reads')
        os.mkdir(reads_folder)
        for i in range(0, len(self._tool_inputs['FASTQ']), 2):
            forward_reads = self._tool_inputs['FASTQ'][i]
            reverse_reads = self._tool_inputs['FASTQ'][i+1]
            sample_name = self.__get_sample_name(i)
            logging.info("Adding sample '{}' as input".format(sample_name))
            sample_folder = os.path.join(reads_folder, sample_name)
            logging.info("Creating sample folder '{}'".format(sample_folder))
            os.mkdir(sample_folder)
            os.symlink(forward_reads.path, os.path.join(sample_folder, '{}_1.fastq'.format(sample_name)))
            os.symlink(reverse_reads.path, os.path.join(sample_folder, '{}_2.fastq'.format(sample_name)))
        return reads_folder

    def __get_sample_name(self, index):
        """
        Returns the name of the sample based on the read names.
        :param index: Index in the input list
        :return: Sample name
        """
        if 'VAL_Name' in self._tool_inputs:
            return self._tool_inputs['VAL_Name'][index / 2].value.replace(' ', '_')

        filename = self._tool_inputs['FASTQ'][index].path
        m = re.match('.*/(.*)_\d\.[fastq]+$', filename)
        if m is None:
            raise ValueError("Cannot determine sample name from: {}".format(filename))
        return m.group(1)

    def __build_command(self, reads_folder):
        """
        Builds the command.
        :param reads_folder: Reads folder
        :return: None
        """
        command_parts = [
            '. $VIRTUALENV;',
            self._tool_command,
            '-s {}'.format(reads_folder),
            self.__get_reference(),
            ' '.join(self._build_options())
        ]
        if 'TXT' in self._tool_inputs:
            command_parts.insert(2, '-c {}'.format(self._tool_inputs['TXT'][0].path))
        self._command.command = ' '.join(command_parts)

    def __set_output(self):
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

    def __analyze_metrics_file(self):
        """
        Analyzes the metrics file and adds the information to the informs.
        :return: None
        """
        saved_columns = [2] + list(range(5, 14))
        metrics_file = os.path.join(self._folder, 'metrics.tsv')
        if not os.path.isfile(metrics_file):
            raise IOError("No metrics file generated.")
        with open(metrics_file) as handle:
            lines = handle.readlines()
        column_names = lines[0].split('\t')
        for line in lines[1:]:
            print(line)
            line_parts = line.split('\t')
            sample_name = line_parts[0].replace('"', '')
            stats = {column_names[i]: line_parts[i] for i in range(0, len(line_parts)) if i in saved_columns}
            self._informs[sample_name] = stats

    def _check_command_output(self):
        """
        Checks the command output.
        :return: None
        """
        pass
