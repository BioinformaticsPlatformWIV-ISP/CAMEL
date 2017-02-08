import os
import re

from app.error.toolexecutionerror import ToolExecutionError
from app.io.tooliofile import ToolIOFile
from app.tools.velvet.velvetg_output_analyzer import VelvetgOutputAnalyzer
from app.tools.velvet.velvet import Velvet


class VelvetOptimiser(Velvet):
    """
    VelvetOptimiser is a multi-threaded Perl script for automatically optimising the three primary parameter
    options (K, -exp_cov, -cov_cutoff) for the Velvet de novo sequence assembler.
    """
    OUTPUT_FASTA = 'contigs.fa'
    LOG_FILE = 'Log'
    STATS_FILE = 'stats.txt'
    # velvetoptimiser's own logfile, named as a combination of cmd parameter 'log_prefix' and a suffix
    VO_LOG_FILE_SURFIX = '_logfile.txt'

    def __init__(self, camel):
        """
        Initialize VelvetOptimiser
        :param camel: Camel instance
        :return: None
        """
        super(VelvetOptimiser, self).__init__('VelvetOptimiser', '2.2.5', camel)
        self._output_dir = None

    def _execute_tool(self):
        """
        Function to run BWA index
        :return: None
        """
        self._set_input()
        self.__build_command()
        self._execute_command()
        self.__set_output()
        self.__set_inform()

    def __set_output(self):
        """
        Set the outputs for VelvetOptimiser
        :return: None
        """
        self._output_dir = os.path.join(self._folder, self._parameters['output_dir'].value)
        self._tool_outputs['FASTA_Contig'] = [ToolIOFile(os.path.join(self._output_dir, VelvetOptimiser.OUTPUT_FASTA))]

    def __build_command(self):
        """
        Build command to run VelvetOptimiser
        :return: None
        """
        if 'velvethopts' in self._parameters:
            velveth_opts = " ".join([self._input_string, self._parameters['velvethopts'].value])
        else:
            velveth_opts = " ".join([self._input_string])
        self._command.command = "{} {} -f {!r} -o {!r}".format(
            self._tool_command,
            " ".join(self._build_options(excluded_parameters=['velvetgopts', 'velvethopts'])),
            velveth_opts,
            self._parameters['velvetgopts'].value
        )

    @staticmethod
    def __analyze_velvetoptimiser_logfile(log_file):
        """
        Analyze the logfile of VelvetOptimiser to collect assembly information
        :param log_file: velvetg Log output file (with complete path)
        :return: informs, dictionary containing update velvetg run information
        """
        informs = {}
        info_items = (
            'Velvet hash value',
            'Roadmap file size',
            'Assembly score'
        )
        with open(log_file) as f:
            content = f.readlines()
        # skip content till the Final results section
        while True:
            line = content.pop(0)
            if line.strip() == 'Final optimised assembly details:':
                break
        for line in content:
            if re.search(":", line):
                items = line.split(":")
                items[1] = items[1].strip()
                if items[0] in info_items:
                    informs[items[0]] = items[1]
            if re.search("Velvetg parameter string", line):
                res = re.search("-cov_cutoff (\d+\.\d+)", line)
                if res:
                    informs['cov_cutoff'] = res.groups()[0]
        # store hash value in key 'kmer'
        informs['kmer'] = int(informs['Velvet hash value'])
        informs.pop('Velvet hash value')

        return informs

    def _check_command_output(self):
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if re.search('VelvetOptimiser.pl line', self.stderr) is not None:
            raise ToolExecutionError("Command execution failed (stderr: {})".format(self.stderr))

        if re.search('(error|fail)', self.stderr) is not None:
            raise ToolExecutionError("Command execution failed (stderr: {})".format(self.stderr))

    def __set_inform(self):
        """
        Analyze the result of velvetoptimiser to extract statsitics. Three sources are analyzed: analyze velvetg Log
        file for basic assembly statistics; velvetoptimiser Logfile for statistics: kmer, cov_cutoff, and assembly
        score; and stats.txt for general contigs stats
        :return: None
        """
        self.informs['tool_name'] = 'VelvetOptimiser'

        # analyze results utilizing functions in velvetg_output_analyzer module
        self.informs.update(
            VelvetgOutputAnalyzer.analyze_log_file(os.path.join(self._output_dir, VelvetOptimiser.LOG_FILE))
        )

        log_file_path = os.path.join(self._output_dir, self._parameters[
                                     'log_prefix'].value + VelvetOptimiser.VO_LOG_FILE_SURFIX)
        self.informs.update(self.__analyze_velvetoptimiser_logfile(log_file_path))
        # analyze logfile before stats file, as some information extracted from velvet optimiser logfile is needed to
        # analyze the later
        self.informs.update(
            VelvetgOutputAnalyzer.analyze_stats_file(
                self._informs, os.path.join(self._output_dir, VelvetOptimiser.STATS_FILE))
        )
