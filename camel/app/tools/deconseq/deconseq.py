import os
import re
import shutil
from typing import List, Dict

from camel.app.camel import Camel
from camel.app.components.files.fastautils import FastaUtils
from camel.app.components.files.fastqutils import FastqUtils
from camel.app.components.files.fileutils import FileUtils
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.loggers import logger
from camel.app.tools.tool import Tool


class Deconseq(Tool):
    """
    The DeconSeq tool can be used to automatically detect and efficiently remove sequence contaminations from genomic
    and metagenomic datasets. It is easily configurable and provides a user-friendly interface.
    """
    KEEP_INTERMEDIATE_FILES = False

    def __init__(self, camel: Camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('deconseq', '0.4.3', camel)
        self._input_key = None
        self._current_input = None
        self._db = None
        self._db_outputs = {}
        self._commands = []
        self._filehelper = None
        self._cont_file = None

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid:
        - FASTQ or FASTA key is required
        - Only one input file allowed
        - No other input keys are allowed
        :return: None
        """
        for key, value in self._tool_inputs.items():
            if key not in ['FASTQ', 'FASTA']:
                raise InvalidInputSpecificationError('Illegal input key given for DeconSeq, '
                                                     'only FASTQ or FASTA allowed: {!r}'.format(self._tool_inputs))
            if len(value) != 1:
                raise InvalidInputSpecificationError('Illegal number of input files (max = 1) '
                                                     'provided for DeconSeq: {!r}'.format(self._tool_inputs))
        if len(self._tool_inputs.keys()) != 1:
            raise InvalidInputSpecificationError('Too many input keys given for DeconSeq, '
                                                 'only FASTQ or FASTA allowed: {!r}'.format(self._tool_inputs))

    def _execute_tool(self) -> None:
        """
        Runs Deconseq. If the 'sequential' parameter is set, every contamination database will be run
        seperately and informs saved per run. Else, all databases are checked at the same time.
        :return: None
        """
        self.__set_input_key()
        self.__set_filehelper()
        self._cont_file = self.__compose_cont_file_name()
        self.__prepare_informs()
        contamination_dbs = [self._parameters['dbs'].value] if not self._parameters['sequential'].value else self._parameters['dbs'].value.split(',')
        for db in contamination_dbs:
            self.__set_input()
            self._db = db
            self.__build_command()
            self._execute_command()
            self._db_outputs[self._db] = self.__get_output_specs()
            self.__check_or_create_output()
            self.__update_informs()
            if len(contamination_dbs) > 1:
                self.__append_cont_file()
            if self.__results_emtpy():
                self.__update_inform_of_missing_dbs()
                break
        self._tool_outputs.update(self.__get_output_specs())
        self.__move_cont_file()
        self.__update_final_informs()
        self.__cleanup_intermediate_files()

    def __prepare_informs(self) -> None:
        """
        Adds keys to the informs that will be used during the processing of each contamination
        database
        :return: None
        """
        self._informs['deconseq_stats'] = {}
        self._informs['processed_dbs'] = []
        self._informs['removed_reads_counts'] = []
        self._informs['commands'] = []
        self._informs['initial_reads_count'] = self._filehelper.count_reads(self._tool_inputs[self._input_key][0].path)

    def __append_cont_file(self) -> None:
        """
        Appends the reads from the current 'cont' file to a temporary intermediate as this file
        will be overwritten in the next run in case the databases are run sequentially.
        :return: None
        """
        with open(self._db_outputs[self._db][f'{self._input_key}_Cont'][0].path, 'r') as inhandle, open(self._cont_file, 'a') as outhandle:
            for line in inhandle:
                outhandle.write(line)

    def __set_input_key(self) -> None:
        """
        Sets the instance variable self._input_key to either FASTA or FASTQ
        :return: None
        """
        self._input_key = list(self._tool_inputs.keys())[0]

    def __set_filehelper(self) -> None:
        """
        Set the proper file handler based on self._input_key (FASTA or FASTQ)
        :return: None
        """
        if self._input_key == 'FASTA':
            self._filehelper = FastaUtils
        else:
            self._filehelper = FastqUtils

    def __set_input(self) -> None:
        """
        Sets the input for the current run based on whether deconseq is run sequentially
        :return: None
        """
        if self._parameters['sequential']:
            self._current_input = self.__get_sequential_input()
        else:
            self._current_input = self._tool_inputs[self._input_key][0]

    def __get_sequential_input(self) -> str:
        """
        Returns the file name to use for the current run when deconseq is run sequentially.
        When this is the first run, the provided input file is returned. Else, when a retain database
        is given, the 'clean' and 'both' files from the previous run are concatenated as they are
        both valid inputs. When no retain database is given and this is not the first run, the clean
        reads from the previous run are copied to a new location as that file will be overwritten by
        the current run.
        :return: File name of the input to use for this run
        """
        if self._db is None:  # First run of Deconseq sequential
            return self._tool_inputs[self._input_key][0]
        retain_file = self.__compose_retain_file_name()
        if self.__with_retain_db():
            FileUtils.concatenate_files(retain_file, [self._db_outputs[self._db][self._input_key + '_Clean'][0].path,
                                                      self._db_outputs[self._db][self._input_key + '_Both'][0].path])
        else:
            shutil.copyfile(self._db_outputs[self._db][self._input_key + '_Clean'][0].path, retain_file)
        return retain_file

    def __compose_retain_file_name(self) -> str:
        """
        Compose 'retain' file name to hold read retained (cleaned and/or both reads)
        :return: String as retain filename
        """
        basename = self.__get_basename()
        extension = self.__get_extension()
        return basename + "." + self._db + '_retain' + extension

    def __compose_cont_file_name(self) -> str:
        """
        Compose 'contamination' file name to hold reads that were removed by Deconseq
        :return: Contamination filename
        """
        basename = self.__get_basename()
        extension = self.__get_extension()
        return basename + "." + '_contaminations' + extension

    def __get_basename(self) -> str:
        """
        Returns the prefix that will be used in the output.
        :return: String with the prefix used in the output
        """
        infile = os.path.basename(self._tool_inputs[self._input_key][0].path)
        return os.path.join(self._folder, os.path.splitext(infile)[0]) + '.deconseq'

    def __get_extension(self) -> str:
        """
        Returns the extension that the output file will have
        :return: Extension of the output file
        """
        if self._input_key == 'FASTQ':
            return '.fq'
        else:
            return '.fa'

    def __with_retain_db(self) -> bool:
        """
        Returns whether a retain database is given for this instance
        :return: True/False
        """
        return 'dbs_retain' in self._parameters

    def __build_command(self) -> None:
        """
        Concatenates required parameters and options to build the command
        :return: None
        """
        input_string = self.__build_input_string()
        options_string = ' '.join(self._build_options())
        self._command.command = '{} {} {}'.format(self._tool_command, input_string, options_string)

    def __build_input_string(self) -> str:
        """
        Creates the string with the input files and output directories
        :return: String with the input parameters
        """
        command_parts = ['-f {}'.format(self._current_input),
                         '-id {}'.format(os.path.basename(self.__get_basename())),
                         '-out_dir {}'.format(self._folder)]
        return ' '.join(command_parts)

    def _build_options(self, excluded_parameters: List[str] = None, delimiter: str = ' ') -> List[str]:
        """
        Overload the default func '_build_options' to handle 'dbs' specially in order to run on each contamination db separately
        if needed.
        :return: list of formatted parameters
        """
        options_list = super(Deconseq, self)._build_options(excluded_parameters=['dbs', 'sequential'])
        options_list.append(f'-dbs {self._db}')
        return options_list

    def __get_output_specs(self) -> Dict[str, List[ToolIOFile]]:
        """
        Sets the output specifications, i.e. the output key and location of the current output
        :return: output specs (key and files)
        """
        basename = self.__get_basename()
        extension = self.__get_extension()
        outputs = {
            self._input_key + '_Clean': [ToolIOFile(basename + '_clean' + extension)],
            self._input_key + '_Cont': [ToolIOFile(basename + '_cont' + extension)]
        }
        if 'dbs_retain' in self._parameters:
            outputs[self._input_key + '_Both'] = [ToolIOFile(basename + '_both' + extension)]
        return outputs

    def __set_output(self) -> None:
        """
        Sets the name of the output files
        :return: None
        """
        basename = self.__get_basename()
        extension = self.__get_extension()
        self._tool_outputs[self._input_key + '_Clean'] = [ToolIOFile(basename + '_clean' + extension)]
        self._tool_outputs[self._input_key + '_Cont'] = [ToolIOFile(basename + '_cont' + extension)]
        if 'dbs_retain' in self._parameters:
            self._tool_outputs[self._input_key + '_Both'] = [ToolIOFile(basename + '_both' + extension)]

    def __check_or_create_output(self) -> None:
        """
        Check each output file, when missing, generate an empty file. This is needed
        as other functions might be looking for this file
        :return: None
        """
        for key, value in self._db_outputs[self._db].items():
            if not os.path.isfile(value[0].path):
                open(value[0].path, 'w').close()

    def __update_informs(self) -> None:
        """
        Update self._informs with statistics after each run of deconseq
        :return: None
        """
        removed_reads = self._filehelper.count_reads(self._db_outputs[self._db][self._input_key + '_Cont'][0].path)
        self._informs['deconseq_stats'].update({
            self._db: {
                'input_reads_count': self._filehelper.count_reads(str(self._current_input)),
                'removed_reads_count': removed_reads
            }})
        self._informs['processed_dbs'].append(self._db)
        self._informs['removed_reads_counts'].append(removed_reads)
        self._informs['commands'].append(self._informs['_command'])

    def __results_emtpy(self) -> bool:
        """
        Check whether there are still reads left after decontamination on a db
        :return: False if there are reads left in either 'Clean' or 'Both' output, True otherwise.
        """
        output_files = self._db_outputs[self._db]
        for outf_type in output_files.keys():
            res = re.search("_Clean$", outf_type)
            if res and os.path.getsize(output_files[outf_type][0].path) != 0:
                return False
            res = re.search("_Both$", outf_type)
            if res and os.path.getsize(output_files[outf_type][0].path) != 0:
                return False
        return True

    def __update_inform_of_missing_dbs(self) -> None:
        """
        Update informs for dbs not yet handled in case that all reads are removed by a decontamination db
        :return: None
        """
        if self._parameters['sequential']:
            for db in self._parameters['dbs'].value.split(','):
                if db not in self._informs['deconseq_stats']:
                    self._informs['deconseq_stats'].update({
                        db: {'input_reads_count': 0, 'removed_reads_count': 0}
                    })
                    self._informs['processed_dbs'].append(db)
                    self._informs['removed_reads_counts'].append(0)

    def __cleanup_intermediate_files(self) -> None:
        """
        Intermediate file clean up
        :return: None
        """
        if self.KEEP_INTERMEDIATE_FILES:
            outputs = [os.path.basename(f.path) for x in self._db_outputs[self._db].values() for f in x]
            for f in os.listdir(self._folder):
                if f.endswith('.fq') and f not in outputs:
                    logger.debug("Intermedaite file {} cleaned.".format(f))
                    FileUtils.silent_remove(os.path.join(self._folder, f))

    def __update_final_informs(self) -> None:
        """
        Update final results informations
        :return: None
        """
        reads_count = self._filehelper.count_reads(
            self._db_outputs[self._db][self._input_key + '_Clean'][0].path)
        if self.__with_retain_db():
            reads_count += self._filehelper.count_reads(
                self._db_outputs[self._db][self._input_key + '_Both'][0].path)
        self._informs['final_reads_count'] = reads_count
        self._informs['_command'] = '\n'.join(self._informs['commands'])

    def __move_cont_file(self) -> None:
        """
        Move the temporary file with all contaminated sequences to the expected output path
        :return: None
        """
        if os.path.isfile(self._cont_file):
            shutil.move(self._cont_file, self._tool_outputs[f'{self._input_key}_Cont'][0].path)

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError(f"Deconseq execution failed (Exit code: {self._command.returncode})")
        if 'does not exist in config file' in self._command.stderr:
            raise ToolExecutionError(f"Deconseq execution failed, unknown database given: {self._db}")
