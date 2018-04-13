# This script provides some info on how to run tools in camel. It uses the FastQC tool as an example.
import os

from camel.app.camel import Camel
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.fastqc.fastqc import FastQC

# All tools require a Camel instance for their initialisation.
# The camel instance handles the configuration of the logging and the database.
main_directory = os.path.dirname(os.path.dirname(__file__))
camel = Camel(os.path.join(main_directory, 'config/db.yml'),
              os.path.join(main_directory, 'config/logging.yml'))

# You can directly initialize a tool by calling its constructor.
fastqc = FastQC(camel)

##############
# Parameters #
##############
# The tool will automatically load the default parameters from the database.
# You can view them by calling tool.parameter_overview
print(fastqc.parameter_overview)

# Parameter can be added or values can be changed by calling tool.update parameters
fastqc.update_parameters(threads=10)

# Parameters can be removed by setting their value to False
fastqc.update_parameters(quiet=False)
print(fastqc.parameter_overview)

# If you try to update/add/remove a parameter that does not exists, an InvalidParameterError will be raised.
try:
    fastqc.update_parameters(non_existing='10')
except InvalidParameterError:
    pass

# tool.clear_parameters can be used to clear all the parameters (even the default ones)
fastqc.clear_parameters()
print('Parameters: {}'.format(fastqc.parameter_overview))

# Set the parameters for the rest of the example
fastqc.update_parameters(threads=10, extract=True)

##################
# Input / Output #
##################
# Tool inputs and outputs are passed on as ToolIO objects.
# FastQC required one or more FASTQ files as input. Create a ToolIOFile for a fastq file.
fastq_file = ToolIOFile(os.path.join(os.path.dirname(__file__), 'testdata', 'reads_1.fastq'))

# The ToolIOFile object has several properties that can come in handy when working with IO files
print('Path: {}'.format(fastq_file.path))
print('Basename: {}'.format(fastq_file.basename))
print('Extension: {}'.format(fastq_file.file_extension))
print('File size (bytes): {}'.format(fastq_file.size))

# Tool inputs are provided to tools as a dictionary. The keys correspond to types, If a tool accepts multiple
# inputs of the same type, the key can be extended with a short description (e.g. FASTA_Subject & FASTA_Query). The
# values are lists of ToolIO objects.
# FastQC expects one or more FASTQ files as input.
fastqc.add_input_files({'FASTQ': [fastq_file]})

###########
# Running #
###########
# Tools are ran by calling the tool.run() method. The execution folder can be given as argument. The tool checks if
# the provided inputs and parameters are valid before running.
fastqc.run('/data/temp/')

# The tool outputs can be checked by using the tool.tool_outputs property. They are specified in the same manner as
# the tool inputs.
print("Tool outputs: {}".format(fastqc.tool_outputs))

# The informs can be checked by using the tool.informs property.
print("Informs: {}".format(fastqc.informs))

# The stdout and stderr messages are present in the logs, but they can also be explicitly checked.
print("Stdout: {}".format(fastqc.stdout.strip()))
print("Stderr: {}".format(fastqc.stderr.strip()))

# Afterwards we can run FastQC again on other FASTQ file(s).
fastq_file_2 = ToolIOFile(os.path.join(os.path.dirname(__file__), 'testdata', 'reads_2.fastq'))
fastqc.add_input_files({'FASTQ': [fastq_file_2]})
fastqc.run('/data/temp/')
print("Tool outputs: {}".format(fastqc.tool_outputs))
