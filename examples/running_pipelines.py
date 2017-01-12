# This script provides some info on how to run pipelines in camel. The resistance characterization pipeline will
# be used as an example.
import os

from app.camel import Camel
from app.components.html.htmlhelper import HtmlHelper
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.pipeline.pipeline import Pipeline
from resources import CSS_STYLE, YAML_RES_CHAR_FAST


# Initialize Camel
main_directory = os.path.dirname(os.path.dirname(__file__))
camel = Camel(os.path.join(main_directory, 'config/db.yml'),
              os.path.join(main_directory, 'config/logging.yml'))

# The pipeline checks for resistance genes using Blastn and reports the output in a HTML report.
# Pipeline objects have to be initialized with a YAML file and a CAMEL instance.
# There are two optional arguments:
# - db_pipeline_parameters: If True, the pipeline parameters will be loaded from the database
# - db_logging: If True, the initial pipeline input, the step outputs and the job parameters will be logged in the
#               database.
pipeline = Pipeline(YAML_RES_CHAR_FAST, camel, True, False)

# The pipeline steps can be checked by using the pipeline.steps property
print("Pipeline steps:")
for step in pipeline.steps:
    print(step.name)

# Initialize an empty HTML report
report_path = '/data/temp/bebog/report_reschar.html'
if os.path.isfile(report_path):
    os.remove(report_path)
report = HtmlHelper(report_path)
report.initialize('Resistance Characterization Pipeline', CSS_STYLE)
report.close()
report_file = ToolIOFile(report_path)
report_dir = ToolIODirectory(os.path.dirname(report_path))

# Create the other input files and directories
resfinder_dir = ToolIODirectory('/data/blastdb/nucleotide/ResFinder/latest')
contigs = ToolIOFile(os.path.join(os.path.dirname(__file__), 'testdata', 'contigs.fasta'))

# Add all the inputs, they are specified in the same way as tool inputs.
pipeline.set_initial_input({'DIR': [resfinder_dir],
                            'FASTA': [contigs],
                            'DIR_HTML': [report_dir],
                            'HTML': [report_file]})

# Pipeline parameters can also be changed. Parameters are given as a dictionary with key: step name and value a
# dictionary of pipeline parameters in the format parameter_name: value.
pipeline.add_job_options({'Blastn': {'evalue': '150'}})

# Run the pipeline
pipeline.run('/data/temp/')
