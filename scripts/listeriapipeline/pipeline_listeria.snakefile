import sys
sys.path.append('/home/qiafu/Work/camel3/')


import datetime
import os
import shutil
from app.camel import Camel
from app.components.filesystemhelper import FileSystemHelper
from app.components.html.htmlreport import HtmlReport
from app.components.html.htmlreportsection import HtmlReportSection
from app.io.tooliodirectory import ToolIODirectory
from app.io.tooliofile import ToolIOFile
from app.io.tooliovalue import ToolIOValue
from app.snakemake.snakemakeutils import SnakemakeUtils
from app.pipeline.snakestep import SnakeStep


# 1. Get the sub workflows
from resources import CSS_STYLE
from scripts.listeriapipeline.subworkflows import WORKFLOW_INIT_REPORT, WORKFLOW_ASSEMBLY, WORKFLOW_READ_TRIMMING, \
    WORKFLOW_GENE_DETECTION, WORKFLOW_SEQUENCE_TYPING, WORKFLOW_CONTAMINATION_CHECK, WORKFLOW_QUALITY_CHECKS

# 2. Get the working dir from the config
__WORKING_DIR = config.get('working_dir', '.')

# 3. Add the workflows
include:
    WORKFLOW_INIT_REPORT
include:
    WORKFLOW_READ_TRIMMING
include:
    WORKFLOW_ASSEMBLY
include:
    WORKFLOW_CONTAMINATION_CHECK
include:
    WORKFLOW_QUALITY_CHECKS
include:
    WORKFLOW_GENE_DETECTION
include:
    WORKFLOW_SEQUENCE_TYPING

# 4. Set the pipeline version
__PIPELINE_VERSION = "0.1"

# 5. Create a CAMEL instance to run tools
camel = Camel()

# 6. Create the directory to store the output files
if not os.path.isdir(config['output_dir']):
    os.makedirs(config['output_dir'])

# 7. Rule all ensures that the final output is generated.
rule combine_reports:
    input:
        config.get('report'),
        os.path.join(__WORKING_DIR, 'report_read_trimming', 'html.io'),
#        os.path.join(__WORKING_DIR, 'contamination_check', 'html.io'),
        os.path.join(__WORKING_DIR, 'report_assembly', 'html.io'),
#        os.path.join(__WORKING_DIR, 'report_quality_checks', 'html.io'),
#        os.path.join(__WORKING_DIR, 'report_gene_detection', 'html.io') if 'gene_detection' in config else [],
#        os.path.join(__WORKING_DIR, 'report_sequence_typing', 'html.io') if 'sequence_typing' in config else [],
    params:
        output_dir = config['output_dir'],
        version = __PIPELINE_VERSION,
        detection_method = config['detection_method']
    run:
        report = HtmlReport(input[0], params.output_dir)
        report.initialize("Listeria pipeline", CSS_STYLE)
        report.add_pipeline_header("Listeria pipeline")
        section_input = HtmlReportSection('Input')
        section_input.add_table([
            ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
            ['Input files:', ', '.join(['tdo', 'tdo'])],
            ['Pipeline version:', params.version],
            ['Detection method:', params.detection_method]],
            table_attributes=[('class', 'information')])
        report.add_html_object(section_input)

        for pickle in input[1:]:
            section = SnakemakeUtils.load_object(pickle)[0].value
            report.add_html_object(section)
        report.save()


onsuccess:
    print("ONSUCCESS: pipeline run finished.")

onerror:
    print("ONERROR: pipeline fails to finish. log file {!r}".format(log))
    # shutil.copy(log, os.path.join(config.get('output_dir'), 'snake_error.log'))
    # shutil.copy(os.path.join(os.getcwd(), 'camel.log'), os.path.join(config.get('output_dir'), 'snake_camel.log'))
