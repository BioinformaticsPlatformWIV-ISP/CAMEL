import sys
sys.path.append('/home/qiafu/Work/camel3/')


import datetime
import os
import shutil
from app.camel import Camel
from app.components.filesystemhelper import FileSystemHelper
from app.components.html.htmlreport import HtmlReport
from app.components.html.htmlelement import HtmlElement
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
__SUMMARY_DIR = os.path.join(__WORKING_DIR, 'summary_info')

# 3. Add the workflows
include: WORKFLOW_INIT_REPORT
include: WORKFLOW_READ_TRIMMING
include: WORKFLOW_ASSEMBLY
include: WORKFLOW_CONTAMINATION_CHECK
include: WORKFLOW_QUALITY_CHECKS
include: WORKFLOW_GENE_DETECTION
include: WORKFLOW_SEQUENCE_TYPING

# 4. Set the pipeline version
__PIPELINE_VERSION = "0.2"

# 5. Create a CAMEL instance to run tools
camel = Camel()

# 6. Create the directory to store the output files
if not os.path.isdir(config['output_dir']):
    os.makedirs(config['output_dir'])
if not os.path.isdir(__SUMMARY_DIR):
    os.makedirs(__SUMMARY_DIR)

# 7. Final rule to combine sub-reports for the final output
rule combine_reports:
    input:
        READ_TRIMMING_REPORT,
        CONTAMINATION_REPORT,
        ASSEMBLY_REPORT,
        QUALITY_CHECKS_REPORT,
        GENE_DETECTION_REPORT if len(config['gene_detection']) != 0 else [],
        TYPING_REPORT if len(config['sequence_typing']) != 0 else [],
    params:
        report_dir = config.get('report'),
        output_dir = config['output_dir'],
        version = __PIPELINE_VERSION,
        detection_method = config['detection_method']
    run:
        report = HtmlReport(params.report_dir, params.output_dir)
        report.initialize("Listeria pipeline", CSS_STYLE)
        report.add_pipeline_header("Listeria pipeline")
        section_input = HtmlReportSection('Input')
        section_input.add_table([
            ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
            ['Input files:', ', '.join([os.path.basename(f) for f in config['fastq_pe']])],
            ['Pipeline version:', params.version],
            ['Detection method:', params.detection_method]],
            table_attributes=[('class', 'information')])
        report.add_html_object(section_input)

        for pickle in input:
            section = SnakemakeUtils.load_object(pickle)[0].value
            report.add_html_object(section)
        report.save()


onsuccess:
    print("ONSUCCESS: pipeline run finished.")

onerror:
    print("ONERROR: pipeline fails to finish. log file {!r}".format(log))
    GALAXY_ERROR_LOG_DIR = '/scratch/qiafu/listeria_pipeline/Galaxy_runs'
    shutil.copy(log, os.path.join(GALAXY_ERROR_LOG_DIR, os.path.basename(log)))
    shutil.copy(os.path.join(os.getcwd(), 'camel.log'), os.path.join(GALAXY_ERROR_LOG_DIR, os.path.basename(log)+'_camel.log'))
