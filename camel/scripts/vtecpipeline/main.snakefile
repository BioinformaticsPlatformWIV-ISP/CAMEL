import datetime
import os
import shutil
from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.pipeline.snakestep import SnakeStep



# 1. Get the sub workflows
from camel.resources import CSS_STYLE
from camel.scripts.vtecpipeline.subworkflows import WORKFLOW_ASSEMBLY, WORKFLOW_READ_TRIMMING, \
    WORKFLOW_GENE_DETECTION


# 2. Get the working dir from the config
__WORKING_DIR = config.get('working_dir', '.')
__SUMMARY_DIR = os.path.join(__WORKING_DIR, 'summary_info')

# 3. Add the workflows
include: WORKFLOW_READ_TRIMMING
include: WORKFLOW_ASSEMBLY
# include: WORKFLOW_CONTAMINATION_CHECK
# include: WORKFLOW_QUALITY_CHECKS
include: WORKFLOW_GENE_DETECTION
# include: WORKFLOW_SEQUENCE_TYPING
# include: WORKFLOW_SNP_TYPING

# 4. Set the pipeline version
__PIPELINE_VERSION = "0.1"

# 5. Create a CAMEL instance to run tools
camel = Camel()

# 6. Create the directory to store the output files
if not os.path.isdir(config['output_dir']):
    os.makedirs(config['output_dir'])
if not os.path.isdir(__SUMMARY_DIR):
    os.makedirs(__SUMMARY_DIR)

# 7. Combines all reports.
rule combine_reports:
    input:
        os.path.join(__WORKING_DIR, 'report_read_trimming', 'html.io') if not config.get('skip_trimming', False) else [],
        # os.path.join(__WORKING_DIR, 'contamination_check', 'html.io'),
        os.path.join(__WORKING_DIR, 'report_assembly', 'html.io') if not config.get('skip_assembly', False) else [],
        # os.path.join(__WORKING_DIR, 'report_quality_checks', 'html.io'),
        os.path.join(__WORKING_DIR, 'report_gene_detection', 'html.io') if 'gene_detection' in config else [],
        os.path.join(__WORKING_DIR, 'report_sequence_typing', 'html.io') if 'sequence_typing' in config else [],
        os.path.join(__WORKING_DIR, 'snp_typing', 'html.io') if 'snp_typing' in config else []
    params:
        report_path=config.get('report'),
        output_dir=config['output_dir'],
        version=__PIPELINE_VERSION,
        detection_method=config['detection_method']
    run:
        report = HtmlReport(params.report_path, params.output_dir)
        report.initialize("VTEC pipeline", CSS_STYLE)
        report.add_pipeline_header("VTEC pipeline")
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
