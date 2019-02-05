import datetime
import os
import shutil
import logging
import logging.config
import yaml

import sys
sys.path.append('/data/testdir/qiafu/Work/camel3/')

# using snakemake utility func 'srcdir' to get the directory of snakefile
logging_cfg = srcdir("pipeline_logging.yml")
logging.config.dictConfig(yaml.load(open(logging_cfg, 'r')))


from camel.app.camel import Camel
from camel.app.components.filesystemhelper import FileSystemHelper
from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlelement import HtmlElement
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.io.tooliodb import ToolIODb
from camel.app.io.tooliodirectory import ToolIODirectory
from camel.app.io.tooliofile import ToolIOFile
from camel.app.io.tooliovalue import ToolIOValue
from camel.app.snakemake.snakemakeutils import SnakemakeUtils
from camel.app.pipeline.step import Step
from camel.resources import CSS_STYLE, LISTERIA_CITATIONS


# 1. Get the sub workflows
from scripts.listeriapipeline.subworkflows import WORKFLOW_INIT_REPORT, WORKFLOW_ASSEMBLY, WORKFLOW_READ_TRIMMING, \
    WORKFLOW_GENE_DETECTION, WORKFLOW_SEQUENCE_TYPING, WORKFLOW_CONTAMINATION_CHECK, WORKFLOW_QUALITY_CHECKS

# 2. Get the working dir from the config
__WORKING_DIR = config.get('working_dir', '.')
__SUMMARY_DIR = os.path.join(__WORKING_DIR, 'summary_info')
__SPECIES_CONFIRM_ST_DBS = ['species_confirmation', 'MLST-Pasteur']
__OTHER_ST_DBS = [x for x in config.get('sequence_typing', []) if x not in __SPECIES_CONFIRM_ST_DBS]


# 3. Add the workflows
include: WORKFLOW_READ_TRIMMING
include: WORKFLOW_ASSEMBLY
include: WORKFLOW_CONTAMINATION_CHECK
include: WORKFLOW_SEQUENCE_TYPING
include: WORKFLOW_QUALITY_CHECKS
include: WORKFLOW_GENE_DETECTION

# 4. Set the pipeline version
__PIPELINE_VERSION = config.get('pipeline_version')

# 5. Create a CAMEL instance to run tools
camel = Camel(tool_parameter_loc=os.path.join(os.environ['CAMEL_PATH'], 'scripts/listeriapipeline/tool_data/'))

# 6. Create the directory to store the output files
if not os.path.isdir(config['output_dir']):
    os.makedirs(config['output_dir'])
if not os.path.isdir(__SUMMARY_DIR):
    os.makedirs(__SUMMARY_DIR)
REPORT_HTML = config.get('report_html')
REPORT_SUMMARY = config.get('report_summary')


# 7. Final rule to combine sub-reports for the final output
rule all:
    input:
        REPORT_HTML,
        REPORT_SUMMARY

rule combine_summaries:
    input:
        READ_TRIMMING_SUMMARY,
        CONTAMINATION_SUMMARY,
        ASSEMBLY_SUMMARY,
        QUALITY_CHECKS_SUMMARY,
        GENE_DETECTION_DB_SUMMARY.format(db='resfinder') if 'resfinder' in config['gene_detection'] else [],
        GENE_DETECTION_DB_SUMMARY.format(db='card') if 'card' in config['gene_detection'] else [],
        GENE_DETECTION_DB_SUMMARY.format(db='arg_annot') if 'arg_annot' in config['gene_detection'] else [],
        GENE_DETECTION_DB_SUMMARY.format(db='virulencefinder_listeria') if 'virulencefinder_listeria' in config['gene_detection'] else [],
        GENE_DETECTION_DB_SUMMARY.format(db='plasmidfinder_grampositive') if 'plasmidfinder_grampositive' in config['gene_detection'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='species_confirmation') if 'species_confirmation' in config['sequence_typing'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='MLST-Pasteur') if 'MLST-Pasteur' in config['sequence_typing'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='cgMLST') if 'cgMLST' in config['sequence_typing'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='serogroup') if 'serogroup' in config['sequence_typing'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='virulence') if 'virulence' in config['sequence_typing'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='antibiotic_resistance') if 'antibiotic_resistance' in config['sequence_typing'] else [],
        TYPING_SCHEME_SUMMARY.format(scheme='metal_detergent_resistance') if 'metal_detergent_resistance' in config['sequence_typing'] else []
    output:
        REPORT_SUMMARY
    params:
        sample_name = config['sample_name'],
        version = __PIPELINE_VERSION,
        detection_method = config['detection_method']
    run:
        with open(output[0], 'w') as handle:
            handle.write("\n".join([
                'Sample name:\t{}'.format(params.sample_name),
                'Analysis date:\t{}'.format(datetime.datetime.now().strftime('%d/%m/%Y - %X')),
                'Input files:\t{}'.format(', '.join([os.path.basename(f) for f in config['fastq_pe']])),
                'Pipeline version:\t{}'.format(params.version),
                'Detection method:\t{}'.format(params.detection_method),
            ])+"\n")
            for section_summary in input:
                handle.write(open(section_summary, 'r').read())

rule combine_reports:
    input:
        READ_TRIMMING_REPORT,
        CONTAMINATION_REPORT,
        ASSEMBLY_REPORT,
        QUALITY_CHECKS_REPORT,
        SPECIES_CONFIRM_ST_REPORT,
        GENE_DETECTION_REPORT if len(config['gene_detection']) != 0 else [],
        TYPING_REPORT if len(__OTHER_ST_DBS) != 0 else []
    output:
        REPORT_HTML
    params:
        sample_name = config['sample_name'],
        output_dir = config['output_dir'],
        version = __PIPELINE_VERSION,
        detection_method = config['detection_method']
    run:
        report = HtmlReport(output[0], params.output_dir)
        report.initialize("Listeria pipeline", CSS_STYLE)
        report.add_pipeline_header("Listeria pipeline")
        section_input = HtmlReportSection('Input')
        section_input.add_table([
            ['Sample name:', params.sample_name],
            ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
            ['Input files:', ', '.join([os.path.basename(f) for f in config['fastq_pe']])],
            ['Pipeline version:', params.version],
            ['Detection method:', params.detection_method]],
            table_attributes=[('class', 'information')])
        report.add_html_object(section_input)

        for pickle in input:
            section = SnakemakeUtils.load_object(pickle)[0].value
            report.add_html_object(section)

        section = HtmlReportSection('Citations')
        with open(LISTERIA_CITATIONS, 'r', encoding='utf8') as handle:
            section.add_raw(handle.read())
        report.add_html_object(section)
            
        report.save()

onsuccess:
    print("ONSUCCESS: pipeline run finished.")
    section_log = HtmlReportSection('Runnning log')
    log_filepath = 'camel.log'
    section_log.add_file(os.path.join(os.getcwd(), 'camel.log'), log_filepath)
    section_log.add_link_to_file('camel.log', log_filepath)
    section_log.copy_files(config['output_dir'])
    report = HtmlReport(REPORT_HTML, config['output_dir'])
    report.add_html_object(section_log)
    report.close()

onerror:
    print("ONERROR: pipeline fails to finish. log file {!r}".format(log))
    GALAXY_ERROR_LOG_DIR = "/scratch/temp/galaxy_dumps/"
    shutil.copy(log, os.path.join(GALAXY_ERROR_LOG_DIR, os.path.basename(log)))
    shutil.copy(os.path.join(os.getcwd(), 'camel.log'), os.path.join(GALAXY_ERROR_LOG_DIR, os.path.basename(log)+'_camel.log'))
    section_log = HtmlReportSection('Runnning log')
    log_filepath = 'camel.log'
    section_log.add_file(os.path.join(os.getcwd(), 'camel.log'), log_filepath)
    section_log.add_link_to_file('camel.log', log_filepath)
    section_log.copy_files(config['output_dir'])
    report = HtmlReport(REPORT_HTML, config['output_dir'])
    report.add_html_object(section_log)
    report.close()
