import datetime
import os
import shutil

from yattag import Doc

from camel.app.components.filesystemhelper import FileSystemHelper
from camel.resources import CSS_STYLE


class HtmlReporterSnpPhylogeny(object):
    """
    HTML reporter for the SNP phylogeny pipelines.
    """

    def __init__(self, method, output_directory):
        """
        Initializes the reporter.
        :param method: SNP calling method
        :param output_directory: Output directory
        """
        self._doc, self._tag, self._text = Doc().tagtext()
        self._output_directory = output_directory
        self._method = method

    def get_html(self):
        """
        Returns the HTML code.
        :return: HTML code
        """
        return self._doc.getvalue()

    def initialize(self):
        """
        Initializes the HTML report.
        :return: None
        """
        self._doc.asis('<!DOCTYPE HTML>')
        with self._tag('head'):
            with self._tag('title'):
                self._text('SNP Phylogeny - {}'.format(self._method))
            self._doc.stag('meta', charset='UTF-8')
            with open(CSS_STYLE) as css, self._tag('style', type='text/css'):
                self._doc.asis(css.read())

    def add_analysis_info_section(self, nb_of_samples, reference_name):
        """
        Adds the header section.
        :param nb_of_samples: Number of samples
        :param reference_name: Name of the reference
        :return: None
        """
        with self._tag('h1'):
            self._text('SNP Phylogeny - {}'.format(self._method))
        with self._tag('h2'):
            self._text('Analysis info')
        with self._tag('table', ('class', 'information')):
            with self._tag('tr'):
                with self._tag('td'):
                    self._text('Starting time:')
                with self._tag('td'):
                    self._text(datetime.datetime.now().strftime('%d/%m/%Y - %X'))
            with self._tag('tr'):
                with self._tag('td'):
                    self._text('Nb. of samples:')
                with self._tag('td'):
                    self._text(nb_of_samples)
            with self._tag('tr'):
                with self._tag('td'):
                    self._text('Reference:')
                with self._tag('td'):
                    self._text(reference_name)

    _EXCLUDED_PARAMETERS = ['dir_html', 'html', 'sample', 'reference', 'reference_name']

    def add_parameter_section(self, args):
        """
        Adds the parameter section.
        :param args: Command line arguments
        :return: None
        """
        with self._tag('h2'):
            self._text('Parameters')
        with self._tag('table', klass='data'):
            with self._tag('th'):
                self._text('Name')
            with self._tag('th'):
                self._text('Value')
            for arg_name, arg_value in sorted(vars(args).iteritems()):
                if (arg_name in HtmlReporterSnpPhylogeny._EXCLUDED_PARAMETERS) or (arg_value is None):
                    continue
                with self._tag('tr'):
                    with self._tag('td'):
                        self._text(arg_name)
                    with self._tag('td'):
                        self._text(arg_value)

    def add_trimming_section(self, pipelines):
        """
        Adds the trimming section.
        :return: None
        """
        report_folder = os.path.join(self._output_directory, 'fastqc_reports')
        os.mkdir(report_folder)
        with self._tag('h2'):
            self._text('Read Trimming')
        with self._tag('table', ('class', 'data')):
            self.__add_trimming_table_header()
            for sample_name, pipeline in sorted(pipelines.items()):
                self.__add_trimming_table_row(sample_name, pipeline, report_folder)

    def __add_trimming_table_header(self):
        """
        Adds the header of the trimming table.
        :return: None
        """
        with self._tag('tr'):
            for column_name in ('Sample', 'Nb. of reads', 'Both surviving', 'Forward only surviving',
                                'Reverse Only surviving', 'Dropped', 'Report (Forward)', 'Report (Reverse)'):
                with self._tag('th'):
                    self._text(column_name)

    def __add_trimming_table_row(self, sample_name, pipeline, report_folder):
        """
        Adds a row to the trimming table.
        :param sample_name: Sample name
        :param pipeline: Pipeline instance
        :param report_folder: Folder to store the FastQC reports
        :return: None
        """
        with self._tag('tr'):
            with self._tag('td'):
                self._text(sample_name)
            for inform_key in ('paired_reads_in', 'paired_reads_out', 'forward_only_reads', 'reverse_only_reads',
                               'reads_drop'):
                with self._tag('td'):
                    self._text(pipeline.get_step('Read_trimming').informs[inform_key])
            for fastq_report, read in zip(pipeline.get_step('FastQC_post_trimming').outputs['HTML'], [1, 2]):
                with self._tag('td'):
                    report_name = '{}_{}_fastqc.html'.format(sample_name, read)
                    shutil.copy(fastq_report.path, os.path.join(report_folder, report_name))
                    with self._tag('a', href=os.path.join(os.path.basename(report_folder), report_name)):
                        self._text('view')

    def add_output_files_section(self, sample_names, snp_matrix, snp_matrix_size, output_files):
        """
        Adds the output files section.
        :param sample_names: Sample names
        :param snp_matrix: SNP Matrix
        :param snp_matrix_size: Size of the SNP matrix
        :param output_files: Output files
        :return: None
        """
        with self._tag('h2'):
            self._text("Output Files")

        shutil.copy(snp_matrix.path, os.path.join(self._output_directory, 'snp_matrix.fasta'))
        with self._tag('a', href='snp_matrix.fasta'):
            self._text('SNP Matrix (FASTA)')
        with self._tag('p'):
            self._text('Size: {}'.format(snp_matrix_size))

        with self._tag('table', ('class', 'data')):
            output_files_dir = os.path.join(self._output_directory, 'samples')
            os.mkdir(output_files_dir)
            self.__add_output_files_table_header([x[0] for x in output_files])
            for i in range(0, len(sample_names)):
                self.__add_output_files_table_row(output_files_dir, sample_names[i], output_files, i)

    def __add_output_files_table_header(self, column_names):
        """
        Adds the header to the output files table.
        :param column_names: Column names
        :return: None
        """
        with self._tag('tr'):
            for key in ['Sample'] + column_names:
                with self._tag('th'):
                    self._text(key)

    def __add_output_files_table_row(self, output_files_dir, sample_name, output_files, index):
        """
        Adds a row to the output files table.
        :param output_files_dir: Directory to store the output files
        :param sample_name: Sample name
        :param output_files: Output files
        :param index: index
        :return: None
        """
        sample_dir = os.path.join(output_files_dir, FileSystemHelper.make_valid(sample_name))
        os.mkdir(sample_dir)
        with self._tag('tr'):
            with self._tag('td'):
                self._text(sample_name)
            for name, filename_template, files in output_files:
                filename = filename_template.format(FileSystemHelper.make_valid(sample_name))
                shutil.copy(files[index].path, os.path.join(sample_dir, filename))
                with self._tag('td'):
                    with self._tag('a', href=os.path.join('samples', os.path.basename(sample_dir), filename)):
                        self._text('download')

    def add_metrics_section(self, metrics):
        """
        Adds the metric section.
        :return: None
        """
        with self._tag('h2'):
            self._text("Analysis Metrics")
        with self._tag('table', ('class', 'data')):
            with self._tag('tr'):
                for column_name in metrics[0]:
                    with self._tag('th'):
                        self._text(column_name)
            for row in metrics[1:]:
                with self._tag('tr'):
                    for value in row:
                        with self._tag('td'):
                            self._text(value)

    COLUMN_NAMES = {
        'Average_Insert_Size': 'Insert Size (avg)',
        'Average_Pileup_Depth': 'Pileup Depth (avg)',
        'Percent_of_Reads_Mapped': '% Mapped',
        'Phase1_Preserved_SNPs': 'SNPs pres. (P1)',
        'Phase1_SNPs': 'SNPs (P1)',
        'Phase2_Preserved_SNPs': 'SNPs pres. (P2)',
        'Phase2_SNPs': 'SNPs (P2)'
    }

    def add_analysis_metrics_section(self, cfsan):
        """
        Adds the analysis metrics table.
        :param cfsan: CFSAN instance
        :return: None
        """
        with self._tag('h2'):
            self._text("Analysis Metrics")
        with self._tag('table', ('class', 'data')):
            self.__add_analysis_metrics_header()
            sample_names = [x.value for x in cfsan.tool_outputs['Sample_names']]
            for sample_name in sample_names:
                self.__add_analysis_metrics_row(cfsan, sample_name)

    def __add_analysis_metrics_header(self):
        """
        Adds the header of the analysis metrics table.
        :return: None
        """
        with self._tag('tr'):
            with self._tag('th'):
                self._text('Sample')
            for key in sorted(HtmlReporterSnpPhylogeny.COLUMN_NAMES.keys()):
                with self._tag('th'):
                    self._text(HtmlReporterSnpPhylogeny.COLUMN_NAMES[key])

    def __add_analysis_metrics_row(self, cfsan, sample_name):
        """
        Adds a row to the analysis metrics table.
        :param cfsan: CFSAN instance
        :param sample_name: Sample name
        :return: None
        """
        with self._tag('tr'):
            with self._tag('td'):
                self._text(sample_name)
            for column_name in sorted(HtmlReporterSnpPhylogeny.COLUMN_NAMES.keys()):
                with self._tag('td'):
                    self._text(cfsan.informs[sample_name][column_name])

    def add_model_selection_section(self, model_selection):
        """
        Adds the model selection section.
        :param model_selection: Model selection instance
        :return: None
        """
        with self._tag('h2'):
            self._text("Model Selection")
        with self._tag('p'):
            self._text("Selected model: {}".format(model_selection.informs['model_full']))
            self._doc.stag('br')
            self._text("Rates among sites: {}".format(model_selection.informs['rates_among_sites_full']))
        shutil.copy(model_selection.tool_outputs['CSV'][0].path,
                    os.path.join(self._output_directory, 'model_selection_overview.csv'))
        with self._tag('a', href='model_selection_overview.csv'):
            self._text('Model selection overview (CSV)')
        self._doc.stag('br')

    def add_tree_building_section(self, tree_building, rendered_tree, replications):
        """
        Adds the tree building section.
        :param tree_building: Tree building instance
        :param rendered_tree: Rendered image of the tree
        :param replications: Number of bootstrap replications
        :return: None
        """
        with self._tag('h2'):
            self._text('Tree Building')
        tree_directory = os.path.join(self._output_directory, 'trees')
        os.mkdir(tree_directory)
        shutil.copy(tree_building.tool_outputs['NWK'][0].path, os.path.join(tree_directory, 'phylo.nwk'))
        with self._tag('a', href=os.path.join('trees', 'phylo.nwk')):
            self._text("Phylogenetic tree (Newick)")

        filename = 'tree.png'
        shutil.copy(rendered_tree.path, os.path.join(tree_directory, filename))
        with self._tag('figure'):
            self._doc.stag('img', heigth='512', width='512', src=os.path.join('trees', filename), border='1')
            with self._tag('figcaption'):
                self._text('Phylogentic Tree (With bootstrap support values, {} replications)'.format(replications))

    def add_filtering_section(self, sample_names, trimming_info):
        """
        Adds the SNP filtering section. 
        :param sample_names: Sample names
        :param trimming_info: Trimming information
        :return: None
        """
        with self._tag('h2'):
            self._text('SNP Filtering')
        with self._tag('p'):
            self._text('This table show the number of SNPs that passed each of the filtering steps.')
        with self._tag('table', klass='data'):
            # Header
            with self._tag('tr'):
                for column_name in ['Sample'] + trimming_info.keys():
                    with self._tag('th'):
                        self._text(column_name)
            # Rows
            for i in range(0, len(sample_names)):
                with self._tag('tr'):
                    with self._tag('td'):
                        self._text(sample_names[i])
                    for key in trimming_info.keys():
                        with self._tag('td'):
                            self._text(trimming_info[key][i])

    def add_error_message(self, message):
        """
        Adds an error message.
        :param message: Message
        :return: None
        """
        with self._tag('div', klass='alert-box error'):
            with self._tag('span'):
                self._text('error: ')
            self._text(message)
