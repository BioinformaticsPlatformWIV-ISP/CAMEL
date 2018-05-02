from app.command.command import Command
from app.components.html.htmlreport import HtmlReport
from app.snakemake.snakemakeutils import SnakemakeUtils
from resources import CSS_STYLE


class SnakemakeTestUtils(object):
    """
    This class contains utility functions to test Snakemake workflows.
    """

    @staticmethod
    def run_snakemake(snakefile: str, config_file_path: str, working_dir: str,
                      target_output_file: str, threads: int=8) -> None:
        """
        Runs Snakemake with the given config file.
        :param snakefile: Snakefile to execute
        :param config_file_path: Config file path
        :param target_output_file: Target output file of the workflow
        :param working_dir: Working directory
        :param threads: Number of threads
        :return: None
        """
        command = Command('snakemake --cores {} --snakefile {} --configfile {} {}'.format(
            threads, snakefile, config_file_path, target_output_file))
        command.run_command(working_dir)
        if command.returncode != 0:
            raise ValueError("Snakemake execution failed")

    @staticmethod
    def save_report(report_path: str, output_report_section_pickle: str) -> None:
        """
        Saves the HTML report.
        :param output_report_section_pickle: Path to the pickled output report section
        :param report_path: Report path
        :return: None
        """
        report = HtmlReport(report_path)
        report.initialize('Test report', CSS_STYLE)
        section = SnakemakeUtils.load_object(output_report_section_pickle)[0].value
        report.add_html_object(section)
        report.save()
