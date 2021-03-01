from pathlib import Path
from typing import Any, Dict, List, Optional

from camel.app.components.html.htmlreport import HtmlReport
from camel.app.components.html.htmlreportsection import HtmlReportSection
from camel.app.snakemake.snakepipelineutils import SnakePipelineUtils
from camel.resources import CSS_STYLE


def save_report_section(section: HtmlReportSection, title: str, output_path: Path,
                        informs_commands: Optional[List[Dict[str, Any]]] = None) -> None:
    """
    Saves a report section output of a workflow to file.
    :param section: Section to save
    :param title: Section title
    :param output_path: Output path
    :param informs_commands: Informs for the commands
    :return: None
    """
    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)
    report = HtmlReport(str(output_path), str(output_path.parent))
    report.initialize(title, CSS_STYLE)
    report.add_pipeline_header(title)
    report.add_html_object(section)
    section.copy_files(report.output_dir)
    section_commands = SnakePipelineUtils.create_commands_section(informs_commands, output_path.parent)
    report.add_html_object(section_commands)
    report.save()
