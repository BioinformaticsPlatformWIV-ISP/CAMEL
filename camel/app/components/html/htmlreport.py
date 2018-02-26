import logging
import os
import shutil

from camel.app.components.html.htmlbase import HtmlBase
from camel.app.components.html.htmlelement import HtmlElement
from camel.resources import LOGO_WIV


class HtmlReport(HtmlBase):
    """
    This class represents an HTML report.
    """

    def __init__(self, filename, output_dir=None):
        """
        Initializes the report.
        :param filename: Filename
        :param output_dir: Output directory
        """
        super(HtmlReport, self).__init__()
        self._filename = filename
        self._output_dir = output_dir

    @property
    def output_dir(self):
        """
        Returns the output directory.
        :return: Output directory
        """
        return self._output_dir

    def save(self):
        """
        Saves the report.
        :return: None
        """
        logging.info("Saving report '{}'".format(os.path.basename(self._filename)))
        if self._filename is None:
            raise ValueError("Report with filename 'None' cannot be saved")
        with open(self._filename, 'w', encoding='utf-8') as handle:
            self.add_raw('<!DOCTYPE HTML>')
            handle.write(self.to_html())

    def initialize(self, title, css_style=None):
        """
        Initializes an HTML report.
        :param title: Report title
        :param css_style: CSS style
        :return: None
        """
        logging.info("Initializing report")
        with self._doc.tag('head'):
            with self._doc.tag('title'):
                self._doc.text(title)
            self._doc.stag('meta', charset='UTF-8')
            if css_style is not None:
                with open(css_style) as css, self._doc.tag('style', type="text/css"):
                    self._doc.asis(css.read())

    def add_pipeline_header(self, pipeline_name):
        """
        Adds the pipeline header to the report.
        :param pipeline_name: Name of the pipeline
        :return: None
        """
        if self._output_dir is None:
            raise ValueError("Can't add the pipeline header without an output directory")
        with self.get_tag('div', [('class', 'header')]):
            with self.get_tag('div', [('id', 'header_title')]):
                self._doc.stag('img', src='logo-wiv-isp.png', alt='WIV-ISP Belgium', id='header_logo')
                self._doc.text("{} Report".format(pipeline_name))
        shutil.copy(LOGO_WIV, self._output_dir)

    def to_html(self):
        """
        Returns the report HTML code.
        :return: HTML code.
        """
        parent = HtmlElement('html')
        parent.add_raw(self._doc.getvalue())
        return parent.to_html()
