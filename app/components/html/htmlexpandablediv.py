from typing import List, Tuple

from yattag import Doc

from app.components.html.htmlelement import HtmlElement
from app.components.html.htmlreport import HtmlReport
from resources import CSS_STYLE


class HtmlExpandableDiv(HtmlElement):
    """
    This class contains a div that can be shown / hidden by clicking.
    """

    SCRIPT_TEMPLATE = """
    $(document).ready(function(){{
        $('#content-{id_}').hide()
    }});
    $('input:checkbox#toggle-{id_}').change(
        function() {{
            if ($(this).is(':checked')) {{
                $('#content-{id_}').show()
            }} else {{
                $('#content-{id_}').hide()
            }}
        }});        
    """

    def __init__(self, id_: str, label: str, attributes: List[Tuple[str, str]]=None):
        """
        Initializes an expandable div.
        :param id_: Id of the element that is hidden / shown
        :param label: Label that is used
        :param attributes: Attributes
        """
        super().__init__('div', attributes=attributes)
        self._id = id_
        self._label = label

    # noinspection PyArgumentList
    def to_html(self) -> str:
        """
        Converts this element to HTML code. A novel Doc() instance is created in order to nest the content of this
        elements Doc() inside the tag associated with this HtmlExpandableDiv.
        :return: HTML code
        """
        doc, tag, text = Doc().tagtext()
        with tag('div'):
            with tag('p'):
                with tag('b'):
                    text('Show / hide {}'.format(self._label))
            with tag('label', klass='switch'):
                doc.stag('input', type='checkbox', id='toggle-{}'.format(self._id))
                with tag('span', klass='slider round'):
                    text('')
            with tag('div', id='content-{}'.format(self._id)):
                doc.asis(self._doc.getvalue())
            with tag('script', type='text/javascript'):
                    text(HtmlExpandableDiv.SCRIPT_TEMPLATE.format(id_=self._id))
        return doc.getvalue()


if __name__ == '__main__':

    report = HtmlReport('/scratch/bebog/working/divreport.html')
    report.initialize('Test report', CSS_STYLE)
    report.add_html_object(HtmlElement('script', '', attributes=[('type', 'text/javascript'),
                                                                 ('src', 'jquery-3.2.1.min.js')]))
    report.add_paragraph('This is a test file!')
    report.add_html_object(HtmlElement('script', HtmlExpandableDiv.SCRIPT_TEMPLATE, [('type', 'text/javascript')]))
    report.add_html_object(HtmlElement('p', ' Click to show / hide', [('id', 'toggle'), ('class', 'flip')]))
    div_content = HtmlElement('div', attributes=[('id', 'content')])
    div_content.add_paragraph('This is the div content')
    report.add_html_object(div_content)
    t = report.get_tag('test')
    print(type(t))
    report.save()

"""
script =
$(document).ready(function(){
$("#spoligo_toggle").click(function() {
$("#spoligo_spacers").toggle(0, function() {});
});
$("#spoligo_spacers").hide()
});
\"""
self._section.add_html_object(
HtmlElement('script', '', attributes=[('type', 'text/javascript'), ('src', 'jquery-3.2.1.min.js')]))
self._section.add_html_object(HtmlElement('script', script, [('type', 'text/javascript')]))
        
"""