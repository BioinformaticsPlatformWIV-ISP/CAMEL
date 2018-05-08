from camel.app.components.html.htmlelement import HtmlElement


class TsvExporter(object):
    """
    Exports data in tab separated values (TSV) format.
    """

    @staticmethod
    def export(output_data, header, filename, drop_columns=()):
        """
        Exports output data in TSV format.
        :param output_data: Output data
        :param filename: Filename of the output file
        :param header: Header
        :param drop_columns: Columns not included in the exported file
        :return: None
        """
        with open(filename, 'w') as output_file:
            if header:
                full_header = []
                for i in range(0, len(header)):
                    if i not in drop_columns:
                        full_header.append(header[i].strip())
                output_file.write('\t'.join(full_header))
                output_file.write('\n')
            for row in output_data:
                row_data = []
                for i in range(0, len(row)):
                    if i in drop_columns:
                        continue
                    if isinstance(row[i], str):
                        row_data.append(row[i])
                    if isinstance(row[i], HtmlElement):
                        row_data.append(row[i].text)
                output_file.write('\t'.join(row_data))
                output_file.write('\n')
