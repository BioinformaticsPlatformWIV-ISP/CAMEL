from app.command.command import Command


class SVGConvert(object):
    """
    Converts SVG images to other formats.
    """

    @staticmethod
    def convert_svg(svg_file, output, size=None):
        """
        Converts a SVG image. The output format is determined based on the extension of the output file.
        :param svg_file: SVG input file
        :param output: Output file
        :param size: Image size
        :return: None
        """
        convert_command = Command('convert {} {}'.format(svg_file, output))
        if size:
            convert_command.command += ' -resize {}'.format(size)
        convert_command.run_command('.')
