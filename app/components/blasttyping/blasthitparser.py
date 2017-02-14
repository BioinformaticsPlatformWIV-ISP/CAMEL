import logging

from app.components.blasttyping.blasthit import BlastHit


class BlastHitParser(object):
    """
    Parser a tabular blast output file and returns it as BlastHit objects.
    """

    @staticmethod
    def parse(filename):
        """
        Parses the given Blast output file.
        :param filename: filename
        :return: Parsed output
        """
        blast_hits = []
        logging.info("Parsing blast hits from '{}'".format(filename))
        with open(filename) as input_file:
            for line in input_file.readlines():
                blast_hits.append(BlastHit.convert_blast_output_line(line.strip()))
        logging.info("{} blast hits parsed".format(len(blast_hits)))
        return blast_hits
