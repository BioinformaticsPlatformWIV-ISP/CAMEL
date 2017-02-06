import logging

import os
import re

from app.components.blasttyping.blasthitfiltering import BlastHitFiltering
from app.components.blasttyping.blasthitparser import BlastHitParser
from app.io.tooliofile import ToolIOFile
from app.tools.tool import Tool


class BestHitExtractor(Tool):
    """
    Tool that extracts the best hit from BLAST output. The BLAST output must be generated using BLASTN or BLASTX with
    the assembled genome or contigs as query and the database sequences as subject/db. The following output format has
    to be used:
    "6 pident sseqid sseq slen qseqid qstart qend".

    The tools also filters hits based on:
    - percent identity (min_percent_identity)
    - percent of the target gene covered by the alignment (min_percent_covered)

    A dictionary with each locus as a key is generated as informs. For each locus the following keys are reported:
    - allele_id: Allele id
    - hit: BlastHit object
    - hit_type: 'perfect', 'no_hits_detected', 'multiple_hits_detected', 'imperfect_identity' or 'imperfect_short'
    - type: 'DNA' or 'peptide'
    - url: PubMLST url of the allele

    A tabular output file is also generated and reported added to the tool_outputs (TSV).
    """

    NO_HIT = 'no_hit'
    MULTIPLE_HITS = 'multi_hit'
    PERFECT_HIT = 'perfect'
    IMPERFECT_IDENTITY_HIT = 'imperfect_identity'
    IMPERFECT_LENGTH_HIT = 'imperfect_short'

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(BestHitExtractor, self).__init__('Sequence Typing: Best Hit Extractor', '0.1', camel)

    def _execute_tool(self):
        """
        Runs this tool.
        :return: None
        """
        detected_alleles = []
        if 'TSV_Blastn' in self._tool_inputs:
            for tsv_file, fasta_file in zip(self._tool_inputs['TSV_Blastn'], self._tool_inputs['FASTA_Nucl']):
                detected_alleles.append(self.__handle_tabular_output(tsv_file, fasta_file))
        if 'TSV_Blastx' in self._tool_inputs:
            for tsv_file, fasta_file in zip(self._tool_inputs['TSV_Blastx'], self._tool_inputs['FASTA_Prot']):
                detected_alleles.append(self.__handle_tabular_output(tsv_file, fasta_file))
        self._tool_outputs['TSV'] = [ToolIOFile(
            self.__create_output_file(detected_alleles, self.__get_output_filename()))]

    def _check_input(self):
        """
        Checks if the required input files were specified.
        :return: None
        """
        if not any(key in self._tool_inputs for key in ('FASTA_Nucl', 'FASTA_Prot')):
            raise ValueError("No input FASTA file found.")
        if not any(key in self._tool_inputs for key in ('TSV_Blastn', 'TSV_Blastx')):
            raise ValueError("No TSV input file found.")
        super(BestHitExtractor, self)._check_input()

    def __handle_tabular_output(self, tsv_file, fasta_file):
        """
        Handles tabular BLAST output.
        :param tsv_file: TSV file
        :param fasta_file: FASTA file
        :return: detected allele, information
        """
        logging.info("Processing hits on '{}'".format(fasta_file.basename))
        hits = BlastHitParser.parse(tsv_file.path)
        hits = BlastHitFiltering.filter_percent_identity(hits, float(self._parameters['min_percent_identity'].value))
        hits = BlastHitFiltering.filter_coverage(hits, float(self._parameters['min_percent_covered'].value))

        if len(hits) == 0:
            detected_allele = BestHitExtractor.NO_HIT
        else:
            detected_allele = self.__get_detected_allele(hits)

        allele_info = self._input_informs['locus_set_info'][fasta_file.path]
        self._informs[allele_info['name']] = self.__get_allele_informs(detected_allele, allele_info)
        return [detected_allele, allele_info]

    @staticmethod
    def __get_detected_allele(hits):
        """
        Return the detected allele from the given set of hits.
        :param hits: Hits
        :return: Best hits
        """
        best_hits = BlastHitFiltering.detect_best_hits(hits)
        if len(best_hits) == 1:
            return best_hits[0]
        elif len(best_hits) > 0:
            return BestHitExtractor.MULTIPLE_HITS
        else:
            return BestHitExtractor.NO_HIT

    def __get_allele_informs(self, detected_allele, allele_info):
        """
        Returns the informs for the given allele.
        :return: Informs
        """
        allele_informs = {}
        if detected_allele == BestHitExtractor.NO_HIT or detected_allele == BestHitExtractor.MULTIPLE_HITS:
            allele_informs['url'] = None
            allele_informs['allele_id'] = '-'
        else:
            detected_allele_id = self.__get_allele_id(detected_allele.database_gene,
                                                      allele_info['allele_id_regex'])
            allele_informs['url'] = self.__get_allele_information_url(allele_info['url'], detected_allele_id)
            allele_informs['allele_id'] = detected_allele_id
        allele_informs['hit_type'] = BestHitExtractor.__get_hit_type(detected_allele)
        allele_informs['hit'] = detected_allele
        allele_informs['type'] = allele_info['type']
        return allele_informs

    def __create_output_file(self, detected_alleles, filename):
        """
        Creates the output file.
        :param detected_alleles: Detected alleles
        :param filename: Output filename
        :return: None
        """
        output_path = os.path.join(self._folder, filename)
        with open(output_path, 'w') as handle:
            handle.write('\t'.join(['Locus', 'Allele', '% Identity', 'HSP length/Locus length', 'Type', 'Hit']))
            handle.write('\n')
            for hit, info in detected_alleles:
                if hit == BestHitExtractor.MULTIPLE_HITS:
                    allele_id = '?'
                    percent_identity = '-'
                    length = '-'
                elif hit == BestHitExtractor.NO_HIT:
                    allele_id = '-'
                    percent_identity = '-'
                    length = '-'
                else:
                    allele_id = self.__get_allele_id(hit.database_gene, info['allele_id_regex'])
                    percent_identity = str(hit.percent_identity)
                    length = '{}/{}'.format(hit.alignment_length, hit.database_gene_length)
                handle.write('\t'.join([info['name'], allele_id, percent_identity, length, info['type'],
                                        BestHitExtractor.__get_hit_type(hit)]))
                handle.write('\n')
        return output_path

    def __get_output_filename(self):
        """
        Returns the name for the output file.
        :return: Name
        """
        if 'LocusSetMetadata' in self._input_informs:
            locus_set_name = self._input_informs['locus_set_info']['name']
            return 'best_hit_detection-{}.tsv'.format(locus_set_name.lower().replace(' ', '_'))
        else:
            return 'best_hit_detection.tsv'

    def __get_allele_id(self, allele_name, regex):
        """
        Returns the allele identifier.
        :param allele_name: Full name of the allele
        :param regex: Regular expression of the allele id
        :return: Allele identifier
        """
        if not regex:
            regex = self._parameters['default_allele_id_regex'].value
        m = re.findall(BestHitExtractor.__cleanup_regex(regex), allele_name)
        if not len(m) == 1:
            raise StandardError("Cannot determine allele identifier for '{}' (RE: {})".format(allele_name, regex))
        return m[0]

    @staticmethod
    def __cleanup_regex(regex):
        """
        Cleans a regular expression, PubMLST regexps have to be adapted to be used.
        :param regex: Regular expression
        :return: Cleaned up regular expression
        """
        characters_to_replace = [('^', ''), ('(', '['), (')', ']')]
        cleaned_regex = regex
        for char, repl in characters_to_replace:
            cleaned_regex = cleaned_regex.replace(char, repl)
        return cleaned_regex

    def __get_allele_information_url(self, locus_url, allele_id):
        """
        Returns the link to the allele information for PubMLST.
        :param locus_url: URL of the locus
        :param allele_id: id of the detected allele
        :return: URL of the allele information
        """
        if allele_id == '-' or locus_url is None:
            return None
        locus_id = locus_url.split('/')[-1]
        database = self._parameters['database'].value
        allele_info_url = 'http://pubmlst.org/perl/bigsdb/bigsdb.pl?db={}&page=alleleInfo&locus={}&allele_id={}'.format(
            database, locus_id, BestHitExtractor.__clean_allele_id(allele_id))
        return allele_info_url

    @staticmethod
    def __clean_allele_id(allele_id):
        """
        Removes mismatch (*) or uncertainty (?) indicators from the allele id.
        :param allele_id: Allele identifier
        :return: Allele id
        """
        return allele_id.replace('*', '').replace('?', '')

    @staticmethod
    def __get_hit_type(hit):
        """
        Returns the type of hit.
        :param hit: BLAST hit
        :return: Type of hit ('Perfect', 'No hit', 'Imperfect short', 'Imperfect Identity')
        """
        if hit == BestHitExtractor.NO_HIT:
            return BestHitExtractor.NO_HIT
        if hit == BestHitExtractor.MULTIPLE_HITS:
            return BestHitExtractor.MULTIPLE_HITS
        elif hit.database_gene_length == hit.alignment_length and hit.percent_identity == 100.0:
            return BestHitExtractor.PERFECT_HIT
        elif hit.database_gene_length == hit.alignment_length:
            return BestHitExtractor.IMPERFECT_IDENTITY_HIT
        else:
            return BestHitExtractor.IMPERFECT_LENGTH_HIT
