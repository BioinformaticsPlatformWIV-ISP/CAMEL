import logging
import os
import random
from typing import Union, Dict

from camel.app.components.blasthit.influenzablastnasnparser import InfluenzaBlastnAsnParser
from camel.app.components.files.fastautils import FastaUtils
from camel.app.components.segmenttyping.segmenttypingreads import SegmentTypingReads
from camel.app.components.segmenttyping.segmenttypingcontigs import SegmentTypingContigs
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.invalidparametererror import InvalidParameterError
from camel.app.tools.tool import Tool
from camel.app.io.tooliofile import ToolIOFile
from camel.app.components.seqid.seqidparser import SeqIDParser


class GenomeTyping(Tool):

    """
    Class that performs genome typing. For Influenza, this can be on several segments.
    """

    def __init__(self, camel):
        """
        Initialize tool
        :param camel: Camel instance
        :return: None
        """
        super().__init__('Genome typing', '0.1', camel)
        self._segment_informs = {}
        self._random_seed = None
        self._failed = False

    def _execute_tool(self) -> None:
        """
        Runs the segment typing if the BLAST run before did not fail (empty output).
        :return: None
        """
        if not self._failed:
            self._set_random_seed()
            self._run_sequence_typing()
            self._quality_check()
            self._set_output()
            if 'influenza_a' in self._parameters:
                self._extract_influenza_a_subtype()

    def _check_parameters(self) -> None:
        """
        Checks whether the given parameters for the tool are correct.
        :return: None
        """
        super(GenomeTyping, self)._check_parameters()
        if self._parameters['multi_segment'].value:
            if 'genome_segments' not in self._parameters:
                raise ValueError(f'Multi segments set to True but no genome segments provided!')

    def _check_input(self) -> None:
        """
        Checks whether the given inputs are valid. If the BLAST input is an empty file, the
        failed variable is set to True
        :return: None
        """
        super(GenomeTyping, self)._check_input()
        if 'ASN' not in self._tool_inputs:
            raise InvalidInputSpecificationError(f'Required input key ASN missing from tool inputs: {self._tool_inputs}')
        if 'DB_BLAST' not in self._tool_inputs:
            raise InvalidInputSpecificationError(f'Required input key BLAST_DB missing from tool inputs: {self._tool_inputs}')
        if os.path.getsize(self._tool_inputs['ASN'][0].path) == 0:
            self._failed = True

    def update_parameters(self, **kwargs: Union[str, int, None, bool, Dict[str, Union[str, int, None, bool]]]) -> None:
        """
        Updates the tool parameters but sets the multi_segment parameter to a boolean instead of string
        :param kwargs: Arguments
        :return: None
        """
        super(GenomeTyping, self).update_parameters(**kwargs)
        if self._parameters['multi_segment'].value.lower() == 'true':
            self._parameters['multi_segment'].value = True
        elif self._parameters['multi_segment'].value.lower() == 'false':
            self._parameters['multi_segment'].value = False
        else:
            raise InvalidParameterError(f'Invalid value given for multi_segment parameter: {self._parameters["multi_segment"].value}')
        if 'genome_segments' in self._parameters:
            self._parameters['genome_segments'].value = self._parameters['genome_segments'].value.strip().split(',')

    def _run_sequence_typing(self) -> None:
        """
        Runs the sequence typing by parsing the BLASTn results, grouping the hits per segment and then
        setting the informs for each segment.
        :return: None
        """
        blast_parser = InfluenzaBlastnAsnParser(self._tool_inputs['ASN'][0], self._parameters['multi_segment'].value,
                                                self._parameters['seqIDParser_type'].value, self._parameters['genometyping_method'].value)
        blast_parser.group_hits_per_segment()
        self._informs['segment_informs'] = {}
        segment_typer_class = SegmentTypingReads if self._parameters['genometyping_method'].value == 'blast' else SegmentTypingContigs
        if self._parameters['multi_segment'].value:
            for segment in self._parameters['genome_segments'].value:
                segment_hits = blast_parser.get_segment_hits(segment)
                if segment_hits:
                    segment_typer = segment_typer_class(segment, segment_hits, self._random_seed)
                    self._informs['segment_informs'][segment] = segment_typer.stats
            self._informs['expected_segments'] = self._parameters['genome_segments'].value
        else:
            segment_hits = blast_parser.get_segment_hits('single_segment')
            if segment_hits:
                segment_typer = segment_typer_class('single_segment', segment_hits, self._random_seed)
                self._informs['segment_informs']['single_segment'] = segment_typer.stats
            self._informs['expected_segments'] = ['Single segment']

    def _set_random_seed(self) -> None:
        """
        Sets the random seed that need so be used for tie breaking.
        :return: None
        """
        if 'random_seed' in self._parameters:
            self._random_seed = self._parameters['random_seed'].value
        else:
            self._random_seed = random.randint(1, 10000000)
            logging.info(f'Random seed for SegmentTyping not provided, setting seed to: {self._random_seed}')

    def _quality_check(self) -> None:
        """
        Sets the quality check informs. This includes the segments that were or were not covered and
        the percentage of segments that were found
        :return: None
        """
        self._informs['segment_coverage'] = {'segment_covered': [],
                                             'segment_missing': []}
        if self._parameters['multi_segment'].value:
            segments_found = self._informs['segment_informs'].keys()
            for segment in self._parameters['genome_segments'].value:
                if segment not in segments_found:
                    self._informs['segment_coverage']['segment_missing'].append(segment)
                else:
                    self._informs['segment_coverage']['segment_covered'].append(segment)
            self._informs['segment_coverage']['coverage'] = len(segments_found) * 100.0 / len(self._parameters['genome_segments'].value)
        else:
            if 'single_segment' in self._informs['segment_informs']:
                self._informs['segment_coverage']['segment_covered'].append('Single segment')
            else:
                self._informs['segment_coverage']['segment_missing'].append('Single segment')

    def _set_failed_informs(self) -> None:
        """
        Sets the failed informs so that it can be read later by other tools
        :return: None
        """
        self._informs['failed'] = self._failed

    def _obtain_reference_genome(self) -> str:
        """
        Obtain reference genomes from segment typing results.
        :return: reference genome fasta file
        """
        refseqs = self._retrieve_cluster_sequences()
        reference_segments = [refseqs[segment_inform['refseqid']] for _, segment_inform in self._informs['segment_informs'].items()]
        reference_fasta = os.path.join(self._folder, 'genome_reference_segments.fasta')
        FastaUtils.write(reference_segments, reference_fasta)
        return reference_fasta

    def _retrieve_cluster_sequences(self) -> str:
        """
        Load cluster representative sequences from the genome typing database
        :return: None
        """
        return FastaUtils.read_as_dict(self._tool_inputs['REF_FASTA'][0].path)

    def _set_output(self) -> None:
        """
        Sets the tool outputs
        :return: None
        """
        print(self._informs['segment_informs'])
        self._tool_outputs['FASTA'] = [ToolIOFile(self._obtain_reference_genome())]

    def _extract_influenza_a_subtype(self) -> None:
        """
        Extracts the subtype for Influenza A and sets the informs accordingly.
        :return: None
        """
        parts = self._extract_segment_subtypes()
        if parts['HA'] is None and parts['NA'] is None:
            self._informs['hana_subtyping'] = {'failure_message': 'Influenza A HANA subtyping failed to identify both the HA and NA subtype',
                                               'subtype': 'Unknown', 'ha': 'Unknown', 'na': 'unknown'}
        elif parts['HA'] is None:
            self._informs['hana_subtyping'] = {'failure_message': 'Influenza A HANA subtyping failed to identify the HA subtype',
                                               'subtype': 'Unknown', 'ha': 'Unknown', 'na': parts['NA']}
        elif parts['NA'] is None:
            self._informs['hana_subtyping'] = {'failure_message': 'Influenza A HANA subtyping failed to identify the NA subtype',
                                               'subtype': 'Unknown', 'ha': parts['HA'], 'na': 'Unknown'}
        else:
            subtype = f'{parts["HA"]}{parts["NA"]}'
            self._informs['hana_subtyping'] = {'subtype': subtype, 'ha': parts['HA'], 'na': parts['NA'], 'failure_message': None}

    def _extract_segment_subtypes(self) -> Dict[str, Union[str, None]]:
        """
        Extracts the subtype from the best reference for the HA and NA segment and returns
        a dictionary with that information.
        :return: Dictionary with subtypes of the HA and NA best reference
        """
        parts = {'HA': None, 'NA': None}
        for segment in ['HA', 'NA']:
            if segment in self._informs['segment_coverage']['segment_covered']:
                subtype = SeqIDParser(self._informs['segment_informs'][segment]['refseqid'], self._parameters['seqIDParser_type'].value).subtype
                logging.debug(f"Refseqid: {self._informs['segment_informs'][segment]['refseqid']} -- parser_type: {self._parameters['seqIDParser_type'].value}")
                if segment == 'HA':
                    parts['HA'] = subtype[:2]
                elif segment == 'NA':
                    parts['NA'] = subtype[-2:]
        return parts
