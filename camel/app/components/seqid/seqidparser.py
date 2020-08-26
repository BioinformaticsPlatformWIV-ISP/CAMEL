import logging
import re
from typing import Optional


class SeqIDParser(object):

    def __init__(self, seqid: str, parser_type: str, add_gbi: bool = True):
        """
        Initializes the parser object
        :param seqid: Sequence ID that needs to be parsed
        :param parser_type: Format of the sequence ID
        :param add_gbi: Add the GBI field to a newly composed ID (for Influenza)
        """
        # Parser types:
        #   - single = no segment information
        #   - default = strain and segment divided by '|'
        #   - avian = e.g. A-DUCK-Czech_Republic-1-2011-H1N1-NA
        #   - cyril = e.g. A-Alabama-05-2010-H3N2-HA
        #   - inf_a = e.g. A-ruddy-Delaware-34-1993-H2N1(CY015141)-PB1
        #   - inf_b = e.g. B-Wellington-7-2008(CY149972)-PB2
        #   - inf_c = e.g. C-Yamagata-64(LC123489)-P3
        self._type = parser_type
        self._seqid = seqid
        # XXXXXX-H#N#(YYYYYY)-NA, e.g. A-swine-Argentina-CIP051-C05-M1-5-2012-H1N2(KR863420)-HA
        # or XXXXX-H#N#-NA, e.g. A-Alabama-05-2010-H3N2-HA
        self._inf_a_regex = r'(\S*)-([H0-9]+N[0-9]+[\w?]*)(\(*[\w]*\)*)-([HNAEMPBS123]+)$'
        # XXXXXXXXX-YYYYYY(ZZZZZZ)-NA, e.g. B-Wellington-7-2008(CY149972)-PB2
        # or XXXXXX-YYYYYY-NA, e.g. B-North_Dakota-02-2017-PB2
        self._inf_bc_regex = r'(\S*-[\w]*)(\(*[\w]*\)*)-([HNAEMPBS123]+)$'
        self._re_match = self._get_re_match()
        self._add_gbi = add_gbi

    @property
    def seqid(self) -> str:
        """
        Returns the sequence ID with '-extracted' removed in case of an avian sample
        or a sample coming from Cyril.
        :return: Sequence ID
        """
        if self._type in {'avian', 'cyril'}:
            return self._seqid.replace('-extracted', '')
        return self._seqid

    @property
    def strain(self) -> str:
        """
        Returns the strain part of the sequence ID for instance, for A-ruddy-Delaware-34-1993-H2N1(CY015141)-PB1
        A-ruddy-Delaware-34-1993-H2N1 would be returned.
        :return:
        """
        if self._type == 'single':
            return self._seqid
        elif self._type == 'default':
            return self._seqid.split('|')[0]
        elif self._type in {'avian', 'cyril', 'inf_a'}:
            return f'{self._re_match.group(1)}-{self._re_match.group(2)}'
        elif self._type in {'inf_b', 'inf_c'}:
            return self._re_match.group(1)

    @property
    def segment(self) -> Optional[str]:
        """
        Returns the segment contained in the sequence ID if it is present. For instance,
        A-ruddy-Delaware-34-1993-H2N1(CY015141)-PB1 will return PB1
        :return: Segment if present else None
        """
        if self._type == 'single':
            return None
        elif self._type == 'default':
            return self._seqid.split('|')[1]
        elif self._type in {'avian', 'cyril', 'inf_a'}:
            return self._re_match.group(4)
        elif self._type in {'inf_b', 'inf_c'}:
            return self._re_match.group(3)

    @property
    def subtype(self) -> Optional[str]:
        """
        Returns the subtype in case of an Influenza A sequence ID, else None.
        For instance, A-ruddy-Delaware-34-1993-H2N1(CY015141)-PB1 will return H2N1
        :return: Subtype for Influenza A if present
        """
        if self._type in {'avian', 'cyril', 'inf_a'}:
            return self._re_match.group(2)

    @property
    def gbi(self) -> Optional[str]:
        """
        Returns the GBI for Influenza sequence IDs. For instance, A-ruddy-Delaware-34-1993-H2N1(CY015141)-PB1
        will return (CY015141)
        :return:
        """
        if self._type in {'avian', 'cyril', 'inf_a'}:
            return self._re_match.group(3)
        elif self._type in {'inf_b', 'inf_c'}:
            return self._re_match.group(2)

    @property
    def composed_id(self) -> Optional[str]:
        """
        Creates a newly composed ID from the given sequence ID. The ID can contain the GBI depending on whether
        it is present in the original sequence ID and if it is requested to be added. When the subtype is missing
        from the original sequence ID, it is assumed that the subtype information was present in the strain name.
        :return: Newly composed ID
        """
        if self._type in {'avian', 'cyril', 'inf_a', 'inf_b', 'inf_c'} and self.subtype is None:
            logging.warning(f'Composing new sequence ID: Influenza subtype information is missing from original seqid, '
                            f'assuming subtype information in strain name: {self.strain}.')
        if self._type == 'default':
            return f'{self.strain}|{self.segment}'
        elif self._type in {'avian', 'cyril', 'inf_a', 'inf_b', 'inf_c'}:
            if self.subtype is None:
                return f'{self.strain}{self._get_gbi_string()}-{self.segment}'
            else:
                return f'{self.strain}-{self.subtype}{self._get_gbi_string()}-{self.segment}'

    def _get_gbi_string(self) -> str:
        """
        Return the string to be used for the GID when composing a new ID. Returns an empty string in case of
        avian or Cyril samples as they do not have/need this information.
        :return: String with GBI information
        """
        if self._type in {'avian', 'cyril'} or not self._add_gbi:
            return ''
        elif self._add_gbi:
            return f'{self.gbi}' if self.gbi else f'(GBIDna)'

    def _get_re_match(self) -> re.Match:
        """
        Returns the regex match object that can be used to extract the different subparts from the sequence ID.
        :return: Regex match object
        """
        if self._type in {'avian', 'cyril', 'inf_a'}:
            return re.match(self._inf_a_regex, self.seqid)
        elif self._type in {'inf_b', 'inf_c'}:
            return re.match(self._inf_bc_regex, self.seqid)
