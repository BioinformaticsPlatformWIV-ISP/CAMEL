import os
from typing import List

from camel.app.command.command import Command
from camel.app.components.blasthit.blastnhit import BlastnHit


class BlastnAsnParser(object):

    """
    Class to parse blastn output from an asn file and return hits
    """

    def __init__(self, blastn_file: str, columns: List[str] = None, seq_columns: bool = False, exclude_tax_columns: bool = False,
                 folder: str = None):
        """
        Initialize the object
        :param blastn_file: Location of the BLASTn file
        :param columns: Optional list of columns that need to be parsed
        :param seq_columns: Include sequences in the BLASTn hit objects
        :param exclude_tax_columns: Should columns that require a taxonomy db be excluded
        :param folder: Optional folder to run the commands in
        """
        self._column_types = {'qseqid': 'str', 'qgi': 'int', 'qacc': 'str', 'qaccver': 'str', 'qlen': 'int',
                              'sseqid': 'str', 'sallseqid': 'str', 'sgi': 'int', 'sallgi': 'int', 'sacc': 'str',
                              'saccver': 'str', 'sallacc': 'str', 'slen': 'int', 'qstart': 'int', 'qend': 'int',
                              'sstart': 'int', 'send': 'int', 'qseq': 'str', 'sseq': 'str', 'evalue': 'float',
                              'bitscore': 'float', 'score': 'int', 'length': 'int', 'pident': 'float',
                              'nident': 'int', 'mismatch': 'int', 'positive': 'int', 'gapopen': 'int',
                              'gaps': 'int', 'ppos': 'float', 'frames': 'str', 'qframe': 'int', 'sframe': 'int',
                              'btop': 'str', 'staxids': 'str', 'sscinames': 'str', 'scomnames': 'str',
                              'sblastnames': 'str', 'sskingdoms': 'str', 'stitle': 'str', 'sstrand': 'str',
                              'salltitles': 'str', 'qcovs': 'int', 'qcovhsp': 'int'}
        self._all_columns = ('qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi',
                             'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq',
                             'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch',
                             'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop',
                             'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'sstrand',
                             'salltitles', 'qcovs', 'qcovhsp')  # Fixed order
        self._blastn_file = blastn_file
        self._columns = columns if columns else self._get_columns(seq_columns, exclude_tax_columns)
        self._hits = []
        self._folder = folder if folder else os.getcwd()
        self._parse()

    @property
    def hits(self) -> List[BlastnHit]:
        return self._hits

    @property
    def columns(self) -> List[str]:
        return self._columns

    def _get_columns(self, seq_columns: bool, exclude_tax_columns: bool) -> List[str]:
        """
        Returns the possible output columns from the BLASTn file. If sequences are
        not requested, these columns are removed.
        :param seq_columns: Include columns with sequences
        :param exclude_tax_columns: Should columns that require a taxonomy db be excluded
        :return: List of column names
        """
        columns = list(self._all_columns)
        if not seq_columns:
            columns.remove('sseq')
            columns.remove('qseq')
        if exclude_tax_columns:
            for col in ['sscinames', 'sskingdoms', 'scomnames', 'sblastnames']:
                try:
                    columns.remove(col)
                except ValueError:
                    pass
        return columns

    def _parse(self) -> None:
        """
        Parses the BLASTn ASN file and adds BlastnHit objects to the _hits object variable.
        :return: None
        """
        format_columns = ' '.join(self._columns)
        cmd = Command(f'module load blast; blast_formatter -archive {self._blastn_file} -outfmt "6 {format_columns}" > formatted_blast_output.txt')
        cmd.run_command(self._folder)
        with open(os.path.join(self._folder, 'formatted_blast_output.txt')) as handle:
            for line in handle:
                if line:
                    self._hits.append(BlastnHit(**dict(zip(self._columns, line.split('\t')))))
