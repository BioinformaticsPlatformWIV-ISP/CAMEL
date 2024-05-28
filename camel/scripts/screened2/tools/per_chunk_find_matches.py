import math
from itertools import product
from pathlib import Path
from typing import List, Tuple
import logging

import Bio.Data.IUPACData as Bdi
import pandas as pd
from Bio import SeqIO
from fuzzysearch import find_near_matches

from camel.app.camel import Camel
from camel.app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from camel.app.error.toolexecutionerror import ToolExecutionError
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.tool import Tool


class FindMatches(Tool):
    """
    This function processes FASTA chunks, extracts primer sequences, and performs matching.
    """

    def __init__(self, camel: Camel) -> None:
        """
        Initialize FindMatches.
        :param camel: Camel instance
        :return: None
        """
        super().__init__('FindMatch', '0.1', camel)
        self._param = None

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        self._retrieve_parameters()
        self.__per_chunk_find_matches(self._parameters['output'].value)
        self.__set_output(Path(self._parameters['output'].value))

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'FASTA' not in self._tool_inputs:
            raise InvalidInputSpecificationError('FASTA input is required.')
        super()._check_input()

    def _retrieve_parameters(self) -> None:
        """
        Retrieve required parameters.
        :return: None
        """
        self.primer = str(self._parameters['primer'].value)
        self.fasta_primer_name = str(self._parameters['fasta_primer_name'].value)
        self.end_mismatch = int(self._parameters['end_mismatch'].value)
        self.perc_mismatch = float(self._parameters['perc_mismatch'].value)

    @staticmethod
    def extend_ambiguous_dna(primer_sequence: str) -> List:
        """
        Function has as input a string from the primer sequence and transforms it into a list. If degenerate nucleotides
        (anything other than A, C, T, G) are found, the primer sequence is split in a list all the possible primer
        sequences.
        :param primer_sequence: Input primer sequence as string
        :return: List of all possible sequences given the ambiguous DNA input.
        """
        logging.info('Transform ambiguous nucleotides and give all possibilities of the primer sequence')
        ambiguous_nucleotide_mappings = Bdi.ambiguous_dna_values
        primer_seq_ext = []
        for primer_seq in product(*[ambiguous_nucleotide_mappings[j] for j in primer_sequence]):
            primer_seq_ext += ["".join(primer_seq)]
        return primer_seq_ext

    @staticmethod
    def find_near_matches_with_last_nucleotides(sequence: SeqIO.SeqRecord, prim: str, primer_name: str,
                                                mismatch_end_input: int, max_sub: int) -> List[Tuple[int, str]]:
        """
        Function that checks for matches allowing a specific number of mismatches: The fuzzysearch function
        "find_near_matches" is used to match the primer sequence against one specific FASTA sequence allowing a certain
        number of mismatches between the two. If there cannot be a certain number of mismatches at the end of the
        primer sequence, then it is checked that whether there is perfect match between the primer sequence and the
        matched part of the FASTA sequence.
        :param sequence: Input sequence.
        :param prim: Primer sequence.
        :param primer_name: The name of the primer.
        :param mismatch_end_input: The allowed mismatch at the end of the sequence.
        :param max_sub: Maximum number of substitutions allowed.
        :return: List of tuples containing index and matched sequence.
        """
        result_multiple_prim = []
        if mismatch_end_input == 0:
            result_multiple_prim = find_near_matches(str(prim), str(sequence.seq).upper(),
                                                     max_deletions=max_sub, max_insertions=max_sub,
                                                     max_substitutions=max_sub, max_l_dist=max_sub)
        else:
            result_prim = find_near_matches(str(prim), str(sequence.seq).upper(), max_deletions=max_sub,
                                            max_insertions=max_sub, max_substitutions=max_sub, max_l_dist=max_sub)
            if result_prim:
                for matches in result_prim:
                    matched_seq = matches.matched
                    if "RV" in primer_name:
                        last_nucleotides_seq = matched_seq[:mismatch_end_input]
                        last_nucleotides_prim = prim[:mismatch_end_input]
                    else:
                        last_nucleotides_seq = matched_seq[-mismatch_end_input:]
                        last_nucleotides_prim = prim[-mismatch_end_input:]

                    if last_nucleotides_prim == last_nucleotides_seq:
                        result_multiple_prim.append(matches)

        return result_multiple_prim

    def find_presence_matches(self, sequence: SeqIO.SeqRecord, primer_seq: list, primer_name: str,
                              mismatch_end_input: int, max_perc_mismatch_input: float) -> pd.DataFrame:
        """
        Function that checks for if the sequences match. It calculates the number of maximum mismatches based on the
        percentage that is given and the size of the primer sequence. For each primer sequence the function
        find_near_matches_with_last_nucleotides is used to get the results for one FASTA sequence if there are matches
        with the information between which positions, the distance from the primer sequence and the FASTA sequence that
        matched.
        :param sequence: The input sequence.
        :param primer_seq: List of primer sequences.
        :param primer_name: The name of the primer.
        :param mismatch_end_input: The allowed mismatch at the end of the sequence.
        :param max_perc_mismatch_input: Maximum percentage of allowed mismatches.
        :return: DataFrame with match information.
        """
        max_sub = math.floor(max_perc_mismatch_input * len(primer_seq[0]))
        near_matches = []

        for prim in primer_seq:
            result_prim = self.find_near_matches_with_last_nucleotides(sequence, prim, primer_name, mismatch_end_input,
                                                                       max_sub)
            near_matches.append(result_prim)

        near_matches = [m for m in near_matches if m is not None]
        df_near_matches = pd.DataFrame([(m.start, m.end, m.dist, m.matched) for match_list in near_matches for m in
                                        match_list], columns=['start', 'end', 'dist', 'matched'])
        #Add a column with the id to the dataframe df_near_matches
        df_near_matches['id'] = [sequence.description] * df_near_matches.shape[0]
        #Drop duplicate columns
        df_near_matches.drop_duplicates(inplace=True)
        #Check that the number of mismatches (dist) are indeed below the number of allowed mismatches
        df_near_matches = df_near_matches[df_near_matches['dist'] <= max_sub]
        #Remove NA from dataframe
        df_near_matches.dropna(inplace=True)
        return df_near_matches

    def __per_chunk_find_matches(self, outfile_name: Path) -> None:
        """
        This function processes a fasta file by extracting primer sequences and match them to the sequences of the
        input FASTA file: First it will make a list of the primer sequence(s) that need to be checked and if degenerate
        nucleotides are included in the primer sequence, this list will contain all possible sequences.
        The function "find_presence_matches" is applied on each sequence from the FASTA file resulting in a dataframe
        with id (of fasta sequence), start (position within the sequence where match starts), stop (position within the
        sequence where match stops), dist (the number of mismatches between the primer sequence and the FASTA sequence)
        and matched (the sequence from the FASTA sequence that matched). If there is match found, it's name is added to
        the list fasta_id. If there is no match found, an empty record is made (based on difference between all id's
        and the ones that were found). After conversion to panda dataframe it is saved as an Excel.
        :param outfile_name: CSV file paths
        :return: None
        """
        logging.info('Check the matches of the primers against the FASTA file')
        primer_seq_ext = self.extend_ambiguous_dna(self.primer)

        # Initialize list to store match information as dictionaries
        count_primer_records = []
        fasta_id = []

        with open(str(self._tool_inputs['FASTA'][0])) as handle:
            for seq in SeqIO.parse(handle, 'fasta'):
                # Find matches of primer sequence in the current sequence
                nb_seq = self.find_presence_matches(seq, primer_seq_ext, self.fasta_primer_name,
                                                    self.end_mismatch, self.perc_mismatch)

                # Accumulate match information as dictionaries in the list
                for match in nb_seq.itertuples(index=False):
                    count_primer_records.append({
                        'id': seq.description,
                        'start': match.start,
                        'end': match.end,
                        'dist': match.dist,
                        'matched': match.matched
                    })
                fasta_id.append(seq.description)

        # Identify missing FASTA IDs in the results to add empty rows to the dataframe
        missing_fasta_id = set(fasta_id) - set(record['id'] for record in count_primer_records)
        for f in missing_fasta_id:
            count_primer_records.append({
                'id': f,
                'start': None,
                'end': None,
                'dist': None,
                'matched': None
            })

        # Create DataFrame from the accumulated list of dictionaries
        count_primer = pd.DataFrame(count_primer_records)

        # Save the results to a CSV file
        count_primer.to_csv(outfile_name, index=False)

    def __set_output(self, outfile_path: Path) -> None:
        """
        Sets the tool output.
        :return: None
        """
        self._tool_outputs['CSV'] = [ToolIOFile(outfile_path)]

    def _check_command_output(self) -> None:
        """
        Checks if the command was executed successfully.
        :return: None
        """
        if self._command.returncode != 0:
            raise ToolExecutionError("Command execution failed (Exit code: {})".format(self._command.returncode))
