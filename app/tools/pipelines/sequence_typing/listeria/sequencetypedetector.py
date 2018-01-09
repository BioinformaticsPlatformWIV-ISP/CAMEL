import logging

from app.components.filesystemhelper import FileSystemHelper
from app.error.invalidinputspecificationerror import InvalidInputSpecificationError
from app.tools.tool import Tool


class SequenceTypeDetector(Tool):

    """
    Tool that manages MLST schemes. Also reports scheme metadata information in the informs.
    """

    SYMBOL_NO_ST = 'ND'

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(SequenceTypeDetector, self).__init__('Typing: Sequence Type Detector', '0.1', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        # Get the sequence type
        profiles_file = self._tool_inputs['TSV'][0]
        allele_ids = {h.value.locus: h.value.allele_id if h.value.allele_id != '-' else '0' for h in self._tool_inputs['VAL_Hits']}
        logging.info("Detected allele ids: {}".format(allele_ids))
        sequence_type = self.__get_sequence_type(profiles_file, allele_ids)
        self._informs['sequence_type'] = sequence_type
        logging.info("Detected sequence type: {}".format(sequence_type))

        # Retrieve metadata
        column_indices = self.__get_metadata_columns(profiles_file, allele_ids.keys())
        if sequence_type != SequenceTypeDetector.SYMBOL_NO_ST:
            line = self.__get_sequence_type_line(profiles_file, sequence_type, 0)
            # Check for Neo type (which is not contained in profiles_file)
            if line is not None:
                line_parts = line.split('\t')
                metadata = {column: line_parts[index].strip() for index, column in column_indices.items()}
                logging.info("Metadata for sequence type: {!r}".format(metadata))
            else:
                metadata = {column: '-' for index, column in column_indices.items()}
                logging.info("Neo typing found for sequence type: {!r}".format(metadata))
        else:
            metadata = {column: '-' for index, column in column_indices.items()}
        self._informs['metadata'] = metadata

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise InvalidInputSpecificationError("Sequence type profiles are required")
        if 'VAL_Hits' not in self._tool_inputs:
            raise InvalidInputSpecificationError("No hits input found")
        super(SequenceTypeDetector, self)._check_input()

    @staticmethod
    def __get_gene_indices(gene_names, profiles_file):
        """
        Returns the gene indices from the given MLST profiles file.
        :param gene_names: Names of the genes
        :param profiles_file: MLST profiles file
        :return: Dictionary of gene name : index
        """
        gene_indices = {}
        with open(profiles_file.path) as profiles:
            header = profiles.readline()
            genes_in_header = [FileSystemHelper.make_valid(x) for x in header.strip().split('\t')]
            for gene_name in gene_names:
                try:
                    gene_indices[gene_name] = genes_in_header.index(gene_name)
                except ValueError:
                    raise Exception("Gene {} not found in '{}'".format(gene_name, profiles_file))
        return gene_indices

    def __get_sequence_type(self, profiles_file, gene_alleles):
        """
        Returns the sequence type corresponding to the given gene alleles.
        :param profiles_file: File containing the ST profiles
        :param gene_alleles: Alleles for the genes.
        :return: Sequence type and metadata
        """
        gene_indices = self.__get_gene_indices(list(gene_alleles.keys()), profiles_file)
        with open(profiles_file.path) as profiles:
            content = profiles.readlines()
            for line in content[1:]:
                logging.debug(f'profiles line: {line.strip()}')
                match = True
                for gene_name in gene_alleles:
                    if match:
                        st_allele_id = line.strip().split('\t')[gene_indices[gene_name]]
                        detected_allele_id = gene_alleles[gene_name]
                        logging.debug(f'gene_name {gene_name}: st_allele_id {st_allele_id}, detected_allele_id {detected_allele_id}')
                        if st_allele_id not in ('N') and st_allele_id != detected_allele_id:
                        #if st_allele_id not in ('N', '0') and st_allele_id != detected_allele_id:
                            match = False
                            break
                if match:
                    return line.split('\t')[0]
        return SequenceTypeDetector.SYMBOL_NO_ST

    @staticmethod
    def __get_sequence_type_line(scheme_file, sequence_type, index=0):
        """
        Returns the line corresponding to the given sequence type.
        :param scheme_file: Scheme file
        :param sequence_type: Sequence type
        :param index: Index of the sequence type
        :return: Line containing the sequence type
        """
        clean_sequence_type = sequence_type.replace('?', '').replace('*', '')
        with open(scheme_file.path) as input_mlst:
            lines = input_mlst.readlines()
            for line in lines:
                st = line.split('\t')[index]
                if st == clean_sequence_type:
                    return line

    @staticmethod
    def __get_metadata_columns(mlst_scheme_file, gene_names):
        """
        Returns the metadata columns from the profiles file.
        :param mlst_scheme_file: Profiles file
        :param gene_names: Names of the genes (are not included in metadata)
        :return: Metadata columns
        """
        with open(mlst_scheme_file.path) as handle:
            header = handle.readline().strip().split('\t')

        metadata_columns = {}
        for i in range(0, len(header)):
            if header[i] not in gene_names:
                metadata_columns[i] = header[i].strip().replace('_', ' ')
        return metadata_columns
