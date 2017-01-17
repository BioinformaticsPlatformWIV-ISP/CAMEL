import logging

from app.tools.tool import Tool


class SequenceTypeDetector(Tool):
    """
    Tool that manages MLST schemes. Also reports scheme metadata information in the informs.
    """

    def __init__(self, camel):
        """
        Initialize this tool.
        :param camel: Camel instance
        :return: None
        """
        super(SequenceTypeDetector, self).__init__('Sequence Typing: Sequence Type Detector', '0.1', camel)

    def _execute_tool(self):
        """
        Executes this tool.
        :return: None
        """
        hits = self._input_informs['hits']
        allele_ids = {locus: hits[locus]['allele_id'] for locus in hits}
        sequence_type = self.__get_sequence_type(self._tool_inputs['TSV'][0], allele_ids)
        self._informs['sequence_type'] = sequence_type
        logging.info("Detected sequence type: {}".format(sequence_type))

    def _check_input(self):
        """
        Checks whether the input is correct.
        :return: None
        """
        if 'TSV' not in self._tool_inputs:
            raise ValueError("No sequence type definitions 'TSV' input found.")
        if 'hits' not in self._input_informs:
            raise ValueError("No hits info found.")
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
            for gene_name in gene_names:
                try:
                    gene_indices[gene_name] = header.strip().split('\t').index(gene_name)
                except ValueError:
                    raise StandardError("Gene {} not found in '{}'".format(gene_name, profiles_file))
        return gene_indices

    def __get_sequence_type(self, profiles_file, gene_alleles):
        """
        Returns the sequence type corresponding to the given gene alleles.
        :param profiles_file: File containing the ST profiles
        :param gene_alleles: Alleles for the genes.
        :return: Sequence type and metadata
        """
        gene_indices = self.__get_gene_indices(gene_alleles.keys(), profiles_file)
        with open(profiles_file.path) as profiles:
            content = profiles.readlines()
            for line in content[1:]:
                match = True
                for gene_name in gene_alleles:
                    if match:
                        st_allele_id = line.split('\t')[gene_indices[gene_name]]
                        detected_allele_id = gene_alleles[gene_name]
                        if not st_allele_id == detected_allele_id:
                            match = False
                            break
                if match:
                    return line.split('\t')[0]
        return 'ND'
