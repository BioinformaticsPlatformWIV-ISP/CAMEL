import datetime
import json
import logging
from typing import List, Dict

import humanize
import os
from Bio import SeqIO

from camel.app.camel import Camel
from camel.app.io.tooliofile import ToolIOFile
from camel.app.tools.blast.makeblastdb import MakeBlastDb
from camel.app.tools.bowtie2.bowtie2index import Bowtie2Index
from camel.app.tools.cdhit.cdhitest import Cluster, CDHitEst
from camel.app.tools.samtools.samtoolsfastaindex import SamtoolsFastaIndex


class DBHelper(object):
    """
    This helper class is used to construct gene detection databases.
    """

    def __init__(self, db_name: str, working_dir: str) -> None:
        """
        Initializes the update helper.
        :param db_name: Database name
        :param working_dir: Working directory
        """
        self._db_name = db_name
        self._working_dir = working_dir

    def get_working_subdir(self, name: str) -> str:
        """
        Returns the path to the given sub directory.
        The directory is created if it does not exist yet.
        :param name: Directory name
        :return: Path
        """
        path = os.path.join(self._working_dir, name)
        if not os.path.isdir(path):
            logging.debug(f"Creating directory: '{path}'")
            os.makedirs(path)
        return path

    def get_clusters_form_fasta(self, fasta_file: str, clustering_cutoff: float) -> List[Cluster]:
        """
        Returns the clusters of similar sequences from the given FASTA file.
        :param fasta_file: Input FASTA file
        :param clustering_cutoff: Clustering cutoff (0.0 - 1.0)
        :return: List of clusters
        """
        cdhit = CDHitEst(Camel.get_instance())
        cdhit.update_parameters(identitiy_threshold=str(clustering_cutoff / 100))
        cdhit.add_input_files({'FASTA': [ToolIOFile(fasta_file)]})
        cdhit.run(self.get_working_subdir('clustering'))
        return cdhit.informs['clusters']

    def index_samtools_faidx(self, fasta_file: str, working_dir: str) -> None:
        """
        Indexes the given FASTA file with samtools faidx.
        :param fasta_file: Input FASTA file
        :param working_dir: Working directory
        :return: None
        """
        samtools_faindex = SamtoolsFastaIndex(Camel.get_instance())
        samtools_faindex.add_input_files({'FASTA': [ToolIOFile(fasta_file)]})
        samtools_faindex.run(working_dir)

    def index_bowtie2(self, fasta_file: str, working_dir: str) -> None:
        """
        Creates a bowtie2 index for the given FASTA file.
        :param fasta_file: Input FASTA file
        :param working_dir: Working directory
        :return: None
        """
        bowtie2_index = Bowtie2Index(Camel.get_instance())
        bowtie2_index.add_input_files({'FASTA_REF': [ToolIOFile(fasta_file)]})
        bowtie2_index.run(working_dir)

    def index_blast(self, fasta_file: str, working_dir: str) -> None:
        """
        Indexes the given FASTA file with makeblastdb.
        :param fasta_file: Input FASTA file
        :param working_dir: Working directory
        :return: None
        """
        makeblastdb = MakeBlastDb(Camel.get_instance())
        makeblastdb.add_input_files({'FASTA': [ToolIOFile(fasta_file)]})
        makeblastdb.run(working_dir)

    def convert_fasta_headers_to_seq(self, input_file: str, output_file: str) -> Dict[str, str]:
        """
        Creates a fasta file where all ids are replaced by seq_{number}.
        This ensures that CD-HIT can work properly.
        :param input_file: Input FASTA file
        :param output_file: Output FASTA file
        :return: Mapping of original headers to novel headers
        """
        seq_ids = {}
        output_seqs = []
        with open(input_file) as handle_in:
            for seq in SeqIO.parse(handle_in, 'fasta'):
                new_name = 'seq_{}'.format(len(seq_ids))
                seq_ids[new_name] = seq.description
                seq.id = new_name
                seq.description = ''
                output_seqs.append(seq)
        with open(output_file, 'w') as handle_out:
            SeqIO.write(output_seqs, handle_out, 'fasta')
        return seq_ids

    @staticmethod
    def create_srst2_fasta(input_fasta: str, output_fasta: str, clusters: List[Cluster]) -> None:
        """
        Creates a FASTA file compatible with SRST2.
        The format is: >[clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]
        :param input_fasta: FASTA containing the renamed sequences
        :param output_fasta: Output FASTA file
        :param clusters: Sequence clusters
        :return: Path to generated FASTA file
        """
        seq_record_by_id = {}
        with open(input_fasta) as handle_in:
            for seq in SeqIO.parse(handle_in, 'fasta'):
                seq_record_by_id[seq.id] = seq
                seq.description = ''

        output_seqs = []
        for i in range(0, len(clusters)):
            for sequence_name in clusters[i].seq_ids:
                seq = seq_record_by_id[sequence_name]
                full_name = seq.id
                seq.id = '__'.join([str(i), clusters[i].name, full_name, full_name])
                output_seqs.append(seq)

        with open(output_fasta, 'w') as handle_out:
            SeqIO.write(output_seqs, handle_out, 'fasta')

    def export_metadata(self, name: str, dir_output: str) -> None:
        """
        Exports the database metadata.
        :param name: Database name
        :param dir_output: Output directory
        :return: None
        """
        metadata_file = os.path.join(dir_output, 'db_metadata.txt')
        logging.info(f'Exporting metadata: {metadata_file}')
        metadata = {'name': name.lower(), 'title': name, 'last_updated': datetime.date.today().strftime("%d-%m-%Y")}
        with open(metadata_file, 'w') as handle:
            json.dump(metadata, handle, indent=4, sort_keys=True)
        logging.info(f"Metadata exported: {metadata_file}")

    def standardize_fasta_headers(self, input_fasta: str) -> str:
        """
        Reformat the headers of the input FASTA file.
        :param input_fasta: Input FASTA file
        :return: Reformatted FASTA file
        """
        dir_reformat = os.path.join(self._working_dir, 'reformat')
        if not os.path.isdir(dir_reformat):
            os.makedirs(dir_reformat)
        output_path = os.path.join(dir_reformat, f'{os.path.splitext(os.path.basename(input_fasta))[0].lower()}.fasta')
        with open(input_fasta) as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
        for s in seqs:
            data = {'allele': s.id}
            s.description = json.dumps(data)
        with open(output_path, 'w') as handle:
            SeqIO.write(seqs, handle, 'fasta')
        logging.info(f"Reformatted FASTA file created ({humanize.naturalsize(os.path.getsize(output_path))}): "
                     f"{output_path}")
        return output_path

    def export_mapping(self, mapping: Dict[str, str], output_directory: str) -> None:
        """
        Exports the mapping of the novel headers to the original headers.
        :param mapping: Mapping
        :param output_directory: Output directory
        :return: None
        """
        with open(os.path.join(output_directory, 'mapping.txt'), 'w') as handle:
            json.dump(mapping, handle, indent=4, sort_keys=True)
        logging.info(f"Metadata exported: {output_directory}")
