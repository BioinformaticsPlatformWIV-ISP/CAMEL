class BlastHit(object):
    """
    Class that represents a Blast hit.
    The blast output file has to be generated with the following output format:
    -outfmt "6 pident sseqid sseq slen qseqid qstart qend"
    """

    def __init__(self, percent_identity, database_gene, database_gene_aligned_sequence, database_gene_length, query,
                 query_start, query_end):
        """
        Initializes a Blast hit.
        :param percent_identity: Percent identity
        :param database_gene: Gene in the database that was used as subject
        :param database_gene_aligned_sequence: Aligned sequence of the gene
        :param database_gene_length: Length of the full gene
        :param query: Query sequence id
        :param query_start: Location of the alignment start in the query
        :param query_end: Location of the alignment end in the query
        """
        self.percent_identity = float(percent_identity)
        self.database_gene = database_gene
        self.database_gene_aligned_sequence = database_gene_aligned_sequence
        self.database_gene_length = int(database_gene_length)
        self.query = query
        self.query_start = min(int(query_start), int(query_end))
        self.query_end = max(int(query_start), int(query_end))
        self.nb_of_gaps = database_gene_aligned_sequence.count('-')

    @property
    def type(self):
        """
        Returns the hit type.
        :return: Hit type
        """
        if self.database_gene_length == self.alignment_length and self.percent_identity == 100.0:
            return 'Perfect'
        if self.database_gene_length == self.alignment_length:
            return 'Imperfect identity'
        else:
            return 'Imperfect short'

    @property
    def color(self):
        """
        Returns the color for this hit.
        :return: Color
        """
        if self.type == 'Perfect':
            return 'green'
        elif self.type == 'Imperfect identity':
            return 'lightgreen'
        elif self.type == 'Imperfect short':
            return 'grey'
        raise ValueError("Unknown type: {}".format(self.type))

    @property
    def alignment_length(self):
        """
        Returns the alignment length.
        :return: Alignment length
        """
        return len(self.database_gene_aligned_sequence)

    @property
    def database_gene_coverage(self):
        """
        Returns the percentage of the database gene that is covered by the alignment.
        :return: % query covered
        """
        return 100 * float(self.alignment_length) / self.database_gene_length

    def __str__(self):
        """
        Returns the string representation of a Blast hit.
        :return: String representation
        """
        return '{} ({}/{})\t{}\t(%ID: {}, %QCOV: {})\t{}:{}'.format(
            self.database_gene, self.alignment_length, self.database_gene_length, self.query, self.percent_identity,
            self.database_gene_coverage, self.query_start, self.query_end)

    @staticmethod
    def convert_blast_output_line(line):
        """
        Converts a blast output line to a BlastHit object.
        :param line: Output line
        :return: BlastHit
        """
        output_values = line.split('\t')
        if not len(output_values) == 7:
            raise ValueError("Invalid Blast output line: {}".format(line))

        hit = BlastHit(output_values[0], output_values[1], output_values[2], output_values[3], output_values[4],
                       output_values[5], output_values[6])
        return hit
