from app.tools.blastloop.blastloop import BlastLoop


class BlastnLoop(BlastLoop):
    """
    Tool that loops over input files to run blastn multiple times.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        :param camel: Camel instance
        """
        super(BlastnLoop, self).__init__('blastn (looping)', '2.2.30', camel)
