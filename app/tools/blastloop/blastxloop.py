from app.tools.blastloop.blastloop import BlastLoop


class BlastxLoop(BlastLoop):
    """
    Tool that loops over input files to run blastx multiple times.
    """

    def __init__(self, camel):
        """
        Initializes this tool.
        """
        super(BlastxLoop, self).__init__('blastx (looping)', '2.6.0', camel)
