from app.tools.picard.picard import Picard


class AddOrReplaceReadGroups(Picard):

    """
    Class for Picard AddOrReplaceReadGroups function
    """
    DEFAULT_SAMPLE_NAME = 'sampleA'

    def __init__(self, camel):
        """
        Initialize a picard tool
        :param camel: Camel instance
        :return: None
        """
        super(AddOrReplaceReadGroups, self).__init__('Picard AddOrReplaceReadGroups', '2.8.3', camel)

        self._function_name = 'AddOrReplaceReadGroups'

    def _set_input(self):
        """
        Set input for the picard function
        :return: None
        """
        if 'SAMPLE_NAME' in self._tool_inputs:
            sample_name = self._tool_inputs['SAMPLE_NAME'][0].value
        else:
            sample_name = AddOrReplaceReadGroups.DEFAULT_SAMPLE_NAME
        self._input_string += " RGSM={0} RGID={0}".format(sample_name)
