class Parameter(object):
    """
    Represents a tool parameter.
    """

    def __init__(self, name, option, value):
        """
        Initializes a parameter.
        """
        self._name = name
        self._option = option
        self._value = value

    @property
    def name(self):
        """
        Returns the parameter name.
        :return: Parameter name
        """
        return self._name

    @property
    def option(self):
        """
        Returns the parameter option.
        :return: Option
        """
        return self._option

    @property
    def value(self):
        """
        Returns the parameter value.
        :return: parameter value
        """
        return self._value

    @value.setter
    def value(self, value):
        """
        Changes the value of this parameter.
        :param value: New value
        :return: None
        """
        self._value = value

    def __str__(self):
        """
        Retruns the parameter in string form.
        :return: Parameter string
        """
        if self._value is not None:
            return '{} {}'.format(self._option, self._value)
        else:
            return self._option
