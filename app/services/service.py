class Service(object):
    """
    Parent class for services that interact with the database.
    """

    def __init__(self, connection):
        """
        Initializes a connection object for transactions (update and insert)
        Initializes a cursor to perform queries and iterate through objects
        :param connection: Connection
        :return: None
        """
        self.db_connection = connection
