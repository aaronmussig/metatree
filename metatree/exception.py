class MetaTreeException(Exception):
    """Base class for all MetaTree exceptions."""

    def __init__(self, message):
        Exception.__init__(self, message)


class MetaTreeExit(MetaTreeException):
    """Raised when a controlled exit is requested."""

    def __init__(self, message=''):
        super(MetaTreeExit, self).__init__(message)
