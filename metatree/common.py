import logging
import os
import sys


def make_sure_path_exists(path):
    """Create directory if it does not exist.

    Parameters
    ----------
    path : str
        The path to the directory which should be created.

    Returns
    -------
    bool
        True if the path exists.

    Raises
    ------
    OSError
        If an error was encountered while creating the directory.
    """
    if not path:
        # lack of a path qualifier is acceptable as this
        # simply specifies the current directory
        return True
    elif os.path.isdir(path):
        return True

    try:
        os.makedirs(path)
        return True
    except OSError:
        logger = logging.getLogger('timestamp')
        logger.error('Specified path could not be created: ' + path)
        raise


def is_executable(fpath):
    """Check if file is executable.

    This is a Python implementation of the linux
    command 'which'.

    Parameters
    ----------
    fpath : str
        Path to file.

    Returns
    -------
    boolean
        True if executable, else False.
    """

    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """Return path to program.

    This is a Python implementation of the linux
    command 'which'.

    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Parameters
    ----------
    program : str
        Name of executable for program.

    Returns
    -------
    str
        Path to executable, or None if it isn't on the path.
    """

    fpath, _fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_executable(exe_file):
                return exe_file

    return None


def check_on_path(program, exit_on_fail=True):
    """Check if program is on the system path.

    Parameters
    ----------
    program : str
        Name of executable for program.
    exit_on_fail : boolean
        Exit program with error code -1 if program in not on path.

    Returns
    -------
    boolean
        True if program is on path, else False.
    """

    if which(program):
        return True

    if exit_on_fail:
        print('%s is not on the system path.' % program)
        sys.exit(1)

    return False


def rgba_hex(r, g, b, a=1):
    r, g, b = [int(x*a) for x in (r, g, b)]
    return '#%02x%02x%02x' % (r, g, b)


def hex_rgb(v):
    h = v.lstrip('#')
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))
