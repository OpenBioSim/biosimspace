from .._Utils import _try_import


import os as _os

_yaml = _try_import("yaml")


# Set the default node directory.
_node_dir = _os.path.dirname(__file__) + "/_nodes"

__all__ = ["list", "help", "run", "setNodeDirectory", "getNodeDirectory"]


def list():
    """Return a list of the available nodes."""
    from glob import glob as _glob

    # Glob all Python scripts in the _nodes directory.
    nodes = _glob("%s/*.py" % _node_dir)

    # Strip the extension.
    nodes = [_os.path.basename(x).split(".py")[0] for x in nodes]

    return nodes


def help(name):
    """
    Print the help message for the named node.

    Parameters
    ----------

    name : str
        The name of the node.
    """
    from .. import _Utils
    import subprocess as _subprocess
    from sire.legacy import Base as _SireBase

    if not isinstance(name, str):
        raise TypeError("'name' must be of type 'str'.")

    # Apped the node directory name.
    full_name = _node_dir + "/" + name

    # Make sure the node exists.
    if not _os.path.isfile(full_name):
        if not _os.path.isfile(full_name + ".py"):
            raise ValueError(
                "Cannot find node: '%s'. " % name
                + "Run 'Node.list()' to see available nodes!"
            )
        else:
            full_name += ".py"

    # Create the command.
    command = "%s/python %s --help" % (_SireBase.getBinDir(), full_name)

    # Run the node as a subprocess.
    proc = _subprocess.run(
        _Utils.command_split(command), shell=False, text=True, stdout=_subprocess.PIPE
    )

    # Print the standard output, decoded as UTF-8.
    print(proc.stdout)


def run(name, args={}, work_dir=None):
    """
    Run a node.

    Parameters
    ----------

    name : str
        The name of the node.

    args : dict
        A dictionary of arguments to be passed to the node.

    work_dir : str, optional
        The working directory in which to run the node. If not specified,
        the current working directory is used. Note that inputs should
        use absolute paths if this is set.

    Returns
    -------

    output : dict
        A dictionary containing the output of the node.
    """
    from .. import _Utils
    import subprocess as _subprocess
    from sire.legacy import Base as _SireBase

    # Validate the input.

    if not isinstance(name, str):
        raise TypeError("'name' must be of type 'str'.")

    if not isinstance(args, dict):
        raise TypeError("'args' must be of type 'dict'.")

    if work_dir is not None:
        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
    else:
        work_dir = _os.getcwd()

    # Apped the node directory name.
    full_name = _node_dir + "/" + name

    # Make sure the node exists.
    if not _os.path.isfile(full_name):
        if not _os.path.isfile(full_name + ".py"):
            raise ValueError(
                "Cannot find node: '%s'. " % name
                + "in directory '%s'. " % _node_dir
                + "Run 'Node.list()' to see available nodes!"
            )
        else:
            full_name += ".py"

    with _Utils.cd(work_dir):
        # Write a YAML configuration file for the BioSimSpace node.
        if len(args) > 0:
            with open("input.yaml", "w") as file:
                _yaml.dump(args, file, default_flow_style=False)

            # Create the command.
            command = "%s/python %s --config input.yaml" % (
                _SireBase.getBinDir(),
                full_name,
            )

        # No arguments.
        else:
            command = "%s/python %s" % (_SireBase.getBinDir(), full_name)

        # Run the node as a subprocess.
        proc = _subprocess.run(
            _Utils.command_split(command),
            shell=False,
            text=True,
            stderr=_subprocess.PIPE,
        )

        if proc.returncode == 0:
            # Read the output YAML file into a dictionary.
            with open("output.yaml", "r") as file:
                output = _yaml.safe_load(file)

            # Delete the redundant YAML files.
            _os.remove("input.yaml")
            _os.remove("output.yaml")

            return output

        else:
            # Print the standard error, decoded as UTF-8.
            print(proc.stderr)


def setNodeDirectory(dir):
    """
    Set the directory of the node library.

    Parameters
    ----------

    dir : str
        The path to the node library.
    """

    if not _os.path.isdir(dir):
        raise IOError("Node directory '%s' doesn't exist!" % dir)

    # Use the absolute path.
    dir = _os.path.abspath(dir)

    global _node_dir
    _node_dir = dir


def getNodeDirectory():
    """
    Get the directory of the node library.

    Returns
    -------

    dir : str
        The path to the node library.
    """
    return _node_dir
