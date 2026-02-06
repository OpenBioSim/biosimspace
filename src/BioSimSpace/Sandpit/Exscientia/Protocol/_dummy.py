__all__ = ["Dummy"]
from ._protocol import Protocol as _Protocol


class Dummy(_Protocol):
    """A class for storing Dummy protocols."""

    def __init__(self):
        super().__init__()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return f"<BioSimSpace.Protocol.Dummy>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return f"BioSimSpace.Protocol.Dummy()"
