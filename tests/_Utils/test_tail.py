import os

from BioSimSpace._Utils import Tail


def test_tail_incremental(tmp_path):
    """Only lines appended since the previous read are returned."""
    path = str(tmp_path / "log")

    with open(path, "w") as f:
        f.write("a\nb\n")

    tail = Tail(path)
    assert list(tail) == ["a\n", "b\n"]
    assert list(tail) == []

    with open(path, "a") as f:
        f.write("c\n")

    assert list(tail) == ["c\n"]


def test_tail_partial_line(tmp_path):
    """A partial trailing line is held back until it is complete."""
    path = str(tmp_path / "log")

    with open(path, "w") as f:
        f.write("a\nb")

    tail = Tail(path)
    assert list(tail) == ["a\n"]

    with open(path, "a") as f:
        f.write("c\n")

    assert list(tail) == ["bc\n"]


def test_tail_truncated(tmp_path):
    """Reading restarts from the beginning if the file shrinks."""
    path = str(tmp_path / "log")

    with open(path, "w") as f:
        f.write("a\nb\nc\n")

    tail = Tail(path)
    assert len(list(tail)) == 3

    with open(path, "w") as f:
        f.write("d\n")

    assert list(tail) == ["d\n"]


def test_tail_replaced(tmp_path):
    """Reading restarts from the beginning if the file is replaced."""
    path = str(tmp_path / "log")

    with open(path, "w") as f:
        f.write("a\nb\n")

    tail = Tail(path)
    assert len(list(tail)) == 2

    # Rename over the original so the new file has a different inode.
    with open(path + ".new", "w") as f:
        f.write("c\nd\n")
    os.replace(path + ".new", path)

    assert list(tail) == ["c\n", "d\n"]


def test_tail_missing(tmp_path):
    """A missing file yields nothing and is picked up once created."""
    path = str(tmp_path / "log")

    tail = Tail(path)
    assert list(tail) == []

    with open(path, "w") as f:
        f.write("a\n")

    assert list(tail) == ["a\n"]


def test_tail_crlf(tmp_path):
    """Windows line endings are normalised."""
    path = str(tmp_path / "log")

    with open(path, "wb") as f:
        f.write(b"a\r\nb\r\n")

    assert list(Tail(path)) == ["a\n", "b\n"]
