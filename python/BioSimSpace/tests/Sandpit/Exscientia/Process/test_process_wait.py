from unittest.mock import MagicMock

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Process._process import Process


def test_max_time():
    process = MagicMock()
    Process.wait(process, max_time=1)
    process._process.wait.assert_called_once_with(60000)


def test_None_inactivity_timeout():
    process = MagicMock()
    Process.wait(process, max_time=None, inactivity_timeout=None)
    process._process.wait.assert_called_once()


def test_inactivity_timeout_no_getLastTime():
    process = MagicMock()
    process._getLastTime.return_value = None
    Process.wait(process, max_time=None, inactivity_timeout=BSS.Units.Time.nanosecond)
    process._process.wait.assert_called_once()


def test_hang(tmp_path):
    process = MagicMock()
    process.workDir.return_value = str(tmp_path)
    process._name = "test"
    # Using TEST_HANG_COUNTER to mimic simulation progress
    global TEST_HANG_COUNTER
    TEST_HANG_COUNTER = 0
    process.isRunning.return_value = True

    def _getLastTime():
        global TEST_HANG_COUNTER
        TEST_HANG_COUNTER += 1
        # Mimic simulation hang after 10 calls
        return min(TEST_HANG_COUNTER, 10)

    process._getLastTime = _getLastTime

    def mock_kill():
        # Mock kill to stop the simulation
        process.isRunning.return_value = False

    process.kill.side_effect = mock_kill

    Process.wait(process, max_time=None, inactivity_timeout=BSS.Units.Time.nanosecond)

    assert process._process.wait.call_count == 10
    process.kill.assert_called_once()

    with open(f"{tmp_path}/test.out", "r") as f:
        assert f.read() == "Process Hung. Killed."
