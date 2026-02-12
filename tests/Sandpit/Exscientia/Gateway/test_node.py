import subprocess
import sys

import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from tests.conftest import root_fp

# Store the name of the test script.
script_name = f"{root_fp}/Gateway/node.py"

# Store the name of the python interpreter.
exe = sys.executable

# Create a list of valid command-line arguments.
args = [
    "--bool=True",
    "--int=42",
    "--float=3.14",
    '--string="hello"',
    "--file=tests/input/ala.crd",
    "--fileset=tests/input/ala.crd,tests/input/ala.top",
    '--temperature="298 kelvin"',
    '--time="100 nanoseconds"',
    '--length="10 angstroms"',
    '--area="256 nanometers**2"',
    '--volume="1024 picometers cubed"',
    '--angle="3.14 radians"',
    '--charge="-1 electron charge"',
    '--energy="-1000 kcal/mol"',
    '--pressure="1 atmosphere"',
]


def test_correct_args():
    """Test that correct input from the command line is validated."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


def test_optional_args():
    """Test that the node works when optional arguments aren't set."""

    # Copy the args list.
    local_args = args

    # Delete the optional boolean argument.
    args.remove("--bool=True")

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(local_args)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


# Parameterise with all non-optional arguments.
@pytest.mark.parametrize("value", args[1:])
def test_missing_args(value):
    """Test that the node fails when required arguments are missing."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Strip the argument from the string.
    command = command.replace(value, "")

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command failed.
    assert proc.returncode != 0


@pytest.mark.skip(reason="Intermittently failing during CI only.")
def test_invalid_args():
    """Test that the node fails when required arguments are of the wrong type."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Replace the argument types.
    command = command.replace("42", "hello")
    command = command.replace("3.14", "world")

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command failed.
    assert proc.returncode != 0


def test_multi_args():
    """Test that multi valued arguments can be handled correctly."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0

    # Now strip one of the files from the file-set string.
    command = command.replace(" tests/input/ala.top", "")

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


def test_missing_file():
    """Test that the node fails when a file requirement doesn't exist."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Replace with a missing file.
    invalid_command = command.replace("file=tests/input/ala.crd", "file=missing.txt")

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(invalid_command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command failed.
    assert proc.returncode != 0

    # Replace with a missing file.
    invalid_command = command.replace("tests/input/ala.top", "file=missing.txt")

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(invalid_command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command failed.
    assert proc.returncode != 0


@pytest.mark.parametrize(
    "value",
    [
        "100 K",
        "23 k",
        "1.02e3 K",
        "45 KelVIn",
        " -50C",
        "-1.6e2 c",
        "  1024celsius",
        "22 f",
        "-200 fahrenheit",
    ],
)
def test_temperature(value):
    """Test that different format temperature strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("298 kelvin", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "1.5 hours",
        "0.3hr",
        "500 NanoSECONd",
        "   56ps",
        "15 mins",
        "6.7e3 milliseC onds",
        "1000 ms",
        "2 Fs",
        "3.6 DaYs",
    ],
)
def test_time(value):
    """Test that different format time strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("100 nanoseconds", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "3 meters",
        "1e-9 mEte R",
        "        31 c M",
        " 0.2e3 naNometer",
        "15a",
        "33 Angstrom",
        "  3.6e4 pm",
        "22 milli   MEters",
        "    73 a",
    ],
)
def test_length(value):
    """Test that different format length strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("10 angstroms", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "3 meters squared",
        "1e-9 mEte R ^2",
        "        31 n M2",
        " 0.2e3 naNometer square",
        "15a^2",
        "33 Angstroms squaRed",
        "  3.6e4 pm2",
        "22 MEters      ^2",
        "    73 a2",
    ],
)
def test_area(value):
    """Test that different format area strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("256 nanometers squared", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "3 meters cubed",
        "1e-9 mEte R ^3",
        "        31 n M3",
        " 0.2e3 naNometer cube",
        "15a^3",
        "33 Angstroms cuBed",
        "  3.6e4 pm3",
        "22 MEters      ^3",
        "    73 a3",
    ],
)
def test_volume(value):
    """Test that different format volume strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("1024 picometers cubed", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "1.5 Ele cTron C HArge",
        "-  3.6 cou Lomb S",
        "        -2 es",
        " 0.2es",
        "-1.7e3 e charge",
        "   2.7e-4 coU lS",
        " -1003 Cs",
        "   15e4 eleC",
        " -2.3e-2 E  lEc T",
    ],
)
def test_charge(value):
    """Test that different format charge strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("-1 electron charge", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "1.5 R a DI an S",
        "-  3.6 degr EES",
        "        -2 r",
        " 0.2deG",
        "-1.7e3 dE g",
        "   2.7e-4 Ra DS",
        " -1003 Ds",
        "   15e4 RAD",
        " -2.3e-2 r  ad IAN",
    ],
)
def test_angle(value):
    """Test that different format angle strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("3.14 radians", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        "3 kilo CAlo  ries per   MOle",
        "1e-9 kcal   /mol",
        "        31 KJ/m ole",
        " 0.2e3 KT",
        "15Kilo Joules /mol",
        "33 Kcalories / mol",
        "  3.6e4 kjOul es /  mole",
        "22       K    t",
        "    73 k  caLories / mol",
    ],
)
def test_energy(value):
    """Test that different format energy strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("1000 kcal/mol", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0


@pytest.mark.parametrize(
    "value",
    [
        " 3.6 Atm oS pheres",
        "1e4 baR",
        "        1 atmospheric Pres",
        " 1e-3 b   Ars",
        " 25e4 atm pressure",
    ],
)
def test_pressure(value):
    """Test that different format pressure strings are supported."""

    # Generate the shell command.
    command = "%s %s " % (exe, script_name) + " ".join(args)

    # Modify the string.
    command = command.replace("1 atmosphere", value)

    # Run the command.
    proc = subprocess.run(
        BSS._Utils.command_split(command), shell=False, stdout=subprocess.PIPE
    )

    # Make sure the command completed successfully.
    assert proc.returncode == 0
