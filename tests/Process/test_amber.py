from collections import OrderedDict

import math
import pytest
import shutil
import socket

import BioSimSpace as BSS

from tests.conftest import url, has_amber

# Store the allowed restraints.
restraints = BSS.Protocol._position_restraint_mixin._PositionRestraintMixin.restraints()


@pytest.fixture(scope="module")
def rna_system():
    """An RNA system for re-use."""
    return BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), ["rna_6e1s.rst7", "rna_6e1s.prm7"])
    )


@pytest.fixture(scope="module")
def large_protein_system():
    """A large protein system for re-use."""
    return BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), ["complex_vac0.prm7", "complex_vac0.rst7"])
    )


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_minimise(system, restraint):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100, restraint=restraint)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_equilibrate(system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_heat(system):
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(0, "kelvin"),
        temperature_end=BSS.Types.Temperature(300, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_cool(system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(300, "kelvin"),
        temperature_end=BSS.Types.Temperature(0, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_production(system, restraint):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol, check_data=True)


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_args(system):
    """Test setting an manipulation of command-line args."""

    # Create a default minimisation protocol. This doesn't matter since
    # we're going to clear the default arguments anyway.
    protocol = BSS.Protocol.Minimisation()

    # Create the process object.
    process = BSS.Process.Amber(system, protocol, name="test")

    # Clear the existing arguments.
    process.clearArgs()

    # Now add some arguments. Firstly one-by-one, using a mixture of
    # arguments and flags.
    process.setArg("-a", "A")  # Regular argument.
    process.setArg("-b", "B")  # Regular argument.
    process.setArg("-c", True)  # Boolean flag.
    process.setArg("-d", "D")  # Regular argument.
    process.setArg("-e", True)  # Boolean flag.
    process.setArg("-f", 6)  # Argument value is an integer.

    # Get the arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 6

    # Make sure the string is correct.
    assert len(arg_string_list) == 10
    assert arg_string == "-a A -b B -c -d D -e -f 6"

    # Turn off one of the flags.
    process.setArg("-c", False)

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the same number of arguments.
    assert len(args) == 6

    # Make sure the new string is correct. (The "False" flag should be missing.)
    assert len(arg_string_list) == 9
    assert arg_string == "-a A -b B -d D -e -f 6"

    # Create a new dictionary of extra arguments. This could be a regular
    # dictionary, but we use an OrderedDict for testing purposes.
    extra_args = OrderedDict()

    # Populate the arguments.
    extra_args["-g"] = True
    extra_args["-h"] = "H"
    extra_args["-i"] = False
    extra_args["-k"] = "K"

    # Add the arguments.
    process.addArgs(extra_args)

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 10

    # Make sure the new string is correct.
    assert len(arg_string_list) == 14
    assert arg_string == "-a A -b B -d D -e -f 6 -g -h H -k K"

    # Now we'll delete an argument.
    process.deleteArg("-d")

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 9

    # Make sure the new string is correct.
    assert len(arg_string_list) == 12
    assert arg_string == "-a A -b B -e -f 6 -g -h H -k K"

    # Now test insertion of additional arguments.
    process.insertArg("-x", "X", 0)  # Insert at beginning.
    process.insertArg("-y", True, 4)  # Insert at middle.
    process.insertArg("-z", "Z", 11)  # Insert at end.

    # Get the updated arguments and the string representation.
    args = process.getArgs()
    arg_string = process.getArgString()
    arg_string_list = process.getArgStringList()

    # Make sure there is the correct number of arguments.
    assert len(args) == 12

    # Make sure the new string is correct.
    assert len(arg_string_list) == 17
    assert arg_string == "-x X -a A -b B -y -e -f 6 -g -h H -k K -z Z"


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_backbone_restraint_mask_protein(large_protein_system):
    """
    Test that the amber backbone restraint mask is correct for a protein system.
    We need a large protein system otherwise the logic we want to test will be
    skipped, and individual atoms will be specified in the config.
    """

    # Create an equilibration protocol with backbone restraints.
    protocol = BSS.Protocol.Equilibration(restraint="backbone")

    # Create the process object.
    process = BSS.Process.Amber(large_protein_system, protocol, name="test")

    # Check that the correct restraint mask is in the config.
    config = process.getConfig()
    assert '   restraintmask="@N,CA,C,O",' in config


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_backbone_restraint_mask_rna(rna_system):
    """
    Test that the amber backbone restraint mask is correct for an RNA system.
    """

    # Create an equilibration protocol with backbone restraints.
    protocol = BSS.Protocol.Equilibration(restraint="backbone")

    # Create the process object.
    process = BSS.Process.Amber(rna_system, protocol, name="test")

    # Check that the correct restraint mask is in the config.
    config = process.getConfig()
    assert "   restraintmask=\"@P,C5',C3',O3',O5'\"," in config


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_perturbable_restraint(perturbable_system):
    """Test a free energy perturbation protocol."""

    # Create a short minimisation prototocol with a restraint.
    protocol = BSS.Protocol.Minimisation(steps=100, restraint="heavy")

    # Run the process, check that it finished without error, and returns a system.
    run_process(perturbable_system, protocol)


def run_process(system, protocol, check_data=False):
    """Helper function to run various simulation protocols."""

    # Initialise the AMBER process.
    process = BSS.Process.Amber(system, protocol, name="test")

    # Start the AMBER simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None

    # Make sure the correct amount of data is generated.
    if check_data:
        # Get the config from the process.
        config = process.getConfig()

        # Parse the config to get the report frequency and total number of steps.
        freq = None
        nsteps = None
        for line in config:
            if "ntpr" in line:
                freq = int(line.split("=")[1].split(",")[0])
            elif "nstlim" in line:
                nsteps = int(line.split("=")[1].split(",")[0])

        if freq and nsteps:
            # Work out the number of records. (Add one since the zero step is recorded.)
            nrec = int(nsteps / freq) + 1

            # Get the record data.
            data = process.getRecords()

            for k, v in data.items():
                assert len(v) == nrec


@pytest.mark.skipif(
    has_amber is False, reason="Requires AMBER and pyarrow to be installed."
)
@pytest.mark.parametrize(
    "protocol",
    [
        BSS.Protocol.FreeEnergy(temperature=298 * BSS.Units.Temperature.kelvin),
        BSS.Protocol.FreeEnergyMinimisation(),
    ],
)
def test_parse_fep_output(perturbable_system, protocol):
    """Make sure that we can correctly parse AMBER FEP output."""

    from sire.legacy.Base import findExe

    # Copy the system.
    system_copy = perturbable_system.copy()

    # Use the first instance of sander in the path so that we can
    # test without pmemd.
    exe = findExe("sander").absoluteFilePath()
    process = BSS.Process.Amber(system_copy, protocol, exe=exe)

    # Assign the path to the output file.
    if isinstance(protocol, BSS.Protocol.FreeEnergy):
        out_file = "tests/output/amber_fep.out"
    else:
        out_file = "tests/output/amber_fep_min.out"

    # Copy the existing output file into the working directory.
    shutil.copyfile(out_file, process.workDir() + "/amber.out")

    # Update the stdout record dictionaries.
    process.stdout(0)

    # Get back the records for each region and soft-core part.
    records_ti0 = process.getRecords(region=0)
    records_sc0 = process.getRecords(region=0, soft_core=True)
    records_ti1 = process.getRecords(region=1)
    records_sc1 = process.getRecords(region=1, soft_core=True)

    # Make sure NSTEP is present.
    assert "NSTEP" in records_ti0

    # Get the number of records.
    num_records = len(records_ti0["NSTEP"])

    # Now make sure that the records for the two TI regions contain the
    # same number of values.
    for v0, v1 in zip(records_ti0.values(), records_ti1.values()):
        assert len(v0) == len(v1) == num_records

    # Now check that are records for the soft-core parts contain the correct
    # number of values.
    for v in records_sc0.values():
        assert len(v) == num_records
    for k, v in records_sc1.items():
        assert len(v) == num_records
    if isinstance(protocol, BSS.Protocol.FreeEnergy):
        assert len(records_sc0) == len(records_sc1)
    else:
        assert len(records_sc0) == 0
        assert len(records_sc1) != 0


@pytest.mark.skipif(
    socket.gethostname() != "porridge",
    reason="Local test requiring pmemd installation.",
)
def test_pmemd(system):
    """Single-point energy tests for pmemd."""

    # Path to the pmemd conda environment bin directory.
    bin_dir = "/home/lester/.conda/envs/pmemd/bin"

    # Single-point minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=1)

    # First perform single-point comparisons in solvent.

    # Compute the single-point energy using sander.
    process = BSS.Process.Amber(system, protocol)
    process.start()
    process.wait()
    assert not process.isError()
    nrg_sander = process.getTotalEnergy().value()

    # Compute the single-point energy using pmemd.
    process = BSS.Process.Amber(system, protocol, exe=f"{bin_dir}/pmemd")
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd = process.getTotalEnergy().value()

    # Compute the single-point energy using pmemd.cuda.
    process = BSS.Process.Amber(system, protocol, exe=f"{bin_dir}/pmemd.cuda")
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd_cuda = process.getTotalEnergy().value()

    # Check that the energies are the same.
    assert math.isclose(nrg_sander, nrg_pmemd, rel_tol=1e-4)
    assert math.isclose(nrg_sander, nrg_pmemd_cuda, rel_tol=1e-4)

    # Now perform single-point comparisons in vacuum.

    vac_system = system[0].toSystem()

    # Compute the single-point energy using sander.
    process = BSS.Process.Amber(vac_system, protocol)
    process.start()
    process.wait()
    assert not process.isError()
    nrg_sander = process.getTotalEnergy().value()

    # Compute the single-point energy using pmemd.
    process = BSS.Process.Amber(vac_system, protocol, exe=f"{bin_dir}/pmemd")
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd = process.getTotalEnergy().value()

    # Compute the single-point energy using pmemd.cuda.
    process = BSS.Process.Amber(vac_system, protocol, exe=f"{bin_dir}/pmemd.cuda")
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd_cuda = process.getTotalEnergy().value()

    # Check that the energies are the same.
    assert math.isclose(nrg_sander, nrg_pmemd, rel_tol=1e-4)
    assert math.isclose(nrg_sander, nrg_pmemd_cuda, rel_tol=1e-4)


@pytest.mark.skipif(
    socket.gethostname() != "porridge",
    reason="Local test requiring pmemd installation.",
)
def test_pmemd_fep(solvated_perturbable_system):
    """Single-point FEP energy tests for pmemd."""

    # Path to the pmemd conda environment bin directory.
    bin_dir = "/home/lester/.conda/envs/pmemd/bin"

    # Single-point minimisation protocol.
    protocol = BSS.Protocol.FreeEnergyMinimisation(steps=1)

    # First perform single-point comparisons in solvent.

    # Compute the single-point energy using pmemd.
    process = BSS.Process.Amber(
        solvated_perturbable_system, protocol, exe=f"{bin_dir}/pmemd"
    )
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd = process.getTotalEnergy().value()

    # Compute the single-point energy using pmemd.cuda.
    process = BSS.Process.Amber(
        solvated_perturbable_system, protocol, exe=f"{bin_dir}/pmemd.cuda"
    )
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd_cuda = process.getTotalEnergy().value()

    # Check that the energies are the same.
    assert math.isclose(nrg_pmemd, nrg_pmemd_cuda, rel_tol=1e-4)

    # Now perform single-point comparisons in vacuum.

    vac_system = solvated_perturbable_system[0].toSystem()

    # Compute the single-point energy using pmemd.
    process = BSS.Process.Amber(
        vac_system, protocol, exe=f"{bin_dir}/pmemd", extra_options={"gti_bat_sc": 2}
    )
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd = process.getTotalEnergy().value()

    # Compute the single-point energy using pmemd.cuda.
    process = BSS.Process.Amber(
        vac_system,
        protocol,
        exe=f"{bin_dir}/pmemd.cuda",
        extra_options={"gti_bat_sc": 2},
    )
    process.start()
    process.wait()
    assert not process.isError()
    nrg_pmemd_cuda = process.getTotalEnergy().value()

    # Check that the energies are the same.
    assert math.isclose(nrg_pmemd, nrg_pmemd_cuda, rel_tol=1e-3)
