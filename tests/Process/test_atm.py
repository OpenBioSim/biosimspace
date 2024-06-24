import pytest

import BioSimSpace as BSS


def test_ATM_min(TEMOA_hostguest):
    # First get a system with data
    system, data = TEMOA_hostguest
    # Generate a minimisation protocol
    prot_min = BSS.Protocol.AToMMinimisation(data=data, steps=100)

    run_process(system, prot_min)
    del system, data


def test_ATM_equil_noATMForce(TEMOA_hostguest):
    # First get a system with data
    system, data = TEMOA_hostguest
    # Generate an equilibration protocol
    prot_equil = BSS.Protocol.AToMEquilibration(
        data=data, runtime="0.2 ps", use_atm_force=False
    )

    run_process(system, prot_equil)
    del system, data


def test_ATM_equil_withATMForce(TEMOA_hostguest):
    # First get a system with data
    system, data = TEMOA_hostguest
    # Generate an equilibration protocol
    prot_equil = BSS.Protocol.AToMEquilibration(
        data=data, runtime="0.2 ps", use_atm_force=True
    )

    run_process(system, prot_equil)
    del system, data


def test_ATM_anneal(TEMOA_hostguest):
    # First get a system with data
    system, data = TEMOA_hostguest
    # Generate an annealing protocol
    prot_anneal = BSS.Protocol.AToMAnnealing(data=data, runtime="0.2 ps")

    run_process(system, prot_anneal)
    del system, data


def test_custom_ATM_anneal(TEMOA_hostguest):
    # First get a system with data
    system, data = TEMOA_hostguest
    # now test passing a valid dictionary
    annealing_dict = {
        "lambda1_start": 0.0,
        "lambda1_end": 0.5,
        "lambda2_start": 0.0,
        "lambda2_end": 0.5,
        "alpha_start": 0.0,
        "alpha_end": 0.5,
        "uh_start": 0.0,
        "uh_end": 0.5,
        "W0_start": 0.0,
        "W0_end": 0.5,
    }
    protocol = BSS.Protocol.AToMAnnealing(
        data=data,
        anneal_values=annealing_dict,
        anneal_numcycles=10,
        runtime="0.2 ps",
    )
    run_process(system, protocol)


def test_ATM_production(TEMOA_hostguest):
    # First get a system with data
    system, data = TEMOA_hostguest
    # Generate a production protocol
    prot_prod = BSS.Protocol.AToMProduction(data=data, runtime="0.2 ps")

    run_process(system, prot_prod)

    # now test "MBAR" analysis method
    prot_prod = BSS.Protocol.AToMProduction(
        data=data, runtime="0.2 ps", analysis_method="MBAR"
    )
    run_process(system, prot_prod)

    # finally, test the "both" analysis method
    prot_prod = BSS.Protocol.AToMProduction(
        data=data, runtime="0.2 ps", analysis_method="both"
    )
    run_process(system, prot_prod)


def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the OpenMM  process.
    process = BSS.Process.OpenMM(system, protocol, name="test")

    # Start the OpenMM simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None
