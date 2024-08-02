import pytest
import BioSimSpace as BSS


def test_atm_minimisation(TEMOA_hostguest):
    # We will use this as a test for all of the parent class inputs

    # First need to test that both forms of the `data` input work
    system, data = TEMOA_hostguest
    BSS.Protocol.AToMMinimisation(data=data)
    BSS.Protocol.AToMMinimisation(system=system)
    # Now test the optional inputs, first using biosimspace Units
    protocol_units = BSS.Protocol.AToMMinimisation(
        data=data,
        core_alignment=False,
        CMCM_restraint=False,
        restraint="all",
        force_constant=1.0 * (BSS.Units.Energy.kcal_per_mol / BSS.Units.Area.angstrom2),
        pos_rest_width=0.1 * BSS.Units.Length.angstrom,
        align_kf_sep=1.0 * BSS.Units.Energy.kcal_per_mol / BSS.Units.Area.angstrom2,
        align_k_theta=1.0 * BSS.Units.Energy.kcal_per_mol,
        align_k_psi=1.0 * BSS.Units.Energy.kcal_per_mol,
        SC_umax=10.0 * BSS.Units.Energy.kcal_per_mol,
        SC_u0=1.0 * BSS.Units.Energy.kcal_per_mol,
        SC_a=0.01,
        cm_kf=1.0 * BSS.Units.Energy.kcal_per_mol / BSS.Units.Area.angstrom2,
        cm_tol=1.0 * BSS.Units.Length.angstrom,
    )

    # Now test parsing options as floats
    protocol_floats = BSS.Protocol.AToMMinimisation(
        data=data,
        force_constant=1.0,
        pos_rest_width=0.1,
        align_kf_sep=1.0,
        align_k_theta=1.0,
        align_k_psi=1.0,
        SC_umax=10.0,
        SC_u0=1.0,
        SC_a=0.01,
        cm_kf=1.0,
        cm_tol=1.0,
    )

    # Finally try parsing strings
    protocol_strings = BSS.Protocol.AToMMinimisation(
        data=data,
        force_constant="1.0 kcal mol^-1 angstrom^-2",
        pos_rest_width="0.1 angstrom",
        align_kf_sep="1.0 kcal mol^-1 angstrom^-2",
        align_k_theta="1.0 kcal mol^-1",
        align_k_psi="1.0 kcal mol^-1",
        SC_umax="10.0 kcal mol^-1",
        SC_u0="1.0 kcal mol^-1",
        SC_a=0.01,
        cm_kf="1.0 kcal mol^-1 angstrom^-2",
        cm_tol="1.0 angstrom",
    )
    # using getters, check that all protocols have the same values
    # (skip force constant and as it is not atm exclusive)
    assert (
        protocol_units._getPosRestWidth()
        == protocol_floats._getPosRestWidth()
        == protocol_strings._getPosRestWidth()
    )
    assert (
        protocol_units._getAlignKfSep()
        == protocol_floats._getAlignKfSep()
        == protocol_strings._getAlignKfSep()
    )
    assert (
        protocol_units._getAlignKTheta()
        == protocol_floats._getAlignKTheta()
        == protocol_strings._getAlignKTheta()
    )
    assert (
        protocol_units._getAlignKPsi()
        == protocol_floats._getAlignKPsi()
        == protocol_strings._getAlignKPsi()
    )
    assert (
        protocol_units._getSCUmax()
        == protocol_floats._getSCUmax()
        == protocol_strings._getSCUmax()
    )
    assert (
        protocol_units._getSCU0()
        == protocol_floats._getSCU0()
        == protocol_strings._getSCU0()
    )
    assert (
        protocol_units._getSCa()
        == protocol_floats._getSCa()
        == protocol_strings._getSCa()
    )
    assert (
        protocol_units._getCMKf()
        == protocol_floats._getCMKf()
        == protocol_strings._getCMKf()
    )
    assert (
        protocol_units._getCMTol()
        == protocol_floats._getCMTol()
        == protocol_strings._getCMTol()
    )


def test_atm_equilibration(TEMOA_hostguest):
    # Testing equilibration-specific inputs
    system, data = TEMOA_hostguest

    protocol_units = BSS.Protocol.AToMEquilibration(
        data=data,
        timestep=1 * BSS.Units.Time.femtosecond,
        runtime=0.1 * BSS.Units.Time.nanosecond,
        temperature_start=200 * BSS.Units.Temperature.kelvin,
        temperature_end=300 * BSS.Units.Temperature.kelvin,
        pressure=0.99 * BSS.Units.Pressure.atm,
        thermostat_time_constant=1.5 * BSS.Units.Time.picosecond,
        report_interval=1000,
        restart_interval=1001,
        use_atm_force=True,
        direction=-1,
        lambda1=0.1,
        lambda2=0.2,
        alpha=0.1 * BSS.Units.Energy.kcal_per_mol,
        uh=0.1 * BSS.Units.Energy.kcal_per_mol,
        W0=0.1 * BSS.Units.Energy.kcal_per_mol,
    )

    # test setting alpha,uh and w0 as floats
    protocol_floats = BSS.Protocol.AToMEquilibration(
        data=data,
        alpha=0.1,
        uh=0.1,
        W0=0.1,
    )

    # test setting alpha,uh and w0 as strings
    protocol_strings = BSS.Protocol.AToMEquilibration(
        data=data,
        timestep="1 fs",
        runtime="0.1 ns",
        temperature_start="200 K",
        temperature_end="300 K",
        pressure="0.99 atm",
        thermostat_time_constant="1.5 ps",
        alpha="0.1 kcal mol^-1",
        uh="0.1 kcal mol^-1",
        W0="0.1 kcal mol^-1",
    )

    # Check that all protocols have the same values
    assert protocol_units.getTimeStep() == protocol_strings.getTimeStep()
    assert protocol_units.getRunTime() == protocol_strings.getRunTime()
    assert (
        protocol_units._getStartTemperature() == protocol_strings._getStartTemperature()
    )
    assert protocol_units._getEndTemperature() == protocol_strings._getEndTemperature()
    assert protocol_units._getPressure() == protocol_strings._getPressure()
    assert (
        protocol_units._getThermostatTimeConstant()
        == protocol_strings._getThermostatTimeConstant()
    )
    assert (
        protocol_units._getAlpha()
        == protocol_floats._getAlpha()
        == protocol_strings._getAlpha()
    )
    assert (
        protocol_units._getUh() == protocol_floats._getUh() == protocol_strings._getUh()
    )
    assert (
        protocol_units._getW0() == protocol_floats._getW0() == protocol_strings._getW0()
    )


def test_atm_annealing(TEMOA_hostguest):
    # Testing annealing-specific inputs
    system, data = TEMOA_hostguest

    # first test passing an invalid key in the annealing dictionary
    annealing_dict = {
        "lambda1_start": 0.0,
        "lambda1_end": 0.5,
        "lambda2_start": 0.0,
        "lambda2_end": 0.5,
        "invalid_key": 0.0,
    }
    with pytest.raises(ValueError):
        BSS.Protocol.AToMAnnealing(data=data, anneal_values=annealing_dict)

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
        alpha="0.1 kcal mol^-1",
        uh="0.1 kcal mol^-1",
        W0="0.1 kcal mol^-1",
    )


def test_atm_production(TEMOA_hostguest):
    # Testing production-specific inputs
    system, data = TEMOA_hostguest

    # fist create a production protocol with num_lambda=6
    protocol = BSS.Protocol.AToMProduction(
        data=data,
        num_lambda=6,
    )
    # get values for direction, lambda1, lambda2, alpha, uh and W0
    assert len(protocol._getDirection()) == 6
    assert len(protocol._getLambda1()) == 6
    assert len(protocol._getLambda2()) == 6
    assert len(protocol._getAlpha()) == 6
    assert len(protocol._getUh()) == 6
    assert len(protocol._getW0()) == 6

    # Define custom values for direction that are not valid and check that an error is raised
    d = [1, 2, 3, 4, 5, 6]
    with pytest.raises(ValueError):
        protocol._setDirection(d)

    # Define custom values for lambda1 that are not valid and check that an error is raised
    l1 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    with pytest.raises(ValueError):
        protocol._setLambda1(l1)

    # Define custom values for lambda2 that are not valid and check that an error is raised
    l2 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    with pytest.raises(ValueError):
        protocol._setLambda2(l2)

    # check that a list of strings, ints, floats and BSS units can be parsed for alpha,uh and w0
    list_of_units = [
        "0.1 kcal mol^-1",
        0.1,
        0.1 * BSS.Units.Energy.kcal_per_mol,
        0.1 * BSS.Units.Energy.kcal_per_mol,
        0.1,
        1,
    ]

    protocol = BSS.Protocol.AToMProduction(
        data=data,
        num_lambda=6,
        alpha=list_of_units,
        uh=list_of_units,
        W0=list_of_units,
    )

    end_product = [0.1 * BSS.Units.Energy.kcal_per_mol] * 5
    end_product.append(1 * BSS.Units.Energy.kcal_per_mol)
    assert protocol._getAlpha() == end_product
    assert protocol._getUh() == end_product
    assert protocol._getW0() == end_product

    # now check that all of the allowed analysis options can be set
    protocol = BSS.Protocol.AToMProduction(
        data=data,
        num_lambda=6,
        analysis_method="UWHAM",
    )

    protocol = BSS.Protocol.AToMProduction(
        data=data,
        num_lambda=6,
        analysis_method="MBAR",
    )

    protocol = BSS.Protocol.AToMProduction(
        data=data,
        num_lambda=6,
        analysis_method="both",
    )