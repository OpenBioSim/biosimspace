import math
import pytest
import requests
import tarfile
import tempfile
import json

import BioSimSpace as BSS


def test_makeSystem(TEMOA_host, TEMOA_lig1, TEMOA_lig2):

    atm_generator = BSS.FreeEnergy.AToM(
        protein=TEMOA_host, ligand1=TEMOA_lig1, ligand2=TEMOA_lig2
    )
    # check that an error is thrown in the rigid core atoms are not given to prepare
    with pytest.raises(TypeError):
        atm_system, atm_data = atm_generator.prepare()

    rigid_core = [1, 2, 3]

    atm_system, atm_data = atm_generator.prepare(
        ligand1_rigid_core=rigid_core, ligand2_rigid_core=rigid_core
    )

    # Check that the system contains an atm data property
    data_from_system = json.loads(atm_system._sire_object.property("atom_data").value())
    to_ignore = ["displacement"]
    # check that atm_data and data_from_system are the same, ignoring anything in to_ignore
    assert all(
        [
            data_from_system[key] == atm_data[key]
            for key in atm_data
            if key not in to_ignore
        ]
    )

    # check that data[ligand1_rigid_core] and data[ligand2_rigid_core] are the same as the input
    assert data_from_system["ligand1_rigid_core"] == rigid_core
    assert data_from_system["ligand2_rigid_core"] == rigid_core

    # get the coordinates of the ligands
    lig1_coords = atm_system[atm_data["ligand1_index"]]._sire_object.coordinates()
    lig2_coords = atm_system[atm_data["ligand2_index"]]._sire_object.coordinates()
    # make sure the displacement is correct for the default value of 20A
    assert pytest.approx((lig2_coords - lig1_coords).length().value(), rel=1) == 20.0

    vector = BSS.Types.Vector(10.0, 10.0, 10.0)

    system_withvec, data_withvec = atm_generator.prepare(
        ligand1_rigid_core=rigid_core,
        ligand2_rigid_core=rigid_core,
        displacement=vector,
    )

    data_from_system = json.loads(
        system_withvec._sire_object.property("atom_data").value()
    )
    assert pytest.approx(data_from_system["displacement"], rel=1e-3) == [
        vector.x(),
        vector.y(),
        vector.z(),
    ]
    lig1_coords = system_withvec[
        data_withvec["ligand1_index"]
    ]._sire_object.coordinates()
    lig2_coords = system_withvec[
        data_withvec["ligand2_index"]
    ]._sire_object.coordinates()

    d = lig2_coords - lig1_coords
    assert pytest.approx(d.x().value(), 1) == vector.x()
    assert pytest.approx(d.y().value(), 1) == vector.y()
    assert pytest.approx(d.z().value(), 1) == vector.z()

    # make a new atm_generator and check the parsing of a full system
    atm_generator = BSS.FreeEnergy.AToM(system=atm_system)
