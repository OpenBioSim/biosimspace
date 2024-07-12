import math
import pytest
import requests
import tarfile
import tempfile
import json
import pandas as pd

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


def test_single_point_energies(TEMOA_host, TEMOA_lig1, TEMOA_lig2):
    # Tests the single point energies of the
    # Mirroring inputs for G. lab code
    lig1_cm_atoms_absolute = [
        196,
        197,
        198,
        199,
        200,
        201,
        202,
        203,
        204,
        205,
        206,
        207,
        208,
        209,
        210,
        211,
        212,
        213,
        214,
        215,
        216,
    ]
    lig2_cm_atoms_absolute = [
        217,
        218,
        219,
        220,
        221,
        222,
        223,
        224,
        225,
        226,
        227,
        228,
        229,
        230,
        231,
        232,
        233,
    ]
    lig1_cm_rel = [x - 196 for x in lig1_cm_atoms_absolute]
    lig2_cm_rel = [x - 217 for x in lig2_cm_atoms_absolute]
    prot_cm_atoms = [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
        51,
        52,
        53,
        54,
        55,
        56,
        57,
        58,
        59,
        60,
        61,
        62,
        63,
        64,
        65,
        66,
        67,
        68,
        69,
        70,
        71,
        72,
        73,
        74,
        75,
        76,
        77,
        78,
        79,
        80,
        81,
        82,
        83,
        84,
        85,
        86,
        87,
        88,
        89,
        90,
        91,
        92,
        93,
        94,
        95,
        96,
        97,
        98,
        99,
        100,
        101,
        102,
        103,
        104,
        105,
        106,
        107,
        108,
        109,
        110,
        111,
        112,
        113,
        114,
        115,
        116,
        117,
        118,
        119,
        120,
        121,
        122,
        123,
        124,
        125,
        126,
        127,
        128,
        129,
        130,
        131,
        132,
        133,
        134,
        135,
        136,
        137,
        138,
        139,
        140,
        141,
        142,
        143,
        144,
        145,
        146,
        147,
        148,
        149,
        150,
        151,
        152,
        153,
        154,
        155,
        156,
        157,
        158,
        159,
        160,
        161,
        162,
        163,
        164,
        165,
        166,
        167,
        168,
        169,
        170,
        171,
        172,
        173,
        174,
        175,
        176,
        177,
        178,
        179,
        180,
        181,
        182,
        183,
        184,
        185,
        186,
        187,
        188,
        189,
        190,
        191,
        192,
        193,
        194,
        195,
    ]
    atm_generator = BSS.FreeEnergy.AToM(
        protein=TEMOA_host, ligand1=TEMOA_lig1, ligand2=TEMOA_lig2
    )
    system, data = atm_generator.prepare(
        displacement=[22, 22, 22],
        ligand1_rigid_core=[8, 6, 4],
        ligand2_rigid_core=[3, 5, 1],
        ligand1_com_atoms=lig1_cm_rel,
        ligand2_com_atoms=lig2_cm_rel,
        protein_com_atoms=prot_cm_atoms,
    )

    pos_rst_atoms = [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        34,
        35,
        36,
        37,
        38,
        39,
    ]

    production_atm = BSS.Protocol.AToMProduction(
        system=system,
        CMCM_restraint=True,
        cm_kf=25.0,
        cm_tol=5.0,
        restraint=pos_rst_atoms,
        pos_rest_width=0.5,
        force_constant=25.0,
        align_k_psi=10.0,
        align_k_theta=10.0,
        align_kf_sep=2.5,
        runtime="100 ps",
        num_lambda=22,
        SC_umax=100.0,
        SC_a=0.0625,
        SC_u0=50.0,
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        production = BSS.Process.OpenMM(
            system,
            production_atm,
            platform="CUDA",
            setup_only=True,
            work_dir=tmpdirname,
            **{"_is_testing": True},
        )
        production.start()
        production.wait()
        # now get the file containing single points
        df = pd.read_csv(tmpdirname + "/energies_singlepoint.csv")
        ens = df.to_dict()

        # Here we are specifically verifying the energies of the ATMForce
        ens_GL = {
            0.0: 2847.6,
            0.0476: 2849.1,
            0.0952: 2850.7,
            0.1429: 2852.2,
            0.1905: 2853.7,
            0.2381: 2855.3,
            0.2857: 2856.8,
            0.3333: 2858.4,
            0.381: 2859.9,
            0.4286: 2861.5,
            0.4762: 2863.0,
        }
        # Need to add an offset due to treatment of 1-4 forces in GL code
        offset = 803.3
        # now check that the energies are the same
        for lam, en in ens_GL.items():
            assert pytest.approx(ens[str(lam)][0], rel=1) == en + offset

        # Now check the rest of the forces
        df_nonlam = pd.read_csv(tmpdirname + "/non_lambda_forces.csv")
        ens_nonlam = df_nonlam.to_dict()

        ens_GL_nolam = {
            "com": 0.0,
            "distance": 0.02983,
            "angle": 0.0072010,
            "dihedral": 6.234e-07,
            "position_restraint": 0.0,
        }
        print(ens_nonlam)
        for key, en in ens_GL_nolam.items():
            assert pytest.approx(ens_nonlam[key][0], rel=1e-3) == en
