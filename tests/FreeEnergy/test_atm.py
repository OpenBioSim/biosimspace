import math
import pytest
import requests
import tarfile
import tempfile
import json
import pandas as pd
import os

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


def test_run(TEMOA_hostguest):
    system, _ = TEMOA_hostguest
    production_atm = BSS.Protocol.AToMProduction(
        system=system,
        com_distance_restraint=True,
        runtime="2 fs",
        report_interval=1,
        restart_interval=1,
        num_lambda=2,
        analysis_method="UWHAM",
    )
    production_atm2 = BSS.Protocol.AToMProduction(
        system=system,
        com_distance_restraint=True,
        runtime="4 fs",
        report_interval=1,
        restart_interval=1,
        num_lambda=2,
        analysis_method="UWHAM",
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        production = BSS.FreeEnergy.AToM.run(
            system, production_atm, work_dir=tmpdirname
        )
        production.wait()
        # read openmm.csv and make sure it has a single row
        df = pd.read_csv(os.path.join(tmpdirname, "lambda_0.0000/openmm.csv"))
        assert len(df) == 1

        production2 = BSS.FreeEnergy.AToM.run(
            system, production_atm2, work_dir=tmpdirname
        )
        production2.wait()
        df = pd.read_csv(os.path.join(tmpdirname, "lambda_0.0000/openmm.csv"))
        assert len(df) == 2


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
        com_distance_restraint=True,
        com_k=25.0,
        com_restraint_width=5.0,
        restraint=pos_rst_atoms,
        positional_restraint_width=0.5,
        force_constant=25.0,
        align_k_psi=10.0,
        align_k_theta=10.0,
        align_k_distance=2.5,
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
            platform="CPU",
            setup_only=True,
            work_dir=tmpdirname,
            **{"_is_testing": True},
        )
        production.start()
        production.wait()

        assert not production.isError()
        # now get the file containing single points
        df = pd.read_csv(os.path.join(tmpdirname, "energies_singlepoint.csv"))
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
        df_nonlam = pd.read_csv(os.path.join(tmpdirname, "non_lambda_forces.csv"))
        ens_nonlam = df_nonlam.to_dict()

        ens_GL_nolam = {
            "com": 0.0,
            "distance": 0.02983,
            "angle": 0.0072010,
            "dihedral": 3.55355e-13,
            "position_restraint": 0.0,
        }
        for key, en in ens_GL_nolam.items():
            assert pytest.approx(ens_nonlam[key][0], rel=1e-3) == en


def test_UWHAM():
    import numpy as np

    # To try and ensure parity with the Gallachio lab code
    # we will test each individual element of the UWHAM calculation
    potential = -69702.79
    e_pert = 70.48908
    beta = 1.678238963
    lambda1 = 0.0
    lambda2 = 0.0
    alpha = 0.0
    u0 = 0.0
    w0 = 0.0

    n_pot = 116977.9

    from BioSimSpace.FreeEnergy._ddg import _npot_fcn

    npot = _npot_fcn(
        e0=potential,
        epert=e_pert,
        bet=beta,
        lam1=lambda1,
        lam2=lambda2,
        alpha=alpha,
        u0=u0,
        w0=w0,
    )

    assert pytest.approx(npot, rel=1e-3) == n_pot

    # Now testing agreement with known values from
    # UWHAM-R analysis
    ln_q_array = np.array(
        [
            [
                117251.85638785,
                117147.70578372,
                117259.71235395,
                117184.35793014,
                116934.45115,
                117405.64541825,
                116930.39936544,
                117131.36660758,
                117072.35871073,
                117041.11910054,
                117166.97160247,
            ],
            [
                117246.36847836,
                117141.83269751,
                117254.16656958,
                117181.69670683,
                116933.05433974,
                117404.90353356,
                116930.17474063,
                117130.18191686,
                117071.63726339,
                117041.08918386,
                117167.13191866,
            ],
            [
                117240.88056888,
                117135.95961129,
                117248.62078521,
                117179.03548353,
                116931.65752948,
                117404.16164887,
                116929.95011581,
                117128.99722614,
                117070.91581606,
                117041.05926718,
                117167.29223485,
            ],
            [
                117235.3926594,
                117130.08652508,
                117243.07500085,
                117176.37426023,
                116930.26071922,
                117403.41976418,
                116929.725491,
                117127.81253542,
                117070.19436872,
                117041.0293505,
                117167.45255104,
            ],
            [
                117229.90474991,
                117124.21343886,
                117237.52921648,
                117173.71303693,
                116928.86390896,
                117402.67787949,
                116929.50086619,
                117126.6278447,
                117069.47292138,
                117040.99943382,
                117167.61286723,
            ],
            [
                117224.41684043,
                117118.34035265,
                117231.98343212,
                117171.05181363,
                116927.4670987,
                117401.9359948,
                116929.27624138,
                117125.44315397,
                117068.75147404,
                117040.96951714,
                117167.77318343,
            ],
            [
                117218.92893094,
                117112.46726644,
                117226.43764775,
                117168.39059033,
                116926.07028844,
                117401.1941101,
                116929.05161657,
                117124.25846325,
                117068.0300267,
                117040.93960046,
                117167.93349962,
            ],
            [
                117213.44102146,
                117106.59418022,
                117220.89186339,
                117165.72936702,
                116924.67347817,
                117400.45222541,
                116928.82699176,
                117123.07377253,
                117067.30857936,
                117040.90968378,
                117168.09381581,
            ],
            [
                117207.95311198,
                117100.72109401,
                117215.34607902,
                117163.06814372,
                116923.27666791,
                117399.71034072,
                116928.60236695,
                117121.88908181,
                117066.58713202,
                117040.8797671,
                117168.254132,
            ],
            [
                117202.46520249,
                117094.84800779,
                117209.80029465,
                117160.40692042,
                116921.87985765,
                117398.96845603,
                116928.37774214,
                117120.70439109,
                117065.86568468,
                117040.84985042,
                117168.41444819,
            ],
            [
                117196.97729301,
                117088.97492158,
                117204.25451029,
                117157.74569712,
                116920.48304739,
                117398.22657134,
                116928.15311733,
                117119.51970037,
                117065.14423734,
                117040.81993374,
                117168.57476438,
            ],
        ]
    )

    n_samples = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    known_answer = 12.42262
    known_error = 1.422156

    from BioSimSpace.FreeEnergy._ddg import _estimate_f_i

    f_i, d_i, weights = _estimate_f_i(ln_q_array, n_samples)
    ddg = f_i[-1] - f_i[0]
    ddg = ddg / beta
    ddg_error = np.sqrt(d_i[-1] + d_i[0]) / beta

    assert pytest.approx(ddg, rel=1e-3) == known_answer
    assert pytest.approx(ddg_error, rel=1e-3) == known_error
