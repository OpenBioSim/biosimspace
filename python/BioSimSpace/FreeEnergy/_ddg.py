######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>, Matthew Burman <matthew@openbiosim.org>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
######################################################################


# Alchemical transfer analysis methods. UWHAM implementation adapted from
# both the `femto` and `AToM-openmm` packages.
__all__ = ["analyse_UWHAM", "analyse_MBAR"]

import pandas as _pd

import numpy as _numpy
import scipy.optimize as _optimize
import scipy.special as _special
import functools as _functools
import pathlib as _pathlib
import os as _os
import warnings as _warnings


def _compute_weights(ln_z, ln_q, factor):
    q_ij = _numpy.exp(ln_q - ln_z)
    return q_ij / (factor * q_ij).sum(axis=-1, keepdims=True)


def _compute_kappa_hessian(ln_z, ln_q, factor, n):
    ln_z = _numpy.insert(ln_z, 0, 0.0)

    w = (factor * _compute_weights(ln_z, ln_q, factor))[:, 1:]
    return -w.T @ w / n + _numpy.diag(w.sum(axis=0) / n)


def _compute_kappa(ln_z, ln_q, factor, n):
    ln_z = _numpy.insert(ln_z, 0, 0.0)

    ln_q_ij_sum = _special.logsumexp(a=ln_q - ln_z, b=factor, axis=1)
    kappa = ln_q_ij_sum.sum() / n + (factor * ln_z).sum()

    w = factor * _compute_weights(ln_z, ln_q, factor)
    grad = -w[:, 1:].sum(axis=0) / n + factor[1:]

    return kappa, grad


def _compute_variance(ln_z, w, factor, n):
    o = w.T @ w / n

    b = o * factor - _numpy.eye(len(ln_z))
    b = b[1:, 1:]

    b_inv_a = -o + o[0, :]
    b_inv_a = b_inv_a[1:, 1:]

    var_matrix = (b_inv_a @ _numpy.linalg.inv(b.T)) / n
    return _numpy.insert(_numpy.diag(var_matrix), 0, 0.0)


def _bias_fcn(epert, lam1, lam2, alpha, u0, w0):
    """
    This is for the bias ilogistic potential
    (lambda2-lambda1) ln[1+exp(-alpha (u-u0))]/alpha + lambda2 u + w0
    """
    ebias1 = _numpy.zeros_like(epert)
    if alpha > 0:
        ee = 1 + _numpy.exp(-alpha * (epert - u0))
        ebias1 = (lam2 - lam1) * _numpy.log(ee) / alpha
    return ebias1 + lam2 * epert + w0


def _npot_fcn(e0, epert, bet, lam1, lam2, alpha, u0, w0):
    # This is the negative reduced energy
    # -beta*(U0+bias)
    return -bet * (e0 + _bias_fcn(epert, lam1, lam2, alpha, u0, w0))


def _estimate_f_i(ln_q, n_k):
    """Estimates the free energies of a set of *sampled* states.


    Args:
        n_k: The number of samples at state ``k``.
        ln_q: array of netgative potentials with ``shape=(n_states,n_samples)``.

    Returns:
        The estimated reduced free energies and their estimated variance.
    """
    n_k = _numpy.array(n_k)

    ln_q = _numpy.array(ln_q).T

    n_samples, n_states = ln_q.shape

    if n_states != len(n_k):
        raise RuntimeError(
            "The number of states do not match: %d != %d" % (n_states, len(n_k))
        )
    if n_samples != n_k.sum():
        raise RuntimeError(
            "The number of samples do not match: %d != %d" % (n_samples, n_k.sum())
        )

    ln_z = _numpy.zeros(len(n_k) - 1)  # ln_z_0 is always fixed at 0.0
    ln_q -= ln_q[:, :1]

    n = n_k.sum()
    factor = n_k / n

    result = _optimize.minimize(
        _functools.partial(_compute_kappa, ln_q=ln_q, n=n, factor=factor),
        ln_z,
        method="trust-ncg",
        jac=True,
        hess=_functools.partial(_compute_kappa_hessian, ln_q=ln_q, n=n, factor=factor),
    )

    if not result.success:
        raise RuntimeError("The UWHAM minimization failed to converge.")

    f_i = _numpy.insert(-result.x, 0, 0.0)
    ln_z = _numpy.insert(result.x, 0, 0.0)

    weights = _compute_weights(ln_z, ln_q, factor)

    if not _numpy.allclose(weights.sum(axis=0) / n, 1.0, atol=1e-2):
        w = weights.sum(axis=0) / n
        _warnings.warn(f"The UWHAM weights do not sum to 1.0 ({w})")

    df_i = _compute_variance(ln_z, weights, factor, n)

    return f_i, df_i, weights / n


def _sort_folders(work_dir):
    """Sorts folder names by lambda value, ensuring they are read correctly.

    Parameters
    ----------
    work_dir : str
        The directory containing the simulation data.

    Returns
    -------
    folders : dict
        A dictionary of folder names and their corresponding lambda values.
    """
    folders = {}
    for folder in _pathlib.Path(work_dir).iterdir():
        if folder.is_dir() and folder.name.startswith("lambda_"):
            try:
                lambda_val = float(folder.name.split("_")[-1])
            except ValueError:
                continue
            folders[lambda_val] = folder
    return {k: v for k, v in sorted(folders.items())}


def _get_inflection_indices(folders):
    # Find folders at which 'direction' goes from 1 to -1
    # This is the point at which the direction of the lambda windows changes
    # NOTE: this assumes that the folders are correctly sorted

    # check that the keys are sorted
    keys = list(folders.keys())
    if keys != sorted(keys):
        raise ValueError(f"Folders are not sorted correctly. {keys} != {sorted(keys)}")

    directions = []
    for folder in folders.values():
        df = _pd.read_csv(folder / "openmm.csv")
        direction = df["direction"].values[0]
        directions.append(direction)

    # get the indices at which the direction changes
    for i in range(len(directions) - 1):
        if directions[i] != directions[i + 1]:
            inflection_indices = (i, i + 1)
            break

    return inflection_indices


def analyse_UWHAM(work_dir, ignore_lower, ignore_upper, inflection_indices=None):
    """
    Analyse the output of BioSimSpace AToM simulations.

    Parameters
    ----------
    work_dir : str
        The directory containing the simulation data.
    ignore_lower : int
        The number of rows to ignore at the start of each file.
    ignore_upper : int
        The number of rows to ignore at the end of each file.
    inflection_indices : tuple, optional
        The point at which 'direction' changes.
        Should be (last index of direction 1, first index of direction 2).
        If not provided not provided, will be implied from files.

    Returns
    -------
    ddg_total : :class:`BioSimSpace.Types.Energy`
        The free energy.
    ddg_total_error : :class:`BioSimSpace.Types.Energy`
        The error in the free energy.

    """
    # NOTE: This code is not designed to work with repex
    # It always assumes that each window is at the same temperature
    dataframes = []
    slices = {}
    total_states = 0
    total_samples = 0
    folders = _sort_folders(work_dir)
    if inflection_indices is None:
        inflection_indices = _get_inflection_indices(folders)
    for folder in folders.values():
        df = _pd.read_csv(folder / "openmm.csv")
        # drop the first `ignore_lower` rows of each df
        if ignore_upper is not None:
            df = df.iloc[ignore_lower:ignore_upper]
        else:
            df = df.iloc[ignore_lower:]
        # Beta values, assuming that energies are in kj/mol
        df["beta"] = 1 / (0.001986209 * df["temperature"])
        total_states += 1
        total_samples += len(df)
        for sub_df in df.groupby("window"):
            # get value of window
            window = sub_df[0]
            # check if window is in slices
            if window not in slices:
                slices[window] = []
            # append the dataframe to the list of dataframes for that window
            # now get the tuple 'sub_df' and convert it to a dataframe
            s = sub_df[1]
            slices[window].append(s)

    # now combine all dataframes in each slice
    for window in slices:
        # get the dataframes for the current window
        dfs = slices[window]
        # combine the dataframes
        combined_df = _pd.concat(dfs)
        dataframes.append(combined_df)

    # sort 'dataframes' based on 'window'
    dataframes = sorted(dataframes, key=lambda x: x["window"].values[0])

    pots = []
    pert_es = []
    n_samples = []
    # check that all dataframes are the same length, throw a warning if they are not
    for df in dataframes:
        n_samples.append(len(df))
        e0 = df["pot_en"].values
        pert_e = df["pert_en"].values
        pots.append(e0)
        pert_es.append(pert_e)

    # Should only matter in cases where states are at different temps,
    # leaving here for debugging and parity with GL code
    for be in range(len(n_samples)):
        pots[be] = pots[be] - _bias_fcn(
            pert_es[be],
            lam1=dataframes[be]["lambda1"].values[0],
            lam2=dataframes[be]["lambda2"].values[0],
            alpha=dataframes[be]["alpha"].values[0],
            u0=dataframes[be]["uh"].values[0],
            w0=dataframes[be]["w0"].values[0],
        )
    # We will assume that the point at which leg1 and leg2 are split is halfway through
    n_samples_first_half = n_samples[: inflection_indices[0] + 1]
    pots_first_half = _numpy.concatenate(pots[: inflection_indices[0] + 1])
    pert_es_first_half = _numpy.concatenate(pert_es[: inflection_indices[0] + 1])
    ln_q = _numpy.zeros((inflection_indices[0] + 1, len(pots_first_half)))
    sid = 0

    for be in range(len(n_samples_first_half)):
        lnq = _npot_fcn(
            e0=pots_first_half,
            epert=pert_es_first_half,
            bet=dataframes[be]["beta"].values[0],
            lam1=dataframes[be]["lambda1"].values[0],
            lam2=dataframes[be]["lambda2"].values[0],
            alpha=dataframes[be]["alpha"].values[0],
            u0=dataframes[be]["uh"].values[0],
            w0=dataframes[be]["w0"].values[0],
        )
        ln_q[sid] = lnq
        sid += 1
    f_i, d_i, weights = _estimate_f_i(ln_q, n_samples_first_half)
    ddg = f_i[-1] - f_i[0]
    ddg1 = ddg / dataframes[0]["beta"].values[0]
    # print(f"Forward leg: {ddg1}")
    ddg_error_1 = _numpy.sqrt(d_i[-1] + d_i[0]) / dataframes[0]["beta"].values[0]

    n_samples_second_half = n_samples[inflection_indices[1] :]
    pots_second_half = _numpy.concatenate(pots[inflection_indices[1] :])
    pert_es_second_half = _numpy.concatenate(pert_es[inflection_indices[1] :])
    ln_q = _numpy.zeros((total_states - inflection_indices[1], len(pots_second_half)))
    sid = 0

    # note the order of (be, te)
    for be in range(len(n_samples_second_half)):
        lnq = _npot_fcn(
            e0=pots_second_half,
            epert=pert_es_second_half,
            bet=dataframes[be]["beta"].values[0],
            lam1=dataframes[be]["lambda1"].values[0],
            lam2=dataframes[be]["lambda2"].values[0],
            alpha=dataframes[be]["alpha"].values[0],
            u0=dataframes[be]["uh"].values[0],
            w0=dataframes[be]["w0"].values[0],
        )
        ln_q[sid] = lnq
        sid += 1
    f_i, d_i, weights = _estimate_f_i(ln_q, n_samples_second_half)
    ddg = f_i[-1] - f_i[0]
    ddg2 = ddg / dataframes[0]["beta"].values[0]
    # print(f"Reverse leg: {ddg2}")
    ddg_error_2 = _numpy.sqrt(d_i[-1] + d_i[0]) / dataframes[0]["beta"].values[0]

    ddg_total = ddg1 - ddg2
    ddg_total_error = _numpy.sqrt(ddg_error_1**2 + ddg_error_2**2)
    from BioSimSpace.Units import Energy as _Energy

    ddg_total = ddg_total * _Energy.kcal_per_mol
    ddg_total_error = ddg_total_error * _Energy.kcal_per_mol

    return ddg_total, ddg_total_error


def analyse_MBAR(work_dir):
    """
    Analyse the MBAR-compatible outputs.
    Adapted version of BioSimSpace _analyse_internal function
    """
    from ._relative import Relative as _Relative
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol
    from .. import Units as _Units

    try:
        from alchemlyb.estimators import AutoMBAR as _AutoMBAR
    except ImportError:
        from alchemlyb.estimators import MBAR as _AutoMBAR

    if not isinstance(work_dir, str):
        raise TypeError("work_dir must be a string")
    if not _os.path.isdir(work_dir):
        raise ValueError("work_dir must be a valid directory")

    glob_path = _pathlib.Path(work_dir)
    files = sorted(glob_path.glob("**/energies*.csv"))

    # Slightly more complicated than a standard FE calculation
    # the key complication comes from the need to split the forward and reverse legs
    # instead of being inherently separate as in a standard FE calculation, they
    # are dictated by 'direction'. This means that the energy arrays need
    # to be re-numbered in to separate forward and reverse legs.

    # need to make sure that all lambdas were run at the same temp
    temps = []
    dataframes_forward = []
    dataframes_backward = []
    for file in files:
        # read the csv to a dataframe
        df = _pd.read_csv(file)
        # read the temperature column and make sure all values in it are equal
        temp = df["temperature"].unique()
        if len(temp) != 1:
            raise ValueError(f"Temperature column in {file} is not uniform")
        # check if the last column in the dataframe is full of NaNs
        if df.iloc[:, -1:].isnull().values.all():
            reverse = False
        else:
            reverse = True
        temps.append(temp[0])
        # now drop the temperature column
        df = df.drop(columns=["temperature"])
        # remove columns with NaN values
        df = df.dropna(axis=1)
        # we will need to match the fep-lambda value to the correct new value
        # first get fep-lambda, should be the same value for all entries in the 'fep-lambda' column
        fep_lambda = df["fep-lambda"].unique()
        if len(fep_lambda) != 1:
            raise ValueError(f"fep-lambda column in {file} is not uniform")
        # find all columns whose titles are only numbers
        cols = []
        num_lams = 0
        for col in df.columns:
            try:
                val = float(col)
                if val == fep_lambda:
                    index_fep_lambda = num_lams
                num_lams += 1
            except ValueError:
                cols.append(col)
        new_lambdas = list(_numpy.linspace(0, 1, num_lams))
        new_fep_lambda = new_lambdas[index_fep_lambda]
        new_cols = cols + new_lambdas
        # rename the columns
        df.columns = new_cols
        # now replace all values in the fep-lambda column with the new value
        df["fep-lambda"] = new_fep_lambda
        df.set_index(cols, inplace=True)
        if reverse:
            dataframes_backward.append(df)
        else:
            dataframes_forward.append(df)

            # check that all temperatures are the same
    if len(set(temps)) != 1:
        raise ValueError("All temperatures must be the same")
    data_forward = _Relative._preprocess_data(dataframes_forward, "MBAR")
    data_backward = _Relative._preprocess_data(dataframes_backward, "MBAR")
    print("\n\n\n\n\n")
    print(type(data_forward))
    data_forward.attrs = {
        "temperature": temps[0],
        "energy_unit": "kJ/mol",
    }
    data_backward.attrs = {
        "temperature": temps[0],
        "energy_unit": "kJ/mol",
    }
    try:
        alchem_forward = _AutoMBAR().fit(data_forward)
    except ValueError as e:
        raise ValueError(f"Error in fitting forward leg of MBAR calculation: {e}")

    try:
        alchem_backward = _AutoMBAR().fit(data_backward)
    except ValueError as e:
        raise ValueError(f"Error in fitting backward leg of MBAR calculation: {e}")

    alchem_forward.delta_f_.attrs = {
        "temperature": temps[0],
        "energy_unit": "kJ/mol",
    }
    delta_f_for = _to_kcalmol(alchem_forward.delta_f_)
    alchem_forward.d_delta_f_.attrs = {
        "temperature": temps[0],
        "energy_unit": "kJ/mol",
    }
    d_delta_f_for = _to_kcalmol(alchem_forward.d_delta_f_)
    data_forward_final = []
    for lamb in new_lambdas:
        x = new_lambdas.index(lamb)
        mbar_value = delta_f_for.iloc[0, x]
        mbar_error = d_delta_f_for.iloc[0, x]

        data_forward_final.append(
            (
                lamb,
                (mbar_value) * _Units.Energy.kcal_per_mol,
                (mbar_error) * _Units.Energy.kcal_per_mol,
            )
        )
    alchem_backward.delta_f_.attrs = {
        "temperature": temps[0],
        "energy_unit": "kJ/mol",
    }
    delta_f_back = _to_kcalmol(alchem_backward.delta_f_)
    alchem_backward.d_delta_f_.attrs = {
        "temperature": temps[0],
        "energy_unit": "kJ/mol",
    }
    d_delta_f_back = _to_kcalmol(alchem_backward.d_delta_f_)
    data_backward_final = []
    for lamb in new_lambdas:
        x = new_lambdas.index(lamb)
        mbar_value = delta_f_back.iloc[0, x]
        mbar_error = d_delta_f_back.iloc[0, x]

        data_backward_final.append(
            (
                lamb,
                (mbar_value) * _Units.Energy.kcal_per_mol,
                (mbar_error) * _Units.Energy.kcal_per_mol,
            )
        )

    return data_forward_final, data_backward_final
