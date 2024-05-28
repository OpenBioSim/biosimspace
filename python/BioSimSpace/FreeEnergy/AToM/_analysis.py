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

__all__ = ["analyse"]

import pandas as pd

import numpy
import scipy.optimize
import scipy.special
import pathlib
import pandas as pd
import functools
import pymbar.timeseries


def _compute_weights(ln_z, ln_q, factor):
    q_ij = numpy.exp(ln_q - ln_z)
    return q_ij / (factor * q_ij).sum(axis=-1, keepdims=True)


def _compute_kappa_hessian(ln_z, ln_q, factor, n):
    ln_z = numpy.insert(ln_z, 0, 0.0)

    w = (factor * _compute_weights(ln_z, ln_q, factor))[:, 1:]
    return -w.T @ w / n + numpy.diag(w.sum(axis=0) / n)


def _compute_kappa(ln_z, ln_q, factor, n):
    ln_z = numpy.insert(ln_z, 0, 0.0)

    ln_q_ij_sum = scipy.special.logsumexp(a=ln_q - ln_z, b=factor, axis=1)
    kappa = ln_q_ij_sum.sum() / n + (factor * ln_z).sum()

    w = factor * _compute_weights(ln_z, ln_q, factor)
    grad = -w[:, 1:].sum(axis=0) / n + factor[1:]

    return kappa, grad


def _compute_variance(ln_z, w, factor, n):
    o = w.T @ w / n

    b = o * factor - numpy.eye(len(ln_z))
    b = b[1:, 1:]

    b_inv_a = -o + o[0, :]
    b_inv_a = b_inv_a[1:, 1:]

    var_matrix = (b_inv_a @ numpy.linalg.inv(b.T)) / n
    return numpy.insert(numpy.diag(var_matrix), 0, 0.0)


def _bias_fcn(epert, lam1, lam2, alpha, u0, w0):
    """
    This is for the bias ilogistic potential
    (lambda2-lambda1) ln[1+exp(-alpha (u-u0))]/alpha + lambda2 u + w0
    """
    ebias1 = numpy.zeros_like(epert)
    if alpha > 0:
        ee = 1 + numpy.exp(-alpha * (epert - u0))
        ebias1 = (lam2 - lam1) * numpy.log(ee) / alpha
    return ebias1 + lam2 * epert + w0


def npot_fcn(e0, epert, bet, lam1, lam2, alpha, u0, w0):
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
    n_k = numpy.array(n_k)

    ln_q = numpy.array(ln_q).T

    n_samples, n_states = ln_q.shape

    if n_states != len(n_k):
        raise RuntimeError(
            "The number of states do not match: %d != %d" % (n_states, len(n_k))
        )
    if n_samples != n_k.sum():
        raise RuntimeError(
            "The number of samples do not match: %d != %d" % (n_samples, n_k.sum())
        )

    ln_z = numpy.zeros(len(n_k) - 1)  # ln_z_0 is always fixed at 0.0
    ln_q -= ln_q[:, :1]

    n = n_k.sum()
    factor = n_k / n

    result = scipy.optimize.minimize(
        functools.partial(_compute_kappa, ln_q=ln_q, n=n, factor=factor),
        ln_z,
        method="trust-ncg",
        jac=True,
        hess=functools.partial(_compute_kappa_hessian, ln_q=ln_q, n=n, factor=factor),
    )

    if not result.success:
        raise RuntimeError("The UWHAM minimization failed to converge.")

    f_i = numpy.insert(-result.x, 0, 0.0)
    ln_z = numpy.insert(result.x, 0, 0.0)

    weights = _compute_weights(ln_z, ln_q, factor)

    if not numpy.allclose(weights.sum(axis=0) / n, 1.0, atol=1e-3):
        raise RuntimeError("The UWHAM weights do not sum to 1.0")

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
    for folder in pathlib.Path(work_dir).iterdir():
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
        df = pd.read_csv(folder / "openmm.csv")
        direction = df["direction"].values[0]
        directions.append(direction)

    # get the indices at which the direction changes
    for i in range(len(directions) - 1):
        if directions[i] != directions[i + 1]:
            inflection_indices = (i, i + 1)
            break

    return inflection_indices


def analyse(work_dir, inflection_indices=None):
    """
    Analyse the output of BioSimSpace AToM simulations.

    Parameters
    ----------
    work_dir : str
        The directory containing the simulation data.
    inflection_indices : tuple, optional
        The point at which 'direction' changes.
        Should be (last index of direction 1, first index of direction 2).
        If not provided not provided, will be implied from files.

    """
    dataframes = []
    slices = {}
    total_states = 0
    total_samples = 0
    folders = _sort_folders(work_dir)
    if inflection_indices is None:
        inflection_indices = _get_inflection_indices(folders)
    for folder in folders.values():
        df = pd.read_csv(folder / "openmm.csv")
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
        combined_df = pd.concat(dfs)
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

    # We will assume that the point at which leg1 and leg2 are split is halfway through
    n_samples_first_half = n_samples[: inflection_indices[0] + 1]
    pots_first_half = numpy.concatenate(pots[: inflection_indices[0] + 1])
    pert_es_first_half = numpy.concatenate(pert_es[: inflection_indices[0] + 1])
    ln_q = numpy.zeros((inflection_indices[0] + 1, len(pots_first_half)))
    sid = 0

    for be in range(len(n_samples_first_half)):
        lnq = npot_fcn(
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

    print(ln_q)
    f_i, d_i, weights = _estimate_f_i(ln_q, n_samples_first_half)
    ddg = f_i[-1] - f_i[0]
    ddg1 = ddg / dataframes[0]["beta"].values[0]
    ddg_error_1 = numpy.sqrt(d_i[-1] + d_i[0]) / dataframes[0]["beta"].values[0]

    n_samples_second_half = n_samples[inflection_indices[1] :]
    pots_second_half = numpy.concatenate(pots[inflection_indices[1] :])
    pert_es_second_half = numpy.concatenate(pert_es[inflection_indices[1] :])
    ln_q = numpy.zeros((total_states - inflection_indices[1], len(pots_second_half)))
    sid = 0

    # note the order of (be, te)
    for be in range(len(n_samples_second_half)):
        lnq = npot_fcn(
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
    ddg_error_2 = numpy.sqrt(d_i[-1] + d_i[0]) / dataframes[0]["beta"].values[0]

    ddg_total = ddg1 - ddg2
    ddg_total_error = numpy.sqrt(ddg_error_1**2 + ddg_error_2**2)
    return ddg_total, ddg_total_error
