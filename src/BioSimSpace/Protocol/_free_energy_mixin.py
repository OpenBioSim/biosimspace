######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
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
#####################################################################

"""A mixin for alchemical free energy protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["_FreeEnergyMixin"]


class _FreeEnergyMixin:
    """A mixin for alchemical free energy protocols."""

    def __init__(
        self,
        lam=0.0,
        lam_vals=None,
        min_lam=0.0,
        max_lam=1.0,
        num_lam=11,
        perturbation_type="full",
    ):
        """
        Constructor.

        Parameters
        ----------

        lam : float
            The perturbation parameter: [0.0, 1.0]

        lam_vals : [float]
            The list of lambda parameters.

        min_lam : float
            The minimum lambda value.

        max_lam : float
            The maximum lambda value.

        num_lam : int
            The number of lambda values.

        perturbation_type : str
            The type of perturbation to perform. Options are:
             "full" : A full perturbation of all terms (default option).
             "discharge_soft" : Perturb all discharging soft atom charge terms (i.e. value->0.0).
             "vanish_soft" : Perturb all vanishing soft atom LJ terms (i.e. value->0.0).
             "flip" : Perturb all hard atom terms as well as bonds/angles.
             "grow_soft" : Perturb all growing soft atom LJ terms (i.e. 0.0->value).
             "charge_soft" : Perturb all charging soft atom LJ terms (i.e. 0.0->value).

             Currently perturubation_type != "full" is only supported by
             BioSimSpace.Process.Somd.
        """

        # Call the base class constructor.
        super().__init__()

        # Validate and set the lambda values.
        self.setLambdaValues(lam, lam_vals, min_lam, max_lam, num_lam)

        # Set the perturbation type. Default is "full", i.e. onestep protocol.
        self.setPerturbationType(perturbation_type)

    def _get_parm(self):
        """Return a string representation of the parameters."""
        return (
            f"lam={self._lambda}, "
            f"lam_vals={self._lambda_vals}, "
            f"perturbation_type='{self._perturbation_type}'"
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        return f"<BioSimSpace.Protocol._FreeEnergyMixin: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return f"BioSimSpace.Protocol._FreeEnergyMixin({self._get_parm()})"

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, _FreeEnergyMixin):
            return False

        return (
            self._lambda == other._lambda
            and self._lambda_vals == other._lambda_vals
            and self._perturbation_type == other._perturbation_type
        )

    def getPerturbationType(self):
        """
        Get the perturbation type.

        Returns
        -------

        perturbation_type : str
            The perturbation type.
        """
        return self._perturbation_type

    def setPerturbationType(self, perturbation_type):
        """
        Set the perturbation type.

        Parameters
        ----------

        perturbation_type : str
            The type of perturbation to perform. Options are:
             "full" : A full perturbation of all terms (default option).
             "discharge_soft" : Perturb all discharging soft atom charge terms (i.e. value->0.0).
             "vanish_soft" : Perturb all vanishing soft atom LJ terms (i.e. value->0.0).
             "flip" : Perturb all hard atom terms as well as bonds/angles.
             "grow_soft" : Perturb all growing soft atom LJ terms (i.e. 0.0->value).
             "charge_soft" : Perturb all charging soft atom LJ terms (i.e. 0.0->value).
        """
        if not isinstance(perturbation_type, str):
            raise TypeError("'perturbation_type' must be of type 'str'")

        # Convert to lower case and strip whitespace.
        perturbation_type = perturbation_type.lower().replace(" ", "")

        allowed_perturbation_types = [
            "full",
            "discharge_soft",
            "vanish_soft",
            "flip",
            "grow_soft",
            "charge_soft",
        ]

        if perturbation_type not in allowed_perturbation_types:
            raise ValueError(
                f"'perturbation_type' must be one of: {allowed_perturbation_types}"
            )

        self._perturbation_type = perturbation_type

    def getLambda(self):
        """
        Get the value of the perturbation parameter.

        Returns
        -------

        lam : float
            The value of the perturbation parameter.
        """
        return self._lambda

    def getLambdaIndex(self):
        """
        Get the index of the lambda value within the lambda array.

        Returns
        -------

        index : int
            The index of the lambda value in teh lambda array.
        """
        return self._lambda_vals.index(self._lambda)

    def getLambdaValues(self):
        """
        Get the list of lambda values.

        Returns
        -------

        lam_vals : [float]
            The list of lambda values.
        """
        return self._lambda_vals

    def setLambdaValues(
        self, lam, lam_vals=None, min_lam=None, max_lam=None, num_lam=None
    ):
        """
        Set the list of lambda values.

        Parameters
        ----------

        lam : float
            The perturbation parameter: [0.0, 1.0]

        lam_vals : [float]
            A list of lambda values.

        min_lam : float
            The minimum lambda value.

        max_lam : float
            The maximum lambda value.

        num_lam : int
            The number of lambda values.
        """

        # Convert int to float.
        if type(lam) is int:
            lam = float(lam)

        # Validate the lambda parameter.
        if not isinstance(lam, float):
            raise TypeError("'lam' must be of type 'float'.")

        self._lambda = lam

        # A list of lambda values takes precedence.
        if lam_vals is not None:
            # Make sure list (or tuple) has been passed.
            if not isinstance(lam_vals, (list, tuple)):
                raise TypeError("'lam_vals' must be of type 'list'.")

            # Make sure all lambda values are of type 'float'.
            if not all(isinstance(x, float) for x in lam_vals):
                raise TypeError("'lam_vals' must contain values of type 'float'.")

            # Make sure all values are in range [0.0, 1.0]
            for lam in lam_vals:
                if lam < 0:
                    raise ValueError(
                        "'lam_vals' must contain values in range [0.0, 1.0]"
                    )
                elif lam > 1:
                    raise ValueError(
                        "'lam_vals' must contain values in range [0.0, 1.0]"
                    )

            # Sort the values.
            lam_vals.sort()

            # Make sure there are no duplicates.
            if len(set(lam_vals)) < len(lam_vals):
                raise ValueError("'lam_vals' cannot contain duplicate values!")

            # Make sure the lambda value is in the list.
            if not lam in lam_vals:
                raise ValueError("'lam' is not a member of the 'lam_vals' list!")

            # Set the values.
            self._lambda_vals = lam_vals

        else:
            # Convert int to float.

            if type(min_lam) is int:
                min_lam = float(min_lam)

            if type(max_lam) is int:
                max_lam = float(max_lam)

            # Validate type.

            if not isinstance(min_lam, float):
                raise TypeError("'min_lam' must be of type 'float'.")

            if not isinstance(max_lam, float):
                raise TypeError("'max_lam' must be of type 'float'.")

            if not type(num_lam) is int:
                raise TypeError("'num_lam' must be of type 'int'.")

            # Validate values.

            if min_lam < 0 or min_lam > 1:
                raise ValueError("'min_lam' must be in range [0.0, 1.0]")

            if max_lam < 0 or max_lam > 1:
                raise ValueError("'max_lam' must be in range [0.0, 1.0]")

            if num_lam <= 2:
                raise ValueError("'num_lam' must be greater than 2!")

            if min_lam >= max_lam:
                raise ValueError("'min_lam' must be less than 'max_lam'!")

            # Set values.

            self._lambda_vals = []
            step = (max_lam - min_lam) / (num_lam - 1)

            for x in range(0, num_lam):
                self._lambda_vals.append(round(min_lam + (x * step), 5))

            # Make sure the lambda value is in the list.
            if not self._lambda in self._lambda_vals:
                raise ValueError(
                    "'lam' (%5.4f) is not a member of the 'lam_vals' list: %s"
                    % (self._lambda, self._lambda_vals)
                )
