__all__ = ["_FreeEnergyMixin"]


from ._protocol import Protocol as _Protocol


class _FreeEnergyMixin(_Protocol):
    """A mixin for storing free energy protocols."""

    def __init__(
        self,
        lam=0.0,
        lam_vals=None,
        min_lam=0.0,
        max_lam=1.0,
        num_lam=11,
        perturbation_type="full",
    ):
        """Constructor.

        Parameters
        ----------

        lam : float or pandas.Series
            The perturbation parameter: [0.0, 1.0]

        lam_vals : [float] or pandas.DataFrame
            A list of lambda values.

        min_lam : float or pandas.Series
            The minimum lambda value.

        max_lam : float or pandas.Series
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
            "restraint" : Perturb the receptor-ligand restraint strength by linearly
                        scaling the force constants (0.0->value).
            "release_restraint" : Used with multiple distance restraints to release all
                                  restraints other than the "permanent" one when the ligand
                                  is fully decoupled. Note that lambda = 0.0 is the fully
                                  released state, and lambda = 1.0 is the fully restrained
                                  state (i.e. 0.0 -> value).

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
            f"lam={self._lambda.to_string()}, "
            f"lam_vals={self._lambda_vals.to_string()}"
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol._FreeEnergyMixin: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "BioSimSpace.Protocol.Custom"
        else:
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
            "restraint" : Perturb the receptor-ligand restraint strength by linearly
                        scaling the force constants (0.0->value).
            "release_restraint" : Used with multiple distance restraints to release all
                                  restraints other than the "permanent" one when the ligand
                                  is fully decoupled. Note that lambda = 0.0 is the fully
                                  released state, and lambda = 1.0 is the fully restrained
                                  state (i.e. 0.0 -> value).

        """
        if type(perturbation_type) is not str:
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
            "restraint",
            "release_restraint",
        ]

        if perturbation_type not in allowed_perturbation_types:
            raise ValueError(
                f"'perturbation_type' must be one of: {allowed_perturbation_types}"
            )

        self._perturbation_type = perturbation_type

    @staticmethod
    def _check_column_name(df):
        """Check if the dataframe or series has the right column name."""
        import pandas as _pd
        import warnings

        permitted_names = [
            "fep",
            "bonded",
            "coul",
            "vdw",
            "restraint",
            "mass",
            "temperature",
        ]
        if isinstance(df, _pd.Series):
            for name in df.index:
                if not name in permitted_names:
                    warnings.warn(
                        f"{name} not in the list of permitted names, "
                        f"so may not be supported by the MD Engine "
                        f'({" ".join(permitted_names)}).'
                    )
        elif isinstance(df, _pd.DataFrame):
            for name in df.columns:
                if not name in permitted_names:
                    warnings.warn(
                        f"{name} not in the list of permitted names, "
                        f"so may not be supported by the MD Engine "
                        f'({" ".join(permitted_names)}).'
                    )

    def getLambda(self, type="float"):
        """
        Get the value of the perturbation parameter.

        Parameters
        ----------

        type : string
            The type of lambda to be returned. ('float', 'series')

        Returns
        -------

        lam : float or pd.Series
            The value of the perturbation parameter.
        """
        import warnings

        if type.lower() == "float":
            if len(self._lambda) == 1:
                return float(self._lambda.iloc[0])
            else:
                warnings.warn(
                    f"The {self._lambda} has more than one value, return as pd.Series."
                )
                return self._lambda
        elif type.lower() == "series":
            return self._lambda
        else:
            raise ValueError(f"{type} could only be 'float' or 'series'.")

    def getLambdaIndex(self):
        """
        Get the index of the lambda value within the lambda array.

        Returns
        -------

        index : int
            The index of the lambda value in teh lambda array.
        """
        return self._lambda_vals.index[
            (self._lambda_vals == self._lambda).all(axis=1)
        ].tolist()[0]

    def getLambdaValues(self, type="list"):
        """
        Get the lambda values.

        Parameters
        ----------

        type : string
            The type of lambda values to be returned. ('list', 'dataframe')

        Returns
        -------

        lam_vals : [float, ] or pd.DataFrame
            The lambda values.
        """
        import warnings

        if type.lower() == "list":
            lambda_list = self._lambda_vals.values.tolist()
            if len(lambda_list[0]) == 1:
                # Ensure that return [float, ] when lambda_list only has one
                # column
                return [lam[0] for lam in lambda_list]
            else:
                warnings.warn(
                    f"The {self._lambda_vals} has more than one column, return as pd.DataFrame."
                )
                return self._lambda_vals
        elif type.lower() == "dataframe":
            return self._lambda_vals
        else:
            raise ValueError(f"{type} could only be 'list' or 'dataframe'.")

    def setLambda(self, lam):
        """
        Set the current lambda.

        Parameters
        ----------

        lam : float or pandas.Series
            The perturbation parameter: [0.0, 1.0]
        """
        import pandas as _pd

        if not isinstance(lam, _pd.Series):
            # For pandas < 1.4, TypeError won't be raised if the type cannot
            # be converted to float
            lam = _pd.Series(data={"fep": lam}, dtype=float)
        else:
            self._check_column_name(lam)

        # Make sure the lambda value is in the list.
        if not (lam == self._lambda_vals).all(axis=1).any():
            raise ValueError("'lam' is not a member of the 'lam_vals' list!")

        self._lambda = lam

    @staticmethod
    def checkLambdaValues(lam_vals, min_lam=None, max_lam=None, num_lam=None):
        """
        Sanity check of the list of lambda values.

        Parameters
        ----------

        lam_vals : [float] or pandas.DataFrame
            A list of lambda values.

        min_lam : float or pandas.Series
            The minimum lambda value.

        max_lam : float or pandas.Series
            The maximum lambda value.

        num_lam : int
            The number of lambda values.

        Returns
        -------
        lam_vals : pandas.DataFrame
            The pd.DataFrame representing the checked lambda values.
        """
        import pandas as _pd
        import numpy as _np

        # A list of lambda values takes precedence.
        if lam_vals is not None:
            if not isinstance(lam_vals, _pd.DataFrame):
                # For pandas < 1.4, TypeError won't be raised if the type
                # cannot be converted to float
                lam_vals = _pd.DataFrame(data={"fep": lam_vals}, dtype=float)
            else:
                _FreeEnergyMixin._check_column_name(lam_vals)

            # Make sure all values are in range [0.0, 1.0]
            if not ((lam_vals <= 1).all(axis=None) and (lam_vals >= 0).all(axis=None)):
                raise ValueError("'lam_vals' must contain values in range [0.0, 1.0]")

            # Sort the values.
            lam_vals.sort_values(lam_vals.columns.values.tolist(), inplace=True)

            # Make sure there are no duplicates.
            if lam_vals.duplicated().any():
                raise ValueError("'lam_vals' cannot contain duplicate values!")

            # Update the index to the sorted rows
            lam_vals.reset_index(drop=True, inplace=True)

        else:
            # Convert int to float.
            if not isinstance(min_lam, _pd.Series):
                # For pandas < 1.4, TypeError won't be raised if the type cannot
                # be converted to float
                min_lam = _pd.Series(data={"fep": min_lam}, dtype=float)
            else:
                _FreeEnergyMixin._check_column_name(lam_vals)

            if not isinstance(max_lam, _pd.Series):
                # For pandas < 1.4, TypeError won't be raised if the type cannot
                # be converted to float
                max_lam = _pd.Series(data={"fep": max_lam}, dtype=float)
            else:
                _FreeEnergyMixin._check_column_name(lam_vals)

            if not type(num_lam) is int:
                raise TypeError("'num_lam' must be of type 'int'.")

            # Validate values.
            if (min_lam < 0).any() or (min_lam > 1).any():
                raise ValueError("'min_lam' must be in range [0.0, 1.0]")

            if (max_lam < 0).any() or (max_lam > 1).any():
                raise ValueError("'max_lam' must be in range [0.0, 1.0]")

            if num_lam <= 2:
                raise ValueError("'num_lam' must be greater than 2!")

            if not (min_lam.index == max_lam.index).all():
                raise ValueError("'min_lam' and 'max_lam' must have the same index!")

            if (min_lam >= max_lam).any():
                raise ValueError("'min_lam' must be less than 'max_lam'!")

            # Set values.
            lambda_vals = _pd.DataFrame(
                data=_np.linspace(min_lam, max_lam, num_lam), columns=min_lam.index
            )

            lam_vals = lambda_vals.round(5)
        return lam_vals

    def setLambdaValues(
        self, lam, lam_vals=None, min_lam=None, max_lam=None, num_lam=None
    ):
        """
        Set the list of lambda values.

        Parameters
        ----------

        lam : float or pandas.Series
            The perturbation parameter: [0.0, 1.0]

        lam_vals : [float] or pandas.DataFrame
            A list of lambda values.

        min_lam : float or pandas.Series
            The minimum lambda value.

        max_lam : float or pandas.Series
            The maximum lambda value.

        num_lam : int
            The number of lambda values.
        """
        self._lambda_vals = self.checkLambdaValues(lam_vals, min_lam, max_lam, num_lam)
        self.setLambda(lam)
