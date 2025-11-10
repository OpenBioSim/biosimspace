#####################################################################
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

"""
Functionality for defining and validating BioSimSpace input and output requirements.
Author: Lester Hedges <lester.hedges@gmail.com>.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "Boolean",
    "Integer",
    "Float",
    "String",
    "File",
    "FileSet",
    "Length",
    "Area",
    "Volume",
    "Angle",
    "Charge",
    "Energy",
    "Mass",
    "Pressure",
    "Temperature",
    "Time",
]


class Requirement:
    """Base class for BioSimSpace Node requirements."""

    # Set the argparse argument type.
    _arg_type = None

    # Default to single arguments.
    _is_multi = False

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
        optional=False,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default :
            The default value.

        unit : str
            The unit.

        minimum :
            The minimum allowed value.

        maximum :
            The maximum allowed value.

        allowed : list
            A list of allowed values.

        optional : bool
            Whether the requirement is optional.
        """

        # Don't allow user to create an instance of this base class.
        if type(self) is Requirement:
            raise Exception("<Requirement> must be subclassed.")

        # Required keyword arguments.

        if help is None:
            raise ValueError("Missing 'help' keyword argument!")
        elif not isinstance(help, str):
            raise TypeError("'help' keyword argument must be of type 'str'")
        else:
            self._help = help

        # Optional keywords arguments.

        if not isinstance(optional, bool):
            raise TypeError("'optional' keyword argument must be of type 'bool'")

        # Set defaults
        self._value = None
        self._default = None
        self._unit = unit
        self._min = None
        self._max = None
        self._allowed = None
        self._is_optional = optional

        # Set the minimum value.
        if minimum is not None:
            self._min = self._validate(minimum)

        # Set the maximum value.
        if maximum is not None:
            self._max = self._validate(maximum)

        # Min and max.
        if minimum is not None and maximum is not None:
            if self._max < self._min:
                raise ValueError(
                    "The maximum value '%s' is less than the minimum '%s'"
                    % (maximum, minimum)
                )

        # Convert tuple to list.
        if isinstance(allowed, tuple):
            allowed = list(allowed)

        # Set the allowed values.
        if allowed is not None:
            if not isinstance(allowed, list):
                allowed = [allowed]
            self._allowed = [self._validate(x) for x in allowed]

            # Conflicting requirements.
            if self._min is not None or self._max is not None:
                raise ValueError(
                    "Conflicting requirement: cannot have allowed values and min/max."
                )

        # Set the default value.
        if default is not None:
            # Make sure the default is the right type.
            self._default = self._validate(default)
            # Make sure the default satisfies the constraints.
            self._validate_default()
            # Flag the the requirement is optional.
            self._is_optional = True

    def setValue(self, value, name=None):
        """
        Validate and set the value.

        Parameters
        ----------

        value :
            The value of the input requirement.

        name : str
            The name of the requirement.
        """

        if value is None and not self._is_optional:
            if isinstance(name, str):
                raise ValueError("Value is unset for requirement '%s'!" % name)
            else:
                raise ValueError("Value is unset!")

        # Validate the value.
        value = self._validate(value)

        # Now check value against any constraints.

        # Minimum.
        if self._min is not None and value < self._min:
            raise ValueError(
                "The value (%s) is less than the allowed "
                "minimum (%s)" % (value, self._min)
            )

        # Maximum.
        if self._max is not None and value > self._max:
            raise ValueError(
                "The value (%s) is greater than the allowed "
                "maximum (%s)" % (value, self._max)
            )

        # Allowed values.
        if self._allowed is not None and value not in self._allowed:
            # For String requirements, strip whitespace and ignore case.
            if isinstance(self, String):
                new_value = value.replace(" ", "").upper()
                allowed = [x.replace(" ", "").upper() for x in self._allowed]

                # If we find a match, then set to the unmodified allowed value
                # at the matching index.
                if new_value in allowed:
                    value = self._allowed[allowed.index(new_value)]
                else:
                    raise ValueError(
                        "The value (%s) is not in the list of allowed values: "
                        "%s" % (value, str(self._allowed))
                    )
            else:
                raise ValueError(
                    "The value (%s) is not in the list of allowed values: "
                    "%s" % (value, str(self._allowed))
                )

        # All is okay. Set the value.
        self._value = value

    def getValue(self):
        """Return the value."""
        return self._value

    def getDefault(self):
        """Return the default value."""
        return self._default

    def getUnit(self):
        """
        Return the unit.

        Returns
        -------

        unit : str
            The unit associated with the requirement.
        """
        return self._unit

    def getHelp(self):
        """
        Return the documentation string.

        Returns
        -------

        help : str
            The help string.
        """
        return self._help

    def isMulti(self):
        """
        Whether the requirement has multiple values.

        Returns
        -------

        is_multi : bool
            Whether the requirement has multiple values.
        """
        return self._is_multi

    def isOptional(self):
        """
        Whether the requirement is optional.

        Returns
        -------

        is_optional : bool
            Whether the requirement is optional.
        """
        return self._is_optional

    def getArgType(self):
        """
        The command-line argument type.

        Returns
        -------

        arg_type : bool, int, float, str
           The command-line argument type.
        """
        return self._arg_type

    def getMin(self):
        """Return the minimum allowed value."""
        return self._min

    def getMax(self):
        """Return the maximum allowed value."""
        return self._max

    def getAllowedValues(self):
        """
        Return the allowed values.

        Returns
        -------

        allowed : list
            The list of allowed values that the requirement can take.
        """
        return self._allowed

    def _validate_default(self):
        """Validate the default value."""

        if self._min is not None and self._default < self._min:
            raise ValueError(
                "The default '%s' is less than the minimum "
                "allowed value '%s'" % (self._default, self._min)
            )

        if self._max is not None and self._default > self._max:
            raise ValueError(
                "The default '%s' is greater than the maximum "
                "allowed value '%s'" % (self._default, self._max)
            )

        if self._allowed is not None and self._default not in self._allowed:
            raise ValueError(
                "The default '%s' is not one of the allowed values %s"
                % (self._default, str(self._allowed))
            )


class Boolean(Requirement):
    """A boolean requirement.

    Example
    -------

    Create a boolean flag with a default of False.

    >>> import BioSimSpace as BSS
    >>> flag = BSS.Gateway.Boolean(help="A boolean flag", default=False)
    """

    # Set the argparse argument type.
    _arg_type = bool

    def __init__(self, help=None, default=None):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : bool
            The default value.
        """

        # Call the base class constructor.
        super().__init__(help=help, default=default)

    def _validate(self, value):
        """Validate the value."""

        if isinstance(value, bool):
            return value
        else:
            raise TypeError("The value should be of type 'bool'.")


class Integer(Requirement):
    """An integer requirement.

    Examples
    --------

    Create an integer requirement with an allowed range and no default.

    >>> import BioSimSpace as BSS
    >>> my_int = BSS.Gateway.Integer(help="An integer requirement.", minimum=0, maximum=10)

    Create an integer requirement with a given set of allowed values.

    >>> import BioSimSpace as BSS
    >>> my_int = BSS.Gateway.Integer(help="An integer requirement.", allowed=[1,2,3,4,5])

    Create an integer requirement with a maximum value of 10 and default of 3.

    >>> import BioSimSpace as BSS
    >>> my_int = BSS.Gateway.Integer(help="An integer requirement.", default=3, maximum=10)
    """

    # Set the argparse argument type.
    _arg_type = int

    def __init__(
        self, help=None, default=None, minimum=None, maximum=None, allowed=None
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : int
            The default value.

        minimum : int
            The minimum allowed value.

        maximum : int
            The maximum allowed value.

        allowed : [int]
            A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if type(value) is int:
            return value
        else:
            raise TypeError("The value should be of type 'int'.")


class Float(Requirement):
    """A floating point requirement.

    Examples
    --------

    Create a float requirement with an allowed range and no default.

    >>> import BioSimSpace as BSS
    >>> my_float = BSS.Gateway.Float(help="A float requirement.", minimum=-13.2, maximum=27.3)

    Create a float requirement with a given set of allowed values.

    >>> import BioSimSpace as BSS
    >>> my_float = BSS.Gateway.Float(help="A float requirement.", allowed=[1.0,3.5,12.8])

    Create a float requirement with a maximum value of 57.3 and default of 18.2.

    >>> import BioSimSpace as BSS
    >>> my_float = BSS.Gateway.Float(help="A float requirement.", default=18.2, maximum=57.3)
    """

    # Set the argparse argument type.
    _arg_type = float

    def __init__(
        self, help=None, default=None, minimum=None, maximum=None, allowed=None
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : float
            The default value.

        minimum : float
            The minimum allowed value.

        maximum : float
            The maximum allowed value.

        allowed : [float]
            A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if isinstance(value, float):
            return value
        elif type(value) is int:
            return float(value)
        else:
            raise TypeError("The value should be of type 'float' or 'int'.")


class String(Requirement):
    """A string requirement.

    Examples
    --------

    Create a string requirement with a default value.

    >>> import BioSimSpace as BSS
    >>> my_string = BSS.Gateway.String(help="A string requirement.", default="dog")

    Create a string requirement with a list of allowed values and a default of "cat".

    >>> import BioSimSpace as BSS
    >>> my_string = BSS.Gateway.String(help="A string requirement.", allowed=["cat", "dog", "fish"], default="cat")
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, default=None, allowed=None):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : str
            The default value.

        allowed : [str]
            A list of allowed values.
        """

        # Call the base class constructor.
        super().__init__(help=help, default=default, allowed=allowed)

    def _validate(self, value):
        """Validate that the value is of the correct type."""

        if isinstance(value, str):
            return value
        else:
            raise TypeError("The value should be of type 'str'")


class File(Requirement):
    """A file requirement.

    Example
    -------

    Create an optional file requirement.

    >>> import BioSimSpace as BSS
    >>> my_file = BSS.Gateway.File(help="A file requirement.", optional=True)
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(self, help=None, optional=False):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        optional : bool
            Whether the file is optional.
        """

        # Call the base class constructor.
        super().__init__(help=help, optional=optional)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        import os as _os

        # Handle optional requirement.
        if self._is_optional and value is None:
            return None

        # Check the type.
        if isinstance(value, str):
            file = _unarchive(value)
            if file is None:
                file = value
            elif isinstance(file, (list, tuple)):
                raise ValueError(
                    "The archive contains multiple files: use a FileSet instead!"
                )
        else:
            raise TypeError("The value should be of type 'str'")

        # Make sure the file exists.
        if not _os.path.isfile(file):
            raise IOError("File doesn't exist: '%s'" % file)
        else:
            return file


class FileSet(Requirement):
    """A file requirement.

    Example
    -------

    Create a file set requirement.

    >>> import BioSimSpace as BSS
    >>> my_files = BSS.Gateway.FileSet(help="A file set requirement.")
    """

    # Set the argparse argument type.
    _arg_type = str

    # Multiple files can be passed.
    _is_multi = True

    def __init__(self, help=None, optional=False):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        optional : bool
            Whether the file set is optional.
        """

        # Call the base class constructor.
        super().__init__(help=help, optional=optional)

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : [str]
            A list of the files associated with this requirement.
        """
        if self._value is None:
            return None
        else:
            return self._value.copy()

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        import re as _re
        import os as _os

        # Handle optional requirement.
        if self._is_optional and value is None:
            return None

        # Handle single strings.
        if isinstance(value, str):
            value = [value]

        # The user can pass a list of compressed files so we need to keep
        # track of the names of the uncompressed files.
        uncompressed_files = []

        # We should receive a list of strings.
        if isinstance(value, (list, tuple)):
            # A single file was passed.
            if len(value) == 1:
                # Remove whitespace and split on commas.
                value = _re.sub(r"\s+", "", value[0]).split(",")

            # Loop over all strings.
            for file in value:
                # Check the types.
                if not isinstance(file, str):
                    raise TypeError("The value should be of type 'str'")

                # Check whether this is an archive.
                files = _unarchive(file)

                if files is not None:
                    if isinstance(files, list):
                        uncompressed_files += files
                    else:
                        uncompressed_files.append(files)
                else:
                    # Make sure the file exists.
                    if not _os.path.isfile(file):
                        raise IOError("File doesn't exist: '%s'" % file)

        if len(uncompressed_files) > 0:
            return uncompressed_files
        else:
            return value


class Length(Requirement):
    """A length requirement.

    Examples
    --------

    Create a length requirement with a default of 10 Angstrom.

    >>> import BioSimSpace as BSS
    >>> my_length = BSS.Gateway.Length(help="A length requirement", default=10, unit="angstrom")

    The same, but explicitly passing a :class:`Length <BioSimSpace.Types.Length>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_length = BSS.Gateway.Length(
    ...     help="A length requirement",
    ...     default=10*BSS.Units.Length.angstrom
    ... )

    Create a length requirement with a default of 10 Angstrom and a maximum
    of 50 nanometers. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_length = BSS.Gateway.Length(
    ...     help="A length requirement",
    ...     default=10*BSS.Units.Length.angstrom,
    ...     maximum=50*BSS.Units.Length.nanometer
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Length <BioSimSpace.Types.Length>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Length <BioSimSpace.Types.Length>`
            The minimum allowed value.

        maximum : :class:`Length <BioSimSpace.Types.Length>`
            The maximum allowed value.

        allowed : [:class:`Length <BioSimSpace.Types.Length>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            length = _Types.Length("1 %s" % unit)
            self._unit = length.unit()
            self._print_unit = length._print_format[length.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Length <BioSimSpace.Types.Length>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Length):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "length")

            if unit is None:
                return _Types.Length(value, self._unit)
            else:
                return _Types.Length(value, unit)._convert_to(self._unit)


class Area(Requirement):
    """An area requirement.

    Examples
    --------

    Create an area requirement with a default of 10 square Angstrom.

    >>> import BioSimSpace as BSS
    >>> my_area = BSS.Gateway.Area(help="An area requirement", default=10, unit="angstrom2")

    The same, but explicitly passing a :class:`Area <BioSimSpace.Types.Area>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_area = BSS.Gateway.Area(help="An area requirement", default=10*BSS.Units.Area.angstrom2)

    Create an area requirement with a default of 10 square Angstrom and a maximum
    of 50 square nanometers. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_area = BSS.Gateway.Area(
    ...     help="An area requirement",
    ...     default=100*BSS.Units.Area.angstrom2,
    ...     maximum=50*BSS.Units.Area.nanometer2
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Area <BioSimSpace.Types.Area>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Area <BioSimSpace.Types.Area>`
            The minimum allowed value.

        maximum : :class:`Area <BioSimSpace.Types.Area>`
            The maximum allowed value.

        allowed : [:class:`Area <BioSimSpace.Types.Area>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            area = _Types.Area("1 %s" % unit)
            self._unit = area.unit()
            self._print_unit = area._print_format[area.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Area <BioSimSpace.Types.Area>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Area):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "area")

            if unit is None:
                return _Types.Area(value, self._unit)
            else:
                return _Types.Area(value, unit)._convert_to(self._unit)


class Volume(Requirement):
    """A volume requirement.

    Examples
    --------

    Create a volume requirement with a default of 10 cubed Angstrom.

    >>> import BioSimSpace as BSS
    >>> my_volume = BSS.Gateway.Volume(help="A volume requirement", default=10, unit="angstrom3")

    The same, but explicitly passing a :class:`Volume <BioSimSpace.Types.Volume>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_volume = BSS.Gateway.Volume(help="A volume requirement", default=10*BSS.Units.Volume.angstrom3)

    Create a volume requirement with a default of 10 cubed Angstrom and a maximum
    of 50 cubed nanometers. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_volume = BSS.Gateway.Volume(
    ...     help="A volume requirement",
    ...     default=10*BSS.Units.Volume.angstrom3,
    ...     maximum=50*BSS.Units.Volume.nanometer3
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Volume <BioSimSpace.Types.Volume>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Volume <BioSimSpace.Types.Volume>`
            The minimum allowed value.

        maximum : :class:`Volume <BioSimSpace.Types.Volume>`
            The maximum allowed value.

        allowed : [:class:`Volume <BioSimSpace.Types.Volume>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            volume = _Types.Volume("1 %s" % unit)
            self._unit = volume.unit()
            self._print_unit = volume._print_format[volume.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Volume <BioSimSpace.Types.Volume>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Volume):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "volume")

            if unit is None:
                return _Types.Volume(value, self._unit)
            else:
                return _Types.Volume(value, unit)._convert_to(self._unit)


class Angle(Requirement):
    """An angle requirement.

    Examples
    --------

    Create an angle requirement with a default of 3.14 radians.

    >>> import BioSimSpace as BSS
    >>> my_angle = BSS.Gateway.Angle(help="An angle requirement", default=3.14, unit="radian")

    The same, but explicitly passing a :class:`Angle <BioSimSpace.Types.Angle>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_angle = BSS.Gateway.Angle(help="An angle requirement", default=3.14*BSS.Units.Angle.radian)

    Create an angle requirement with a default of 3.14 radian and a maximum
    of 360 degrees. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_angle = BSS.Gateway.Angle(
    ...     help="An angle requirement",
    ...     default=3.14*BSS.Units.Angle.radian,
    ...     maximum=360*BSS.Units.Angle.degree
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Angle <BioSimSpace.Types.Angle>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Angle <BioSimSpace.Types.Angle>`
            The minimum allowed value.

        maximum : :class:`Angle <BioSimSpace.Types.Angle>`
            The maximum allowed value.

        allowed : [:class:`Angle <BioSimSpace.Types.Angle>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            angle = _Types.Angle("1 %s" % unit)
            self._unit = angle.unit()
            self._print_unit = angle._print_format[angle.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Angle <BioSimSpace.Types.Angle>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Angle):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "angle")

            if unit is None:
                return _Types.Angle(value, self._unit)
            else:
                return _Types.Angle(value, unit)._convert_to(self._unit)


class Charge(Requirement):
    """A charge requirement.

    Examples
    --------

    Create a charge requirement with a default of 3 electron charge.

    >>> import BioSimSpace as BSS
    >>> my_charge = BSS.Gateway.Charge(help="A charge requirement", default=3, unit="electron charge")

    The same, but explicitly passing a :class:`Charge <BioSimSpace.Types.Charge>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_charge = BSS.Gateway.Charge(help="A charge requirement", default=3*BSS.Units.Charge.electron_charge)

    Create a charge requirement with a default of 3 electron charge and a
    maximum of -10 Coulomb. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_charge = BSS.Gateway.Charge(
    ...     help="A charge requirement",
    ...     default=3*BSS.Units.Charge.electron_charge,
    ...     maximum=10*BSS.Units.Charge.coulomb
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Charge <BioSimSpace.Types.Charge>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Charge <BioSimSpace.Types.Charge>`
            The minimum allowed value.

        maximum : :class:`Charge <BioSimSpace.Types.Charge>`
            The maximum allowed value.

        allowed : [:class:`Charge <BioSimSpace.Types.Charge>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            nrg = _Types.Charge("1 %s" % unit)
            self._unit = nrg.unit()
            self._print_unit = nrg._print_format[nrg.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Charge <BioSimSpace.Types.Charge>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Charge):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "energy")

            if unit is None:
                return _Types.Charge(value, self._unit)
            else:
                return _Types.Charge(value, unit)._convert_to(self._unit)


class Energy(Requirement):
    """An energy requirement.

    Examples
    --------

    Create an energy requirement with a default of 3 kcal per mol.

    >>> import BioSimSpace as BSS
    >>> my_energy = BSS.Gateway.Energy(help="An energy requirement", default=3, unit="kcal per mol")

    The same, but explicitly passing a :class:`Energy <BioSimSpace.Types.Energy>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_energy = BSS.Gateway.Energy(help="An energy requirement", default=3*BSS.Units.Energy.kcal_per_mol)

    Create an energy requirement with a default of 3 kcal per mol and a
    maximum of 50 kJ per mol. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_energy = BSS.Gateway.Energy(
    ...     help="An energy requirement",
    ...     default=3*BSS.Units.Energy.kcal_per_mol,
    ...     maximum=50*BSS.Units.Energy.kj_per_mol
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Energy <BioSimSpace.Types.Energy>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Energy <BioSimSpace.Types.Energy>`
            The minimum allowed value.

        maximum : :class:`Energy <BioSimSpace.Types.Energy>`
            The maximum allowed value.

        allowed : [:class:`Energy <BioSimSpace.Types.Energy>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            nrg = _Types.Energy("1 %s" % unit)
            self._unit = nrg.unit()
            self._print_unit = nrg._print_format[nrg.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Energy <BioSimSpace.Types.Energy>`
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Energy):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "energy")

            if unit is None:
                return _Types.Energy(value, self._unit)
            else:
                return _Types.Energy(value, unit)._convert_to(self._unit)


class Mass(Requirement):
    """A mass requirement.

    Examples
    --------

    Create a mass requirement with a default of 10 grams.

    >>> import BioSimSpace as BSS
    >>> my_mass = BSS.Gateway.Mass(help="A mass requirement with a default of 10 grams", default=10, unit="gram")

    The same, but explicitly passing a :class:`Mass <BioSimSpace.Types.Mass>`

    >>> import BioSimSpace as BSS
    >>> my_mass = BSS.Gateway.Mass(help="A mass requirement with a default of 10 grams", default=10*BSS.Units.Mass.gram)

    Create a mass requirement with a default of 10 grams and a maximum of 10 kilograms.
    Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_mass = BSS.Gateway.Mass(
    ...     help="A mass requirement with a default of 10 grams",
    ...     default=10*BSS.Units.Mass.gram,
    ...     maximum=10*BSS.Units.Mass.kilogram
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Mass <BioSimSpace.Types.Mass>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Mass <BioSimSpace.Types.Mass>`
            The minimum allowed value.

        maximum : :class:`Mass <BioSimSpace.Types.Mass>`
            The maximum allowed value.

        allowed : [:class:`Mass <BioSimSpace.Types.Mass>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            mass = _Types.Mass("1 %s" % unit)
            self._unit = mass.unit()
            self._print_unit = mass._print_format[mass.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Mass <BioSimSpace.Types.Mass>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Mass):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "pressure")

            if unit is None:
                return _Types.Mass(value, self._unit)
            else:
                return _Types.Mass(value, unit)._convert_to(self._unit)


class Pressure(Requirement):
    """A pressure requirement.

    Examples
    --------

    Create a pressure requirement with a default of 1 atmosphere.

    >>> import BioSimSpace as BSS
    >>> my_pressure = BSS.Gateway.Pressure(help="A pressure requirement", default=1, unit="atm")

    The same, but explicitly passing a :class:`Pressure <BioSimSpace.Types.Pressure>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_pressure = BSS.Gateway.Pressure(help="A pressure requirement", default=BSS.Units.Pressure.atm)

    Create a pressure requirement with a default of 1 atomosphere and a
    maximum of 10 bar. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_pressure = BSS.Gateway.Pressure(
    ...     help="A pressure requirement",
    ...     default=BSS.Units.Pressure.atm,
    ...     maximum=10*BSS.Units.Pressure.bar
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The minimum allowed value.

        maximum : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The maximum allowed value.

        allowed : [:class:`Pressure <BioSimSpace.Types.Pressure>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            press = _Types.Pressure("1 %s" % unit)
            self._unit = press.unit()
            self._print_unit = press._print_format[press.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Pressure):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "pressure")

            if unit is None:
                return _Types.Pressure(value, self._unit)
            else:
                return _Types.Pressure(value, unit)._convert_to(self._unit)


class Temperature(Requirement):
    """A temperature requirement.

    Examples
    --------

    Create a temperature requirement with a default of 300 kelvin.

    >>> import BioSimSpace as BSS
    >>> my_temperature = BSS.Gateway.Temperature(help="A temperature requirement", default=300, unit="kelvin")

    The same, but explicitly passing a :class:`Temperature <BioSimSpace.Types.Temperature>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_temperature = BSS.Gateway.Temperature(help="A temperature requirement", default=300*BSS.Units.Temperature.kelvin)

    Create a temperature requirement with a default of 300 Kelvin and a
    maximum of 100 Celsius. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_temperature = BSS.Gateway.Temperature(
    ...     help="A temperature requirement",
    ...     default=300*BSS.Units.Temperature.kelvin,
    ...     maximum=100*BSS.Units.Temperature.celsius
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The minimum allowed value.

        maximum : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The maximum allowed value.

        allowed : [:class:`Temperature <BioSimSpace.Types.Temperature>`]
            A list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            temp = _Types.Temperature("1 %s" % unit)
            self._unit = temp.unit()
            self._print_unit = temp._print_format[temp.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The value of the requirement.
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Temperature):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "temperature")

            if unit is None:
                return _Types.Temperature(value, self._unit)
            else:
                return _Types.Temperature(value, unit)._convert_to(self._unit)


class Time(Requirement):
    """A time requirement.

    Examples
    --------

    Create a time requirement with a default of 35 minutes.

    >>> import BioSimSpace as BSS
    >>> my_time = BSS.Gateway.Time(help="A time requirement", default=35, unit="minutes")

    The same, but explicitly passing a :class:`Time <BioSimSpace.Types.Time>`
    for the default.

    >>> import BioSimSpace as BSS
    >>> my_time = BSS.Gateway.Time(help="A time requirement", default=35*BSS.Units.Time.minute)

    Create a time requirement with a default of 35 minutes and a maximum
    of 5 hours. Note that the unit is taken from the default value.

    >>> import BioSimSpace as BSS
    >>> my_time = BSS.Gateway.Time(
    ...     help="A time requirement",
    ...     default=35*BSS.Units.Time.minute,
    ...     maximum=5*BSS.Units.Time.hour
    ... )
    """

    # Set the argparse argument type.
    _arg_type = str

    def __init__(
        self,
        help=None,
        default=None,
        unit=None,
        minimum=None,
        maximum=None,
        allowed=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        help : str
            The help string.

        default : :class:`Time <BioSimSpace.Types.Time>`
            The default value.

        unit : str
            The unit.

        minimum : :class:`Time <BioSimSpace.Types.Time>`
            The minimum allowed value.

        maximum : :class:`Time <BioSimSpace.Types.Time>`
            The maximum allowed value.

        allowed : [:class:`Time <BioSimSpace.Types.Time>`]
            The list of allowed values.
        """
        from .. import Types as _Types

        # Validate the unit.
        if unit is not None:
            time = _Types.Time("1 %s" % unit)
            self._unit = time.unit()
            self._print_unit = time._print_format[time.unit()]
        else:
            try:
                self._unit = default.unit()
            except:
                raise ValueError("No unit or default value has been specified!")

        # Call the base class constructor.
        super().__init__(
            help=help,
            default=default,
            unit=self._unit,
            minimum=minimum,
            maximum=maximum,
            allowed=allowed,
        )

    def getValue(self):
        """
        Return the value.

        Returns
        -------

        value : :class:`Time <BioSimSpace.Types.Time>`
        """
        import copy as _copy

        if self._value is None:
            return None
        else:
            return _copy.deepcopy(self._value)

    def _validate(self, value):
        """Validate that the value is of the correct type."""
        from .. import Types as _Types

        if isinstance(value, _Types.Time):
            return value._convert_to(self._unit)

        else:
            # Extract the value and unit from the argument string.
            value, unit = _validate_unit_requirement(value, "time")

            if unit is None:
                return _Types.Time(value, self._unit)
            else:
                return _Types.Time(value, unit)._convert_to(self._unit)


def _validate_unit_requirement(value, unit_type):
    """
    Helper function to validate input requirements with units.

    Parameters
    ----------

    value : str
        The value of the input requirement.

    unit_type : str
        The unit type.

     Returns
     -------

     (value, unit) : tuple
         The value and unit of the requirement.
    """
    import re as _re

    # No unit by default.
    unit = None

    # Float.
    if isinstance(value, float):
        pass

    # Integer.
    elif type(value) is int:
        value = float(value)

    # String.
    elif isinstance(value, str):
        # First try to directly convert to a float.
        try:
            value = float(value)

        # Use a regular expression to extract the value and unit.
        except ValueError:
            # Strip white space from the string.
            value = value.replace(" ", "")

            # Try to match scientific format.
            match = _re.search(r"(\-?\d+\.?\d*e\-?\d+)(.*)", value, _re.IGNORECASE)

            # Try to match decimal format.
            if match is None:
                match = _re.search(r"(\-?\d+\.?\d*)(.*)", value, _re.IGNORECASE)

                # No matches, raise an error.
                if match is None:
                    raise ValueError(
                        "Could not interpret %s: '%s'" % (unit_type, value)
                    )

            # Extract the value and unit.
            value, unit = match.groups()

            # Convert the value to a float.
            value = float(value)

    # Unsupported.
    else:
        raise TypeError(
            "Unsupported value type '%s'. Options are 'float', 'int', or 'str'."
            % type(value)
        )

    return (value, unit)


def _unarchive(name):
    """
    Decompress an archive and return a list of files.

    Parameters
    ----------

    name : str
        The name of the archive (full path).

    Returns
    -------

    files : [ str ]
        A list of file names.
    """
    import tarfile as _tarfile
    import gzip as _gzip
    import os as _os
    import bz2 as _bz2
    import zipfile as _zipfile
    import shutil as _shutil

    # Get the directory name.
    dir = _os.path.dirname(name)

    # If the file compressed file has been passed on the command-line, then
    # we'll extract it to a directory to avoid littering the current workspace.
    if dir == "uploads":
        dir += "/"
    else:
        dir = _os.path.splitext(name)[0] + "/"

    # List of supported tar file formats.
    tarfiles = ["tar.gz", "tar.bz2", "tar"]

    # Check whether this is a tar compressed file.
    for tar_name in tarfiles:
        # Found a match.
        if tar_name in name.lower():
            # The list of decompressed files.
            files = []

            # Decompress the archive.
            with _tarfile.open(name) as tar:
                # We need to call tar.list(), otherwise the tar object will not know
                # about nested directories, i.e. it will appear as if there is a single
                # member.
                print("Decompressing...")
                tar.list(verbose=False)

                # Loop over all of the members and get the file names.
                # If the name has no extension, then we assume that it's a directory.
                for file in tar.members:
                    if _os.path.splitext(file.name)[1] != "":
                        files.append(dir + file.name)

                # Now extract all of the files.
                tar.extractall(dir)

            return files

    # Get the file name and extension.
    file, ext = _os.path.splitext(name)

    # This is a zip file.
    if ext.lower() == ".zip":
        # The list of decompressed files.
        files = None

        with _zipfile.ZipFile(name) as zip:
            files = zip.namelist()
            print("Decompressing...")
            for file in files:
                print(file)
            zip.extractall(dir)

        return [dir + file for file in files]

    # This is a gzip file.
    if ext.lower() == ".gz" or ext == ".gzip":
        with _gzip.open(name, "rb") as f_in:
            with open(file, "wb") as f_out:
                print("Decompressing...\n%s" % name)
                _shutil.copyfileobj(f_in, f_out)

        return file

    # This is a bzip2 file.
    if ext.lower() == ".bz2" or ext == ".bzip2":
        with _bz2.open(name, "rb") as f_in:
            with open(file, "wb") as f_out:
                print("Decompressing...\n%s" % name)
                _shutil.copyfileobj(f_in, f_out)

        return file

    # If we get this far, then this is not a supported archive.
    return None
