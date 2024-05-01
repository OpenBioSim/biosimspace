######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2024
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

"""Functionality for creating BioSimSpace workflow components (nodes)."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Node"]

from .._Utils import _try_import

import configargparse as _argparse
import collections as _collections
import __main__
import os as _os
import shutil as _shutil
import sys as _sys
import textwrap as _textwrap
import warnings as _warnings

_yaml = _try_import("yaml")

from .. import _is_notebook
from .. import setVerbose

# Enable Jupyter widgets.
if _is_notebook:
    from IPython.display import FileLink as _FileLink

    import ipywidgets as _widgets
    import zipfile as _zipfile

from ..Types._type import Type as _Type

from ._requirements import Area as _Area
from ._requirements import Boolean as _Boolean
from ._requirements import File as _File
from ._requirements import FileSet as _FileSet
from ._requirements import Float as _Float
from ._requirements import Charge as _Charge
from ._requirements import Energy as _Energy
from ._requirements import Integer as _Integer
from ._requirements import Length as _Length
from ._requirements import Pressure as _Pressure
from ._requirements import Requirement as _Requirement
from ._requirements import String as _String
from ._requirements import Temperature as _Temperature
from ._requirements import Time as _Time
from ._requirements import Volume as _Volume

# Float types (including those with units).
_float_types = [
    _Float,
    _Charge,
    _Energy,
    _Pressure,
    _Length,
    _Area,
    _Volume,
    _Temperature,
    _Time,
]

# Unit types.
_unit_types = [
    _Charge,
    _Energy,
    _Pressure,
    _Length,
    _Area,
    _Volume,
    _Temperature,
    _Time,
]


class Parser(_argparse.ArgumentParser):
    # Pass the message straight through to the exit method.
    def error(self, message):
        return self.exit(status=1, message=message)

    # Print the full help text, then the message, then exit.
    def exit(self, status=0, message=None):
        if message is not None:
            self.print_help()
            print("\nArgument parser failed with the following message:")
            message = "   " + message + "\n"
        return super().exit(status, message)


class CwlAction(_argparse.Action):
    """Helper class to export CWL wrappers from Node metadata."""

    @classmethod
    def bind_node(cls, node):
        """Bind the inputs and outputs of a node to this action."""
        cls.inputs = node._inputs
        cls.outputs = node._outputs

    def __call__(self, parser, namespace, values, option_string=None):
        """Export the CWL wrapper."""

        if values == False:
            parser.exit()
            return

        for value in self.outputs.values():
            # Currently we only support File and FileSet output
            # requirements with CWL.
            output_type = type(value)
            if output_type not in [_File, _FileSet]:
                raise TypeError(
                    "We currently only support File and " "FileSet outputs with CWL."
                )

        # Store the absolute path of the Python interpreter used to run the node.
        exe = _sys.executable

        # Store the absolute path of the node.
        import __main__

        node = _os.path.abspath(__main__.__file__)

        # Create the name of the CWL wrapper.
        cwl_wrapper = __main__.__file__.replace(".py", ".cwl")

        # Write the wrapper.
        with open(cwl_wrapper, "w") as file:
            # Write the header.
            file.write("cwlVersion: v1.0\n")
            file.write("class: CommandLineTool\n")
            file.write(f'baseCommand: ["{exe}", "{node}", "--strict-file-naming"]\n')

            # Write the inputs section.
            file.write("\n")
            file.write("inputs:\n")
            for key, value in self.inputs.items():
                file.write(f"  {key}:\n")

                # Map the requirement to the appropriate CWL type.

                if isinstance(value, _Boolean):
                    cwl_type = "bool"

                elif isinstance(value, _Integer):
                    cwl_type = "int"

                elif isinstance(value, _Float):
                    cwl_type = "float"

                elif isinstance(value, _String):
                    cwl_type = "string"

                elif isinstance(value, _File):
                    cwl_type = "File"

                elif isinstance(value, _FileSet):
                    cwl_type = "array"

                # Use a string for unit-based types since it gives
                # the user greatest flexibility in expressing the input.
                if type(value) in _unit_types:
                    cwl_type = "string"

                # Handle FileSet types separately.
                if cwl_type == "array":
                    file.write("    type:\n")
                    if value.isOptional():
                        file.write('      - "null"\n')
                    file.write("      - type: array\n")
                    file.write("        items: File\n")

                # Handle optional values.
                else:
                    if value.isOptional():
                        cwl_type += "?"
                    file.write(f"    type: {cwl_type}\n")

                # Handle default values.
                default = value.getDefault()
                if default is not None:
                    if type(value) in _unit_types:
                        value = default.value()
                        unit = default.unit()
                        unit = unit.lower()
                        file.write(f"    default: {value} {unit}\n")
                    else:
                        file.write(f"    default: {default}\n")

                # Bind the command-line option name.
                file.write("    inputBinding:\n")
                file.write(f"      prefix: --{key}\n")
                file.write("      separate: true\n")

                file.write("\n")

            # Write the outputs section.
            if len(self.outputs) == 0:
                file.write("outputs: []\n")
            else:
                file.write("outputs:\n")
                for key, value in self.outputs.items():
                    output_type = type(value)
                    file.write(f"  {key}:\n")

                    # Only support File and FileSet for now. This has been
                    # validated at the start of the __call__ method, but we
                    # include an if/elif/else conditional block so that we
                    # can support additional types in future. Note that we
                    # use glob to bind the output, so the prefix used to name
                    # files must match the key used to define the requirement.

                    # File.
                    if output_type is _File:
                        file.write("    type: File\n")
                        file.write("    outputBinding:\n")
                        file.write(f'      glob: "{key}.*"\n')

                    # FileSet.
                    elif output_type is _FileSet:
                        file.write("    type:\n")
                        file.write("      type: array\n")
                        file.write("      items: File\n")
                        file.write("    outputBinding:\n")
                        file.write(f'      glob: "{key}.*"\n')

        # Exit the parser.
        parser.exit()


class Node:
    """A class for interfacing with BioSimSpace nodes.

    Nodes are used to collect and validate user input, document the
    intentions of the workflow component, track and report errors,
    and validate output. Once written, a node can be run from within
    Jupyter, from the command-line, or plugged into a workflow engine,
    such as Knime.

    Example
    -------

    A generic energy minimisation node:

    >>> import BioSimSpace as BSS
    >>> node = BSS.Gateway.Node("Perform energy minimisation")
    >>> node.addAuthor(name="Lester Hedges", email="lester@openbiosim.org", affiliation="OpenBioSim")
    >>> node.setLicence("GPLv3")
    >>> node.addInput("files", BSS.Gateway.FileSet(help="A set of molecular input files."))
    >>> node.addInput("steps", BSS.Gateway.Integer(help="The number of minimisation steps.", minimum=0, maximum=100000, default=10000))
    >>> node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system."))
    >>> node.showControls()
    >>> system = BSS.IO.readMolecules(node.getInput("files"))
    >>> protocol = BSS.Protocol.Minimisation(steps=node.getInput("steps"))
    >>> process = BSS.MD.run(system, protocol)
    >>> node.setOutput("minimised", BSS.IO.saveMolecules("minimised", process.getSystem(block=True), system.fileFormat()))
    >>> node.validate()
    """

    # Whether the node is run from Knime.
    _is_knime = False

    # Whether the node is run from a Jupyter notebook.
    _is_notebook = _is_notebook

    def __init__(self, description, name=None):
        """
        Constructor.

        Parameters
        ----------

        description : str
            A description of the node.

        name : str
            The name of the node.
        """

        if not isinstance(description, str):
            raise TypeError("The 'description' keyword must be of type 'str'.")

        # Set the node name.
        if name is None:
            try:
                self._name = _os.path.basename(__main__.__file__)
            except:
                self._name = None
        else:
            if not isinstance(name, str):
                raise TypeError("The 'name' keyword must be of type 'str'.")
            self._name = name

        # Set the node description string.
        self._description = description

        # Initialise the authors.
        self._authors = None

        # Initialise the license.
        self._license = None

        # Initialise dictionaries for the inputs/outputs.
        self._inputs = _collections.OrderedDict()
        self._outputs = _collections.OrderedDict()

        # A dictionary of Jupyter widgets.
        self._widgets = _collections.OrderedDict()

        # Whether the input is valid.
        self._is_valid_input = False

        # Whether the output has been validated.
        self._is_output_validated = False

        # A list of user error messages.
        self._errors = []

        # Initialise the parser.
        self._parser = None
        self._required = None

        # Initialise the Jupyter input panel.
        self._control_panel = None

        # Strict file naming is off by default.
        self._strict_file_naming = False

        # Running from the command-line.
        if not self._is_knime and not self._is_notebook:
            # Generate the node help description.
            description = self._generate_description()

            # Create the parser.
            self._parser = Parser(
                description=description,
                formatter_class=_argparse.RawTextHelpFormatter,
                add_help=False,
                config_file_parser_class=_argparse.YAMLConfigFileParser,
                add_config_file_help=False,
            )

            # Bind the node inputs and outputs to the CWL action.
            CwlAction.bind_node(self)

            # Add argument groups.
            self._required = self._parser.add_argument_group("Required arguments")
            self._optional = self._parser.add_argument_group("Optional arguments")
            self._optional.add_argument(
                "-h", "--help", action="help", help="Show this help message and exit."
            )
            self._optional.add_argument(
                "-c",
                "--config",
                is_config_file=True,
                help="Path to configuration file.",
            )
            self._optional.add_argument(
                "-v",
                "--verbose",
                type=_str2bool,
                nargs="?",
                const=True,
                default=False,
                help="Print verbose error messages.",
            )
            self._optional.add_argument(
                "--export-cwl",
                action=CwlAction,
                type=_str2bool,
                nargs="?",
                const=True,
                default=False,
                help="Export Common Workflow Language (CWL) wrapper and exit.",
            )
            self._optional.add_argument(
                "--strict-file-naming",
                type=_str2bool,
                nargs="?",
                const=True,
                default=False,
                help="Enforce that the prefix of any file based output matches its name.",
            )

            # Overload the "_check_value" method for more flexible string support.
            # (Ignore whitespace and case insensitive.)
            self._parser._check_value = _check_value

    def __del__(self):
        """Destructor."""

        # Validate the node if the user hasn't already done so.
        if self._is_valid_input:
            if not self._is_output_validated:
                self.validate()

    def addInput(self, name, input):
        """
        Add an input requirement.

        Parameters
        ----------

        name : str
            The name of the input.

        input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
            The input requirement object.
        """

        if not isinstance(name, str):
            raise TypeError("'name' must be of type 'str'.")

        if not isinstance(input, _Requirement):
            raise TypeError("'input' must be of type 'Requirement'.")

        # We already have an input with this name.
        reset = False
        if name in self._inputs:
            if self._is_notebook:
                _warnings.warn("Duplicate input requirement '%s'" % name)
                reset = True
            else:
                raise ValueError("Duplicate input requirement '%s'" % name)

        # Add the input to the dictionary.
        self._inputs[name] = input

        # Create a Knime GUI widget.
        if self._is_knime:
            self._addInputKnime(name, input)

        # Create a Jupyter GUI widget.
        elif self._is_notebook:
            return self._addInputJupyter(name, input, reset)

        # Command-line argparse ArgumentParser.
        else:
            self._addInputCommandLine(name, input)

    def _addInputCommandLine(self, name, input):
        """
        Add an input requirement for the command-line.

        Parameters
        ----------

        name : str
            The name of the input.

        input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
            The input requirement object.
        """

        # Append long-form argument name if not present.
        if len(name) > 2:
            if name[0:2] != "--":
                name = "--" + name
        else:
            name = "--" + name

        if input.isOptional():
            if input.getDefault() is not None:
                if input.isMulti() is not False:
                    self._optional.add_argument(
                        name,
                        type=input.getArgType(),
                        nargs="+",
                        help=self._create_help_string(input),
                        default=input.getDefault(),
                    )
                else:
                    if isinstance(input.getArgType(), bool):
                        self._optional.add_argument(
                            name,
                            type=_str2bool,
                            nargs="?",
                            const=True,
                            default=input.getDefault(),
                            help=self._create_help_string(input),
                        )
                    else:
                        if input.getAllowedValues() is not None:
                            self._optional.add_argument(
                                name,
                                type=input.getArgType(),
                                help=self._create_help_string(input),
                                default=input.getDefault(),
                                choices=input.getAllowedValues(),
                            )
                        else:
                            self._optional.add_argument(
                                name,
                                type=input.getArgType(),
                                help=self._create_help_string(input),
                                default=input.getDefault(),
                            )
            else:
                if input.isMulti() is not False:
                    self._optional.add_argument(
                        name,
                        type=input.getArgType(),
                        nargs="+",
                        help=self._create_help_string(input),
                    )
                else:
                    self._optional.add_argument(
                        name,
                        type=input.getArgType(),
                        help=self._create_help_string(input),
                    )
        else:
            if input.isMulti() is not False:
                self._required.add_argument(
                    name,
                    type=input.getArgType(),
                    nargs="+",
                    help=self._create_help_string(input),
                    required=True,
                )
            else:
                if input.getAllowedValues() is not None:
                    self._required.add_argument(
                        name,
                        type=input.getArgType(),
                        help=self._create_help_string(input),
                        required=True,
                        choices=input.getAllowedValues(),
                    )
                else:
                    if isinstance(input.getArgType(), bool):
                        self._required.add_argument(
                            name,
                            type=_str2bool,
                            nargs="?",
                            const=True,
                            help=self._create_help_string(input),
                        )
                    else:
                        self._required.add_argument(
                            name,
                            type=input.getArgType(),
                            help=self._create_help_string(input),
                            required=True,
                        )

    def _addInputKnime(self, name, input):
        """
        Add an input requirement for Knime.

        Parameters
        ----------

        name : str
            The name of the input.

        input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
            The input requirement object.
        """
        return None

    def _addInputJupyter(self, name, input, reset=False):
        """
        Add an input requirement for Jupyter.

        Parameters
        ----------

        name : str
            The name of the input.

        input : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
            The input requirement object.

        reset : bool
            Whether to reset the widget data.
        """

        # Create a widget button to indicate whether the requirement value
        # has been set.
        button = _widgets.Button(
            tooltip="The input requirement is unset.",
            button_style="warning",
            icon="fa-exclamation-triangle",
            layout=_widgets.Layout(flex="1 1 auto", width="auto"),
            disabled=False,
        )

        # Add a Jupyter widget for each of the supported requirement types.

        # Boolean.
        if isinstance(input, _Boolean):
            # Create a Jupyter toggle button.
            widget = _widgets.ToggleButton(
                value=False,
                description=name,
                tooltip=input.getHelp(),
                button_style="",
                disabled=False,
            )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Get the default value.
            default = input.getDefault()

            if default is not None:
                widget.value = default
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Integer.
        elif isinstance(input, _Integer):
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = _widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False,
                )

            else:
                # Get the range of the input.
                _min = input.getMin()
                _max = input.getMax()

                # Whether the integer is unbounded.
                is_unbounded = True

                if _min is not None:
                    # Set the default.
                    if default is None:
                        default = _min

                    # Bounded integer.
                    if _max is not None:
                        step = int((_max - _min) / 100)
                        # Create an int slider widget.
                        widget = _widgets.IntSlider(
                            value=default,
                            min=_min,
                            max=_max,
                            step=step,
                            description=name,
                            tooltip=input.getHelp(),
                            continuous_update=False,
                            orientation="horizontal",
                            readout=True,
                            readout_format="d",
                            disabled=False,
                        )

                        # Flag that the integer is bounded.
                        is_unbounded = False

                # Unbounded integer.
                if is_unbounded:
                    # Create an integer widget.
                    widget = _widgets.IntText(
                        value=default,
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False,
                    )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is not None:
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Float types (including those with units).
        elif type(input) in _float_types:
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            # Get the value of types with units.
            if isinstance(default, _Type):
                default = default.value()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                    # Get the value of types with units.
                    if isinstance(default, _Type):
                        default = default.value()

                # Create a dropdown for the list of allowed values.
                widget = _widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False,
                )

            else:
                # Get the range of the input.
                _min = input.getMin()
                _max = input.getMax()

                # Get the value of types with units.
                if isinstance(_min, _Type):
                    _min = _min.value()
                if isinstance(_max, _Type):
                    _max = _max.value()

                # Whether the float is unbounded.
                is_unbounded = True

                if _min is not None:
                    # Set the default.
                    if default is None:
                        default = _min

                    # Bounded float.
                    if _max is not None:
                        step = (_max - _min) / 100
                        # Create a float slider widget.
                        widget = _widgets.FloatSlider(
                            value=default,
                            min=_min,
                            max=_max,
                            step=step,
                            description=name,
                            tooltip=input.getHelp(),
                            continuous_update=False,
                            orientation="horizontal",
                            readout=True,
                            readout_format=".2f",
                            disabled=False,
                        )

                        # Flag that the float is bounded.
                        is_unbounded = False

                # Unbounded float.
                if is_unbounded:
                    # Create a float widget.
                    widget = _widgets.FloatText(
                        value=default,
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False,
                    )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is not None:
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # String.
        elif isinstance(input, _String):
            # Get the list of allowed values.
            allowed = input.getAllowedValues()

            # Get the default value.
            default = input.getDefault()

            if allowed is not None:
                # Set the default.
                if default is None:
                    default = allowed[0]

                # Create a dropdown for the list of allowed values.
                widget = _widgets.Dropdown(
                    options=allowed,
                    value=default,
                    description=name,
                    tooltip=input.getHelp(),
                    disabled=False,
                )

            else:
                if default is None:
                    # Create a text widget without a default.
                    widget = _widgets.Text(
                        placeholder="Type something",
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False,
                    )
                else:
                    # Create a text widget.
                    widget = _widgets.Text(
                        value=default,
                        placeholder="Type something",
                        description=name,
                        tooltip=input.getHelp(),
                        disabled=False,
                    )

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Add an attribute to flag whether the widget value has
            # been set by the user.
            if input.getDefault() is not None:
                widget._is_set = True
                widget._button.tooltip = "The input requirement is set."
                widget._button.button_style = "success"
                widget._button.icon = "fa-check"

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_value_change, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # File / File set.
        elif isinstance(input, (_File, _FileSet)):
            # Create a fileupload widget.
            if isinstance(input, _FileSet):
                widget = _widgets.FileUpload(multiple=True)
            else:
                widget = _widgets.FileUpload(multiple=False)

            # Make the widget dynamically resize to the content.
            widget.layout = {"width": "max-content"}

            # Add the 'set' indicator button to the widget.
            widget._button = button

            # Flag that the widget is unset.
            widget._is_set = False

            # Flag that this widget references files.
            widget._is_file = True

            # Store the requirement name.
            widget._name = name

            # Bind the callback function.
            widget.observe(_on_file_upload, names="value")

            # Store the widget.
            self._widgets[name] = widget

        # Unsupported input.
        else:
            raise ValueError("Unsupported requirement type '%s'" % type(input))

    def addOutput(self, name, output):
        """
        Add an output requirement.

        Parameters
        ----------

        name : str
            The name of the output.

        output : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>`
            The output requirement object.
        """

        if not isinstance(name, str):
            raise TypeError("'name' must be of type 'str'.")

        if not isinstance(output, _Requirement):
            raise TypeError("'output' must be of type 'Requirement'.")

        # We already have an output requirement with this name.
        if name in self._outputs:
            _warnings.warn("Duplicate input requirement. Overwriting existing value!")

        # Add the output to the dictionary.
        self._outputs[name] = output

        # Update the parser description.
        if not self._is_notebook:
            self._parser.description = self._generate_description()

    def setOutput(self, name, value):
        """
        Set the value of an output.

        Parameters
        ----------

        name : str
            The name of the output.

        value :
            The value of the output.
        """
        try:
            # Enforce strict naming for all file-based outputs. This ensures
            # that the prefix used matches the requirement name.
            if self._strict_file_naming:
                if isinstance(self._outputs[name], (_File, _FileSet)):
                    is_file = False
                    new_value = []
                    # For convenience, convert single file names to a list with
                    # one entry.
                    if isinstance(value, str):
                        value = [value]
                        is_file = True
                    # Loop over each file.
                    for file in value:
                        # Get the directory name, file prefix, and file extension.
                        basename = _os.path.basename(file)
                        dirname = _os.path.dirname(file)
                        fileprefix = basename.split(".")[0]
                        extension = basename.split(".")[1]

                        # If the file prefix doesn't match the requirement name, then
                        # rename, i.e. move, the file.
                        if fileprefix != name:
                            _warnings.warn(
                                f"Output file prefix '{fileprefix}' "
                                f"doesn't match requirement name '{name}'. "
                                "Renaming file!"
                            )
                            new_name = dirname + f"/{name}.{extension}"
                            _shutil.move(file, new_name)
                            file = new_name

                        # Store the new value of the file name.
                        new_value.append(file)

                    # Convert back into a single entry if this was a File requirement.
                    if is_file:
                        value = new_value[0]
                    else:
                        value = new_value

            self._outputs[name].setValue(value, name=name)
        except KeyError:
            raise

    def getInput(self, name):
        """
        Get the value of the named input.

        Parameters
        ----------

        name : str
            The name of the input requirement.

        Returns
        -------

        input :
            The value of the named input requirement.
        """

        if not isinstance(name, str):
            raise TypeError("The name must be of type 'str'")

        # Validate the inputs.
        self._is_valid_input = self._validateInput()

        try:
            value = self._inputs[name].getValue()
            if isinstance(value, list):
                return value.copy()
            else:
                return value
        except KeyError:
            raise

    def getInputs(self):
        """
        Get all of the input requirements.

        Returns
        -------

        inputs : { str : :class:`Requirement <BioSimSpace.Gateway._requirement.Requirement>` }
            The dictionary of input requirements.
        """

        # Validate the inputs.
        self._is_valid_input = self._validateInput()

        return self._inputs.copy()

    def addError(self, error):
        """
        Add an error message.

        Parameters
        ----------

        error : str
            The error message.
        """

        if not isinstance(error, str):
            raise TypeError("The error message must be of type 'str'")
        else:
            self._errors.append(error)

    def addAuthor(self, name=None, email=None, affiliation=None):
        """
        Add an author for the node.

        Parameters
        ----------

        name : str
            The author's name.

        email : str
            The author's email address.

        affiliation : str
            The author's affiliation.
        """

        if name is None:
            raise ValueError("Missing required 'name' keyword argument.")

        if not isinstance(name, str):
            raise TypeError("'name' must be of type 'str'")

        if email is not None and not isinstance(email, str):
            raise TypeError("'email' must be of type 'str'")

        if affiliation is not None and not isinstance(affiliation, str):
            raise TypeError("'affiliation' must be of type 'str'")

        if self._authors is None:
            self._authors = [{"name": name, "email": email, "affiliation": affiliation}]
        else:
            author = {"name": name, "email": email, "affiliation": affiliation}
            if not author in self._authors:
                self._authors.append(author)

    def getAuthors(self):
        """
        Return the list of authors.

        Returns
        -------

        authors : [dict]
           A list of author dictionaries.
        """
        return self._authors.copy()

    def setLicense(self, license):
        """
        Set the license for the node.

        Parameters
        ----------

        license : str
            The license type.
        """
        if not isinstance(license, str):
            raise TypeError("The license must be of type 'str'")
        else:
            self._license = license

    def getLicense(self):
        """
        Return the license.

        Returns
        -------

        license : str
            The license of the node.
        """
        return self._license

    def showControls(self):
        """
        Show the Jupyter widget GUI to allow the user to enter input.

        Returns
        -------

        controls : ipywidgets.form
           A gui control panel for setting input requirements.
        """

        if not self._is_notebook:
            return

        # Create the layout object.
        layout = _widgets.Layout(
            display="flex", flex_flow="row", justify_content="space-between"
        )

        # Initialise the list of form items.
        requirements = []

        # Initialise the list of input requirement 'set' indicators.
        indicators = []

        # Loop over all of the widgets.
        for name, widget in self._widgets.items():
            # Create the label string.
            string = "%s: %s" % (name, self._inputs[name].getHelp())

            # Add the unit information.
            unit = self._inputs[name].getUnit()
            if unit is not None:
                string += " (%s)" % self._inputs[name]._print_unit

            # Create the widget label.
            label = _widgets.Label(value=string)

            # This is a FileSet requirement with multiple widgets.
            if isinstance(widget, list):
                items = [label] + widget
                indicator = widget[0]._button
            else:
                items = [label, widget]
                indicator = widget._button

            # Create a box for the widget.
            box = _widgets.Box(items, layout=layout)

            # Add the widget box to the list.
            requirements.append(box)

            # Add the indicator to the list.
            indicators.append(indicator)

        # Create the widget form.
        form1 = _widgets.Box(
            requirements,
            layout=_widgets.Layout(
                display="flex",
                flex_flow="column",
                border="solid 2px",
                align_items="stretch",
                width="100%",
            ),
        )

        # Create the indicator form.
        form2 = _widgets.VBox(
            indicators,
            layout=_widgets.Layout(
                display="flex", flex_flow="column", align_items="stretch", width="5%"
            ),
        )

        # Combine the two forms.
        form = _widgets.Box([form1, form2], layout=layout)

        # Store the form.
        self._control_panel = form

        return form

    def _validateInput(self):
        """Validate the parsed inputs."""

        # Knime.
        if self._is_knime:
            pass

        # Jupyter.
        elif self._is_notebook:
            # Loop over the widgets and set the input values.
            for key, widget in self._widgets.items():
                # Use the widget value if it has been set, otherwise, set the value to None.
                # This ensures that the user actually sets a value.

                # File based widget.
                if hasattr(widget, "_is_file"):
                    if widget._is_set:
                        if len(widget._files) == 1:
                            # Single file.
                            value = widget._files[0]
                        else:
                            # File set.
                            value = widget._files
                    else:
                        value = None

                # Non file widget.
                else:
                    if widget._is_set:
                        value = widget.value
                    else:
                        value = None

                self._inputs[key].setValue(value, name=key)

        # Command-line.
        else:
            # Parse the arguments into a dictionary.
            args = vars(
                self._parser.parse_known_args(
                    args=None if _sys.argv[1:] else ["--help"]
                )[0]
            )

            # Now loop over the arguments and set the input values.
            for key, value in args.items():
                if key == "verbose":
                    setVerbose(value)
                elif key == "strict_file_naming":
                    if value is True:
                        self._strict_file_naming = True
                else:
                    if not key in ["config", "export_cwl"]:
                        self._inputs[key].setValue(value, name=key)

    def validate(self, file_prefix="output"):
        """
        Whether the output requirements are satisfied.

        Parameters
        ----------

        file_prefix : str
            The prefix of the output file name.

        Returns
        -------

        output : IPython.lib.display.FileLink, str
            If running interatvely: A link to a zipfile containing the
            validated output, else the name of a YAML file containing
            the node output.
        """

        if not isinstance(file_prefix, str):
            raise TypeError("The 'file_prefix' keyword must be of type 'str'.")

        # Flag that we have validated output.
        self._is_output_validated = True

        # A list of File and FileSet outputs.
        file_outputs = []

        # Check no outputs are None.
        for name, output in self._outputs.items():
            if output.getValue() is None:
                self._errors.append("Missing output for requirement '%s'" % name)
            else:
                if isinstance(output, (_File, _FileSet)):
                    file_outputs.append(output)

        # Node failed.
        if len(self._errors) > 0:
            for error in self._errors:
                print("%s" % error, file=_sys.stderr)

            if self._name is not None:
                raise SystemExit("Node '%s' failed!" % self._name)
            else:
                raise SystemExit("Node failed!")

        # Create a compressed archive containing all file output for the node.
        if self._is_notebook:
            # There are files.
            if len(file_outputs) > 0:
                # Create the archive name.
                zipname = "%s.zip" % file_prefix

                # Append the files to the archive.
                with _zipfile.ZipFile(zipname, "w") as zip:
                    # Loop over all of the file outputs.
                    for output in file_outputs:
                        if isinstance(output, _File):
                            file = output.getValue()
                            zip.write(file, arcname=_os.path.basename(file))
                        else:
                            for file in output.getValue():
                                zip.write(file, arcname=_os.path.basename(file))

                # Create a FileLink to the archive.
                file_link = _FileLink(zipname)

                # Set the download attribute so that JupyterLab doesn't try to open the file.
                file_link.html_link_str = (
                    f"<a href='%s' target='_blank' download='{zipname}'>%s</a>"
                )

                # Return a link to the archive.
                return file_link
        else:
            # Initialise an empty dictionary to store the output data.
            data = {}

            # Create the YAML file name.
            yamlname = "%s.yaml" % file_prefix

            # Populate the dictionary.
            for name, output in self._outputs.items():
                data[name] = output.getValue()

            # Write the outputs to a YAML file.
            with open(yamlname, "w") as file:
                _yaml.dump(data, file, default_flow_style=False)

            return yamlname

    def _create_help_string(self, input):
        """
        Create a nicely formatted argparse help string.

        Parameters
        ----------

        input : :class:`Requirement <BioSimSpace.Gateway.Requirement>`
            The input requirement.

        Returns
        -------

        help : str
            The formatted help string.
        """

        # Initialise the help string.
        help = "\n".join(_textwrap.wrap(input.getHelp(), 56))

        # Get the unit.
        units = input.getUnit()

        # Add the units to the help string.
        if units is not None:
            help += "\n  units=%s" % (units[0] + units[1:].lower())

        # Get the default value.
        default = input.getDefault()

        # Add default value to help string.
        if default is not None:
            help += "\n  default=%s" % default

        # Get the min/max allowed values.
        minimum = input.getMin()
        maximum = input.getMax()

        # Add min/max values to help string.
        if minimum is not None:
            if maximum is not None:
                help += "\n  min=%s, max=%s" % (minimum, maximum)
            else:
                help += "\n  min=%s" % minimum
        else:
            if maximum is not None:
                help += "\n  max=%s" % maximum

        return help

    def _generate_description(self):
        """
        Generate a formatted output section for the argparse help description.

        Returns
        -------

        output : str
            A string listing the output requirements.
        """

        string = "\n".join(_textwrap.wrap(self._description, 80))

        yaml_help = (
            "Args that start with '--' (e.g. --arg) can also be set "
            "in a config file (specified via -c). The config file uses "
            "YAML syntax and must represent a YAML 'mapping' "
            "(for details, see http://learn.getgrav.org/advanced/yaml). "
            "If an arg is specified in more than one place, then "
            "commandline values override config file values which "
            "override defaults."
        )

        string += "\n\n" + "\n".join(_textwrap.wrap(yaml_help, 80))

        # Initialise the output string.
        string += "\n\nOutput:\n"

        # Add documentation for each output.
        for name, output in self._outputs.items():
            num_whitespace = 19 - len(output.__class__.__name__) - len(name)
            s = "  %s: %s " % (name, output.__class__.__name__)
            for i, line in enumerate(_textwrap.wrap(output.getHelp(), 56)):
                if i == 0:
                    for _ in range(0, num_whitespace):
                        s += " "
                else:
                    for _ in range(0, 24):
                        s += " "
                s += line + "\n"
            string += s

        return string


def _on_value_change(change):
    """Helper function to flag that a widget value has been set."""
    change["owner"]._is_set = True
    change["owner"]._button.tooltip = "The input requirement is set."
    change["owner"]._button.button_style = "success"
    change["owner"]._button.icon = "fa-check"


def _on_file_upload(change):
    """Helper function to handle file uploads."""

    # Initialise the widget label.
    label = ""

    # Initialise file counter.
    num_files = 0

    # Clear the list of files.
    change["owner"]._files = []

    # Handle file upload widget API changes.
    if isinstance(change["owner"].value, tuple):
        is_tuple = True
    else:
        is_tuple = False

    # Loop over all uploaded files.
    for filename in change["owner"].value:
        if is_tuple:
            name = filename["name"]
            num_bytes = len(filename["content"])
        else:
            name = filename
            num_bytes = len(change["owner"].value[filename]["content"])

        # Return if there is no data.
        if num_bytes == 0:
            return

        # Separate label with commas.
        if num_files > 0:
            label += ", "

        # Extract the file content.
        if is_tuple:
            content = filename["content"]
        else:
            content = change["owner"].value[filename]["content"]

        # Create the uploads directory if it doesn't already exist.
        if not _os.path.isdir("uploads"):
            _os.makedirs("uploads")

        # Append the upload directory to the file name.
        new_filename = "uploads/%s" % name

        # Has this file already been uploaded?
        if _os.path.isfile(new_filename):
            # We'll append a number to the file name.
            index = 1
            new_filename_append = new_filename + ".%d" % index

            # Keep trying until a unique name is found.
            while _os.path.isfile(new_filename_append):
                index += 1
                new_filename_append = new_filename + ".%d" % index

            # Copy back into the new_filename variable.
            new_filename = new_filename_append

        # Write the file to disk.
        with open(new_filename, "wb") as file:
            file.write(content)

        # Report that the file was uploaded.
        print("Uploaded '{}' ({:.2f} kB)".format(name, num_bytes / 2**10))

        # Truncate the filename string if it is more than 15 characters.
        label += (name[:15] + "...") if len(name) > 15 else name

        # Increment the number of files.
        num_files += 1

        # Store the location of the uploaded file on disk.
        change["owner"]._files.append(new_filename)

    # Flag that the widget value has been set.
    change["owner"]._is_set = True
    change["owner"]._button.tooltip = "The input requirement is set."
    change["owner"]._button.button_style = "success"
    change["owner"]._button.icon = "fa-check"

    # Update the widget description with the name of the uploaded file/files.
    change["owner"].description = label

    # Update the widget counter. For some reason the default shows the total
    # number of files uploaded, rather than the current number of files.
    # This means that the number is incorrect if the user changes the files
    # that are uploaded, e.g. fixing an error, or re-running the same node
    # with different input.
    change["owner"]._counter = len(change["owner"].value)


def _check_value(action, value):
    """Helper function to overload argparse's choice checker."""
    if action.choices is not None and value not in action.choices:
        args = {"value": value, "choices": ", ".join(map(repr, action.choices))}
        msg = _argparse.argparse._(
            "invalid choice: %(value)r (choose from %(choices)s)"
        )

        # If the value is a string, then strip whitespace and try a case insensitive search.
        if isinstance(value, str):
            new_value = value.replace(" ", "").upper()
            choices = [x.replace(" ", "").upper() for x in action.choices]

            # Check whether we now have a match.
            if new_value not in choices:
                raise _argparse.ArgumentError(action, msg % args)
        else:
            raise _argparse.ArgumentError(action, msg % args)


def _str2bool(v):
    """Convert an argument string to a boolean value."""
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise _argparse.ArgumentTypeError("Boolean value expected.")
