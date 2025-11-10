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

"""
A thin wrapper around Sire.Mol.Residue. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Residue"]


from ._sire_wrapper import SireWrapper as _SireWrapper


class Residue(_SireWrapper):
    """A class for storing a residue."""

    def __init__(self, residue):
        """
        Constructor.

        Parameters
        ----------

        residue : Sire.Mol.Residue, :class:`Residue <BioSimSpace._SireWrappers.Residue>`
            A Sire or BioSimSpace Residue object.
        """
        from sire.legacy import Mol as _SireMol

        # Check that the residue is valid.

        # A Sire Residue object.
        if isinstance(residue, _SireMol._Mol.Residue):
            sire_object = residue

        # Another BioSimSpace Residue object.
        elif isinstance(residue, Residue):
            sire_object = residue._sire_object

        # Invalid type.
        else:
            raise TypeError(
                "'residue' must be of type 'Sire.Mol.Residue' "
                "or 'BioSimSpace._SireWrappers.Residue'."
            )

        # Call the base class constructor.
        super().__init__(sire_object)

        # Flag that this object holds multiple atoms.
        self._is_multi_atom = True

        # Store the number of atoms in the residue.
        self._num_atoms = self._sire_object.num_atoms()

        # Store the atom indices in the residue.
        self._atom_idxs = self._sire_object.atom_idxs()

        # Initialise the iterator count.
        self._iter_count = 0

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Residue: name=%r, molecule=%d, index=%d, nAtoms=%d>" % (
            self.name(),
            self.moleculeNumber(),
            self.index(),
            self.nAtoms(),
        )

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Residue: name=%r, molecule=%d, index=%d, nAtoms=%d>" % (
            self.name(),
            self.moleculeNumber(),
            self.index(),
            self.nAtoms(),
        )

    def __contains__(self, other):
        """Return whether other is in self."""
        from ._atom import Atom as _Atom

        if not isinstance(other, _Atom):
            raise TypeError("'other' must be of type 'BioSimSpace._SireWrappers.Atom'.")

        # Return whether the residue contains the atom.
        return self._sire_object.contains(other._sire_object.atom().number())

    def __getitem__(self, key):
        """Get an atom from the residue."""
        from ._atom import Atom as _Atom

        # Slice.
        if isinstance(key, slice):
            # Create a list to hold the atoms.
            atoms = []

            # Iterate over the slice.
            for x in range(*key.indices(self._num_atoms)):
                atoms.append(_Atom(self[x]))

            # Return the list of atoms.
            return atoms

        # Index.
        else:
            try:
                key = int(key)
            except:
                raise TypeError("'key' must be of type 'int'")

            if key < -self._num_atoms or key > self._num_atoms - 1:
                raise IndexError("Residue index is out of range.")

            if key < 0:
                key = key + self._num_atoms

            # Extract and return the corresponding atom.
            return _Atom(self._sire_object.atom(self._atom_idxs[key]))

    def __setitem__(self, key, value):
        """Set an atom in the residue."""
        raise TypeError("'Residue' object does not support assignment.")

    def __iter__(self):
        """An iterator for the object."""
        # Reset the iterator counter and return the object.
        self._iter_count = 0
        return self

    def __next__(self):
        """An iterator for the object."""

        # Stop if we've reached the end of the residue.
        if self._iter_count == self._num_atoms:
            raise StopIteration

        # Extract the next atom in the residue.
        atom = self[self._iter_count]

        # Update the iterator counter.
        self._iter_count += 1

        # Return the atom.
        return atom

    def __len__(self):
        """Return the number of atoms in the residue."""
        return self._num_atoms

    def name(self):
        """
        Return the name of the residue.

        Returns
        -------

        name : str
            The name of the residue.
        """
        return self._sire_object.name().value()

    def index(self):
        """
        Return the index of the residue.

        Returns
        -------

        index : int
            The index of the residue.
        """
        return self._sire_object.index().value()

    def moleculeNumber(self):
        """
        Return the number of the molecule to which this residue belongs.

        Returns
        -------

        number : int
            The number of the molecule to which the residue belongs.
        """
        return self._sire_object.molecule().number().value()

    def coordinates(self, property_map={}):
        """
        Return the coordinates of the atoms in the residue.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        [coordinates] : [class:`Coordinate <BioSimSpace.Types.Coordinate>`]
            The coordinates of the atoms in the residue.
        """
        # Get the "coordinates" property for each atom in the residue.
        try:
            coordinates = []
            for atom in self.getAtoms():
                coordinates.append(atom.coordinates(property_map))
        except:
            return None

        # Return the coordinates.
        return coordinates

    def nAtoms(self):
        """
        Return the number of atoms in the residue.

        Returns
        -------

        num_atoms : int
            The number of atoms in the system.
        """
        return self._num_atoms

    def getAtoms(self):
        """
        Return a list containing all of the atoms in the residue.

        Parameters
        ----------

        Returns
        -------

        atoms : [:class:`Atoms <BioSimSpace._SireWrappers.Atom>`]
            The list of atoms in the residue.
        """
        from ._atom import Atom as _Atom

        atoms = []
        for atom in self._sire_object.atoms():
            atoms.append(_Atom(atom))
        return atoms

    def toMolecule(self):
        """
        Convert a single Residue to a Molecule.

        Returns
        -------

        system : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
        """
        from ._molecule import Molecule as _Molecule
        from sire.legacy import Mol as _SireMol

        return _Molecule(
            _SireMol.PartialMolecule(self._sire_object).extract().molecule()
        )

    def search(self, query, property_map={}):
        """
        Search the residue for atoms and residues.

        Parameters
        ----------

        query : str
            The search query.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        results : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`]
            A list of objects matching the search query.

        Examples
        --------

        Search for all oxygen or hydrogen atoms.

        >>> result = residue.search("element oxygen or element hydrogen")

        Search for atom index 23.

        >>> result = residue.search("atomidx 23")
        """
        from .. import _isVerbose
        from ._search_result import SearchResult as _SearchResult
        from sire.legacy import Mol as _SireMol

        if not isinstance(query, str):
            raise TypeError("'query' must be of type 'str'")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Initialise a list to hold the search results.
        results = []

        try:
            # Query the Sire residue.
            search_result = _SireMol.Select(query)(self._sire_object, property_map)

        except Exception as e:
            msg = "'Invalid search query: %r" % query
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        return _SearchResult(search_result)
