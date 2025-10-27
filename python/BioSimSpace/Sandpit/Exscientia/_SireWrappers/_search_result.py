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
A thin wrapper around Sire.Mol.SelectResult. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["SearchResult"]


class SearchResult:
    """A thin wrapper around Sire.Mol.SelectResult."""

    def __init__(self, select_result):
        """
        Constructor.

        Parameters
        ----------

        select_result : Output of a Sire search
            The Sire result object.
        """

        # Note that the type of a sire search will now be
        # whatever makes most sense for the result, e.g.
        # SelectorMol if this is a list of molecules,
        # SelectorBond if this is a list of bonds etc.

        # Store the Sire select result.
        self._sire_object = select_result.__deepcopy__()

        # Store the number of results.
        self._num_results = len(self._sire_object)

        # Initialise the iterator count.
        self._iter_count = 0

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.SearchResult: nResults=%d>" % self.nResults()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.SearchResult: nResults=%d>" % self.nResults()

    def __len__(self):
        """
        Return the number of results.

        Returns
        -------

        num_results : int
            The number of search results.
        """
        return self._num_results

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object == other._sire_object
        else:
            return False

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object != other._sire_object
        else:
            return False

    def __hash__(self):
        """Hash operator."""
        return hash(self._sire_object)

    def __getitem__(self, key):
        """Get a search result from the container."""
        from ._bond import Bond as _Bond
        from ._atom import Atom as _Atom
        import sire.legacy as _Sire
        from ._molecule import Molecule as _Molecule
        from ._residue import Residue as _Residue

        # Slice.
        if type(key) is slice:
            # Create a list to hold the results.
            results = []

            # Iterate over the slice.
            for x in range(*key.indices(self._num_results)):
                results.append(self[x])

            # Return the results.
            return results

        # Index.
        else:
            try:
                key = int(key)
            except:
                raise TypeError("'key' must be of type 'int'")

            if key < -self._num_results or key > self._num_results - 1:
                raise IndexError("SearchResult index is out of range.")

            if key < 0:
                key = key + self._num_results

            # Extract the result from the Sire object.
            result = self._sire_object[key]

            # Atom.
            if isinstance(result, _Sire.Mol._Mol.Atom):
                return _Atom(result)
            # Residue.
            elif isinstance(result, _Sire.Mol._Mol.Residue):
                return _Residue(result)
            # Molecule.
            elif isinstance(result, _Sire.Mol._Mol.Molecule):
                return _Molecule(result)
            # Bond
            elif isinstance(result, _Sire.MM._MM.Bond):
                return _Bond(result)
            # Unsupported.
            else:
                print(f"WARNING: Type {type(result)} is not supported.")
                return None

    def __setitem__(self, key, value):
        """Set a molecule in the container."""
        raise TypeError("'SearchResult' object does not support assignment.")

    def __iter__(self):
        """An iterator for the object."""
        # Reset the iterator counter and return the object.
        self._iter_count = 0
        return self

    def __next__(self):
        """An iterator for the object."""

        # Stop if we've reached the end of the container.
        if self._iter_count == self._num_results:
            raise StopIteration

        # Extract the next result in the container.
        result = self[self._iter_count]

        # Update the iterator counter.
        self._iter_count += 1

        # Return the result.
        return result

    def copy(self):
        """
        Create a copy of this object.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object.
        """
        return SearchResult(self)

    def atoms(self):
        """
        Return all of the atoms that contain the results of this
        search.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object viewed via its atoms.
        """
        return SearchResult(self._sire_object.atoms())

    def residues(self):
        """
        Return all of the residues that contain the results of this
        search.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object viewed via its residues.
        """
        return SearchResult(self._sire_object.residues())

    def chains(self):
        """
        Return all of the chains that contain the results of this
        search.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object viewed via its chains.
        """
        return SearchResult(self._sire_object.chains())

    def segments(self):
        """
        Return all of the segments that contain the results of this
        search.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object viewed via its segments.
        """
        return SearchResult(self._sire_object.segments())

    def molecules(self):
        """
        Return all of the molecules that contain the results of this
        search.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object viewed via its molecules.
        """
        try:
            return SearchResult(self._sire_object.molecules())
        except Exception:
            # this is a single molecule view
            from sire.mol import SelectorMol

            return SearchResult(SelectorMol(self._sire_object.molecule()))

    def bonds(self):
        """
        Return all of the bonds in this result.

        Returns
        -------

        search_result : :class:`SearchResult <BioSimSpace._SireWrappers.SearchResult>`
            A copy of the object viewed via its bonds.
        """
        try:
            from sire.mm import SelectorMBond

            return SearchResult(SelectorMBond(self._sire_object))
        except Exception:
            from sire.mm import SelectorBond

            return SearchResult(SelectorBond(self._sire_object))

    def nResults(self):
        """
        Return the number of results.

        Returns
        -------

        num_results : int
            The number of search results.
        """
        return self._num_results

    def getResult(self, index):
        """
        Return the result at the given index.

        Parameters
        ----------

        index : int
            The index of the result.

        Returns
        -------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The requested molecule.

        result : :class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                 :class:`Residue <BioSimSpace._SireWrappers.Residue>`, \
                 :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                 :class:`Bond <BioSimSpace._SireWrappers.Bond>`
        """
        return self[index]

    def _getSireObject(self):
        """
        Return the underlying Sire object.

        Returns
        -------

        object : Sire.Mol.SelectResult
        """
        return self._sire_object
