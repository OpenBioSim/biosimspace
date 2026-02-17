import pytest

import BioSimSpace as BSS


def test_bonds():
    system = BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])

    bonds = system.search("bonds")

    assert len(bonds) == 1911

    for bond in bonds:
        atom0 = bond.atom0()
        atom1 = bond.atom1()

        assert atom0 == bond[0]
        assert atom1 == bond[1]

        assert bond.moleculeNumber() == atom0.moleculeNumber()

        c0 = atom0.coordinates()
        c1 = atom1.coordinates()

        d2 = (c0.x() - c1.x()) ** 2 + (c0.y() - c1.y()) ** 2 + (c0.z() - c1.z()) ** 2

        assert d2.value() == pytest.approx(bond.length().value() ** 2)

        # Some of these bonds have zero energy.
        if bond._sire_object.energy().value() == 0:
            import sire

            amber_bond = sire.legacy.MM.AmberBond(
                bond._sire_object.potential(), sire.legacy.CAS.Symbol("r")
            )
            assert amber_bond.r0() == pytest.approx(bond.length().value())
