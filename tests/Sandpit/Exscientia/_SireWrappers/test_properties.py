import BioSimSpace.Sandpit.Exscientia as BSS
from tests.conftest import root_fp


def test_sire_properties():
    s = BSS.IO.readMolecules([f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"])

    m = s[0]._sire_object

    charges = m.property("charge")

    c = charges.to_vector()

    assert len(c) == m.num_atoms()

    masses = m.property("mass")

    v = masses.to_vector()

    assert len(v) == m.num_atoms()

    a = m.property("ambertype")

    a = a.to_vector()

    assert len(a) == m.num_atoms()

    g = m.property("gb_radii")

    g = g.to_vector()

    assert len(g) == m.num_atoms()

    s = m.property("gb_screening")

    s = s.to_vector()

    assert len(s) == m.num_atoms()

    t = m.property("treechain")

    t = t.to_vector()

    assert len(t) == m.num_atoms()


if __name__ == "__main__":
    test_sire_properties()
