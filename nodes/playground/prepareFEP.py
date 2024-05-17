#!/usr/bin/env python
# coding: utf-8

# Author: Julien Michel
#
# email: julien.michel@ed.ac.uk

# # PrepareFEP
# Loads a pair of input files, perform mapping between the first molecule of each input. Write down input files for a SOMD FEP calculation.

# In[ ]:


import os
import zipfile
import BioSimSpace as BSS
from sire.legacy.Mol import AtomIdx


# In[ ]:


def writeLog(ligA, ligB, mapping):
    """Human readable report on atoms used for the mapping."""
    atoms_in_A = list(mapping.keys())
    stream = open("somd.mapping", "w")
    atAdone = []
    atBdone = []
    for atAidx in atoms_in_A:
        atA = ligA._sire_object.select(AtomIdx(atAidx))
        atB = ligB._sire_object.select(AtomIdx(mapping[atAidx]))
        stream.write(
            "%s %s --> %s %s\n" % (atA.index(), atA.name(), atB.index(), atB.name())
        )
        atAdone.append(atA)
        atBdone.append(atB)
    for atom in ligA._sire_object.atoms():
        if atom in atAdone:
            continue
        stream.write("%s %s --> dummy\n" % (atom.index(), atom.name()))
    for atom in ligB._sire_object.atoms():
        if atom in atBdone:
            continue
        stream.write("dummy --> %s %s\n" % (atom.index(), atom.name()))
    stream.close()


# In[ ]:


def loadMapping(mapping_file):
    """Parse a text file that specifies mappings between atomic indices in input1 --> atoms in input2"""
    stream = open(mapping_file, "r")
    buffer = stream.readlines()
    stream.close()
    mapping = {}
    for line in buffer:
        if line.startswith("#"):
            continue
        elems = line.split(",")
        idx1 = int(elems[0])
        idx2 = int(elems[1])
        mapping[idx1] = idx2

    return mapping


# In[ ]:


node = BSS.Gateway.Node(
    "A node to generate input files for a SOMD relative free energy calculation."
)


# In[ ]:


node.addAuthor(
    name="Julien Michel",
    email="julien.michel@ed.ac.uk",
    affiliation="University of Edinburgh",
)
node.addAuthor(
    name="Lester Hedges",
    email="lester.hedges@bristol.ac.uk",
    affiliation="University of Bristol",
)
node.setLicense("GPLv3")


# In[ ]:


node.addInput("input1", BSS.Gateway.FileSet(help="A topology and coordinates file"))
node.addInput("input2", BSS.Gateway.FileSet(help="A topology and coordinates file"))
node.addInput(
    "prematch",
    BSS.Gateway.String(
        help="list of atom indices that are matched between input2 and input1. Syntax is of the format 1-3,4-8,9-11... Ignored if a mapping is provided",
        default="",
    ),
)
node.addInput(
    "mapping",
    BSS.Gateway.File(
        help="csv file that contains atom indices in input1 mapped ot atom indices in input2",
        optional=True,
    ),
)
node.addInput(
    "timeout",
    BSS.Gateway.Time(
        help="The timeout for the maximum common substructure search",
        default=10 * BSS.Units.Time.second,
    ),
)
node.addInput(
    "allow_ring_breaking",
    BSS.Gateway.Boolean(
        help="Whether to allow opening/closing of rings during merge", default=False
    ),
)
node.addInput(
    "allow_ring_size_change",
    BSS.Gateway.Boolean(
        help="Whether to allow ring size changes during merge", default=False
    ),
)
node.addInput(
    "output",
    BSS.Gateway.String(
        help="The root name for the files describing the perturbation input1->input2."
    ),
)
node.addInput(
    "somd2",
    BSS.Gateway.Boolean(
        help="Whether to generate input for SOMD2.",
        default=False,
    ),
)


# In[ ]:


node.addOutput(
    "nodeoutput",
    BSS.Gateway.FileSet(help="SOMD input files for a perturbation of input1->input2."),
)


# In[ ]:


node.showControls()


# In[ ]:


do_mapping = True
custom_mapping = node.getInput("mapping")
# print (custom_mapping)
if custom_mapping is not None:
    do_mapping = False
    mapping = loadMapping(custom_mapping)
    # print (mapping)


# In[ ]:


# Optional input, dictionary of Atom indices that should be matched in the search.
prematch = {}
prematchstring = node.getInput("prematch")
if len(prematchstring) > 0:
    entries = prematchstring.split(",")
    for entry in entries:
        idxA, idxB = entry.split("-")
        prematch[int(idxA)] = int(idxB)
# print (prematch)


# In[ ]:


# Load system 1
system1 = BSS.IO.readMolecules(node.getInput("input1"))


# In[ ]:


# Load system 2
system2 = BSS.IO.readMolecules(node.getInput("input2"))


# In[ ]:


# We assume the molecules to perturb are the first molecules in each system
lig1 = system1[0]
lig2 = system2[0]


# In[ ]:


if do_mapping:
    # Return a maximum of 10 matches, scored by RMSD and sorted from best to worst.
    mappings, scores = BSS.Align.matchAtoms(
        lig1,
        lig2,
        matches=10,
        prematch=prematch,
        return_scores=True,
        scoring_function="RMSDalign",
        timeout=node.getInput("timeout"),
    )
    # We retain the top mapping
    mapping = mappings[0]
    # print (len(mappings))
    # print (mappings)


# In[ ]:


# print (mapping)
# for x in range(0,len(mappings)):
#    print (mappings[x], scores[x])


# In[ ]:


inverted_mapping = dict([[v, k] for k, v in mapping.items()])
# print (inverted_mapping)


# In[ ]:


# Align lig2 to lig1 based on the best mapping (inverted). The molecule is aligned based
# on a root mean squared displacement fit to find the optimal translation vector
# (as opposed to merely taking the difference of centroids).
lig2 = BSS.Align.rmsdAlign(lig2, lig1, inverted_mapping)
# Merge the two ligands based on the mapping.
merged = BSS.Align.merge(
    lig1,
    lig2,
    mapping,
    allow_ring_breaking=node.getInput("allow_ring_breaking"),
    allow_ring_size_change=node.getInput("allow_ring_size_change"),
)
# Create a composite system
system1.removeMolecules(lig1)
system1.addMolecules(merged)
# Make sure the box vectors are in reduced form.
system1.reduceBoxVectors()
system1.rotateBoxVectors()


# In[ ]:


# Log the mapping used
writeLog(lig1, lig2, mapping)
# Are we saving output for SOMD2?
is_somd2 = node.getInput("somd2")
# File root for all output.
root = node.getInput("output")
BSS.IO.saveMolecules(
    "merged_at_lam0.pdb",
    merged,
    "PDB",
    property_map={
        "coordinates": "coordinates0",
        "bond": "bond0",
        "element": "element0",
    },
)
if is_somd2:
    BSS.Stream.save(system1, root)
    stream_file = "%s.bss" % root
else:
    # Generate package specific input
    protocol = BSS.Protocol.FreeEnergy(
        runtime=2 * BSS.Units.Time.femtosecond, num_lam=3
    )
    process = BSS.Process.Somd(system1, protocol)
    process.getOutput()
    with zipfile.ZipFile("somd_output.zip", "r") as zip_hnd:
        zip_hnd.extractall(".")


# In[ ]:


if not is_somd2:
    mergedpdb = "%s.mergeat0.pdb" % root
    pert = "%s.pert" % root
    prm7 = "%s.prm7" % root
    rst7 = "%s.rst7" % root
    mapping_str = "%s.mapping" % root


# In[ ]:


if not is_somd2:
    os.replace("merged_at_lam0.pdb", mergedpdb)
    os.replace("somd.pert", pert)
    os.replace("somd.prm7", prm7)
    os.replace("somd.rst7", rst7)
    os.replace("somd.mapping", mapping_str)
    try:
        os.remove("somd_output.zip")
        os.remove("somd.cfg")
        os.remove("somd.err")
        os.remove("somd.out")
    except Exception:
        pass


# In[ ]:


if is_somd2:
    output = [stream_file]
else:
    output = [mergedpdb, pert, prm7, rst7, mapping_str]
node.setOutput("nodeoutput", output)


# In[ ]:


node.validate()
