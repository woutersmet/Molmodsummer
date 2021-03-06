# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# MolMod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --

from common import BaseTestCase

from molmod.randomize import *
from molmod.molecular_graphs import MolecularGraph
from molmod.units import A
from molmod.transformations import Translation, Rotation

from molmod.io.xyz import XYZFile

import unittest, numpy, os

__all__ = ["RandomizeTestCase"]

# These threshold values are not good for serious applications!!!
# Just for testing the code, they are OK.
nonbond_thresholds = {
    frozenset([1,1]): 0.9*A,
    frozenset([1,6]): 1.4*A,
    frozenset([1,7]): 1.4*A,
    frozenset([1,8]): 1.4*A,
    frozenset([6,6]): 2.2*A,
    frozenset([6,7]): 2.2*A,
    frozenset([6,8]): 2.2*A,
    frozenset([7,7]): 2.2*A,
    frozenset([7,8]): 2.2*A,
    frozenset([8,8]): 2.2*A,
}


class RandomizeTestCase(BaseTestCase):
    def yield_test_molecules(self):
        for filename in ["tpa.xyz", "water.xyz", "thf_single.xyz"]:
            molecule = XYZFile(os.path.join("input", filename)).get_molecule()
            molecule.filename = filename
            graph = MolecularGraph.from_geometry(molecule)
            yield molecule, graph

    def test_randomize_count(self):
        all_counts = {
            "tpa.xyz": (12+4*7, 12, 13*6, 0),
            "water.xyz": (2,0,1,0), # Stretch, Torsion, Bend, DoubleStretch
            "thf_single.xyz": (8, 14, 5*4, 10),
        }
        for molecule, graph in self.yield_test_molecules():
            manipulations = generate_manipulations(graph, molecule)
            randomized_molecule = randomize_molecule(molecule, graph, manipulations, nonbond_thresholds)
            randomized_molecule.write_to_file(os.path.join("output", molecule.filename))
            counts = all_counts.get(molecule.filename)
            if counts is not None:
                for cls, count in zip([RandomStretch, RandomTorsion, RandomBend, RandomDoubleStretch], counts):
                    got = sum(isinstance(mpl, cls) for mpl in manipulations)
                    self.assertEqual(got, count, "%s count problem, should be %i, got %i" % (str(cls), count, got))

    def test_random_dimer(self):
        for molecule1, graph1 in self.yield_test_molecules():
            for molecule2, graph2 in self.yield_test_molecules():
                dimer = random_dimer(molecule1, molecule2, nonbond_thresholds, 0.5*A)
                self.assertEqual(dimer.coordinates.shape, (molecule1.coordinates.shape[0] + molecule2.coordinates.shape[0], 3))
                self.assertEqual(dimer.numbers.shape, (molecule1.numbers.shape[0] + molecule2.numbers.shape[0],))
                dimer.write_to_file(os.path.join("output", "%s_%s" % (molecule1.filename, molecule2.filename)))

    def test_single_manipulation(self):
        for molecule, graph in self.yield_test_molecules():
            manipulations = generate_manipulations(graph, molecule)
            for i in xrange(100):
                randomized_molecule, mol_transformation = single_random_manipulation(molecule, graph, manipulations, nonbond_thresholds)
                randomized_molecule.write_to_file(os.path.join("output", molecule.filename))
                mol_transformation.write_to_file(os.path.join("output", "tmp"))
                check_transformation = MolecularTransformation.read_from_file(os.path.join("output", "tmp"))
                self.assertEqual(mol_transformation.affected_atoms, check_transformation.affected_atoms)
                self.assertArraysAlmostEqual(mol_transformation.transformation.r, check_transformation.transformation.r, 1e-5, doabs=True)
                self.assertArraysAlmostEqual(mol_transformation.transformation.t, check_transformation.transformation.t, 1e-5, doabs=True)

