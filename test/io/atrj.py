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

from molmod.io.atrj import ATRJReader
from molmod.units import ps, kcalmol, angstrom

import numpy, unittest


__all__ = ["ATRJTestCase"]


class ATRJTestCase(BaseTestCase):
    def test_load(self):
        # A) normal
        atrj_reader = ATRJReader("input/bartek.atrj")
        self.assertEqual(atrj_reader.num_atoms, 1293)
        frames = list(atrj_reader)
        self.assertEqual(len(frames), 3)
        # check time
        self.assertAlmostEqual(frames[0].time/ps, 1.0)
        self.assertAlmostEqual(frames[1].time/ps, 2.0)
        self.assertAlmostEqual(frames[2].time/ps, 3.0)
        # check step
        self.assertAlmostEqual(frames[0].step, 1000)
        self.assertAlmostEqual(frames[1].step, 2000)
        self.assertAlmostEqual(frames[2].step, 3000)
        # check total energy
        self.assertAlmostEqual(frames[0].total_energy/kcalmol, 3.4186035768405162e2)
        self.assertAlmostEqual(frames[1].total_energy/kcalmol, 3.3443356630787252e2)
        self.assertAlmostEqual(frames[2].total_energy/kcalmol, 3.3613629561285467e2)
        # check (some of) the coordinates
        self.assertAlmostEqual(frames[0].coordinates[0,0]/angstrom, 1.1953453341349823e1)
        self.assertAlmostEqual(frames[1].coordinates[5,1]/angstrom, 1.2284911431298470e1)
        self.assertAlmostEqual(frames[-1].coordinates[-5,-1]/angstrom, 2.1392983758428979e1)

        # B) sliced
        atrj_reader = ATRJReader("input/bartek.atrj", slice(None, None, 2))
        self.assertEqual(atrj_reader.num_atoms, 1293)
        frames = list(atrj_reader)
        self.assertEqual(len(frames), 2)
        # check time
        self.assertAlmostEqual(frames[0].time/ps, 1.0)
        self.assertAlmostEqual(frames[1].time/ps, 3.0)
