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


from molmod.io.common import slice_match

import numpy


__all__ = ["Error", "DumpReader"]


class Error(Exception):
    pass


class DumpReader(object):
    def __init__(self, filename, units, sub=slice(None)):
        # first read the number of atoms
        self._f = file(filename)
        try:
            while True:
                line = self._f.next()
                if line == "ITEM: NUMBER OF ATOMS\n":
                    break
            try:
                line = self._f.next()
                self.num_atoms = int(line)
            except ValueError:
                raise Error("Could not read the number of atoms. Expected an integer. Got '%s'" % line)
        except StopIteration:
            raise Error("Could not find line 'ITEM: NUMBER OF ATOMS'.")
        self._f.seek(0) # go back to the beginning of the file
        self._f = file(filename)
        self.units = units
        self.sub = sub
        self.counter = 0

    def __del__(self):
        self._f.close()

    def __iter__(self):
        return self

    def next(self):
        # Read one frame, we assume that the current file position is at the
        # line 'ITEM: TIMESTEP' and that this line marks the beginning of a
        # time frame.
        line = self._f.next()
        if line != 'ITEM: TIMESTEP\n':
            raise Error("Expecting line 'ITEM: TIMESTEP' at the beginning of a time frame.")
        try:
            line = self._f.next()
            step = int(line)
        except ValueError:
            raise Error("Could not read the step number. Expected an integer. Got '%s'" % line[:-1])

        # Now we assume that the next section contains (again) the number of
        # atoms.
        line = self._f.next()
        if line != 'ITEM: NUMBER OF ATOMS\n':
            raise Error("Expecting line 'ITEM: NUMBER OF ATOMS'.")
        try:
            line = self._f.next()
            num_atoms = int(line)
        except ValueError:
            raise Error("Could not read the number of atoms. Expected an integer. Got '%s'" % line[:-1])
        if num_atoms != self.num_atoms:
            raise Error("A variable number of atoms is not supported.")

        # The next section contains the box boundaries. We will skip it
        for i in xrange(4):
            self._f.next()

        # The next and last section contains the atom related properties
        line = self._f.next()
        if line != 'ITEM: ATOMS\n':
            raise Error("Expecting line 'ITEM: ATOMS'.")
        fields = [list() for i in xrange(len(self.units))]
        for i in xrange(self.num_atoms):
            line = self._f.next()
            words = line.split()[1:]
            for j in xrange(len(fields)):
                fields[j].append(float(words[j]))
        fields = [step] + [numpy.array(field)*unit for field, unit in zip(fields, self.units)]

        return fields





