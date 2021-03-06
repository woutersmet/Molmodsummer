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


from molmod.units import ps, nm
from molmod.io.common import slice_match

import numpy


__all__ = ["GroReader"]


class Error(Exception):
    pass


class GroReader(object):
    def __init__(self, filename, sub=slice(None)):
        self._f = file(filename)
        self._sub = sub
        self.num_atoms = None
        (time, pos, vel, cell) = self.read_frame()
        self.num_atoms = len(pos)
        self._counter = 0
        self._f.seek(0)

    def get_line(self):
        line = self._f.readline()
        if len(line) == 0:
            raise StopIteration
        return line

    def read_frame(self):
        # Read the first line, ignore the title and try to get the time. The
        # time field is optional.
        line = self.get_line()
        pos = line.rfind("t=")
        if pos >= 0:
            time = float(line[pos+2:])*ps
        else:
            time = 0.0
        # Read the second line, the number of atoms must match with the first
        # frame.
        num_atoms = int(self.get_line())
        if self.num_atoms is not None and self.num_atoms != num_atoms:
            raise ValueError("The number of atoms must be the same over the entire file.")
        # Read the atom lines
        pos = numpy.zeros((num_atoms,3),numpy.float32)
        vel = numpy.zeros((num_atoms,3),numpy.float32)
        for i in xrange(num_atoms):
            words = self.get_line()[22:].split()
            pos[i,0] = float(words[0])
            pos[i,1] = float(words[1])
            pos[i,2] = float(words[2])
            vel[i,0] = float(words[3])
            vel[i,1] = float(words[4])
            vel[i,2] = float(words[5])
        pos *= nm
        vel *= nm/ps
        # Read the cell line
        cell = numpy.zeros((3,3), numpy.float32)
        words = self.get_line().split()
        if len(words) >= 3:
            cell[0,0] = float(words[0])
            cell[1,1] = float(words[1])
            cell[2,2] = float(words[2])
        if len(words) == 9:
            cell[1,0] = float(words[3])
            cell[2,0] = float(words[4])
            cell[0,1] = float(words[5])
            cell[2,1] = float(words[6])
            cell[0,2] = float(words[7])
            cell[1,2] = float(words[8])
        cell *= nm
        return time, pos, vel, cell

    def skip_frame(self):
        line = self.get_line()
        num_atoms = int(self.get_line())
        if self.num_atoms is not None and self.num_atoms != num_atoms:
            raise ValueError("The number of atoms must be the same over the entire file.")
        for i in xrange(num_atoms+1):
            self.get_line()

    def __del__(self):
        self._f.close()

    def __iter__(self):
        return self

    def next(self):
        # skip frames as requested
        while not slice_match(self._sub, self._counter):
            self._counter += 1
            self.skip_frame()

        self._counter += 1
        return self.read_frame()



