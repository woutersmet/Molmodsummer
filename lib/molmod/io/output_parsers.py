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


import os.path


class Error(Exception):
    pass


class OutputParser(object):
    def __init__(self, file_parsers=[]):
        self.clear()
        self.add_parsers(file_parsers)

    def clear(self):
        self.file_parser_groups = {}

    def add_parsers(self, file_parsers):
        for file_parser in file_parsers:
            tag = (file_parser.filename, file_parser.extension)
            file_parser_group = self.file_parser_groups.get(tag)
            if file_parser_group is None:
                file_parser_group = FileParserGroup(file_parser.filename, file_parser.extension)
                self.file_parser_groups[tag] = file_parser_group
            file_parser_group.add_parser(file_parser)

    def sort_groups(self):
        for file_parser_group in self.file_parser_groups.itervalues():
            file_parser_group.depends_on = []
            file_parser_group.works_for = []

        for file_parser_group in self.file_parser_groups.itervalues():
            for file_parser in file_parser_group.items:
                for dependency in file_parser.depends_on:
                    file_parser_group.depends_on.append(dependency.group)
                    #dependency.group.works_for.append(file_parser_group)

        result = self.file_parser_groups.values()
        result.sort(FileParserGroup.compare)
        return result

    def parse(self, prefix):
        sorted_groups = self.sort_groups()

        directory = os.path.dirname(prefix)
        if len(directory) > 0: directory += "/"
        basename = os.path.basename(prefix)

        result = {}
        for file_parser_group in sorted_groups:
            for file_parser in file_parser_group.items:
                file_parser.reset()
            if file_parser_group.extension:
                path = "%s%s%s" % (directory, basename, file_parser_group.filename)
            else:
                path = "%s%s" % (directory, file_parser_group.filename)
            if os.path.isfile(path):
                f = file(path, 'r')
                for line in f:
                    #print line[:-1]
                    for file_parser in file_parser_group.items:
                        file_parser.conditioned_parse(line)
                f.close()
                for file_parser in file_parser_group.items:
                    result[file_parser.label] = file_parser.result()
            else:
                raise Error("File %s not found" % path)

        return result


class FileParserGroup(object):
    def __init__(self, filename, extension):
        self.items = []
        self.depends_on = None
        #self.works_for = None
        self.filename = filename
        self.extension = extension

    def compare(self, other):
        if self == other:
            return 0
        if self in other.depends_on:
            return -1
        if other in self.depends_on:
            return 1
        return 0

    def add_parser(self, file_parser):
        self.items.append(file_parser)
        file_parser.group = self


class FileParser(object):
    extension = None

    def __init__(self, label, condition=None, depends_on=[]):
        self.label = label
        self.condition = condition
        self.depends_on = depends_on
        self.group = None

    def reset(self):
        raise NotImplementedError

    def conditioned_parse(self, line):
        if (self.condition is None) or self.condition():
            self.parse(line)

    def parse(self, line):
        raise NotImplementedError

    def result(self):
        raise NotImplementedError


class MultiLineParser(FileParser):
    def __init__(self, label, activator, deactivator, condition=None, depends_on=[]):
        FileParser.__init__(self, label, condition, depends_on)
        self.activator = activator
        self.deactivator = deactivator

    def reset(self):
        self.active = False

    def parse(self, line):
        if self.active:
            if self.deactivator is not None and self.deactivator.search(line) is not None:
                #print "Deactivated on line: %s" % line[:-1]
                self.active = False
                self.stop_collecting()
            else:
                if self.collect(line) is not None:
                    self.active = False
                    self.stop_collecting()
        elif self.activator is not None and self.activator.search(line) is not None:
            #print "Activated on line:   %s" % line[:-1]
            self.active = True
            self.start_collecting()
        elif self.activator is None and self.deactivator is None:
            self.collect(line)

    def start_collecting(self):
        raise NotImplementedError

    def collect(self, line):
        raise NotImplementedError

    def stop_collecting(self):
        raise NotImplementedError








