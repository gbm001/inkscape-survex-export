#!/usr/bin/env python3
"""
svx_output.py
Python script for exporting survex (.svx) file from Inkscape

Copyright (C) 2015 Patrick B Warren

Email: patrickbwarren@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<http://www.gnu.org/licenses/>.
"""

import os
import sys
import math
from time import strftime
from itertools import combinations
from inkex.paths import Path
from inkex.utils import AbortExtension
import inkex


# Define a (trivial) exception class to catch path errors
class PathError(AbortExtension):
    pass


def sprintd(b):
    "Takes a bearing and returns it as string in 000 format"
    while b < 0:
        b += 360
    b = int(b + 0.5)
    while b >= 360:
        b -= 360
    return f'{b:03}'


def distance(p1, p2):
    "Calculate the distance between two points"
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dl = math.sqrt(dx*dx + dy*dy)
    return dx, dy, dl


def measure(steps):
    "Measure the distance between first two steps"
    return distance(steps[0][1], steps[1][1])


class SurvexOutputExtension(inkex.extensions.OutputExtension):
    def add_arguments(self, pars):
        pars.add_argument('--tab', type=str, dest='tab', default='',
                          help='Dummy argument')

        pars.add_argument('--scale', type=float, dest='scale', default='100.0',
                          help='Length of scale bar (in m)')

        pars.add_argument('--north', type=float, dest='north', default='0.0',
                          help='Bearing for orientation line (in degrees)')

        pars.add_argument('--tol', type=float, dest='tol', default='0.2',
                          help='Tolerance to equate stations (in m)')

        pars.add_argument('--layer', type=str, dest='layer', default='',
                          help='Restrict conversion to a named layer')

        pars.add_argument('--cpaths', type=str, dest='cpaths', default='#ff0000',
                          help='Color of (poly)lines for export (default #ff0000)')

        pars.add_argument('--cnorth', type=str, dest='cnorth', default='#00ff00',
                          help='Color of orientation line (default #00ff00)')

        pars.add_argument('--cscale', type=str, dest='cscale', default='#0000ff',
                          help='Color of scale bar line (default #0000ff)')

    def save(self, stream):

        # Parse the SVG file which is passed as the last command line argument

        # The commented out 'svg.find' and 'svg.findall' statements below show
        # the correct way to to pass namespaces, however it appears they do not
        # work on Windows.  Therefore the actual statements used are kludges
        # which use a named argument in a string format.

        #e.options = e.arg_parser.parse_args()

        # Find out some basic properties

        #e.load_raw()
        svg = self.document.getroot()
        inkex_sodipodi = inkex.NSS['sodipodi']
        inkex_svg = inkex.NSS['svg']
        inkex_inkscape = inkex.NSS['inkscape']

        s = f'{{{inkex_sodipodi}}}docname'

        if s in svg.attrib:
            docname = svg.attrib[s]
        else:
            docname = 'untitled'

        if self.options.layer == '':
            toplevel = os.path.splitext(docname)[0]
        else:
            toplevel = self.options.layer

        # el = svg.find('.//svg:image', namespaces=inkex.NSS)
        el = svg.find(f'.//{{{inkex_svg}}}image')

        s = f'{{{inkex_sodipodi}}}absref'

        if el is not None and s in el.attrib:
            imgfile = os.path.split(el.attrib[s])[1]
        else:
            imgfile = None

        # Find all the (poly)lines in the document

        # list = svg.findall('.//svg:g/svg:path', namespaces=inkex.NSS)
        lines = svg.findall(f'.//{{{inkex_svg}}}g/{{{inkex_svg}}}path')

        # paths is a list of tuples (id, d, stroke, layer)

        paths = []

        for line in lines:
            stroke = dict(inkex.Style.parse_str(line.attrib['style']))['stroke']
            layer = line.getparent().attrib[f'{{{inkex_inkscape}}}label']
            paths.append((line.attrib['id'], line.attrib['d'], stroke, layer))

        # We extract the information from the paths and raise a PathError
        # exception if there is a problem.  This is caught (at the end) and
        # the message (in msg) printed to stderr.

        if not paths:
            msg = 'No paths found at all!'
            raise PathError(msg)

        # Find the orientation line.

        subpaths = [path for path in paths if path[2] == self.options.cnorth]

        if not subpaths:
            msg = f'No orientation line found of color {self.options.cnorth}'
            raise PathError(msg)

        # Construct the unit vector (nx, ny) to point along N, and the unit
        # (ex, ey) to point along E.  We correct for north later.

        steps = Path(subpaths[0][1]).to_arrays()
        dx, dy, dl = measure(steps)
        nx, ny = dx/dl, dy/dl
        ex, ey = -ny, nx

        # Find the scale bar line

        subpaths = [path for path in paths if path[2] == self.options.cscale]

        if not subpaths:
            msg = f'No scale bar line found of color {self.options.cscale}'
            raise PathError(msg)

        # Calculate the scale factor

        steps = Path(subpaths[0][1]).to_arrays()
        scalelen = measure(steps)[2]
        scalefac = self.options.scale / scalelen

        # Find the exportable (poly)lines

        paths = [path for path in paths if path[2] == self.options.cpaths]

        if not paths:
            msg = f'No exportable lines found of color {self.options.cpaths}'
            raise PathError(msg)

        if self.options.layer != '':
            paths = [path for path in paths if path[3] == self.options.layer]
            if not paths:
                msg = (f'No exportable lines found of color {self.options.cpaths} '
                       + f'in layer {self.options.layer}')
                raise PathError(msg)

        # Now build the survex traverses.  Keep track of stations and
        # absolute positions to identify equates and exports.

        # Stations is a list of tuples of (x, y, traverse_name, station_id)
        # Traverses is a list of tuples of (traverse_name, legs), where
        # Legs is a list of tuples of (from_id, to_id, tape, compass)

        stations = []
        traverses = []

        for path in paths:
            legs = []
            steps = Path(path[1]).to_arrays()
            for i, step in enumerate(steps):
                stations.append((step[1][0], step[1][1], path[0], i))
                if i == 0:
                    continue
                dx, dy, dl = distance(steps[i-1][1], step[1])
                tape = scalefac * dl
                compass = self.options.north + math.degrees(
                    math.atan2(ex*dx+ey*dy, nx*dx+ny*dy))
                legs.append((i-1, i, tape, compass))
            traverses.append((path[0], legs))

        ntraverse = len(traverses)
        nstation = len(stations)

        # Identify the equates.  This is an O(n^2) pairwise comparison and
        # more efficient methods are available but n should not get too large:
        # for a large project it almost always bound to be a good idea to
        # break the survey up into manageable chunks, each of which can be
        # allocated its own survex file.  This can be facilitated by putting
        # different sections into different inkscape layers.

        # Equates is a list of tuples of (station, station, distance)
        # where station is a tuple (traverse_name, station_id)

        equates = []

        for pair in combinations(stations, 2):
            dl = scalefac * distance(pair[0], pair[1])[2]
            if dl < self.options.tol:
                equates.append((pair[0][2:], pair[1][2:], dl))

        # Afficianados will notice this is a job only half done.  What we have
        # generated is an (incomplete) list of equivalence relations between
        # stations.  It may be incomplete because if A is near B, and B is
        # near C, it doesn't always follow that A is near enough C to satisfy
        # the closeness criterion.  What we should really do next is build the
        # set of equivalence classes of stations, then we can generate
        # precisely n-1 *equate directives for each non trivial equivalence
        # class of size n > 1.  In fact, survex allows for mutiple stations to
        # be listed in one *equate line, so we could just generate one *equate
        # directive for each non trivial equivalence class.  However survex
        # doesn't complain if there is redundant information in the equate
        # directives so below we take the easy option of using the list of
        # equivalence relations to generate a 1:1 list of equate directives.
        # Fastidiuous people may wish to tidy this up by hand afterwards.

        # Extract the set of stations required for export from the list of
        # equates.

        # Exports is a *set* of stations where a station is a tuple
        # (traverse_name, station_id)

        exports = set()

        for equate in equates:
            exports.add(equate[0])
            exports.add(equate[1])

        # Exportd is a dictionary to keep track of stations which should be
        # exported from a given traverse.  The key is the traverse name.  The
        # value is a list of stations to export.  If there are no stations to
        # export then the list is empty (rather than there not being a key).

        exportd = dict()
        for traverse in traverses:
            exportd[traverse[0]] = []

        for traverse_name, station_id in exports:
            exportd[traverse_name].append(station_id)

        # If we made it this far we're ready to write the survex file

        print(f'; survex file autogenerated from {docname}')

        if imgfile is not None:
            print(f'; embedded image file name {imgfile}')

        print(f'; generated {strftime("%c")}\n')

        print(f'; SVG orientation: ({nx}, {ny}) is {sprintd(self.options.north)}')
        print(f'; SVG orientation: ({ex}, {ey}) is {sprintd(self.options.north + 90)}')
        print(f'; SVG scale: {scalelen} is {self.options.scale} m, scale factor = {scalefac}')
        print(f'; SVG contained {ntraverse} traverses and {nstation} stations')
        print(f'; tolerance for identifying equates = {self.options.tol} m\n')

        print(f'\n*begin {toplevel}')

        if equates:
            print()
            for equate in equates:
                print(f'*equate {equate[0][0]}.{equate[0][1]}'
                      + f'{equate[1][0]}.{equate[1][1]}'
                      + f'; separation {equate[2]:0.2f} m')

        print('\n*data normal from to tape compass clino')

        for traverse in traverses:
            print(f'\n*begin {traverse[0]}')
            if exportd[traverse[0]]:
                sorted_export_str = [str(x) for x in sorted(exportd[traverse[0]])]
                print(f'*export {" ".join(sorted_export_str)}')
            for leg in traverse[1]:
                print(f'{leg[0]:3} {leg[1]:3} {leg[2]:7.2f} {sprintd(leg[3])}  0')
            print('*end', traverse[0])

        print(f'\n*end {toplevel} \n')
        print('; end of file')

        # End of python script


if __name__ == '__main__':
    e = SurvexOutputExtension()
    e.run()
