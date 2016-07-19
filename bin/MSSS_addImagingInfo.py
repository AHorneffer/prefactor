#!/usr/bin/env python

# addImagingInfo: Python script to add meta info to a CASA image
# Copyright (C) 2012
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id: makebeamtables.cc 18133 2011-05-30 13:53:45Z diepen $
#
# @author Ger van Diepen


import lofar.addImagingInfo as addimg


if __name__ == "__main__":
    import sys
    import getopt
    (opts,args) = getopt.getopt(sys.argv[1:], 't:b:w')
    if len(args) < 5:
        print "Insufficient arguments; run as:"
        print "   addImagingInfo [-t taql] [-b baseline] [-w]"
        print "                  image sourcedb maxbl ms1 ms2 ..."
        print "      taql         optional TaQL selection string used in imager"
        print "                   It always adds: ANTENNA1!=ANTENNA2"
        print "      baseline     optional baseline selection used in imager"
        print "                   (in CASA syntax)"
        print "      -w           Min and max baseline length are given in"
        print "                   wavelengths, not in meters."
        print "                   They are converted to meters using the highest"
        print "                   (for min) and lowest (for max) frequency"
        print "      image        name of the image"
        print "      sourcedb     name of SourceDB containing the sources found"
        print "                   If empty name, no LOFAR_SOURCE subtable is made"
        print "      minbl        minimum baseline length (in m or wvl) used in imager"
        print "                   > 0   always use this value"
        print "                   0     determine from the first MS"
        print "                   < 0   use max(abs(minbl), ms_value)"
        print "      maxbl        same for maximum baseline length"
        print "      ms1 ms2 ...  names of MSs the image is made of"
        print "                   the names can be individual arguments or a"
        print "                   comma-separated list of names (or a mix)"
        print "   E.g.,"
        print "      addImagingInfo myimg mysdb 0. 0. ms1,ms2 ms3 ms4,ms5"
        sys.exit(1)
    taqlStr  = ""
    baseline = ""
    aswvl = False
    for (opt,arg) in opts:
        if opt == '-t':
            taqlStr = arg
        elif opt == '-b':
            baseline = arg
        elif opt == '-w':
            aswvl = True
    msNames = []
    minbl = float(args[2])
    maxbl = float(args[3])
    sourcedb = args[1].strip(' \"')
    for arg in args[4:]:
        msNames.extend (arg.split(","))
    try:
        addimg.addImagingInfo (args[0], msNames, sourcedb, minbl, maxbl,
                               aswvl, taqlStr, baseline);
    except Exception as raised_stuff:
        if str(raised_stuff) == "addImagingInfo already done (keyword ATTRGROUPS already exists)":
            print "addImagingInfo thinks it is already done!"
        else:
            print "Please tell ASTRON to fix the raising of exceptions in \"addImagingInfo.py\", see Issue #9694"
            raise
