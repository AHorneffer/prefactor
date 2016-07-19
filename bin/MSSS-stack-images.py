#! /usr/bin/env python
"""
Script to sort a list of MSs into frequency-bands, and compute additional values needed for initsubtract
"""
import pyrap.tables as pt
import sys, os
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct


def main(input_files, outname):
    """
    Check a list of MS files for missing frequencies

    Parameters
    ----------
    input_files : list or str
        List of input filenames, or string with list, or path to a mapfile
    outname: str
        Name of output file

    Returns
    -------
    result : dict
        Dict with the name of the generated file

    """
    if type(input_files) is str:
        if input_files.startswith('[') and input_files.endswith(']'):
            file_list = [f.strip(' \'\"') for f in input_files.strip('[]').split(',')]
        else:
            map_in = DataMap.load(input_files)
            map_in.iterator = DataMap.SkipIterator
            file_list = []
            for fname in map_in:
                if fname.startswith('[') and fname.endswith(']'):
                    for f in fname.strip('[]').split(','):
                        file_list.append(f.strip(' \'\"'))
                else:
                    file_list.append(fname.strip(' \'\"'))  
    elif type(input_files) is list:
        file_list = [str(f).strip(' \'\"') for f in input_files]
    else:
        raise TypeError('MSSS-stack-images: type of "input_files" unknown!')

    outfile = file_list[(len(file_list)/2)]

    result = {'stacked_image' : outfile }
    return result


