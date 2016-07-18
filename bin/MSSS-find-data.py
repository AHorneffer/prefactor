#! /usr/bin/env python
"""
Script to look for the input MSs for one MSSS field and to generate
the needed mapfiles
"""
import pyrap.tables as pt
import sys, os
import glob
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
from lofarpipe.support.data_map import MultiDataMap
from lofarpipe.support.data_map import MultiDataProduct

def get_direction_deg(mspath):
    # if that doesn't work, then update your casacore installation
    import casacore.tables as pt
    import numpy as np
    pointing_table = pt.table(mspath+'::FIELD', ack=False)
    radec = pointing_table.getcell('PHASE_DIR',0)[0]
    print radec
    RA = radec[0]/np.pi*180.
    DEC = radec[1]/np.pi*180.
    return (RA, DEC)
    

def main(field_name=None, input_directory=None, mapfile_basename=None, mapfile_dir=None,
         uvmin_lowdec=0.1, uvmin_highdec=0., dec_border=35., uvmax_klambda=666.):
    """
    Find the MSs for one MSSS field

    Parameters
    ----------
    field_name : str
      Name of the field, e.g. "H017+39"
    input_directory : str
      Directory where the data can be found. The structure has to be:
      <input_directory>/<field_name>/<obs_number>/BAND[0-7]/L*.MS

    Returns
    -------
    result : dict
        Dict with the names of the generated mapfiles
    """

    if field_name==None or input_directory==None or mapfile_basename==None or mapfile_dir==None:
        print 'MSSS-find-data: field_name, input_directory, mapfile_basename, and mapfile_dir are needed!'
        raise ValueError('MSSS-find-data: field_name, input_directory, mapfile_basename, and mapfile_dir are needed!')

    fieldpath = os.path.join(input_directory,field_name)
    grouped_map = MultiDataMap()
    file_single_map = DataMap([])
    for i in range(8):
        band_pattern = '%s/[0-1]*/BAND%d/L*.MS' % (fieldpath,i) 
        ms_files = glob.glob(band_pattern)
        if len(ms_files) > 0:
            grouped_map.append(MultiDataProduct('localhost', ms_files, False))
            for filename in ms_files:
                file_single_map.append(DataProduct('localhost', filename, False))
    if len(file_single_map) < 2:
        raise ValueError('MSSS-find-data: found less than 2 inputs files for field!')

    (ra, dec) = get_direction_deg(file_single_map[0].file)

    if dec <= dec_border:
        uvmin = uvmin_lowdec
    else:
        uvmin = uvmin_highdec

    grouped_mapname = os.path.join(mapfile_dir, mapfile_basename+'_grouped')
    grouped_map.save(grouped_mapname)
    file_single_mapname = os.path.join(mapfile_dir, mapfile_basename+'_single')
    file_single_map.save(file_single_mapname)

    result = {'groupedmap' : grouped_mapname, 'single_mapfile' : file_single_mapname,
              'RA' : ra , 'DEC' : dec, 'UVmin' : uvmin, 
              'UVmin_lambda' : (uvmin*1000.), 'UVmax_lambda' : (uvmax_klambda*1000.)}
    print "MSSS-find-data.py result:",result
    return result
