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
#from lofarpipe.support.data_map import MultiDataMap
#from lofarpipe.support.data_map import MultiDataProduct

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


class MultiDataProduct(DataProduct):
    """
    Class representing multiple files in a DataProduct.
    """
    def __init__(self, host=None, file=None, skip=True):
        super(MultiDataProduct, self).__init__(host, file, skip)
        if not file:
            self.file = list()
        else:
            self._set_file(file)

    def __repr__(self):
        """Represent an instance as a Python dict"""
        return (
            "{{'host': '{0}', 'file': {1}, 'skip': {2}}}".format(self.host, self.file, str(self.skip))
        )

    def __str__(self):
        """Print an instance as 'host:[filelist]'"""
        return ':'.join((self.host, str(self.file)))

    def _set_file(self, data):
        try:
            # Try parsing as a list
            if isinstance(data, list):
                self.file = data
            if isinstance(data, DataProduct):
                self._from_dataproduct(data)
            if isinstance(data, DataMap):
                self._from_datamap(data)

        except TypeError:
            raise DataProduct("No known method to set a filelist from %s" % str(file))

    def _from_dataproduct(self, prod):
        print 'setting filelist from DataProduct'
        self.host = prod.host
        self.file = prod.file
        self.skip = prod.skip

    def _from_datamap(self, inmap):
        print 'setting filelist from DataMap'
        filelist = {}
        for item in inmap:
            if not item.host in filelist:
                filelist[item.host] = []
            filelist[item.host].append(item.file)
        self.file = filelist['i am']

    def append(self, item):
        self.file.append(item)


class MultiDataMap(DataMap):
    """
    Class representing a specialization of data-map, a collection of data
    products located on the same node, skippable as a set and individually
    """
    @DataMap.data.setter
    def data(self, data):
        if isinstance(data, DataMap):
            mdpdict = {}
            data.iterator = DataMap.SkipIterator
            for item in data:
                if not item.host in mdpdict:
                    mdpdict[item.host] = []
                mdpdict[item.host].append(item.file)
            mdplist = []
            for k, v in mdpdict.iteritems():
                mdplist.append(MultiDataProduct(k, v, False))
            self._set_data(mdplist, dtype=MultiDataProduct)
        elif isinstance(data, MultiDataProduct):
            self._set_data(data, dtype=MultiDataProduct)
        elif not data:
            pass
        else:
            self._set_data(data, dtype=MultiDataProduct)

    def split_list(self, number):
        mdplist = []
        for item in self.data:
            for i in xrange(0, len(item.file), number):
                chunk = item.file[i:i+number]
                mdplist.append(MultiDataProduct(item.host, chunk, item.skip))
        self._set_data(mdplist)

