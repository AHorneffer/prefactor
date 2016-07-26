#! /usr/bin/env python
"""
Script to sort a list of MSs into frequency-bands, and compute additional values needed for initsubtract

Beam class by Rene Breton (2016)
Stacking function by George Heald (2016)

"""
import pyrap.tables as pt
import sys, os
import numpy as np
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct

class Beam(object):
    """
    Create a Gaussian beam object given a semi-major and semi-minor axis, and
    a position angle.

    Parameters
    ----------
    amin : float
        Semi-minor axis
    amaj : float
        Semi-major axis
    theta : float
        Position angle, in degrees, measure from the first axis.

    Examples
    --------
        >>> beam = Beam(0.5, 1., 0.)
        >>> beam.amin
            0.
        >>> beam.amaj
            1.
        >>> beam.theta
            0.
        >>> beam.cov
            array([[ 0.25 ,  0.],
                   [ 0.   ,  1.]])

        >>> beam = Beam(0.5, 1., 90.)
        >>> beam.cov
            array([[ 1. ,  0.],
                   [ 0. ,  0.25]])

        >>> bean = Beam(0.5, 1., 90.)
        >>> beam.cov
            array([[ 0.625,  0.375],
                   [ 0.375,  0.625]])

    """
    def __init__(self, amin, amaj, theta):
        self.amin = amin
        self.amaj = amaj
        self.theta = theta
        return

    @property
    def cov(self):
        theta = self.theta / 180 * np.pi
        amaj = self.amaj
        amin = self.amin
        cov = np.array( [[amin**2*np.cos(theta)**2 + amaj**2*np.sin(theta)**2,
            (amaj**2-amin**2)*np.cos(theta)*np.sin(theta)],
            [(amaj**2-amin**2)*np.cos(theta)*np.sin(theta),
            amin**2*np.sin(theta)**2 + amaj**2*np.cos(theta)**2]] )
        return cov

    @property
    def icov(self):
        return np.linalg.inv(self.cov)

    @property
    def norm(self):
        cov = self.cov
        norm = 1 / ( 2*np.pi * np.sqrt(cov[0,0]*cov[1,1]-2*cov[0,1]) )
        return norm
    
    def Convolve(self, beam):
        """
        Return the convolution of the current beam with the provided beam.

        According to the convolution theorem, the convolution of two Gaussians
        is simply another Gaussian as:
        N(mu1, Sigma1) * N(mu2, Sigma2) = N(mu1+mu2, Sigma1+Sigma2)
        """
        cov = self.cov + beam.cov
        return Beam(*Cov2Polar(cov))

    def Deconvolve(self, beam):
        """
        Return the deconvolution of the current beam by the provided beam.

        According to the convolution theorem, the deconvolution of two
        Gaussians is simply another Gaussian as it is the inverse of the
        convolution:
        N(mu1, Sigma1) * N(mu2, Sigma2) = N(mu1+mu2, Sigma1+Sigma2)
        """
        cov = self.cov - beam.cov
        return Beam(*Cov2Polar(cov))


def Cov2Polar(cov):
    """
    Given a covariance matrix extract the semi-major and semi-minor axes, and the
    position angle.

    Parameters
    ----------
    cov : array (2x2)
        Covariance matrix

    Outputs
    -------
    amin : float
        Semi-minor axis
    amaj : float
        Semi-major axis
    theta : float
        Position angle, in degrees, measure from the first axis.
    """
    theta = 0.5 * np.arctan2(2*cov[0,1], cov[1,1]-cov[0,0]) * 180 / np.pi
    amin = np.sqrt( 0.5 * (cov[0,0] + cov[1,1] - np.sqrt((cov[1,1]-cov[0,0])**2 + 4*cov[0,1]**2) ) )
    amaj = np.sqrt( 0.5 * (cov[0,0] + cov[1,1] + np.sqrt((cov[1,1]-cov[0,0])**2 + 4*cov[0,1]**2) ) )
    return amin, amaj, theta

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def stack_images(input_images,output_file):
    import pyrap.images as pim
    import pyrap.tables as pt  #needed to update image info (beam size)
    import numpy
    import scipy.signal
    
    beam_params = {} # to store beam params for each image for later use
    image_data = {} # store all image data
    image_increments = {} # pixel size of each image
    decval={} # declination of each image (should be the same)
    max_bmaj = 0.
    for inim in input_images:
        print bcolors.OKGREEN+'Opening %s'%(inim)+bcolors.ENDC
        im = pim.image(inim)
        beam_params[inim] = im.imageinfo()['restoringbeam']
        if im.imageinfo()['restoringbeam']['major']['value'] > max_bmaj:
            max_bmaj = im.imageinfo()['restoringbeam']['major']['value']
        # The next line presumes single-frequency, single-Stokes!
        image_data[inim]=im.getdata()[0,0,:,:]
        coords = im.coordinates()
        d=coords.get_coordinate('direction')
        decval[inim]=d.get_referencevalue()[0]
        image_increments[inim]=coords.get_increment()[2][0]*206265.

    print bcolors.OKGREEN+'Maximum beam size is %f arcsec'%(max_bmaj)+bcolors.ENDC
    big_beam_value = numpy.ceil(max_bmaj)+1.
    print bcolors.OKGREEN+'Target beam size will be %f arcsec'%(big_beam_value)+bcolors.ENDC
    # In future, possibly try to reject images with particularly large beam values.

    # Here, convolve data arrays to big common beam
    big_beam = Beam(big_beam_value,big_beam_value,0.)
    for inim in input_images:
        print bcolors.OKGREEN+'Convolving image data from %s'%(inim)+bcolors.ENDC
        bmaj = beam_params[inim]['major']['value']
        bmin = beam_params[inim]['minor']['value']
        bpa = beam_params[inim]['positionangle']['value']
        req = big_beam.Deconvolve(Beam(bmin,bmaj,bpa))
        print 'Starting beam arcsec (maj,min,pa): %.2f, %.2f, %.2f'%(bmaj,bmin,bpa)
        print 'Convolving beam arcsec (maj,min,pa): %.2f, %.2f, %.2f'%(req.amaj,req.amin,req.theta)
        pix = image_increments[inim]
        req_maj_pix = req.amaj/pix
        req_min_pix = req.amin/pix
        b_pix = Beam(req_min_pix,req_maj_pix,req.theta)
        scale_factor_factor = 1.
        print 'Convolving beam pixels (maj,min,pa): %.2f, %.2f, %.2f'%(b_pix.amaj,b_pix.amin,b_pix.theta)
        xr = numpy.ceil(req.amaj/pix)*5. # convolving beam defined out to +/- 5 sigma
        x,y = numpy.mgrid[-xr:xr+1,-xr:xr+1]
        r = numpy.array([x,y]).T
        g_a = b_pix.norm/scale_factor_factor * numpy.exp( -0.5 * (r * (b_pix.icov*r[...,None,:]).sum(-1)).sum(-1).T )
        image_data[inim] = scipy.signal.fftconvolve(image_data[inim],g_a,mode='same')

    # Here, obtain noise (proxy) for inverse-variance weighting
    print bcolors.OKGREEN+'Computing noise proxy for each convolved image'+bcolors.ENDC
    image_weights = {}
    sumweight = 0.
    for inim in input_images:
        #image_weights[inim]=1.
        noise_proxy = numpy.median(numpy.absolute(image_data[inim]-numpy.median(image_data[inim])))
        print 'Noise proxy for %s is %f'%(inim,noise_proxy)
        image_weights[inim]=1./noise_proxy**2
        sumweight += image_weights[inim]

    # Now stack
    stack_data = numpy.zeros(image_data[input_images[0]].shape)
    for inim in input_images:
        stack_data += image_weights[inim]*image_data[inim]/sumweight

    # Write out to disk
    print bcolors.OKGREEN+'Writing image to %s'%(output_file)+bcolors.ENDC
    im.saveas(output_file)
    im = pim.image(output_file)
    stack_data_out = numpy.expand_dims(numpy.expand_dims(stack_data,axis=0),axis=0)
    im.putdata(stack_data_out)
    im = 0

    # Update beam size in image info struct
    print bcolors.OKGREEN+'Updating image info (beam size)'+bcolors.ENDC
    t = pt.table(output_file,readonly=False,ack=False)
    iminfo = t.getkeyword('imageinfo')
    iminfo['restoringbeam']['major']['value'] = big_beam_value
    iminfo['restoringbeam']['minor']['value'] = big_beam_value
    iminfo['restoringbeam']['positionangle']['value'] = 0.
    t.putkeyword('imageinfo',iminfo)
    t.close()



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

    stack_images(file_list, outname)

    result = {'stacked_image' : outname }
    return result


