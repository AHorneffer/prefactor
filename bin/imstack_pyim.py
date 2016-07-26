#!/usr/bin/env python

"""

IMSTACK_PYIM.PY

A script to stack images in casa format, for use in the MSSS pipeline.

Original version by G. Heald, completed 26 July 2016.

"""

import glob
import pyrap.images as pim
import pyrap.tables as pt  #needed to update image info (beam size)
import argparse
import numpy
import scipy.signal
import Beam # module from Rene Breton

"""
The "Beam" module is used here to get the parameters needed
to convolve all images to the same beam size. From Rene's docs:

    band0 = Beam(130.528, 155.099, 101.34)
    band3 = Beam(122.359, 142.744, 106.75)
    conv = Beam(157.0, 157.0, 0.)
    req_band0 = conv.Deconvolve(band0)
    req_band3 = conv.Deconvolve(band3)

Also to create a beam image

    x,y = np.mgrid[-3:3:101j,-3:3:101j]
    r = np.array([x,y]).T
    a = Beam(0.5, 1., 20.)
    ## Calculating the Gaussian beam maps
    g_a = np.exp( -0.5 * (r * (a.icov*r[...,None,:]).sum(-1)).sum(-1).T )

"""

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def main(args):

	input_images = glob.glob(args.indir + '/' + args.inpat)
	assert(len(input_images)<=8 and len(input_images)>0)
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
	big_beam = Beam.Beam(big_beam_value,big_beam_value,0.)
	for inim in input_images:
		print bcolors.OKGREEN+'Convolving image data from %s'%(inim)+bcolors.ENDC
		bmaj = beam_params[inim]['major']['value']
		bmin = beam_params[inim]['minor']['value']
		bpa = beam_params[inim]['positionangle']['value']
		req = big_beam.Deconvolve(Beam.Beam(bmin,bmaj,bpa))
		print 'Starting beam arcsec (maj,min,pa): %.2f, %.2f, %.2f'%(bmaj,bmin,bpa)
		print 'Convolving beam arcsec (maj,min,pa): %.2f, %.2f, %.2f'%(req.amaj,req.amin,req.theta)
		pix = image_increments[inim]
		req_maj_pix = req.amaj/pix
		req_min_pix = req.amin/pix
		b_pix = Beam.Beam(req_min_pix,req_maj_pix,req.theta)
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
	print bcolors.OKGREEN+'Writing image to %s'%(args.outim)+bcolors.ENDC
	im.saveas(args.outim)
	im = pim.image(args.outim)
	stack_data_out = numpy.expand_dims(numpy.expand_dims(stack_data,axis=0),axis=0)
	im.putdata(stack_data_out)
	im = 0

	# Update beam size in image info struct
	print bcolors.OKGREEN+'Updating image info (beam size)'+bcolors.ENDC
	t = pt.table(args.outim,readonly=False,ack=False)
	iminfo = t.getkeyword('imageinfo')
	iminfo['restoringbeam']['major']['value'] = big_beam_value
	iminfo['restoringbeam']['minor']['value'] = big_beam_value
	iminfo['restoringbeam']['positionangle']['value'] = 0.
	t.putkeyword('imageinfo',iminfo)
	t.close()

ap = argparse.ArgumentParser()
ap.add_argument('--indir',help='Directory for input images [.]',default='.')
ap.add_argument('--outim',help='Output image name [stacked.img]',default='stacked.img')
ap.add_argument('--inpat',help='Input image name pattern [*.restored.corr]',default='*.restored.corr')
args = ap.parse_args()
if __name__=='__main__':
	main(args)

