#copied from /home/pratley/analyze_im.py
#modified by Georgi Kokotanekov 27.09.2014 to compute the noise in a crentral rectangular region with sides half of the dimensions of the image, i.e. a region of area 1/4 of the image
#before it used a narrow stripe-like region from y=0 to y_max and the central 1/8 of x-size of the image, i.e. using 1/8 of the area of the image.

#! /usr/bin/env python

from pyrap import images as pim
import numpy as np
import sys
#import warnings
import matplotlib

matplotlib.use('Agg')

import pylab
#from scipy.optimize import curve_fit



#f = 16
#f = 4

#warnings.filterwarnings("ignore",category=DeprecationWarning)

def myfit(x,y, fn):
	# Find the width of a Gaussian distribution by computing the second moment of 
	# the data (y) on a given axis (x)
	
	w = np.sqrt(abs(sum(x**2*y)/sum(y)))
	
	# Or something like this. Probably a Gauss plus a power-law with a lower cutoff is necessary
	#func = lambda x, a, b: np.exp(0.5*x**2/a**2) + x**(-1.*b)	 

	#[popt, pvar] = curve_fit(func, x, y)
	# Need a newer version of scipy for this to work...

	#a = popt[0]
	#b = popt[1]
	pylab.ioff()
	pylab.semilogy(x, y, '+', label="Data")
	#pylab.semilogy(x, np.exp(-0.5*x**2/a**2) + x**(-1.*b), label="Noise fit")
	pylab.semilogy(x, np.exp(-0.5*x**2/w**2), label="Noise fit")
	pylab.title('Normalized pixel distribution')
	pylab.ylabel('Rel. Num. of Pixels')
	pylab.xlabel('Pixel brightness (Jy/beam)')
	pylab.legend(loc=3)
	pylab.savefig(fn)
	#pylab.show()
	pylab.close()
	return w



def computeInoise(Id, fn, hrange):
	
	Ih = np.histogram(Id, bins=100, range=hrange) # 0 = values, 1 = bin edges
	Ix = Ih[1][:-1] + 0.5*(Ih[1][1] - Ih[1][0])
	Iv = Ih[0]/float(max(Ih[0]))    
	Inoise = myfit(Ix, Iv, fn+'_histI.png')
	
	##GDK
	#Ih = np.histogram(Id, bins=100)
	#Ix = Ih[1][:-1] + 0.5*(Ih[1][1] - Ih[1][0])
	#Iv = Ih[0]#/float(max(Ih[0]))
	#Inoise = myfit(Ix, Iv, fn+'_histI.png')
	
	##	print 'Ih[1][:-1]',Ih[1][:-1]
	##	print 'Ih[1][1]',Ih[1][1]
	##	print 'Ih[1][0]',Ih[1][0]
	#	print 'Ix', Ix
	#	print 'Iv', Iv
	
	return Inoise


# image file name
#try:
#	fn = sys.argv[1]
#except IndexError:
#	raise Exception("You must specify an image name.")
def obtainRMS(fn):
    im = pim.image(fn)
    nfo = im.info()
    d = im.getdata()
    nstokes = d.shape[1]
    nra = d.shape[2]
    ndec = d.shape[3]

    bmaj = nfo['imageinfo']['restoringbeam']['major']['value']
    bmin = nfo['imageinfo']['restoringbeam']['minor']['value']
    barea = 2.*np.pi*bmaj*bmin/(2.3548**2)

    print "Restoring beam:"
    print "		BMaj      = %.3f %s"%(bmaj, nfo['imageinfo']['restoringbeam']['major']['unit'])
    print "		BMin      = %.3f %s"%(bmin, nfo['imageinfo']['restoringbeam']['minor']['unit'])
    print "		BPA       =  %.3f %s"%(nfo['imageinfo']['restoringbeam']['positionangle']['value'], nfo['imageinfo']['restoringbeam']['positionangle']['unit'])
    print "		Beam Area = %.3f sq. %s"%(barea, nfo['imageinfo']['restoringbeam']['major']['unit'])

    #print 'nra', nra
    
    f = 4

    #Id = d[0,0, (nra/2 - nra/f):(nra/2 + nra/f)].flatten() #default
    #Id = d[0,0, (nra/2 - nra/f):(nra/2 + nra/f), (ndec/2 - ndec/4):(ndec/2 + ndec/4)].flatten()
    Id = d[0,0, (nra/2 - nra/f):(nra/2 + nra/f), (ndec/2 - ndec/f):(ndec/2 + ndec/f)].flatten()  #central square #888
    
    
    #print 'Id', len(Id)
    
    
    #Id_center = d[0,0, (nra/2 - nra/f):(nra/2 + nra/f)].flatten()
    ##Id_center = d[0,0, (nra/2 - nra/f):(nra/2 + nra/f),  (nra/2 - nra/f):(nra/2 + nra/f)]#.flatten()

    ##print 'Id_center len',len(Id_center) 
    ##print 'Id_center ',Id_center 
    
    ##Id_out = d[0,0, (nra - 2*nra/f):nra].flatten()#right side
    #Id_out = d[0,0, 0:2*nra/f].flatten()#left side
    
    ##print 'Id_out ',len(Id_out)
    #Id_all = d[0,0, 0:nra].flatten()
    ##print 'Id_all ', Id_all
    

    hrange = (-1,1)    
    Inoise = computeInoise(Id, fn, hrange)
    
    print "Image noise, measured within the inner 1/16th of the area: "
    print "		Stokes I noise = %.3f %s"%(Inoise, nfo['unit'])

















    if nstokes==4:
	    Qd = d[0,1, (nra/2 - nra/f):(nra/2 + nra/f)].flatten()
	    Ud = d[0,2, (nra/2 - nra/f):(nra/2 + nra/f)].flatten()
	    Vd = d[0,3, (nra/2 - nra/f):(nra/2 + nra/f)].flatten()

    if nstokes==4:
	    hrange = (-0.1, 0.1)
	    Qh = np.histogram(Qd, bins=100,range=hrange) # 0 = values, 1 = left bin edges
	    Qx = Qh[1][:-1] + 0.5*(Qh[1][1] - Qh[1][0])
	    Qv = Qh[0]/float(max(Qh[0]))
	    Uh = np.histogram(Ud, bins=100, range=hrange) # 0 = values, 1 = left bin edges
	    Ux = Uh[1][:-1] + 0.5*(Uh[1][1] - Uh[1][0])
	    Uv = Uh[0]/float(max(Uh[0]))
	    Vh = np.histogram(Vd, bins=100, range=hrange) # 0 = values, 1 = left bin edges
	    Vx = Vh[1][:-1] + 0.5*(Vh[1][1] - Vh[1][0])
	    Vv = Vh[0]/float(max(Vh[0]))

	    Qnoise = myfit(Qx, Qv, fn+'_histQ.png')
	    Unoise = myfit(Ux, Uv, fn+'_histU.png')
	    Vnoise = myfit(Vx, Vv, fn+'_histV.png')
	    print "		Stokes Q noise = %.3f %s"%(Qnoise, nfo['unit']) 
	    print "		Stokes U noise = %.3f %s"%(Unoise, nfo['unit']) 
	    print "		Stokes V noise = %.3f %s"%(Vnoise, nfo['unit'])
    else:
	    Qnoise=0
	    Unoise=0
	    Vnoise=0
	    
	    
    return Inoise, Qnoise, Unoise, Vnoise
#obtainRMS(fn)


def main(imagename):
	Inoise = obtainRMS(imagename)[0] 
	threshold = str(2.0*Inoise)+'Jy'
	return {'threshold_Jy' : threshold }
