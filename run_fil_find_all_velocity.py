'''
Python script to run the filfind algorithm on all velocity slices of the GALFA
cube

LL2016
'''

# imports
import sys
from astropy.io import fits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import astropy.units as u #currently not used        

fil_finder_dir = '/Users/larryli/Documents/' \
                 'CC/16-17/research/GALFA_filfind/fil_finder'
sys.path.append(fil_finder_dir)

import filfind_class as filfind


# import cube
img, hdr = fits.getdata("../data/usharpbg30.fits", header=True)

num_v_channels = img.shape[0]
print "There are %d velocity channels" % num_v_channels

mask_objs = []

for x in xrange(0, num_v_channels):
	vslice = img[x, 300:1100, 1800:]
	print "\n\nworking on vslice %d" % x
	fils = filfind.fil_finder_2D(vslice, header=hdr, beamwidth=10.0, glob_thresh=20,
	                             distance=100, flatten_thresh=95, standard_width=1,
	                             size_thresh=1000)
	this_mask_obj = fils.create_mask(verbose=False, regrid=False, border_masking=True, save_png=True, run_name=str(x), output_mask_objs=True)
	# note size_thresh, adapt_thresh, smooth_size, fill_hole_size can all be set by args
	
	mask_objs.append(this_mask_obj)


	def get_mask_objs():
		pass



