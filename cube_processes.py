'''
contains all the functions related to processing a cube and parsing data

'''

import scipy.ndimage
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord as coord


def index_to_radec(xs, ys, hdr, verbose=True):
    '''
    turns arrays of indecies to ra, dec based on header

    Arguments:
        xs {[type]} -- [description]
        ys {[type]} -- [description]
        hdr {[type]} -- [description]
    '''

    if type(xs) != np.ndarray:
        xs = np.array(xs)
    if type(ys) != np.ndarray:
        ys = np.array(ys)

    if xs.size != ys.size:
        if verbose:
            print "\t x & y array sizes different in conversion to RA-DEC"

    ras = (xs - hdr['CRPIX1']) * hdr['CDELT1'] - hdr['CRVAL1']
    decs = (ys - hdr['CRPIX2']) * hdr['CDELT2'] - hdr['CRVAL2']

    return ras, decs


def mask_lb(data, hdr, b_cutoff, toNaN=False, l_cutoff=None):
    pass


def umask_and_save(data, hdr, save_dir, file_name, radius=None):
    '''
    unsharp masking a data slice and saving it
    '''
    save_path = save_dir + file_name.rsplit('.', 1)[0] + '_umask.fits'
    if radius is None:
        umask_data = umask(data)
    else:
        umask_data = umask(data, radius)
    fits.writeto(save_path, umask_data, header=hdr)

    return umask_data


def radecs_to_lb(ras, decs):
    '''
    Transformation between lists of ras, decs, to ls, bs. Assumes ra, dec in degrees
    Conforms to astropy 0.4.3
    taken from https://github.com/seclark/FITSHandling/commit/f04a6e54c6624741e4f3077ba8ba96af620871ac

    for lb masks
    '''
    obj = coord.SkyCoord(ras, decs, unit="deg", frame="icrs")
    obj = obj.galactic

    ls = obj.l.degree
    bs = obj.b.degree

    return ls, bs


def circ_kern(diameter):
    '''
    Performs a circle-cut of given diameter on inkernel.
    Outkernel is 0 anywhere outside the window.
    taken from https://github.com/seclark/RHT/blob/master/rht.py

    for umask step
    '''
    assert diameter % 2
    r = diameter // 2   # int(np.floor(diameter / 2))
    mnvals = np.indices((diameter, diameter)) - r
    rads = np.hypot(mnvals[0], mnvals[1])
    return np.less_equal(rads, r).astype(np.int)


def umask(data, radius=15, smr_mask=None):
    '''
    unsharp masking of a data slice
    taken from https://github.com/seclark/RHT/blob/master/rht.py
    ^conversion to binary data step skipped

    radius is set to 15 by default for GALFA data: diameter of 30 arcmin
    '''
    assert data.ndim == 2

    kernel = circ_kern(2 * radius + 1)
    outdata = scipy.ndimage.filters.correlate(data, kernel)

    # Correlation is the same as convolution here because kernel is symmetric
    # Our convolution has scaled outdata by sum(kernel), so we will divide out these weights.
    kernweight = np.sum(kernel)
    subtr_data = data - outdata / kernweight

    # No values < 0
    subtr_data[np.where(subtr_data < 0.0)] = 0

    # set NaNs to 0??
    # subtr_data[np.where(subtr_data == np.NaN)] = 0

    if smr_mask is None:
        return subtr_data
    else:
        return np.logical_and(smr_mask, subtr_data)
