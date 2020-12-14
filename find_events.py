#!/usr/bin/env python

# Post-processing script that finds events and particle tracks in
# Geant4 data. Reads generic Geant4 FITS pixel lists generated by
# pre-processing script 'convet_to_fits_*.py'. Input is a list of raw
# pixel FITS files, including the relative path if needed. The
# invocation should be of the form:
#
#   find_events.py path/rawpix_<runid1>.fits [ path/rawpix_<runid2>.fits ]
#
# where <runid> is a sequence of numbers indicating a unique Geant4
# run. The script will process each FITS file in turn, sorting
# primaries into radomly selected frames of specified length, and then
# finding and characterizing events and particle tracks ('blobs'). It
# will write out:
#
#   path/geant_events_<runid>_frames.fits
#       Summary information for each frame.
#   path/geant_events_<runid>_pix.fits
#       List of all pixels with signal. Basically a copy of the
#       input but with frame and blob info added.
#   path/geant_events_<runid>_evt.fits
#       List of events, which are identified as 3x3 islands around
#       a local maximum. 
#   path/geant_events_<runid>_blobs.fits
#       List of particle tracks, which are contiguous islands of 
#       5 or more pixels (set by npixthresh), or any number of pixels
#       including at least one pixel above the MIP threshold.
#
# EDM Mon Dec 14 11:04:27 EST 2020
# First generic version that works. Adapted from 'doit.pl' and
# 'doit.py' scripts that run on the Open University Geant4 data.

import numpy as np
from skimage import measure
from astropy.io import fits
from astropy.table import Table
import sys, glob, os
import datetime
import re


num_args = len(sys.argv) - 1    
if num_args < 1:
    sys.exit(f"Usage: {sys.argv[0]} <path/to/raw pixel FITS file> ...")

# These things are Geant4-source-specific, so should be
# written to header of rawpix file. Probably.
#sphere_radius = 20. # radius of boundary sphere in cm
#num_protons_per_run = 16000 # number of proton primaries 
#sphere_radius = 70. # radius of boundary sphere in cm
#num_protons_per_run = 1.e6 # number of proton primaries 


# define defaults
texp_frame = .005   # frame exposure time in sec
evtth = .1          # 100 eV for WFI Geant4 simulations
splitth = evtth    # same as evtth for WFI Geant4 simulations
npixthresh = 5      # minimum number of pixels in a blob
mipthresh = 15.     # minimum ionizing particle threshold in keV
clip_energy = 22.   # maximum pixel value reported by WFI, in keV
skip_writes = -1    # writes FITS images for every skip_writes primary; 
                    # set to -1 to turn off writes
evperchan = 1000.   # why not? PHA here is really PI
mipthresh_chan = mipthresh * 1000. / evperchan  # minimum ionizing particle threshold in PHA units
spec_maxkev = 100.
numchans = int((spec_maxkev*1000.) / evperchan)
gain_intercept = 0. # use this for Geant4 data
gain_slope = 1000.  # use this for Geant4 data (pixel PH units are keV)
#gain_intercepts = (0, 0, 0, 0)         # in ADU
#gain_slopes = (1000, 1000, 1000, 1000) # in eV/ADU
# rate and frame defaults
proton_flux = 4.1 * 7./5.   # protons/s/cm2; 7/5 accounts for alphas, etc.
detector_area = 4. * (130.e-4 * 512.)**2 # detector area in cm2, 
                                              # assuming 130 um pixels

"""
Comments...
"""

epicpn_pattern_table = np.array([
     0, 13,  3, 13, 13, 13, 13, 13,  4, 13,  8, 12, 13, 13, 13, 13,
     2, 13,  7, 13, 13, 13, 11, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
     1, 13, 13, 13, 13, 13, 13, 13,  5, 13, 13, 13, 13, 13, 13, 13,
     6, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13,  9, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    10, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
    13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13 
])

# dictionary of particle types
# these are the same codes as the PDG in MPE data
# and are used for both primary type and particle
# (primary or secondary) that deposits energy;
# can be bit-wise added for multiple secondaries
# depositing in the same pixel for the same primary
ptypes = { 'proton' : 1,
             'e-' : 2,
             'gamma' : 4,
             'e+' : 8,
             'pi-' : 16,
             'pi+' : 32,
             'neutron' : 64,
             'alpha' : 128 }

# initialize rng
rng = np.random.RandomState(1234)

# temporary variables
x, y = 0, 0

# 0-511; indexed by detector number 0-3
x_offset = 513
y_offset = 513
imgsize = 1027
actmin = 0
actmax = 1026
xydep_min = -513
xydep_max = 513

def match(regex: str, string: str):
    return re.compile(regex).search(string)
def indexND(array, ind, value=None):
    if value is None:
        return array[ind[1],ind[0]]     
    array[ind[1],ind[0]] = value        
def wfits(data, fitsfile, ow=True, hdr=None):
    if isinstance(data, dict):
        pass
    elif isinstance(data, tuple):
        t = Table(data[1], names = data[0])
        hdu = fits.table_to_hdu(t)
        if hdr:
            for key, val in hdr.items():
                hdu.header[key] = val
    else:
        hdu = fits.PrimaryHDU()
        hdu.data = data
    hdu.writeto(fitsfile, overwrite = ow)
def which(condition):
    if condition.ndim != 1:
        condition = condition.flat
    return np.where(condition)[0]
def whichND(condition):
    return np.array(np.where(condition)[::-1])
def which_both(condition):
    return which(condition), which(condition == False)

date = datetime.datetime.now().astimezone()
date_string = date.strftime("%a %d %b %Y %I:%M:%S %p %Z").rstrip()
print("############################################################")
print(f"### Started {sys.argv[0]} on {date_string}")

#######################################
# Main loop.
# Step through input raw pixel FITS files.
# For each one, create a frame per primary,
# find blobs and MIPs, and then find events in that frame.

for filename in sys.argv[1:] :

    if not re.search('.*rawpix_[0-9]+\.fits$', filename) : 
        print(f'Error: file {filename} does not look like a rawpix FITS')
        continue
    if not os.path.isfile(filename) :
        print(f'Error reading file {filename}')
        continue

    infile = filename
    path = os.path.dirname(filename)
    this_run = int( re.sub(r'.*rawpix_([0-9]+)\.fits$', r'\1', infile) )
    print(f"### Processing {infile}.")

    # read in raw pixel energy deposition FITS file
    # some of this will be written out the pixel list (which will
    # also have frame info)
    rawpix = Table.read(infile)
    xdep = rawpix['ACTX']
    ydep = rawpix['ACTY']
    endep = rawpix['EDEP']
    rundep = rawpix['RUNID']
    detectordep = rawpix['DETID']
    eiddep = rawpix['PRIMID']
    ptypedep = rawpix['PDG']
    framedep = np.zeros(xdep.size, dtype=int)
    piddep = np.zeros(xdep.size, dtype=int)
    blobid = np.zeros(xdep.size, dtype=int)

    # get header info
    sphere_radius = rawpix.meta['SPH_RAD']
    num_protons_per_run = rawpix.meta['NPRI_GEN']

    # get number of pixels with signal; this is the maximum number
    # we'll use in allocating pixel, event, and blob arrays
    numpix_with_signal = xdep.size

    # get unique primids and number of primaries that interact
    uniq_eid = np.unique(eiddep)
    numprimaries = uniq_eid.size

    # total exposure time in sec for this run
    texp_run = num_protons_per_run/3.14159/proton_flux/(sphere_radius**2)
    # average num. primaries per frame
    mu_per_frame = num_protons_per_run * texp_frame / texp_run 

    print(f"### There are {num_protons_per_run} primaries in this run.")
    print(f"### The run exposure time is {texp_run} sec.")
    print(f"### The frame exposure time is {texp_frame} sec")
    print(f"### for an expected mean of {mu_per_frame} primaries per frame.")

    # initialize event piddles, which will be written out or used later
    evt_runid = np.zeros(numpix_with_signal, dtype=np.uint64)
    evt_detectorid = np.zeros(numpix_with_signal, dtype=np.int32)
    evt_primid = np.zeros(numpix_with_signal, dtype=np.int32)
    evt_actx = np.zeros(numpix_with_signal, dtype=int)
    evt_acty = np.zeros(numpix_with_signal, dtype=int)
    evt_phas = np.zeros((numpix_with_signal,25), dtype=float)
    evt_pha = np.zeros(numpix_with_signal, dtype=float)
    evt_ptype = np.zeros(numpix_with_signal, dtype=np.uint8)
    evt_energy = np.zeros(numpix_with_signal, dtype = float)     # in keV
    evt_evttype = np.zeros(numpix_with_signal, dtype=np.int16)
    evt_blobdist = np.zeros(numpix_with_signal, dtype=float)
    evt_mipdist = np.zeros(numpix_with_signal, dtype=float)
    evt_pattern = np.zeros(numpix_with_signal, dtype=np.uint8)
    evt_vfaint = np.zeros(numpix_with_signal, dtype=np.uint8)
    # assign frames and times to each primary
    # to start, assume mean of one primary per second
    evt_time = np.zeros(numpix_with_signal, dtype=float)
    evt_frame = np.zeros(numpix_with_signal, dtype=np.int32)
    pix_time = np.zeros(numpix_with_signal, dtype=float)
    pix_frame = np.zeros(numpix_with_signal, dtype=int)

    # initialize structures to hold the secondary particle columns
    # piddles for the numeric columns; these are enough for now
    #run = np.zeros(0, dtype=int)        # Geant4 run (* in *_detector[0123])
    #detector = np.zeros(0, dtype=int)   # WFI detector (? in *_detector?)
    #eid = np.zeros(0, dtype=int)        # primary ID
    #particleid = np.zeros(0, dtype=int) # interacting particle ID
    #parentid = np.zeros(0, dtype=int)   # don't really need probably

    # initialize piddles to hold blob-specific things to go in FITS table
    blob_frame = np.zeros(numpix_with_signal, dtype=int)
    blob_blobid = np.zeros(numpix_with_signal, dtype=int)
    blob_cenx = np.zeros(numpix_with_signal, dtype=float)
    blob_ceny = np.zeros(numpix_with_signal, dtype=float)
    blob_cenxcl = np.zeros(numpix_with_signal, dtype=float)
    blob_cenycl = np.zeros(numpix_with_signal, dtype=float)
    blob_npix = np.zeros(numpix_with_signal, dtype=int)
    blob_energy = np.zeros(numpix_with_signal, dtype=float)
    blob_energycl = np.zeros(numpix_with_signal, dtype=float)

    # initialize things for the running frames which we will
    # randomly populate
    # frame settings
    # we know there are $num_protons_per_run, so generate enough
    # random frames to hold them
    framedist = rng.poisson(mu_per_frame, int(2*num_protons_per_run/mu_per_frame))

    cumframedist = framedist.cumsum()
    # get the total number of frames needed to capture all the primaries; 
    # will write this to FITS header so we can combine runs
    numtotframes = which(cumframedist >= num_protons_per_run)[0] + 1
    # this is wrong, because it will remove the last bin which we need
    # it's also unnecessary
    #cumframedist = cumframedist[cumframedist <= num_protons_per_run]

    # initialize piddles to hold frame-specific things to go in FITS table
    frame_frame = np.zeros(numtotframes, dtype=int)
    frame_time = np.zeros(numtotframes, dtype=float)
    frame_runid = np.zeros(numtotframes, dtype=int)
    frame_npix = np.zeros(numtotframes, dtype=int)
    frame_npixmip = np.zeros(numtotframes, dtype=int)
    frame_nevt = np.zeros(numtotframes, dtype=int)
    frame_nevtgood = np.zeros(numtotframes, dtype=int)
    frame_nevt27 = np.zeros(numtotframes, dtype=int)
    frame_nblob = np.zeros(numtotframes, dtype=int)
    frame_nprim = np.zeros(numtotframes, dtype=int)

    # running variables
    numevents = 0
    numtotblobs = 0

    # loop through unique primaries to sort them into frames
    # also figure out the number of primaries with signal in 
    # multiple quadrants, just as a diagnostic
    num_in_different_quadrants = [0, 0, 0, 0]
    """alternate
    for primary in uniq_eid:
    """
    for i in range(numprimaries):
        primary = uniq_eid[i]
        indx = which(eiddep==primary)
        #print("#####################")
        #print(f"Doing primary {primary}")
        #print(f"indx: {indx}")
        #print(f"eiddep: {eiddep[indx]}")
        #print(f"OLD framedep {framedep[indx]}")
        # assign each primary to a frame
        # first get the frame ID (indexed starts at 0)
        #print(cumframedist[-1])
        frame = which(cumframedist >= primary)[0]
        #print(f"THIS IS FRAME {frame}")
        # then set the primary and pixel piddles
        framedep[indx] = frame
        #print(f"NEW framedep framedep[indx]")
        num_quadrants = np.unique(detectordep[indx]).size
        num_in_different_quadrants[num_quadrants-1] += 1

    print(f"### Run {this_run}: number of primaries in 1 2 3 4 quads: {num_in_different_quadrants}")
    # min and max X,Y values

    # figure out the unique frames that are populated by 
    # primaries that interacted
    frame_frame = np.unique(framedep)
    numframes = frame_frame.size
    # now we can allocate frame piddles since we know how many there are
    frame_runid = np.full(numframes, this_run)
    frame_time = np.zeros(numframes, dtype=float)
    frame_npix = np.zeros(numframes)
    frame_npixmip = np.zeros(numframes)
    frame_nevt = np.zeros(numframes)
    frame_nevtgood = np.zeros(numframes)
    frame_nevt27 = np.zeros(numframes)
    frame_nblob = np.zeros(numframes)
    frame_nprim = np.zeros(numframes)

    pct_interact = 100. * numprimaries / num_protons_per_run
    print(f"### Run {this_run}: generated {numtotframes} total frames,")
    print(f"### of which {numframes} frames with the {numprimaries}")
    print(f"### interacting primaries will be written.")
    print(f"### {num_protons_per_run} total primaries were simulated.")
    print(f"### {numprimaries} or {pct_interact}% of these produced a WFI interaction.")

    # loop through frames and make a raw frame for each

    for i in range(numframes):
        # set the frame ID
        this_frame = frame_frame[i]

        # set the frame time
        frame_time[i] = this_frame * texp_frame

        # keep us updated
        if i % 100 == 0:
            print(f"### Run {this_run}: done {i} of {numframes} frames")
            pass

        # are we writing out?
        if skip_writes > 0 and i % skip_writes == 0:
            writeit = 1
        else:
            writeit = 0

        #############################################
        # make a raw frame
        #x = np.zeros(0, float)
        #y = np.zeros(0, float)
        #en = np.zeros(0, float)
        #pixptype = np.zeros(0, float)
        pixel_indx = which(framedep==this_frame)
        x = np.copy(xdep[pixel_indx])
        y = np.copy(ydep[pixel_indx])
        en = np.copy(endep[pixel_indx])
        pixptype = np.copy(ptypedep[pixel_indx])
        pixprimid = np.copy(eiddep[pixel_indx])
        # we want to populate this so don't sever it
        """comments..."""
        xoff = np.amin(x) + x_offset - 2
        yoff = np.amin(y) + y_offset - 2
        x -= np.amin(x)
        x += 2
        y -= np.amin(y)
        y += 2
        xdim = np.amax(x) + 3
        ydim = np.amax(y) + 3

        #line 542

        # img is the energy (pulseheight) image
        # ptypeimg is an image encoding the particle type responsible
        # primidimg is an image encoding the primary responsible
        # for each pixel
        img = np.zeros((ydim, xdim), float) 
        ptypeimg = np.zeros((ydim, xdim), float)
        primidimg = np.zeros((ydim, xdim), float)
        img_xvals = np.fromfunction(lambda x,y: y, (ydim, xdim), dtype = int)
        img_yvals = np.fromfunction(lambda x,y: x, (ydim, xdim), dtype = int)
        #for j in range(en.size):
        #    img[x[j], y[j]] += en[j]
        #    ptypeimg[x[j], y[j]] += en[j]
        # better way to do the above mapping to an image without a loop
        # indexND wants a 2xN piddle, where N is the number of pixels
        coos = np.array((x,y))

        indexND(img, coos, en)
        indexND(ptypeimg, coos, pixptype)
        indexND(primidimg, coos, pixprimid)

        #print(pixprimid, x, y, en, sep='\n')
        #print(img)

        # add to some frame piddles
        frame_npix[i] = pixel_indx.size
        frame_npixmip[i] = which(en >= mipthresh).size
        frame_nprim[i] = np.unique(pixprimid).size

        # write out the frame as FITS image
        if writeit:
            fitsfile = os.path.join(path, f'run{this_run}_frame{this_frame}_img.fits')
            print(f"### Writing raw image {fitsfile}.")
            wfits(img, fitsfile)

        #############################################
        # find blobs

        # segment image into blobs
        # original IDL code
        # blobRegions=label_region(phimg * (phimg GE evtthresh), /ulong)
        # NLR PDL 'cc8compt' is equivalentish to IDL 'label_regions'
        blobimg = measure.label(img > evtth, connectivity=2)

        # the first blob is 1, not 0 (which indicates no blob)
        # blobadjust decrements blob IDs when we chuck one, so the
        # IDs are continuous from 1
        blobadjust = 0
        for j in range(1, np.amax(blobimg) + 1):
            indx = which(blobimg == j)
            #indx2d = whichND(blobimg == j) dsn
            # set the blob to zeros and skip it if there are too
            # few elements; if it contains a MIP, it's a good blob
            # irregardless
            if (indx.size < npixthresh and np.amax(img.flat[indx]) < mipthresh):
                blobimg.flat[indx] = 0
                blobadjust += 1
                continue
            blobimg.flat[indx] -= blobadjust
            this_blobindx = numtotblobs + j - blobadjust - 1

            # this is the running blobid which we need to add to
            # blob piddles
            blob_blobid[this_blobindx] = numtotblobs + j - blobadjust
            blob_frame[this_blobindx] = this_frame
            blob_npix[this_blobindx] = indx.size
            # calculate unclipped blob centroid and summed energy
            tmp_en = img.flat[indx]
            tmp_toten = np.sum(tmp_en)
            tmp_wtd_x = img_xvals.flat[indx] * tmp_en
            tmp_wtd_y = img_yvals.flat[indx] * tmp_en
            blob_cenx[this_blobindx] = xoff + tmp_wtd_x.sum() / tmp_toten
            blob_ceny[this_blobindx] = yoff + tmp_wtd_y.sum() / tmp_toten
            blob_energy[this_blobindx] = tmp_toten
            # calculate clipped blob centroid and summed energy
            tmp_en = np.clip(tmp_en, None, clip_energy)
            tmp_toten = tmp_en.sum()
            tmp_wtd_x = img_xvals.flat[indx] * tmp_en
            tmp_wtd_y = img_yvals.flat[indx] * tmp_en
            blob_cenxcl[this_blobindx] = xoff + tmp_wtd_x.sum() / tmp_toten
            blob_cenycl[this_blobindx] = yoff + tmp_wtd_y.sum() / tmp_toten
            blob_energycl[this_blobindx] = tmp_toten

        # record number of blobs in this frame
        frame_nblob[i] = np.amax(blobimg)

        # if we found some blobs, change their IDs so they reflect the
        # running total for this run, increase that running blob number 
        blobimg[blobimg > 0] += numtotblobs
        numtotblobs = np.amax(blobimg) if np.amax(blobimg) > 0 else numtotblobs

        # mark each pixel in a blob
        blobid[pixel_indx] = indexND(blobimg, coos)

        # reset the blobimg to ones where there are blobs, zeros
        # otherwise; this way all the blobs have distance 0 automagically
        blobimg[blobimg > 0] = 1

        #############################################
        # find MIPs
        # set up some piddles
        #mip_xy = np.zeros((2, 0), float);

        # segment image into mips
        mipimg = np.copy(img)
        mipimg[mipimg < mipthresh] = 0

        # reset the mipimg to ones where there are mips, zeros
        # otherwise; this way all the mips have distance 0 automagically
        mipimg[mipimg > 0]  = 1

        #############################################
        # find events
        # based on Bev's 'findev.pro' rev3.1, 2011-01-26
        # get indices and number of pixels with PH above event threshold
        # need both 1D and 2D indices, latter to ease local max finding
        indx_thcross = which(img > evtth)
        indx2d_thcross = whichND(img > evtth)
        num_thcross = indx_thcross.size

        #line 676  
        if num_thcross > 0:
            # get piddles containing X,Y coords of threshold crossings
            evtx = indx2d_thcross[0]
            evty = indx2d_thcross[1]
            """coments"""
            tmp_evtx = evtx + xoff
            tmp_evty = evty + yoff
            indx_border, indx_notborder = which_both(
                (tmp_evtx<=1) | ((tmp_evtx>=510) & (tmp_evtx<=516)) | (tmp_evtx>=1025) | 
                (tmp_evty<=1) | ((tmp_evty>=510) & (tmp_evty<=516)) | (tmp_evty>=1025) | 
                (img[evty,evtx] < 0)
            )
            num_notborder = indx_notborder.size

            if num_notborder > 0:
                evtx = evtx[indx_notborder]
                evty = evty[indx_notborder]
                indx_thcross = indx_thcross[indx_notborder]
                # find local maxima
                # make a copy of the threshold crossing piddle, which we will
                # subsitute with the largest neighbor value and then compare to the original
                localmax = np.copy(img)
                """comments..."""

                for xi,yi in [(-1,0),(-1,1),(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1)]:             
                    indx = which(img[evty + yi,evtx + xi] > localmax[evty,evtx])
                    localmax.flat[indx_thcross[indx]] = img[evty + yi, evtx + xi].flat[indx]

                # finally compare the original central pixel pulseheight
                # with the maximum neighbor; if the former is greater than
                # or equal to the latter, it's a local maximum
                indx_localmax = which(img.flat[indx_thcross] >= localmax.flat[indx_thcross])
                num_localmax = indx_localmax.size

                #line 748
                if (num_localmax > 0):
                    evtx = evtx[indx_localmax]
                    evty = evty[indx_localmax]
                    indx_thcross = indx_thcross[indx_localmax]
                    # get 3x3 PHAS piddle in correct
                    # 1D order (central pixel first)
                    this_phas = np.zeros((num_localmax, 25), float)

                    phas_order = iter(range(1, 25))
                    this_phas[:,0] = img[evty, evtx]  
                    for yi in range(-1,2):
                        for xi in range(-1,2):
                            if (xi, yi) != (0,0):
                                this_phas[:,next(phas_order)] = img[evty + yi,evtx + xi]              

                    # get outer 5x5 PHAS piddle
                    # first pad the image with zeros all around
                    tmp_img = np.pad(img, 1, mode='constant')
                    # now evtx and evty are too low by one, so correct for that in index
                    for yi in range(-1,4):
                        for xi in range(-1,4):
                            if not(-1 < xi < 3 and -1 < yi < 3):
                                this_phas[:, next(phas_order)] = tmp_img[evty + yi,evtx + xi]

                    # set the particle type for all these events
                    # based on whatever ptype produced the local max
                    this_ptype = np.copy(ptypeimg[evty,evtx])
                    this_primid = np.copy(primidimg[evty,evtx])
                else:
                    continue
            else:
                continue
        else:
            continue

        ## get X,Y coords of any pixels in a blob (don't matter which blob)
        indx_blobs = whichND(blobimg > 0).T

        ## get X,Y coords of any pixels in a mip (don't matter which mip)
        indx_mips = whichND(mipimg > 0).T

        numevents_thisframe = num_localmax
        #print(f"### Found {numevents_thisframe} events, processing.")

        # process the detected events
        # append events to running piddles
        indx = np.arange(numevents,numevents+numevents_thisframe)
        evt_runid[indx] = this_run
        #evt_detectorid[indx] = 0
        evt_primid[indx] = this_primid
        evt_frame[indx] = this_frame
        evt_actx[indx] = evtx
        evt_acty[indx] = evty
        evt_phas[indx,:] = this_phas
        #evt_pha[indx] = np.zeros(numevents_thisframe)
        #evt_ptype[indx] = np.array(this_ptype, np.uint8)
        #evt_evttype[indx] = np.zeros(numevents_thisframe, np.int16)
        #evt_energy[indx] = np.zeros(numevents_thisframe)
        #evt_blobdist[indx] = np.zeros(numevents_thisframe)
        #evt_mipdist[indx] = np.zeros(numevents_thisframe)
        #evt_pattern[indx] = np.zeros(numevents_thisframe, np.uint8)
        #evt_vfaint[indx] = np.zeros(numevents_thisframe, np.uint8)

        # step through all events to determine
        # EPIC-pn pattern, summed PHA, VFAINT flag, minimum distance to a
        # blob and mip.
        for j in range(numevents, numevents+numevents_thisframe):
            # below is already done in event finding and pasted above
            # get X,Y of center pixel
            x = evt_actx[j]
            y = evt_acty[j]
            # below is deprecated, since we've assumed this in event finding
            ## eliminate events on edges so we can use 5x5
            #next if ($x<2 or $x>1021 or $y<2 or $y>1021);

            # get ACIS flt grade; "which" returns the indices of
            # non-central pixels greater than or equal to event threshold,
            # these are used as bits to raise 2 to the power, and summed
            # (this can probably be removed from the loop if I'm smart)
            indx = which(evt_phas[j,1:9] >= splitth)
            fltgrade = pow(2,indx).sum()
            # convert to EPIC-pn pattern (from look-up table)
            evt_pattern[j] = epicpn_pattern_table[fltgrade]

            # sum 3x3 pixels over split threshold to get PHA
            # this is an ACIS PHA, not Suzaku
            # (this can probably be removed from the loop if I'm smart)
            evt_pha[j] = evt_phas[j][evt_phas[j] >= splitth].sum()
            # apply gain correction for this node
            # get the gain parameters from the node (0-3)
            #print(f"{evt_phas[j]} {gain_intercept} {gain_slope} {evt_pha[j]}")
            evt_pha[j] = (evt_pha[j] - gain_intercept) * gain_slope / evperchan

            # convert pha to energy
            # (this can probably be removed from the loop if I'm smart)
            evt_energy[j] = evt_pha[j] * evperchan / 1000
            #print(f"{(evt_phas[j])} {gain_intercept} {gain_slope} {evt_pha[j]} {evt_energy[j]}")

            # perform simple VFAINT filtering; also update the EPIC-pn pattern
            # of doubles, triples, and quads based on it
            # get outer 16 pixels and if any are above split flag it

            noutermost = which(evt_phas[j, 9:25] > splitth).size #maybe don't need which
            if noutermost > 0:
                evt_vfaint[j] = 1
                # EDM Fri Jan  4 11:24:20 EST 2019 
                # Reinstated the 5x5 filtering on PATTERN.
                # EDM Thu May 31 13:46:28 EDT 2018 
                # change to remove 5x5 criterion on EPIC-pn patterns
                # for doubles, triples, quads
                if evt_pattern[j] > 0: 
                    evt_pattern[j] = 13

            # get minimum distance to a blob
            # first find delta X,Y from this event to list of blob pixels
            delta_blobs = indx_blobs - np.array([x,y])
            #print(delta_blobs)
            # square that, sum it, square root it to get the distance
            if delta_blobs.shape[0] > 0:
                evt_blobdist[j] = np.amin(np.sqrt(np.square(delta_blobs).sum(axis = 1)))
            # unless there aren't any blobs, in which case set blobdist to -1
            else:
                evt_blobdist[j] = -1
            #print(f"{evt_blobdist[j]}")

            # get minimum distance to a mip
            # first find delta X,Y from this event to list of mip pixels
            delta_mips = indx_mips - np.array([x,y])
            #print $delta_mips;
            # square that, sum it, square root it to get the distance
            if delta_mips.shape[0] > 0:

                evt_mipdist[j] = np.amin(np.sqrt(np.square(delta_mips).sum(axis = 1)))
            # unless there aren't any mips, in which case set mipdist to -1
            else:
                evt_mipdist[j] = -1

            #print(f"{evt_mipdist[j]}")

            # we really want ACTX,Y in real WFI coords in the event list, 
            # so fix that here; $x and $y are children of $actx and $acty
            # for this event, which is why this works
            evt_actx[j] += xoff
            evt_acty[j] += yoff

            # add info to frame piddles
            frame_nevt[i] += 1
            if evt_energy[j] < mipthresh:
                frame_nevtgood[i] += 1
            if evt_energy[j] >= 2 and evt_energy[j] < 7:
                frame_nevt27[i] += 1


        # increment the number of events and hit the next frame
        numevents += numevents_thisframe
    # done loop through frames

    # segregate the reds and greens as defined now by EPIC-pn pattern
    # reds are bad patterns only
    indx_goodpatterns, indx_badpatterns = which_both(evt_pattern < 13)
    # cyans are singles and doubles, which are "best" of the good
    indx_goodbest, indx_badbest = which_both(evt_pattern < 5);
    # combine indices for the filters we want via intersection
    indx_reds = indx_badpatterns
    indx_greens = np.intersect1d(indx_goodpatterns, indx_badbest)
    indx_cyans = np.intersect1d(indx_goodpatterns, indx_goodbest)
    # determine (event) count rate in cts/cm2/s/keV in the important band
    indx_goodenergy = which((evt_energy > 2) & (evt_energy <= 7))
    flux_goodenergy = indx_goodenergy.size / texp_run / detector_area / 5.

    evt_evttype[indx_reds] = 3
    evt_evttype[indx_greens] = 4
    evt_evttype[indx_cyans] = 6

    hdr = {'NFRAMES': numtotframes, 'NBLOBS': numtotblobs}

    # write out event file
    outevtfile = os.path.join(path, f'geant_events_{this_run}_evt.fits')
    # chop off the unused parts of the arrays
    evt_runid = evt_runid[:numevents]
    evt_detectorid = evt_detectorid[:numevents]
    evt_primid = evt_primid[:numevents]
    evt_actx = evt_actx[:numevents]
    evt_acty = evt_acty[:numevents]
    evt_phas = evt_phas[:numevents]
    evt_pha = evt_pha[:numevents]
    evt_ptype = evt_ptype[:numevents]
    evt_energy = evt_energy[:numevents]
    evt_evttype = evt_evttype[:numevents]
    evt_blobdist = evt_blobdist[:numevents]
    evt_mipdist = evt_mipdist[:numevents]
    evt_pattern = evt_pattern[:numevents]
    evt_vfaint = evt_vfaint[:numevents]
    evt_frame = evt_frame[:numevents]
#    wfits({'RUNID': runid.astype(np.uint16),
#         'DETID': detectorid.astype(np.int32),
#         'PRIMID': primid.astype(np.int32),
#         'FRAME': evt_frame.astype(np.int32),
#         'ACTX': actx.astype(np.uint16),
#         'ACTY': acty.astype(np.uint16),
#         'PHAS': phas.astype(np.float),
#         'PHA': pha.astype(np.float),
#         'PTYPE': ptype.astype(np.uint8),
#         'EVTTYPE': evttype.astype(np.int16),
#         'ENERGY': energy.astype(np.float),
#         'BLOBDIST': blobdist.astype(np.float),
#         'MIPDIST': mipdist.astype(np.float),
#         'PATTERN': pattern.astype(np.uint8),
#         'VFAINT': vfaint.astype(np.uint8)},
#           outevtfile, hdr=hdr)
    wfits(
    (['RUNID', 'DETID', 'PRIMID', 'FRAME', 'ACTX', 'ACTY', 'PHAS', 'PHA', 'PTYPE', 'EVTTYPE', 'ENERGY', 'BLOBDIST', 'MIPDIST', 'PATTERN', 'VFAINT'], [evt_runid.astype(np.int64), evt_detectorid.astype(np.int32), evt_primid.astype(np.int32), evt_frame.astype(np.int32), evt_actx.astype(np.uint16), evt_acty.astype(np.uint16), evt_phas.astype(np.float), evt_pha.astype(np.float), evt_ptype.astype(np.uint8), evt_evttype.astype(np.int16), evt_energy.astype(np.float), evt_blobdist.astype(np.float), evt_mipdist.astype(np.float), evt_pattern.astype(np.uint8), evt_vfaint.astype(np.uint8)]), outevtfile, hdr=hdr)

    # write out pixel list
    print(f"### Writing pixel list for run {this_run}.")
    outpixfile = os.path.join(path, f'geant_events_{this_run}_pix.fits')
    xdep += x_offset
    ydep += y_offset
    indx_sorted = np.argsort(framedep)
    wfits( 
    (['ACTX', 'ACTY', 'BLOBID', 'ENERGY', 'FRAME', 'PRIMID', 'PTYPE', 'RUNID'],
     [xdep[indx_sorted], ydep[indx_sorted], blobid[indx_sorted], endep[indx_sorted], framedep[indx_sorted], eiddep[indx_sorted], ptypedep[indx_sorted], rundep[indx_sorted]]), 
    outpixfile, hdr=hdr)

    # write out frame list
    print(f"### Writing frame list for run {this_run}.")
    outframefile = os.path.join(path, f'geant_events_{this_run}_frames.fits')
    indx_sorted = np.argsort(frame_frame)
    wfits(
    (['FRAME', 'NBLOB', 'NEVT', 'NEVT27', 'NEVTGOOD', 'NPIX', 'NPIXMIP', 'NPRIM', 'RUNID', 'TIME'], 
     [frame_frame[indx_sorted], frame_nblob[indx_sorted], frame_nevt[indx_sorted], frame_nevt27[indx_sorted], frame_nevtgood[indx_sorted], frame_npix[indx_sorted], frame_npixmip[indx_sorted], frame_nprim[indx_sorted], frame_runid[indx_sorted], frame_time[indx_sorted]]),
    outframefile, hdr = hdr)  

    # write out blob list
    print(f"### Writing blob list for run {this_run}.")
    outblobfile = os.path.join(path, f'geant_events_{this_run}_blobs.fits')
    # chop off the unused parts of the arrays
    blob_frame = blob_frame[:numtotblobs]
    blob_blobid = blob_blobid[:numtotblobs]
    blob_cenx = blob_cenx[:numtotblobs]
    blob_ceny = blob_ceny[:numtotblobs]
    blob_cenxcl = blob_cenxcl[:numtotblobs]
    blob_cenycl = blob_cenycl[:numtotblobs]
    blob_npix = blob_npix[:numtotblobs]
    blob_energy = blob_energy[:numtotblobs]
    blob_energycl = blob_energycl[:numtotblobs]
    # sort by blobid
    indx_sorted = np.argsort(blob_blobid)
    wfits(
    (['BLOBID', 'CENX', 'CENXCL', 'CENY', 'CENYCL', 'ENERGY', 'ENERGYCL', 'FRAME', 'NPIX'],
    [blob_blobid[indx_sorted], blob_cenx[indx_sorted], blob_cenxcl[indx_sorted], blob_ceny[indx_sorted], blob_cenycl[indx_sorted], blob_energy[indx_sorted], blob_energycl[indx_sorted], blob_frame[indx_sorted], blob_npix[indx_sorted]]),
    outblobfile, hdr = hdr)

    print(f"### Finished run {this_run}.")


# done loop through runs
