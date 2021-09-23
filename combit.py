#!/usr/bin/env python

# Script to combine the various output FITS files from 'find_events.py'.
# Input argument is a directory path in which there are sets of files of
# the following form:
#
#   geant_events_all_[0-9]+_frames.fits
#   geant_events_all_[0-9]+_pix.fits
#   geant_events_all_[0-9]+_evt.fits
#   geant_events_all_[0-9]+_blobs.fits
# 
# The script loops through the [0-9]+ run IDs and concantenates each of the
# four file types into a single, large FITS file. Some columns are
# incremented so that they monotonically increase without any gaps,
# including FRAME, TIME, and BLOBID, and the header keywords with summary
# information about the data are combined and updated as well. The output
# can be considered a simulated WFI observation obtained from the full
# Geant4 run.
#
# EDM Tue Sep 14 16:43:07 EDT 2021
# Updated to properly increment running PRIMID in events and pixels, as per
# bug reported by Dan Wilkins.
#
# EDM Mon Jul 12 11:57:19 EDT 2021
# Updated to deal with single four-extension input and output FITS files.
# Four extensions are for events, frames, blobs, and pixels.
#
# EDM Tue Dec 22 11:38:33 EST 2020
# Removed PHAS column from combined event list, since astropy.table isn't
# happy with it and it adds a huge amount to file size.
#
# EDM Mon Dec 21 17:05:43 EST 2020
# First Python version, converted from Perl and now allocating arrays
# rather than appending to improve performance.

import numpy as np
from astropy.io import fits
import astropy.table as tbl
import sys, glob, os
import datetime
import re

maxsize = int(2e6)

num_args = len(sys.argv) - 1    
if num_args != 1:
    sys.exit(f"Usage: {sys.argv[0]} <path/to/FITS files/>")

date = datetime.datetime.now().astimezone()
date_string = date.strftime("%a %d %b %Y %I:%M:%S %p %Z").rstrip()
print("############################################################")
print(f"### Started {sys.argv[0]} on {date_string}")

path = sys.argv[1]
if not os.path.isdir(path) :
    print(f'### Error reading directory {path}, exiting.')
    exit()

patt = re.compile('geant_events_[0-9]+.fits')
filenames = [ f for f in os.listdir(path) if patt.match(f) ]
#print(f'{path}')

# initialize some running totals for columns we need
# to increment
numfiles = len(filenames)
numfiles_done = 0
offset_frame = 0
offset_blobid = 0
offset_primid = 0
offset_time = 0.
currow_frames = 0
currow_pix = 0
currow_evt = 0
currow_blobs = 0

# loop through frame FITS files
for fitsfile in filenames :

    # add the path to the frame FITS filename
    fitsfile = os.path.join(path, fitsfile)
    print(f'### Doing file {numfiles_done+1} of {numfiles}, {fitsfile}.')

    # figure out the other FITS filenames and make sure they all exist
    if not os.path.isfile(fitsfile) :
        print(f'### Error reading FITS file {fitsfile}, skipping.')
        continue

    # read FITS file into tables
    hdulist = fits.open(fitsfile) 
    this_frames = tbl.Table.read(hdulist['FRAMES'])
    this_pix = tbl.Table.read(hdulist['PIXELS'])
    this_evt = tbl.Table.read(hdulist['EVENTS'])
    # might need to remove this line eventually
    this_evt.remove_column('PHAS')
    this_blobs = tbl.Table.read(hdulist['BLOBS'])

    # get number of rows
    this_numrows_frames = len(this_frames)
    this_numrows_pix = len(this_pix)
    this_numrows_evt = len(this_evt)
    this_numrows_blobs = len(this_blobs)

    # grab header keywords from frame FITS and assume they're the same in
    # other FITS
    this_sphere_radius = this_frames.meta['SPH_RAD']
    this_numprims_gen = this_frames.meta['NPRI_GEN']
    this_numprims_interact = this_frames.meta['NPRI_INT']
    this_nframes = this_frames.meta['NFRAMES']
    this_nblobs = this_frames.meta['NBLOBS']
    this_texpframe = this_frames.meta['TEXPFRAM']
    this_texp = this_texpframe * this_nframes

    # increment things
    this_frames['FRAME'] += offset_frame
    this_frames['TIME'] += offset_time
    this_pix['FRAME'] += offset_frame
    this_pix['BLOBID'] += offset_blobid
    this_pix['PRIMID'] += offset_primid
    this_evt['FRAME'] += offset_frame
    this_evt['PRIMID'] += offset_primid
    this_blobs['FRAME'] += offset_frame
    this_blobs['BLOBID'] += offset_blobid

    # if this is the first file
    if (numfiles_done == 0) :

        # create the output tables
        # assume each input file has ~ the same number of rows, and allocate 2x that
        # for the output files
        # not sure a better way to grab the dtypes than what I have

        # frames
        numrows = 2 * numfiles * this_numrows_frames
        print(f'### Generating {numrows} rows for frames FITS')
        dtypes = [ this_frames.dtype.base[col] for col in this_frames.colnames ]
        all_frames = tbl.Table(np.zeros([numrows,len(this_frames.colnames)]), names=this_frames.colnames,
                dtype=dtypes)
        all_frames[currow_frames:this_numrows_frames] = this_frames

        # pix
        numrows = 2 * numfiles * this_numrows_pix
        print(f'### Generating {numrows} rows for pix FITS')
        dtypes = [ this_pix.dtype.base[col] for col in this_pix.colnames ]
        all_pix = tbl.Table(np.zeros([numrows,len(this_pix.colnames)]), names=this_pix.colnames,
                dtype=dtypes)
        all_pix[currow_pix:this_numrows_pix] = this_pix

        # evt
        numrows = 2 * numfiles * this_numrows_evt
        print(f'### Generating {numrows} rows for evt FITS')
        dtypes = [ this_evt.dtype.base[col] for col in this_evt.colnames ]
        all_evt = tbl.Table(np.zeros([numrows,len(this_evt.colnames)]), names=this_evt.colnames,
                dtype=dtypes)
        all_evt[currow_evt:this_numrows_evt] = this_evt

        # blobs
        numrows = 2 * numfiles * this_numrows_blobs
        print(f'### Generating {numrows} rows for blobs FITS')
        dtypes = [ this_blobs.dtype.base[col] for col in this_blobs.colnames ]
        all_blobs = tbl.Table(np.zeros([numrows,len(this_blobs.colnames)]), names=this_blobs.colnames,
                dtype=dtypes)
        all_blobs[currow_blobs:this_numrows_blobs] = this_blobs

        # also get global header info
        all_sphere_radius = this_sphere_radius
        all_numprims_gen = this_numprims_gen
        all_numprims_interact = this_numprims_interact
        all_nframes = this_nframes
        all_nblobs = this_nblobs
        all_texpframe = this_texpframe
        all_texptot = this_nframes * this_texpframe

    else :
        # otherwise increment and insert
        all_frames[currow_frames:currow_frames+this_numrows_frames] = this_frames
        all_pix[currow_pix:currow_pix+this_numrows_pix] = this_pix
        all_evt[currow_evt:currow_evt+this_numrows_evt] = this_evt
        all_blobs[currow_blobs:currow_blobs+this_numrows_blobs] = this_blobs

        # also compare header info
        if (this_sphere_radius != all_sphere_radius) :
            print(f'### WARNING: FITS header SPH_RAD {this_sphere_radius} != {all_sphere_radius}.')
        if (this_texpframe != all_texpframe) :
            print(f'### WARNING: FITS header TEXPFRAM {this_texpframe} != {all_texpframe}.')

        # also increment global header info
        all_numprims_gen += this_numprims_gen
        all_numprims_interact += this_numprims_interact
        all_nframes += this_nframes
        all_nblobs += this_nblobs
        all_texptot += this_nframes * this_texpframe

    # increment the incrementers 
    numfiles_done += 1
    offset_frame += this_nframes
    offset_time += this_texp
    offset_blobid += this_nblobs
    offset_primid += this_numprims_gen
    currow_frames += this_numrows_frames
    currow_pix += this_numrows_pix
    currow_evt += this_numrows_evt
    currow_blobs += this_numrows_blobs

# chop off unused parts of the arrays
all_frames = all_frames[0:currow_frames]
all_pix = all_pix[0:currow_pix]
all_evt = all_evt[0:currow_evt]
all_blobs = all_blobs[0:currow_blobs]

# set header/metadata for output FITS files
hdr = fits.Header()
hdr.set('SPH_RAD', all_sphere_radius, 'Radius of Geant4 source sphere in cm')
hdr.set('NPRI_GEN', all_numprims_gen, 'Number of primaries generated')
hdr.set('NPRI_INT', all_numprims_interact, 'Number of primaries producing signal')
hdr.set('NFRAMES', all_nframes, 'Number of frames in this file')
hdr.set('NBLOBS', all_nblobs, 'Number of blobs in this file')
hdr.set('TEXPFRAM', all_texpframe, 'Frame exposure time in sec')
hdr.set('TEXPTOT', all_texptot, 'Total exposure time in sec')

# create an empty primary HDU
empty_primary = fits.PrimaryHDU(header=hdr)
# create HDUs for each table extension
events_hdu = fits.BinTableHDU(all_evt, name='EVENTS', header=hdr)
frames_hdu = fits.BinTableHDU(all_frames, name='FRAMES', header=hdr)
blobs_hdu = fits.BinTableHDU(all_blobs, name='BLOBS', header=hdr)
pix_hdu = fits.BinTableHDU(all_pix, name='PIXELS', header=hdr)

# write out the combined FITS file
outfitsfile = os.path.join(path, 'geant_events_all.fits')
print(f"### Writing combined output FITS file {outfitsfile}.")
hdul = fits.HDUList([empty_primary, events_hdu, frames_hdu, blobs_hdu, pix_hdu])
hdul.writeto(outfitsfile, output_verify='silentfix', overwrite=True, checksum=True)

exit()


