#!/usr/bin/env python

# Pre-processing script that reads Geant4 output from a particular source
# and saves it to a uniform-format FITS file encoding the energy deposited
# in WFI pixels for each primary particle. This is for the MIT-generated
# Geant4 output, which should be sorted into directories for different
# runs.  The directory name is expected as an argument, and it must contain
# sets of files of the following naming convention:
#
# *_EvtLog_*.gdat
# *_StepLog_*.gdat
#
# The EvtLog file contains summary information of each generated primary,
# and the StepLog file contains the energy deposition steps.
#
# EDM Tue Mar  9 13:11:35 EST 2021
# Added SECPARTYPE, which records the particle type of the particle that
# entered the detector and ultimately led to the energy deposition. This
# was mostly added for the MIT version, and so for OU we just copy 
# PARTYPE to the SECPARTYPE column, since (I think) that's what it is 
# anyway. Otherwise the event finder would fail.
#
# EDM Mon Feb  8 14:36:42 EST 2021
# Updated pixel ACTX/ACTY offsets for new OU data (vs. the data provided 
# by Jonathan Keelan in 2018). Also now assumes same 512x512 active region
# for each DEPFET as the MIT version, which comes from 
# the document 'DEPFET_sensitive_area_v7.pdf'.
#
# EDM Thu Jan 14 14:46:40 EST 2021
# Now ignores pixels with no signal.
#
# EDM Fri Jan  8 14:44:29 EST 2021
# First version to read OU format of Geant4 output. Adapted from 
# 'convert_to_fits_mit.py'.  Changed output column EDEP to ENERGY. 
# Probably will break stuff.

import numpy as np
import astropy
import os, sys, glob, re
from astropy.table import Table
from astropy.io import fits
import pandas as pd

num_args = len(sys.argv) - 1    
if num_args < 1:
    sys.exit(f"Usage: {sys.argv[0]} <path/to/Geant4 input file> ...")

# functions
def match(regex: str, string: str):
    return re.compile(regex).search(string)
def wfits(data, fitsfile, ow=True, hdr=None, hdrcomments=None):
    if isinstance(data, dict):
        pass
    elif isinstance(data, tuple):
        t = Table(data[1], names = data[0])
        hdu = fits.table_to_hdu(t)
        if hdr:
            for key, val in hdr.items():
                hdu.header[key] = val
        if hdrcomments:
            for key, val in hdrcomments.items():
                hdu.header.comments[key] = val
    else:
        hdu = fits.PrimaryHDU()
        hdu.data = data
    hdu.writeto(fitsfile, overwrite = ow)

# These things are Geant4-source-specific, so should be
# written to header of rawpix file. Probably.
sphere_radius = 70. # radius of boundary sphere in cm
#numprims_gen = int(1e6) # num. primaries per run
                         # now calculated by rounding off largest primid to one sig digit
array_size = int(1e7)
primary_type = 'proton'

# Updated EDM Mon Feb  8 14:24:12 EST 2021
# Active DEPFET regions are the following (same as MIT):
#       ACTX        ACTY
#   D   59:570      52:563
#   C   59:570      610:1121
#   A   632:1143    52:563
#   B   632:1143    610:1121
# Current OU coords run: X=-592:+582, Y=-807:+397, inclusive.
# This appears to contain 1175x1205 130-µm pixels, as expected
# from 'DEPFET_sensitive_area_v7.pdf'. So I think I just need
# to shift them so LLC is at X,Y=0,0, and then filter out the
# same active DEPFET parts as in the MIT version.
x_offset = 592
y_offset = 807
# Old OU offsets (for J. Keelan data from 2018):
#x_offset = 513
#y_offset = 513
# Not used but left here for safe-keeping:
#imgsize = 1027
#actmin = 0
#actmax = 1026
#xydep_min = -513
#xydep_max = 513

# dictionary of particle types
# these are the same codes as the PDG in MPE data
# and are used for both primary type and particle
# (primary or secondary) that deposits energy;
# can be bit-wise added for multiple secondaries
# depositing in the same pixel for the same primary
ptypes = { 'other' : 0,
           'proton' : 1,
           'electron' : 2,
           'gamma' : 4,
           'e+' : 8,
           'pi-' : 16,
           'pi+' : 32,
           'neutron' : 64,
           'alpha' : 128 }

# loop through input filenames
for filename in sys.argv[1:] :

    if not re.search('.*[0-9]+_detector0(\.gz)*$', filename) : 
        print(f'### Error: file {filename} does not look like an OU Geatn4 file, skipping.')
        continue
    if not os.path.isfile(filename) :
        print(f'### Error reading file {filename}, skipping.')
        continue

    path = os.path.dirname(filename)

    this_runid = int( re.sub('.*?([0-9]+)_detector0(\.gz)*$', r'\1', filename) )
    outfile = os.path.join(path, f'rawpix_{this_runid}.fits')
    print(f'### Converting run {this_runid}.')
    print(f'### Input file {filename}.')
    print(f'### Output file {outfile}.')

    # there should be 4 Geant4 input files, one for each quadrant
    infiles = glob.glob(f"{path}/{this_runid}_detector?")
    print(f'{infiles}')
    if len(infiles) != 4:
        print (f"### Found something other than 4 datafiles for run {this_runid}, skipping.")
        continue

    # allocate output arrays, which are pixel-based
    pix_primid = np.zeros(array_size, dtype=np.uint32)
    pix_detid = np.zeros(array_size, dtype=np.uint8) + 255 # just to make sure this gets set
    pix_actx = np.zeros(array_size, dtype=np.int16)
    pix_acty = np.zeros(array_size, dtype=np.int16)
    pix_energy = np.zeros(array_size, dtype=np.float32)
    pix_partype = np.zeros(array_size, dtype=np.uint16)
    pix_primtype = np.zeros(array_size, dtype=np.uint8) + ptypes[primary_type]
    pix_runid = np.zeros(array_size, dtype=np.uint16) + this_runid

    this_pix = 0

    # loop through the four quadrant data files for this run
    # and combine all four, since single primary can produce signal
    # in multiple quadrants
    for infile in infiles:

        print(f"### Reading {infile}")

        rc = match('[0-9]+_detector([0-9])', infile) #extracts the detector name
        this_detector = int(rc.group(1))

        ptype = {}

        with open(infile, 'r') as IN:
            # step through the input file and accumulate primaries
            for line in IN:
                if match('^\s*#', line): #skip comments
                    continue
                if match('^\s*$', line): #skip blank lines
                    continue
                if not match(',', line): #could be if ',' not in line
                    continue

                fields = line.rstrip().split(',')

                # if the first column is a string, then this is a particle line
                if match('[a-zA-Z]', fields[0]): 
                    # retain the primary for this interaction
                    this_primid = int(float(fields[1]))
                    # particle type is a dict so that the pixel-specific read can pick them up;
                    # they are indexed by the secondary particle ID, and it doesn't matter if the
                    # particle ID is re-used from primary to primary, since this will reset it
                    ptype[int(fields[2])] = fields[0]
                # if the first column is a number, then this is a pixel hit line
                else: 
                    # skip it if it's outside the 512x512 region of a quad
                    this_x, this_y = int(fields[0]), int(fields[1])
                    #if this_x<xydep_min or this_y<xydep_min or this_x>xydep_max or this_y>xydep_max:
                    #    continue 
                    # skip it if the energy deposited is zero
                    this_energy = float(fields[2])
                    if (this_energy <= 0) :
                        continue

                    pix_actx[this_pix] = this_x + x_offset
                    pix_acty[this_pix] = this_y + y_offset
                    pix_energy[this_pix] = this_energy
                    pix_detid[this_pix] = this_detector
                    pix_primid[this_pix] = this_primid
                    # ptype is a dict of particle type strings indexed by the id
                    # ptypes is (constant) dict of my own particle type IDs indexed
                    # by the string (confused yet?)
                    this_tmp = ptype.get(int(fields[3]))
                    #print(f'{this_primid}, {fields[3]}, {this_tmp}')
                    pix_partype[this_pix] = ptypes.get( ptype.get(int(fields[3]), 'other'), 0 )

                    this_pix += 1

    # done primary-by-primary processing

    # eliminate signal outside of the 512-pixel active region of each DEPFET,
    # as defined above
    #       ACTX        ACTY
    #   D   59:570      52:563
    #   C   59:570      610:1121
    #   A   632:1143    52:563
    #   B   632:1143    610:1121
    # remove outer ring
    indx = (pix_actx>=59) & (pix_actx<=1143) & (pix_acty>=52) & (pix_acty<=1121)
    pix_primid = pix_primid[indx]
    pix_detid = pix_detid[indx]
    pix_actx = pix_actx[indx]
    pix_acty = pix_acty[indx]
    pix_energy = pix_energy[indx]
    pix_partype = pix_partype[indx]
    pix_primtype = pix_primtype[indx]
    pix_runid = pix_runid[indx]
    # remove the gaps
    indx = ( (pix_actx<=570) | (pix_actx>=632) ) & ( (pix_acty<=563) | (pix_acty>=610) )
    pix_primid = pix_primid[indx]
    pix_detid = pix_detid[indx]
    pix_actx = pix_actx[indx]
    pix_acty = pix_acty[indx]
    pix_energy = pix_energy[indx]
    pix_partype = pix_partype[indx]
    pix_primtype = pix_primtype[indx]
    pix_runid = pix_runid[indx]

    # lop off unused parts of arrays
    numpix = this_pix
    pix_primid = pix_primid[0:numpix]
    pix_detid = pix_detid[0:numpix]
    pix_actx = pix_actx[0:numpix]
    pix_acty = pix_acty[0:numpix]
    pix_energy = pix_energy[0:numpix]
    pix_partype = pix_partype[0:numpix]
    pix_primtype = pix_primtype[0:numpix]
    pix_runid = pix_runid[0:numpix]

    # theoretically, the number of primaries generated for a run
    # should be the highest primid rounded to one significant digit
    def RoundToSigDig (x, a) :
        return round( x, (-int(np.floor(np.log10(abs(x)))) + (a-1)) )
    numprims_gen = int( RoundToSigDig(pix_primid.max(), 1) )

    # done loop through quadrant data files for this run
    uniq_primid = np.unique(pix_primid)
    numprims_interact = uniq_primid.size
    pct_interact = 100. * numprims_interact / numprims_gen
    print(f"### Run {this_runid}: generated {numprims_gen} primaries.")
    print(f"### Run {this_runid}: found {numprims_interact} primaries ({pct_interact}%) that interacted.")
    print(f"### Run {this_runid}: found {numpix} pixels with deposited energy")

    # add header keywords
    hdr = { 
            'SPH_RAD' : sphere_radius,
            'NPRI_GEN' : numprims_gen,
            'NPRI_INT' : numprims_interact }
    hdrcomments = {
            'SPH_RAD' : 'Radius of Geant4 source sphere in cm',
            'NPRI_GEN' : 'Number of primaries generated',
            'NPRI_INT' : 'Number of primaries producing signal' }

    # make a table and save it to a FITS HDU
    wfits( (['PRIMID', 'DETID', 'ACTX', 'ACTY', 'ENERGY', 'PARTYPE', 'SECPARTYPE', 'PRIMTYPE', 'RUNID'],
        [pix_primid.astype(np.uint32),
        pix_detid.astype(np.uint8),
        pix_actx.astype(np.int16),
        pix_acty.astype(np.int16),
        pix_energy.astype(np.single),
        pix_partype.astype(np.uint8),
        pix_partype.astype(np.uint8),
        pix_primtype.astype(np.uint8),
        pix_runid.astype(np.uint64)]), 
        outfile, hdr=hdr, hdrcomments=hdrcomments)

exit()
