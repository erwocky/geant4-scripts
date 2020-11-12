#!/usr/bin/env python

# Pre-processing script that reads Geant4 output from a particular source
# and saves it to a uniform-format FITS file encoding the energy deposited
# in WFI pixels for each primary particle. This is for the MIT-generated
# Geant4 output.
#
# EDM Wed Nov  4 17:36:39 EST 2020
# First version to read Rick's format of Geant4 output. Adapted from MPE readit.py.

import uproot
import numpy as np
import astropy
import os, sys, glob, re
from astropy.table import Table, vstack
from astropy.io import fits
import pandas as pd

num_args = len(sys.argv) - 1    
if num_args < 1:
    sys.exit(f"Usage: {sys.argv[0]} <path/to/raw pixel FITS file>")

# Geant4 output data should be sorted into directories for different runs.
# The directory name is expected as an argument, and it must contain sets
# of files of the following naming convention:
#
# *_EvtLog_*.gdat
# *_StepLog_*.gdat
#
# The EvtLog file contains summary information of each generated primary,
# and the StepLog file contains the energy deposition steps.


# size of a ionizing particle sub-step sample in µm
substep_size = 5.

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

filename = sys.argv[1]
if not re.search('.*_StepLog_[0-9]+\.gdat$', filename) : 
    raise IOError(f'Error: file {filename} does not look like a StepLog file')
if not os.path.isfile(filename) :
    raise IOError(f'Error reading file {filename}')

path = os.path.dirname(filename)

g4stepfile = filename
g4evtfile = re.sub('StepLog', 'EvtLog', g4stepfile)
runid = int( re.sub(r'.*_([0-9]+)_StepLog_([0-9]+)\.gdat', r'\1\2', g4stepfile) )
outfile = os.path.join(path, f'rawpix_{runid}.fits')
print(f'### Converting run {runid}.')
print(f'### Step file   {g4stepfile}.')
print(f'### Event file  {g4evtfile}.')
print(f'### Output file {outfile}.')

# open the Geant4 files
# below doesn't work if there are string columns
#g4evt = np.loadtxt(g4evtfile, skiprows=2)
#g4step = np.loadtxt(g4stepfile, skiprows=2)
g4evt = pd.read_csv(g4evtfile, sep='\s+')
g4step = pd.read_csv(g4stepfile, sep='\s+')

# sort out Geant4 event table
evt_primid = g4evt['Event']
evt_particletype = g4evt['Particle']
evt_energy = g4evt['Energy']

# sort out Geant4 step table
primid = g4step['Event']
depfet = g4step['Volume']
parent = g4step['Parent']
pid = g4step['ID']
step = g4step['Step']
slen = g4step['SLen']
edep = g4step['Edep']
particletype = g4step['Particle']
prex_mm = g4step['LPre-X']
prey_mm = g4step['LPre-Y']
prez_mm = g4step['LPre-Z']
postx_mm = g4step['LPost-X']
posty_mm = g4step['LPost-Y']
postz_mm = g4step['LPost-Z']
process = g4step['Process']

# number of rows in the step file
numrows = g4step.shape[0]

# total number of primaries generated
numprims_gen = evt_primid.size

# get unique primids and number of primaries that interact
uniq_primid = np.unique(primid)
numprims_interact = uniq_primid.size

# allocate output arrays, which are pixel-based
pix_primid = np.zeros(numrows, dtype=np.uint32)
pix_detid = np.zeros(numrows, dtype=np.uint8) + 8   # just to make sure this gets set
pix_rawx = np.zeros(numrows, dtype=np.uint16)
pix_rawy = np.zeros(numrows, dtype=np.uint16)
pix_actx = np.zeros(numrows, dtype=np.uint16)
pix_acty = np.zeros(numrows, dtype=np.uint16)
pix_edep = np.zeros(numrows, dtype=np.float32)
pix_pdg = np.zeros(numrows, dtype=np.uint16)
pix_primtype = np.zeros(numrows, dtype=np.uint8)
pix_runid = np.zeros(numrows, dtype=np.uint16) + runid

# Convert Rick's LOCAL coords (in mm per depfet) to the full focal plane.
# I'm using Michael Hubbard's sensitive areas from 
# 'sensitive_area/DEPFET_sensitive_area_v7.pdf', with my annotations.
#
#print(f'{prex_mm.max()}, {postx_mm.max()}, {prey_mm.max()}, {posty_mm.max()}')
# 'D' is at origin, so coords don't need to be shifted
indx = (depfet == 'DEPFETA')
#pix_detid[indx] = 3
# 'A' is shifted only in X
indx = (depfet == 'DEPFETA')
#pix_detid[indx] = 0
prex_mm.loc[indx]  = prex_mm.loc[indx]  + 7.64 + 31.36 + 35.3 + 3.7
postx_mm.loc[indx] = postx_mm.loc[indx] + 7.64 + 31.36 + 35.3 + 3.7
# 'C' is shifted only in Y
indx = (depfet == 'DEPFETC')
#pix_detid[indx] = 2
prey_mm.loc[indx]  = prey_mm.loc[indx]  + 6.72 + 31.36 + 35.3 + 2.78
posty_mm.loc[indx] = posty_mm.loc[indx] + 6.72 + 31.36 + 35.3 + 2.78
# 'B' is shifted in both X and Y
indx = (depfet == 'DEPFETB')
#pix_detid[indx] = 1
prex_mm.loc[indx]  = prex_mm.loc[indx]  + 7.64 + 31.36 + 35.3 + 3.7
postx_mm.loc[indx] = postx_mm.loc[indx] + 7.64 + 31.36 + 35.3 + 3.7
prey_mm.loc[indx]  = prey_mm.loc[indx]  + 6.72 + 31.36 + 35.3 + 2.78
posty_mm.loc[indx] = posty_mm.loc[indx] + 6.72 + 31.36 + 35.3 + 2.78
#print(f'{prex_mm.max()}, {postx_mm.max()}, {prey_mm.max()}, {posty_mm.max()}')

prex = np.floor(prex_mm/.13).astype(int)
prey = np.floor(prey_mm/.13).astype(int)
postx = np.floor(postx_mm/.13).astype(int)
posty = np.floor(posty_mm/.13).astype(int)

#print(f'{prex.max()}, {postx.max()}, {prey.max()}, {posty.max()}')

delta_x = postx - prex
delta_y = posty - prey

#pix_primid = pix[b'framenb']
#pix_edep = pix[b'edep']
#pix_localx = pix[b'localx']
#pix_localy = pix[b'localy']
#pix_pdg = pix[b'pdg']
## need these initialized since we populate them before binning
#pix_runid = np.zeros(pix_primid.size, dtype=np.ubyte)
#pix_actx = np.zeros(pix_primid.size, dtype=np.ushort)
#pix_acty = np.zeros(pix_primid.size, dtype=np.ushort)
# things below will be created later, from the binned-up pixel
# data, which could have smaller size
#pix_detid = np.zeros(pix_primid.size, dtype=np.ubyte)
#pix_primtype = np.zeros(pix_primid.size, dtype=np.ubyte) + primtype
#pix_rawx = np.zeros(pix_primid.size, dtype=np.ushort)
#pix_rawy = np.zeros(pix_primid.size, dtype=np.ushort)

## primid starts over every run, so do something complicated to fix that
## first find the boundaries
#diff_primid = pix_primid[1:] - pix_primid[0:-1]
## elements that are negative indicate beginning of a new set, although
## the index is now off by one, so add a zero at the beginning to index
## back into the original pix_primid array
#diff_primid = np.insert(diff_primid, 0, 0)
## now get indices of negative diffs
#indx = np.flatnonzero(diff_primid<0)
## loop through the boundaries (runs) and increment all subsequent
## primids by the number of primaries in the previous run
## also populate the run piddle
#this_runid = 0
#for boundary in indx :
#      pix_primid[boundary:] += run_numprims[this_runid]
#      pix_runid[boundary:] += 1
#      this_runid += 1

# initialize 2D pixel array for summed deposited energy and the secondary 
# particle types ('pdg') responsible; need enough for 0.13mm pixels
# in 78x76.15mm DEPFETS, which is really 600x586
img_edep = np.zeros([1201,1201], dtype=np.float32)
img_pdg = np.zeros([1201,1201], dtype=np.uint16)

# loop through primids
splitstep = 0
this_startrow = 0
for ii in range(numprims_interact) :
    this_primid = uniq_primid[ii]
    this_energy = evt_energy[(evt_primid == this_primid)].values[0]
    this_particletype = evt_particletype[(evt_primid == this_primid)].values[0]
    #print(f'ii = {ii}')
    #print(f'this_primid = {this_primid}')
    #print(f'this_energy = {this_energy}')
    #print(f'this_particletype = {this_particletype}')
    this_primtype = ptypes[this_particletype]
    indx = (primid == this_primid)
    # number of steps in this primary
    numsteps = indx.sum()

    img_edep = img_edep * 0
    img_pdg = img_pdg * 0

    #print(f'### {ii+1} of {numprims_interact}: Primary {this_primid} has {numsteps} steps.')

    for jj in range(numsteps) :

        #print(f'### {jj}: {particletype[indx].iloc[jj]}, {process[indx].iloc[jj]}, X = {prex[indx].iloc[jj]}->{postx[indx].iloc[jj]}, Y = {prey[indx].iloc[jj]}->{posty[indx].iloc[jj]}') 

        # get the secondary particle type
        this_pdg = ptypes[particletype[indx].iloc[jj]]
        this_prex = prex[indx].iloc[jj]
        this_prey = prey[indx].iloc[jj]
        this_postx = postx[indx].iloc[jj]
        this_posty = posty[indx].iloc[jj]
        this_edep = edep[indx].iloc[jj]

        # For any steps that stayed within a single pixel, dump all the deposited enegy in that pixel
        if (this_prex == this_postx) and (this_prey == this_posty) :
            #print(f'### Ionizing within one pixel.')
            img_edep[this_postx,this_posty] = img_edep[this_postx,this_posty] + this_edep
            img_pdg[this_postx,this_posty] = np.bitwise_or(img_pdg[this_postx,this_posty], this_pdg)

        # If the particle was a photon and it just didn't exit the volume, dump all the
        # energy deposited in the post-step pixel
        elif (particletype[indx].iloc[jj] == "gamma") and (process[indx].iloc[jj] != "CoupledTransportation") :
            #print(f'### Photon.')
            img_edep[this_postx,this_posty] = img_edep[this_postx,this_posty] + this_edep
            img_pdg[this_postx,this_posty] = np.bitwise_or(img_pdg[this_postx,this_posty], this_pdg)

        # For ionizing steps that traversed more then 1 pixel, bop along the
        # step in n 5µm substeps and dump 1/nth of the energy in the pixel that lives 
        # at each substep
        elif (this_prex != this_postx or this_prey != this_posty) :
            #print(f'### Ionizing crossing multiple pixels.')
            splitstep = splitstep + 1

            # starting and ending position of the full step, in mm
            pos1 = np.array([ prex_mm[indx].iloc[jj], prey_mm[indx].iloc[jj], prez_mm[indx].iloc[jj] ])
            pos2 = np.array([ postx_mm[indx].iloc[jj], posty_mm[indx].iloc[jj], postz_mm[indx].iloc[jj] ])

            # chop that up into 'substep_size'-µm substeps into which we'll dump charge
            nsubsteps = np.ceil(slen[indx].iloc[jj]/substep_size).astype(int)
            substeps = np.linspace(0,1,nsubsteps)
            pos = np.outer((1-substeps),pos1) + np.outer(substeps,pos2)
            # chop up the charge, too
            subedep = this_edep / nsubsteps

            # convert the substeps into pixels so we know where to dump the charge
            pos_pixel = (np.floor(pos[:,0:2]/0.130)).astype(int)

            # dump it
            for kk in range(nsubsteps) :
                img_edep[pos_pixel[kk,0],pos_pixel[kk,1]] = img_edep[pos_pixel[kk,0],pos_pixel[kk,1]] + subedep
                img_pdg[pos_pixel[kk,0],pos_pixel[kk,1]] = np.bitwise_or(img_pdg[pos_pixel[kk,0],pos_pixel[kk,1]], this_pdg)

    # convert images into pixel lists
    indx = (img_edep > 0)
    this_actx = np.argwhere(indx)[:,0]
    this_acty = np.argwhere(indx)[:,1]
    this_edep = img_edep[indx]
    this_pdg = img_pdg[indx]

    # put it all in the output arrays
    this_endrow = this_startrow + indx.sum()

    # allocate output arrays, which are pixel-based
    pix_primid[this_startrow:this_endrow] = this_primid
    # gotta figure these out
    #pix_detid[this_startrow:this_endrow] = this_detid
    #pix_rawx[this_startrow:this_endrow] = this_rawx
    #pix_rawy[this_startrow:this_endrow] = this_rawy
    pix_actx[this_startrow:this_endrow] = this_actx
    pix_acty[this_startrow:this_endrow] = this_acty
    pix_edep[this_startrow:this_endrow] = this_edep
    pix_pdg[this_startrow:this_endrow] = this_pdg
    pix_primtype[this_startrow:this_endrow] = this_primtype
    # these are already initialized with correct values
    #pix_runid[this_startrow:this_endrow] = this_runid

    this_startrow = this_endrow

# done primary-by-primary processing
# lop off unused parts of arrays
numrows = this_startrow
pix_primid = pix_primid[0:numrows]
pix_detid = pix_detid[0:numrows]
pix_rawx = pix_rawx[0:numrows]
pix_rawy = pix_rawy[0:numrows]
pix_actx = pix_actx[0:numrows]
pix_acty = pix_acty[0:numrows]
pix_edep = pix_edep[0:numrows]
pix_pdg = pix_pdg[0:numrows]
pix_primtype = pix_primtype[0:numrows]
pix_runid = pix_runid[0:numrows]

#    # figure out detector dimensions and coordinates
#    # don't do this because I can't reconcile 'sensarea' with 'localx/y'
#    # min and max of latter divide evenly into 1087 pixels
#    # those of former do not (1086.769 pixels)
#    #detinfo = file[b'BkgdIndetector_tuple_FitsBkgd_run'].arrays()
#    #sensarea_xmin = detinfo[b'sensarea_xmin']
#    #sensarea_xmax = detinfo[b'sensarea_xmax']
#    #sensarea_ymin = detinfo[b'sensarea_ymin']
#    #sensarea_ymax = detinfo[b'sensarea_ymax']
#    #npixx = detinfo[b'npixx']
#    #npixy = detinfo[b'npixy']
#
#    # convert ACT (LDA-specific) to RAW (quad-specific) coords
#    # I'm now using the same coords as the OU version, page 2 of
#    # 'OrientationsE00023277_00_edm.pptx', since Tanja is sending me
#    # LDA-specific locations now.  Gap info still applies:
#    #   In my simulation the gap between neighbouring sensitive areas is
#    #   8.15mm this is larger than the current status.  The current status
#    #   is 4.51mm. So in the new root files you should be able to test the
#    #   effect of different gap sizes.  I choose a sensitive area of
#    #   141.28mm.
#    # A sensitive area of area of 141.28mm is 1086.8 130-um pixels.  If we
#    # set the lower left corner of 1 in that diagram to # (ACTX,ACTY)=(0,0)
#    # then 0-511 will be A, 62.7 pixels (let's say 63 pixels) 512-574 will
#    # be gap, and 575-1086 will be B (if we're going right).  This assumes
#    # the larger gap, which is OK for now.
#
#    gappix = 63
#    imgsize = 512 + 512 + gappix
#    pix_rawx = pix_actx.copy()
#    pix_rawy = pix_acty.copy()
#
#    # set DETID based on ACT coords
#    indx = (pix_actx>=543) & (pix_acty<543)
#    pix_detid[indx] = 0
#    indx = (pix_actx<543) & (pix_acty<543)
#    pix_detid[indx] = 1
#    indx = (pix_actx<543) & (pix_acty>=543)
#    pix_detid[indx] = 2
#    indx = (pix_actx>=543) & (pix_acty>=543)
#    pix_detid[indx] = 3
#
#    # note there is no offset because I have exactly 1087x1087 in ACT
#    # so ACTX=0 is the same pixel as RAWX=0 (and same for Y)
#
#    # quad 0 RAWX = ACTX - 575
#    # quad 0 RAWY = 512 - ACTY
#    indx = (pix_detid==0)
#    pix_rawx[indx] -= 575
#    pix_rawy[indx] = 512 - pix_rawy[indx]
#
#    # quad 1 RAWX = ACTX
#    # quad 1 RAWY = 512 - ACTY
#    indx = (pix_detid==1)
#    #pix_rawx[indx] += 0
#    pix_rawy[indx] = 512 - pix_rawy[indx]
#
#    # quad 2 RAWX = ACTX
#    # quad 2 RAWY = ACTY-575
#    indx = (pix_detid==2)
#    #pix_rawx[indx] += 0
#    pix_rawy[indx] -= 575
#
#    # quad 3 RAWX = ACTX - 575
#    # quad 3 RAWY = ACTY - 575
#    indx = (pix_detid==3)
#    pix_rawx[indx] -= 575
#    pix_rawy[indx] -= 575
#
#    # now eliminate all activated pixels that are in the gap
#    # but just do everything with RAW outside of 0:511 anyway
#    indx = (pix_rawx >= 0) & (pix_rawx <= 511) & (pix_rawy >= 0) & (pix_rawy <= 511)
#    pix_primid = pix_primid[indx]
#    pix_actx = pix_actx[indx]
#    pix_acty = pix_acty[indx]
#    pix_rawx = pix_rawx[indx]
#    pix_rawy = pix_rawy[indx]
#    pix_detid = pix_detid[indx]
#    pix_edep = pix_edep[indx]
#    pix_pdg = pix_pdg[indx]
#    pix_primtype = pix_primtype[indx]
#    pix_runid = pix_runid[indx]
#
#    pix_actx = pix_actx.astype(np.ushort)
#    pix_acty = pix_acty.astype(np.ushort)
#    pix_rawx = pix_rawx.astype(np.ushort)
#    pix_rawy = pix_rawy.astype(np.ushort)

# make a table and save it to a FITS HDU
pix_tab = Table([pix_primid, pix_detid, pix_rawx, pix_rawy, 
    pix_actx, pix_acty, pix_edep, pix_pdg, pix_primtype, pix_runid], 
    names=('PRIMID', 'DETID', 'RAWX', 'RAWY', 'ACTX', 'ACTY', 'EDEP', 'PDG', 'PRIMTYPE', 'RUNID'))

hdu_pix = fits.table_to_hdu(pix_tab)

# put together the FITS HDUs into an HDUList
hdu_primary = fits.PrimaryHDU()
hdulist = fits.HDUList([hdu_primary, hdu_pix])
# add header keywords

# write FITS file
hdulist.writeto(outfile, overwrite=True)

exit()
