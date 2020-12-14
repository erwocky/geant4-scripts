#!/usr/bin/env python

# Pre-processing script that reads Geant4 output from a particular source
# and saves it to a uniform-format FITS file encoding the energy deposited
# in WFI pixels for each primary particle. This is for the MIT-generated
# Geant4 output.
#
# EDM Wed Nov  4 17:36:39 EST 2020
# First version to read Rick's format of Geant4 output. Adapted from MPE
# readit.py.

import numpy as np
import astropy
import os, sys, glob, re
from astropy.table import Table, vstack
from astropy.io import fits
import pandas as pd

num_args = len(sys.argv) - 1    
if num_args < 1:
    sys.exit(f"Usage: {sys.argv[0]} <path/to/Geant4 StepLog file> ...")

# Geant4 output data should be sorted into directories for different runs.
# The directory name is expected as an argument, and it must contain sets
# of files of the following naming convention:
#
# *_EvtLog_*.gdat
# *_StepLog_*.gdat
#
# The EvtLog file contains summary information of each generated primary,
# and the StepLog file contains the energy deposition steps.

# These things are Geant4-source-specific, so should be
# written to header of rawpix file. Probably.
sphere_radius = 20. # radius of boundary sphere in cm
# num_protons_per_run comes from the EvtLog file, since
# it could be different for each MIT run.

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

# loop through input filenames
for filename in sys.argv[1:] :

    if not re.search('.*_StepLog_[0-9]+\.gdat[\.gz]*$', filename) : 
        print(f'### Error: file {filename} does not look like a StepLog file, skipping.')
        continue
    if not os.path.isfile(filename) :
        print(f'### Error reading file {filename}, skipping.')
        continue

    path = os.path.dirname(filename)

    g4stepfile = filename
    g4evtfile = re.sub('StepLog', 'EvtLog', g4stepfile)
    runid = int( re.sub('.*?([0-9]+)_StepLog_([0-9]+)\.gdat(\.gz)*$', r'\1\2', g4stepfile) )
    outfile = os.path.join(path, f'rawpix_{runid}.fits')
    print(f'### Converting run {runid}.')
    print(f'### Step file   {g4stepfile}.')
    print(f'### Event file  {g4evtfile}.')
    print(f'### Output file {outfile}.')

    # open the Geant4 output files
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
    # Layout is (look-down or look-up doesn't matter, methinks):
    #   C B
    #   D A
    # Gotta use pd.loc correctly here in assignment to ensure it's
    # a view and not a copy into the Series. I guess.
    # This should be propagated back to the views I defined above, e.g.
    # g4step['LPre-X'] -> prex_mm
    # g4step['LPost-X'] -> postx_mm
    # g4step['LPre-Y'] -> prey_mm
    # g4step['LPost-Y'] -> posty_mm
    # 'D' is at origin, so coords don't need to be shifted
    # 'A' is shifted only in X
    indx = (depfet == 'DEPFETA')
    g4step.loc[indx, ('LPre-X')] += 7.64 + 31.36 + 35.3 + 3.7
    g4step.loc[indx, ('LPost-X')] += 7.64 + 31.36 + 35.3 + 3.7
    # 'C' is shifted only in Y
    indx = (depfet == 'DEPFETC')
    g4step.loc[indx, ('LPre-Y')] += 6.72 + 31.36 + 35.3 + 2.78
    g4step.loc[indx, ('LPost-Y')] += 6.72 + 31.36 + 35.3 + 2.78
    # 'B' is shifted in both X and Y
    indx = (depfet == 'DEPFETB')
    g4step.loc[indx, ('LPre-X')] += 7.64 + 31.36 + 35.3 + 3.7
    g4step.loc[indx, ('LPost-X')] += 7.64 + 31.36 + 35.3 + 3.7
    g4step.loc[indx, ('LPre-Y')] += 6.72 + 31.36 + 35.3 + 2.78
    g4step.loc[indx, ('LPost-Y')] += 6.72 + 31.36 + 35.3 + 2.78

    # bin the start and end step points into 130µm WFI pixels
    prex = np.floor(prex_mm/.13).astype(int)
    prey = np.floor(prey_mm/.13).astype(int)
    postx = np.floor(postx_mm/.13).astype(int)
    posty = np.floor(posty_mm/.13).astype(int)

    # get the differences to see if the step is fully in a pixel
    delta_x = postx - prex
    delta_y = posty - prey

    # initialize the 2D look-up tables of RAWX and RAWY,
    # defined as DEPFET quadrant pixels 0-511,0-511
    # DETID is 0,1,2,3 for A,B,C,D.
    img_detid = np.full([1201,1201], 8, dtype=np.byte)
    img_detid[:600,:586] = 3
    img_detid[600:,:586] = 0
    img_detid[:600,586:] = 2
    img_detid[600:,586:] = 1

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
        this_primtype = ptypes.get(this_particletype, 0)
        indx = (primid == this_primid)
        # number of steps in this primary
        numsteps = indx.sum()

        img_edep = img_edep * 0
        img_pdg = img_pdg * 0

        #print(f'### {ii+1} of {numprims_interact}: Primary {this_primid} has {numsteps} steps.')

        for jj in range(numsteps) :

            #print(f'### {jj}: {particletype[indx].iloc[jj]}, {process[indx].iloc[jj]}, X = {prex[indx].iloc[jj]}->{postx[indx].iloc[jj]}, Y = {prey[indx].iloc[jj]}->{posty[indx].iloc[jj]}') 

            # get the secondary particle type
            this_pdg = ptypes.get(particletype[indx].iloc[jj], 0)
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
        this_detid = img_detid[indx]

        # put it all in the output arrays
        this_endrow = this_startrow + indx.sum()

        # allocate output arrays, which are pixel-based
        pix_primid[this_startrow:this_endrow] = this_primid
        # gotta figure these out
        pix_detid[this_startrow:this_endrow] = this_detid
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

    # adjust primids so they start at zero and have no gaps
    evt_newprimid = np.arange(evt_primid.size)
    indx_primid = np.zeros(evt_primid.max()+1)
    indx_primid[evt_primid] = evt_newprimid
    pix_newprimid = indx_primid[pix_primid].astype(np.uint32)

    # make a table and save it to a FITS HDU
    pix_tab = Table([pix_newprimid, pix_primid, pix_detid, pix_rawx, pix_rawy, 
        pix_actx, pix_acty, pix_edep, pix_pdg, pix_primtype, pix_runid], 
        names=('PRIMID', 'OPRIMID', 'DETID', 'RAWX', 'RAWY', 'ACTX', 'ACTY', 'EDEP', 'PDG', 'PRIMTYPE', 'RUNID'))
    hdu_pix = fits.table_to_hdu(pix_tab)

    # add header keywords
    hdu_pix.header['SPH_RAD'] = sphere_radius
    hdu_pix.header.comments['SPH_RAD'] = 'Radius of Geant4 source sphere in cm.'
    hdu_pix.header['NPRI_GEN'] = numprims_gen
    hdu_pix.header.comments['NPRI_GEN'] = 'Number of primaries generated.'
    hdu_pix.header['NPRI_INT'] = numprims_interact
    hdu_pix.header.comments['NPRI_INT'] = 'Number of primaries producing signal.'

    # put together the FITS HDUs into an HDUList
    hdu_primary = fits.PrimaryHDU()
    hdulist = fits.HDUList([hdu_primary, hdu_pix])

    # write FITS file
    hdulist.writeto(outfile, overwrite=True)

exit()
