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
# EDM Sun Apr 18 08:55:20 EDT 2021
# Input StepLog file now processed in chunks to reduce memory usage.
# 'chunk_size' is the number of interacting primaries loaded into
# memory and processed at a time, so set this to something huge
# to return to the previous functionality. 10,000 is current default,
# and requires ~ 200 MB to run a single StepLog file, hopefully 
# regardless of size. runid also now a uint64, which is a bug fix.
#
# EDM Mon Apr 12 15:54:56 EDT 2021
# Added a bunch more memory logging, now toggle-able.
# Also changed default sphere size to 50 cm.
#
# EDM Wed Mar 31 11:15:02 EDT 2021
# Added primary particle energy PRIMENERGY to output FTIS pixel table.
#
# EDM Wed Mar 31 10:38:30 EDT 2021
# Added some elapsed time logging for benchmarking. Also added 'flush=True'
# to most of the print statements so they pop out in the supercloud log
# files as they go.
#
# EDM Mon Mar 15 17:31:52 EDT 2021
# Eliminate zero-energy-deposition rows from the step file panda before
# doing anything to them.
#
# EDM Thu Mar 11 08:45:41 EST 2021
# Attempted to reduce memory usage by specifying dtypes in pandas DataFrame
# and converting string columns to integers right away, since that's what
# they become. Also no longer renaming DF columns to new variables, since 
# it's never clear if they are views into the DF or copied as new Series.
# 
# EDM Tue Mar  9 13:11:35 EST 2021
# Added SECPARTYPE, which records the particle type of the particle that
# entered the detector and ultimately led to the energy deposition. This
# can be the primary or an external secondary, and like PARTYPE it is
# bit-wise added so multiple particle types can be recorded (but not
# multiples of the same particle). PARTYPE encodes final particles, and
# in may cases includes electrons ionized by the proton within the 
# detector and immediately dumping their energy.
#
# EDM Mon Feb  8 14:33:30 EST 2021
# Eliminated RAWX and RAWY since they're not used.
#
# EDM Fri Jan  8 13:50:08 EST 2021
# Changed output column EDEP to ENERGY. Probably will break stuff.
#
# EDM Tue Dec 22 09:52:18 EST 2020
# Updated output FITS data types. Implemented 'wfits()' function Brian made
# to simplify FITS writing. 'ptype'->'partype' to eliminate confusion with
# primary type (have to do this downstream too).
#
# EDM Wed Nov  4 17:36:39 EST 2020
# First version to read Rick's format of Geant4 output. Adapted from MPE
# readit.py.

import os, sys, glob, re, datetime
import numpy as np
from astropy.table import Table
from astropy.io import fits
import pandas as pd
import tracemalloc

trace_memory = False

if (trace_memory) :
    tracemalloc.start()

# don't truncate printed arrays
#np.set_printoptions(threshold=sys.maxsize)

num_args = len(sys.argv) - 1    
if num_args < 1:
    sys.exit(f"Usage: {sys.argv[0]} <path/to/Geant4 StepLog file> ...")

# functions
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

# to map a string column to an integer column
def strcol2intcol(strings, trans) :
    if strings in trans :
        return trans[strings]
    return 0

# These things are Geant4-source-specific, so should be
# written to header of rawpix file. Probably.
sphere_radius = 50. # radius of boundary sphere in cm
# num_protons_per_run comes from the EvtLog file, since
# it could be different for each MIT run.

# Number of interacting primaries to read from the StepLog 
# file at a time. Set it to something huge to read
# the whole file (could be bad for memory).
chunk_size = 10000

# size of a ionizing particle sub-step sample in µm
substep_size = 5.

# dictionary of particle types
# these are the same codes as the PDG in MPE data
# and are used for both primary type and particle
# (primary or secondary) that deposits energy;
# can be bit-wise added for multiple secondaries
# depositing in the same pixel for the same primary
ptypes = { 'unknown' : 0,
           'proton' : 1,
           'e-' : 2,
           'gamma' : 4,
           'e+' : 8,
           'pi-' : 16,
           'pi+' : 32,
           'neutron' : 64,
           'alpha' : 128 }

# processes; we only need the ones that affect
# energy deposition
proctypes = { 'unknown' : 0,
              'CoupledTransportation' : 1 }

# quadrants
quadrants = { 'DEPFETA' : 0,
              'DEPFETB' : 1,
              'DEPFETC' : 2,
              'DEPFETD' : 3 }

date_start = datetime.datetime.now().astimezone()
date_string = date_start.strftime("%a %d %b %Y %I:%M:%S %p %Z").rstrip()
print("############################################################", flush=True)
print(f"### Started {sys.argv[0]} on {date_string}", flush=True)
then = date_start

# loop through input filenames
for filename in sys.argv[1:] :

    # make sure filename exists and looks right
    if not os.path.isfile(filename) :
        print(f'### Error reading file {filename}, skipping.', flush=True)
        continue
    if not re.search('.*_StepLog_[0-9]+\.gdat[\.gz]*$', filename) : 
        print(f'### Error: file {filename} does not look like a StepLog file, skipping.', flush=True)
        continue


    path = os.path.dirname(filename)
    g4stepfile = filename
    g4evtfile = re.sub('StepLog', 'EvtLog', g4stepfile)
    runid = int( re.sub('.*?([0-9]+)_StepLog_([0-9]+)\.gdat(\.gz)*$', r'\1\2', g4stepfile) )
    outfile = os.path.join(path, f'rawpix_{runid}.fits')
    print(f'### Converting run {runid}.', flush=True)
    print(f'### Step file   {g4stepfile}.', flush=True)
    print(f'### Event file  {g4evtfile}.', flush=True)
    print(f'### Output file {outfile}.', flush=True)

    if (trace_memory) :
        print(f'[ Memory check ] After init.')
        current, peak = tracemalloc.get_traced_memory()
        print(f"Current memory usage is {current / 10**6}MB.")
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)

    # read the Geant4 EvtLog file
    g4evt = pd.read_csv(g4evtfile, sep='\s+', 
            usecols=['Event', 'Particle', 'Energy'], dtype={'Event':np.uint32, 'Particle':str, 'Energy':np.float32})
    # convert string column(s) to integers to save memory
    g4evt['Particle'] = np.array(list(strcol2intcol(x, ptypes) for x in g4evt['Particle']), dtype=np.uint8)

    # total number of primaries generated
    numprims_gen = g4evt['Event'].size

    if (trace_memory) :
        print(f'[ Memory check ] After EvtLog read.')
        current, peak = tracemalloc.get_traced_memory()
        print(f"Current memory usage is {current / 10**6}MB.")
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')
        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)

    # read the Geant4 EvtLog file 'Event' column to get unique
    # interacting primaries
    g4step_init = pd.read_csv(g4stepfile, sep='\s+', usecols=['Event'], dtype={'Event':np.uint32})

    # number of rows the step file
    numrows = g4step_init.shape[0]

    # get unique primids and number of primaries that interact
    uniq_primid = np.unique(g4step_init['Event'])
    numprims_interact = uniq_primid.size

    print(f'### {numprims_interact} of {numprims_gen} generated primaries interacted.', flush=True)

    # initialize output arrays, which are pixel-based
    # these will be appended to after each chunk
    pix_primid = np.zeros(0, dtype=np.uint32)
    pix_detid = np.zeros(0, dtype=np.uint8)
    pix_actx = np.zeros(0, dtype=np.uint16)
    pix_acty = np.zeros(0, dtype=np.uint16)
    pix_edep = np.zeros(0, dtype=np.float32)
    pix_partype = np.zeros(0, dtype=np.uint16)
    pix_secpartype = np.zeros(0, dtype=np.uint16)
    pix_primtype = np.zeros(0, dtype=np.uint8)
    pix_primenergy = np.zeros(0, dtype=np.float32)
    pix_runid = np.zeros(0, dtype=np.uint64)

    # chunk it
    # this all assumes the event IDs monotonically increase
    chunk_start = 0
    this_skiprows = 0
    counter_primid = 0
    while (chunk_start < numprims_interact) :
        this_chunk = uniq_primid[chunk_start:chunk_start+chunk_size]
        this_nrows = np.isin(g4step_init['Event'], this_chunk).sum()

        print(f'### Processing pimaries chunk {chunk_start}-{chunk_start+chunk_size-1}, {this_nrows} rows.', flush=True)

        # read a chunk of the Geant4 EvtLog file
        # always skip the very first header row and define columns by hand
        g4step = pd.read_csv(g4stepfile, sep='\s+', skiprows=1+this_skiprows, nrows=this_nrows,
                usecols=[0, 1, 2, 3, 5, 6, 7, 14, 15, 16, 17, 18, 19, 20], 
                names=['Event', 'Volume', 'Parent', 'ID', 'SLen', 'Edep', 'Particle', 
                    'LPre-X', 'LPre-Y', 'LPre-Z', 'LPost-X', 'LPost-Y', 'LPost-Z', 'Process'],
                dtype={'Event':np.uint32, 'Volume':str, 'Parent':np.uint32, 'ID':np.uint32, 
                    'SLen':np.float32, 'Edep':np.float32, 'Particle':str, 'LPre-X':np.float32, 
                    'LPre-Y':np.float32, 'LPre-Z':np.float32, 'LPost-X':np.float32, 
                    'LPost-Y':np.float32, 'LPost-Z':np.float32, 'Process':str})
        #g4step = pd.read_csv(g4stepfile, sep='\s+', skiprows=this_skiprows, nrows=this_nrows,
        #        usecols=['Event', 'Volume', 'Parent', 'ID', 'SLen', 'Edep', 'Particle', 
        #            'LPre-X', 'LPre-Y', 'LPre-Z', 'LPost-X', 'LPost-Y', 'LPost-Z', 'Process'], 
        #        dtype={'Event':np.uint32, 'Volume':str, 'Parent':np.uint32, 'ID':np.uint32, 
        #            'SLen':np.float32, 'Edep':np.float32, 'Particle':str, 'LPre-X':np.float32, 
        #            'LPre-Y':np.float32, 'LPre-Z':np.float32, 'LPost-X':np.float32, 
        #            'LPost-Y':np.float32, 'LPost-Z':np.float32, 'Process':str})

        if (trace_memory) :
            print(f'[ Memory check ] After StepLog read.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # convert string column(s) to integers to save memory
        g4step['Particle'] = np.array(list(strcol2intcol(x, ptypes) for x in g4step['Particle']), dtype=np.uint8)
        g4step['Volume'] = np.array(list(strcol2intcol(x, quadrants) for x in g4step['Volume']), dtype=np.uint8)
        g4step['Process'] = np.array(list(strcol2intcol(x, proctypes) for x in g4step['Process']), dtype=np.uint8)

        if (trace_memory) :
            print(f'[ Memory check ] After str->unt8 conversion for four columns.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # delete step file rows that don't deposit energy
        g4step = g4step[(g4step['Edep'] > 0)]

        # number of non-zero deposits in this chunk of the step file
        numsteps = g4step.shape[0]

        if (trace_memory) :
            print(f'[ Memory check ] After removing 0 energy steps.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        print(f'### {numsteps} steps with non-zero energy depositions.', flush=True)

        # allocate output arrays, which are pixel-based
        this_pix_primid = np.zeros(numsteps, dtype=np.uint32)
        this_pix_detid = np.zeros(numsteps, dtype=np.uint8) + 255   # just to make sure this gets set
        this_pix_actx = np.zeros(numsteps, dtype=np.uint16)
        this_pix_acty = np.zeros(numsteps, dtype=np.uint16)
        this_pix_edep = np.zeros(numsteps, dtype=np.float32)
        this_pix_partype = np.zeros(numsteps, dtype=np.uint16)
        this_pix_secpartype = np.zeros(numsteps, dtype=np.uint16)
        this_pix_primtype = np.zeros(numsteps, dtype=np.uint8)
        this_pix_primenergy = np.zeros(numsteps, dtype=np.float32)
        this_pix_runid = np.zeros(numsteps, dtype=np.uint64) + np.uint64(runid)

        if (trace_memory) :
            print(f'[ Memory check ] After output array allocation.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # Convert Rick's LOCAL coords (in mm per depfet) to the full focal plane.
        # DEPFETS are each 78mm wide x 76.15mm high with 0.59mm and 0.37mm gaps.
        # So the full DEPFET area spans X=0.0:156.59, Y=0.0:152.67 in mm, or
        # X=0:1204.5, Y=0:1174.4 in 130-µm pixels. So the array needs to be at 
        # least that big to hold everything. Thus ACT coords are defined as:
        #   ACTX=0:1205, ACTY=0:1205 (1206x1206 pixels)
        # The visible area of each in ACT coords:
        #       LLC_X  URC_X   LLC_Y  URC_Y
        #   D   7.67   74.23   6.76   73.32
        #   C   7.67   74.23   79.3   145.86
        #   A   82.16  148.72  6.76   73.32
        #   B   82.16  148.72  79.3   145.86
        # These have been somewhat jiggled so they fall on (virual) pixel boundaries.
        # And the gaps are rounded to the nearest pixel. These correspond to the foloowing
        # ACT pixels:
        #       ACTX        ACTY
        #   D   59:570      52:563
        #   C   59:570      610:1121
        #   A   632:1143    52:563
        #   B   632:1143    610:1121
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
        indx = (g4step['Volume'] == quadrants['DEPFETA'])
        g4step.loc[indx, ('LPre-X')] += 7.64 + 31.36 + 35.3 + 3.7 + .59
        g4step.loc[indx, ('LPost-X')] += 7.64 + 31.36 + 35.3 + 3.7 + .59
        # 'C' is shifted only in Y
        indx = (g4step['Volume'] == quadrants['DEPFETC'])
        g4step.loc[indx, ('LPre-Y')] += 6.72 + 31.36 + 35.3 + 2.78 + .37
        g4step.loc[indx, ('LPost-Y')] += 6.72 + 31.36 + 35.3 + 2.78 + .37
        # 'B' is shifted in both X and Y
        indx = (g4step['Volume'] == quadrants['DEPFETB'])
        g4step.loc[indx, ('LPre-X')] += 7.64 + 31.36 + 35.3 + 3.7 + .59
        g4step.loc[indx, ('LPost-X')] += 7.64 + 31.36 + 35.3 + 3.7 + .59
        g4step.loc[indx, ('LPre-Y')] += 6.72 + 31.36 + 35.3 + 2.78 + .37
        g4step.loc[indx, ('LPost-Y')] += 6.72 + 31.36 + 35.3 + 2.78 + .37

        if (trace_memory) :
            print(f'[ Memory check ] After coord conversion.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # bin the start and end step points into 130µm WFI pixels
        prex = np.floor(g4step['LPre-X']/.13).astype(np.int16)
        prey = np.floor(g4step['LPre-Y']/.13).astype(np.int16)
        postx = np.floor(g4step['LPost-X']/.13).astype(np.int16)
        posty = np.floor(g4step['LPost-Y']/.13).astype(np.int16)

        if (trace_memory) :
            print(f'[ Memory check ] After coord binning.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # get the differences to see if the step is fully in a pixel
        delta_x = postx - prex
        delta_y = posty - prey

        if (trace_memory) :
            print(f'[ Memory check ] After delta coord calc.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # initialize the 2D look-up tables of DETID # indexed by ACTX and ACTY
        # DETID is 0,1,2,3 for A,B,C,D.
        img_detid = np.full([1205,1205], 8, dtype=np.byte)
        img_detid[:600,:586] = quadrants['DEPFETD']
        img_detid[600:,:586] = quadrants['DEPFETA']
        img_detid[:600,586:] = quadrants['DEPFETC']
        img_detid[600:,586:] = quadrants['DEPFETB']

        # initialize 2D pixel array for summed deposited energy and the secondary 
        # particle types ('partype') responsible; need enough for 0.13mm pixels
        # in 78x76.15mm DEPFETS, which is really 600x586
        # 'partype' encodes all particles that deposited energy
        # 'sectype' encodes all particles produced outside the DEPFET that eventually deposited energy
        img_edep = np.zeros(img_detid.shape, dtype=np.float32)
        img_partype = np.zeros(img_detid.shape, dtype=np.uint8)
        img_secpartype = np.zeros(img_detid.shape, dtype=np.uint8)

        if (trace_memory) :
            print(f'[ Memory check ] After 2D image allocation.')
            current, peak = tracemalloc.get_traced_memory()
            print(f"Current memory usage is {current / 10**6}MB.")
            snapshot = tracemalloc.take_snapshot()
            top_stats = snapshot.statistics('lineno')
            print("[ Top 10 ]")
            for stat in top_stats[:10]:
                print(stat)

        # loop through primids
        splitstep = 0
        this_startrow = 0
        #for ii in range(chunk_size) :
            #this_primid = this_chunk[ii]

        for this_primid in this_chunk :
            this_primtype = g4evt['Particle'][(g4evt['Event'] == this_primid)].values[0]
            this_primenergy = g4evt['Energy'][(g4evt['Event'] == this_primid)].values[0]
            #print(f'this_primid = {this_primid}')
            #print(f'this_energy = {this_energy}')
            #print(f'this_particletype = {this_particletype}')
            #this_primtype = ptypes.get(this_particletype, 0)
            indx = (g4step['Event'] == this_primid)
            # number of steps in this primary
            numsteps = indx.sum()

            # get 'these' steps for this_primary
            these_parent = g4step['Parent'][indx]
            these_pid = g4step['ID'][indx]
            these_partype = g4step['Particle'][indx]
            these_secpartype = these_partype.copy()
            these_process = g4step['Process'][indx]
            these_prex = prex[indx]
            these_prey = prey[indx]
            these_postx = postx[indx]
            these_posty = posty[indx]
            these_prex_mm = g4step['LPre-X'][indx]
            these_prey_mm = g4step['LPre-Y'][indx]
            these_prez_mm = g4step['LPre-Z'][indx]
            these_postx_mm = g4step['LPost-X'][indx]
            these_posty_mm = g4step['LPost-Y'][indx]
            these_postz_mm = g4step['LPost-Z'][indx]
            these_edep = g4step['Edep'][indx]

            # find the unique secondaries; these are PIDs with parents that don't show up in pid
            indx_notsecs = np.isin(these_parent, these_pid)
            indx_secs = np.invert(indx_notsecs)
            num_notsecs = indx_notsecs.sum()
            num_secs = indx_secs.sum()
            #print(f'##### PRIMID {this_primid} ########################################')
            #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            #    print(f'{pd.concat([these_parent, these_pid, these_secparticletype], axis=1)}')

            # an odd thing, but make the secondaries their own parents
            these_parent[indx_secs] = these_pid[indx_secs]
            my_secid = np.unique(these_pid[indx_secs])

            #print(f'### PRIMID {this_primid}')
            #print(f'{num_secs}, {num_notsecs}')
            #if (this_primid == 21591) :
            #    print(f'{these_parent}')

            # now loop through the set of steps from non-secondaries,
            # trace back one step and reassign parent, then iterate
            # actually don't need the iteration because it works in sequence
            num_tofix = num_notsecs
            #print(f'Num. to fix: {num_tofix}')
            #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            #    print(f'{pd.concat([these_parent, these_pid, these_secparticletype], axis=1)}')
            for my_id in np.unique(these_pid[indx_notsecs]) :
                my_old_parent = these_parent[(these_pid==my_id)].iloc[0] 
                my_new_parent = these_parent[(these_pid==my_old_parent)].iloc[0]
                these_parent[(these_pid==my_id)] = my_new_parent
                these_secpartype[(these_pid==my_id)] = these_partype[(these_pid==my_new_parent)].iloc[0]
            indx_notsecs = np.isin(these_parent, my_secid, invert=True)
            num_tofix = indx_notsecs.sum()
            #print(f'Num. to fix: {num_tofix}')
            #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
            #    print(f'{pd.concat([these_parent, these_pid, these_secparticletype], axis=1)}')


            img_edep = img_edep * 0
            img_partype = img_partype * 0
            img_secpartype = img_secpartype * 0

            for jj in range(numsteps) :

                #print(f'### {jj}: {particletype[indx].iloc[jj]}, {process[indx].iloc[jj]}, X = {prex[indx].iloc[jj]}->{postx[indx].iloc[jj]}, Y = {prey[indx].iloc[jj]}->{posty[indx].iloc[jj]}') 

                # get info for 'this' step
                this_partype = these_partype.iloc[jj]
                this_secpartype = these_secpartype.iloc[jj]
                this_process = these_process.iloc[jj]
                this_prex = these_prex.iloc[jj]
                this_prey = these_prey.iloc[jj]
                this_postx = these_postx.iloc[jj]
                this_posty = these_posty.iloc[jj]
                this_prex_mm = these_prex_mm.iloc[jj]
                this_prey_mm = these_prey_mm.iloc[jj]
                this_prez_mm = these_prez_mm.iloc[jj]
                this_postx_mm = these_postx_mm.iloc[jj]
                this_posty_mm = these_posty_mm.iloc[jj]
                this_postz_mm = these_postz_mm.iloc[jj]
                this_edep = these_edep.iloc[jj]

                # For any steps that stayed within a single pixel, dump all the deposited enegy in that pixel
                if (this_prex == this_postx) and (this_prey == this_posty) :
                    #print(f'### Ionizing within one pixel.')
                    img_edep[this_postx,this_posty] = img_edep[this_postx,this_posty] + this_edep
                    img_partype[this_postx,this_posty] = np.bitwise_or(img_partype[this_postx,this_posty], this_partype)
                    img_secpartype[this_postx,this_posty] = np.bitwise_or(img_secpartype[this_postx,this_posty], this_secpartype)

                # If the particle was a photon and it just didn't exit the volume, dump all the
                # energy deposited in the post-step pixel
                elif (this_partype == ptypes['gamma']) and (this_process != proctypes['CoupledTransportation']) :
                    #print(f'### Photon.')
                    img_edep[this_postx,this_posty] = img_edep[this_postx,this_posty] + this_edep
                    img_partype[this_postx,this_posty] = np.bitwise_or(img_partype[this_postx,this_posty], this_partype)
                    img_secpartype[this_postx,this_posty] = np.bitwise_or(img_secpartype[this_postx,this_posty], this_secpartype)

                # For ionizing steps that traversed more then 1 pixel, bop along the
                # step in n 5µm substeps and dump 1/nth of the energy in the pixel that lives 
                # at each substep
                elif (this_prex != this_postx or this_prey != this_posty) :
                    #print(f'### Ionizing crossing multiple pixels.')
                    splitstep = splitstep + 1

                    # starting and ending position of the full step, in mm
                    pos1 = np.array([ this_prex_mm, this_prey_mm, this_prez_mm ])
                    pos2 = np.array([ this_postx_mm, this_posty_mm, this_postz_mm ])

                    # chop that up into 'substep_size'-µm substeps into which we'll dump charge
                    nsubsteps = np.ceil(g4step['SLen'][indx].iloc[jj]/substep_size).astype(int)
                    substeps = np.linspace(0,1,nsubsteps)
                    pos = np.outer((1-substeps),pos1) + np.outer(substeps,pos2)
                    # chop up the charge, too
                    subedep = this_edep / nsubsteps

                    # convert the substeps into pixels so we know where to dump the charge
                    pos_pixel = (np.floor(pos[:,0:2]/0.130)).astype(int)

                    # dump it
                    for kk in range(nsubsteps) :
                        img_edep[pos_pixel[kk,0],pos_pixel[kk,1]] = img_edep[pos_pixel[kk,0],pos_pixel[kk,1]] + subedep
                        img_partype[pos_pixel[kk,0],pos_pixel[kk,1]] = np.bitwise_or(img_partype[pos_pixel[kk,0],pos_pixel[kk,1]], this_partype)
                        img_secpartype[pos_pixel[kk,0],pos_pixel[kk,1]] = np.bitwise_or(img_secpartype[pos_pixel[kk,0],pos_pixel[kk,1]], this_secpartype)

            # convert images into pixel lists
            indx = (img_edep > 0)
            this_actx = np.argwhere(indx)[:,0]
            this_acty = np.argwhere(indx)[:,1]
            this_edep = img_edep[indx]
            this_partype = img_partype[indx]
            this_secpartype = img_secpartype[indx]
            this_detid = img_detid[indx]

            # put it all in the output arrays
            this_endrow = this_startrow + indx.sum()

            # fill in pixel-based output arrays for this primary
            this_pix_primid[this_startrow:this_endrow] = this_primid
            this_pix_detid[this_startrow:this_endrow] = this_detid
            this_pix_actx[this_startrow:this_endrow] = this_actx
            this_pix_acty[this_startrow:this_endrow] = this_acty
            this_pix_edep[this_startrow:this_endrow] = this_edep
            this_pix_partype[this_startrow:this_endrow] = this_partype
            this_pix_secpartype[this_startrow:this_endrow] = this_secpartype
            this_pix_primtype[this_startrow:this_endrow] = this_primtype
            this_pix_primenergy[this_startrow:this_endrow] = this_primenergy
            # these are already initialized with correct values
            #this_pix_runid[this_startrow:this_endrow] = this_runid

            this_startrow = this_endrow

            # keep us updated
            if ((counter_primid+1) % 10 == 0) :
                now = datetime.datetime.now().astimezone()
                print(f"### Run {runid}, chunk {chunk_start}: done {counter_primid+1} of {numprims_interact} primaries, {(now-then).total_seconds()} / {(now-date_start).total_seconds()} sec elapsed.", flush=True)
                then = now
                if (trace_memory) :
                    current, peak = tracemalloc.get_traced_memory()
                    print(f"Current memory usage is {current / 10**6}MB.")
                    snapshot = tracemalloc.take_snapshot()
                    top_stats = snapshot.statistics('lineno')
                    print("[ Top 10 ]")
                    for stat in top_stats[:10]:
                        print(stat)

                pass

            counter_primid += 1

        # done primary-by-primary processing
        # lop off unused parts of arrays and append to existing versions
        numsteps = this_startrow
        pix_primid = np.append(pix_primid, this_pix_primid[0:numsteps])
        pix_detid = np.append(pix_detid, this_pix_detid[0:numsteps])
        pix_actx = np.append(pix_actx, this_pix_actx[0:numsteps])
        pix_acty = np.append(pix_acty, this_pix_acty[0:numsteps])
        pix_edep = np.append(pix_edep, this_pix_edep[0:numsteps])
        pix_partype = np.append(pix_partype, this_pix_partype[0:numsteps])
        pix_secpartype = np.append(pix_secpartype, this_pix_secpartype[0:numsteps])
        pix_primtype = np.append(pix_primtype, this_pix_primtype[0:numsteps])
        pix_primenergy = np.append(pix_primenergy, this_pix_primenergy[0:numsteps])
        pix_runid = np.append(pix_runid, this_pix_runid[0:numsteps])

        chunk_start += chunk_size
        this_skiprows += this_nrows



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
    pix_edep = pix_edep[indx]
    pix_partype = pix_partype[indx]
    pix_secpartype = pix_secpartype[indx]
    pix_primtype = pix_primtype[indx]
    pix_primenergy = pix_primenergy[indx]
    pix_runid = pix_runid[indx]
    # remove the gaps
    indx = ( (pix_actx<=570) | (pix_actx>=632) ) & ( (pix_acty<=563) | (pix_acty>=610) )
    pix_primid = pix_primid[indx]
    pix_detid = pix_detid[indx]
    pix_actx = pix_actx[indx]
    pix_acty = pix_acty[indx]
    pix_edep = pix_edep[indx]
    pix_partype = pix_partype[indx]
    pix_secpartype = pix_secpartype[indx]
    pix_primtype = pix_primtype[indx]
    pix_primenergy = pix_primenergy[indx]
    pix_runid = pix_runid[indx]

    # adjust primids so they start at zero and have no gaps
    evt_newprimid = np.arange(g4evt['Event'].size)
    indx_primid = np.zeros(g4evt['Event'].max()+1)
    indx_primid[g4evt['Event']] = evt_newprimid
    pix_newprimid = indx_primid[pix_primid].astype(np.uint32)
    numprims_interact = np.unique(pix_newprimid).size

    # fix the num. of interacting primaries to only those in the sensitive area

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
    wfits( (['PRIMID', 'OPRIMID', 'DETID', 'ACTX', 'ACTY',
        'ENERGY', 'PARTYPE', 'SECPARTYPE', 'PRIMTYPE', 'PRIMENERGY', 'RUNID'],
        [pix_newprimid.astype(np.uint32),
        pix_primid.astype(np.uint32),
        pix_detid.astype(np.uint8),
        pix_actx.astype(np.uint16),
        pix_acty.astype(np.uint16),
        pix_edep.astype(np.single),
        pix_partype.astype(np.uint8),
        pix_secpartype.astype(np.uint8),
        pix_primtype.astype(np.uint8),
        pix_primenergy.astype(np.single),
        pix_runid.astype(np.uint64)]), 
        outfile, hdr=hdr, hdrcomments=hdrcomments)


if (trace_memory) :
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()

date_stop = datetime.datetime.now().astimezone()
elapsed = (date_stop-date_start).total_seconds()
date_string = date_stop.strftime("%a %d %b %Y %I:%M:%S %p %Z").rstrip()
print(f"### Finished {sys.argv[0]} on {date_string}, took {elapsed} sec elapsed.", flush=True)
print("############################################################", flush=True)

exit()
