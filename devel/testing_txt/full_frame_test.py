import numpy as np
from skimage import measure
from astropy.io import fits
from astropy.table import Table
import sys, glob, os
import datetime
import re


num_args = len(sys.argv) - 1    
if num_args == 1:
    start_run, stop_run = [int(sys.argv[1])]*2
    run_str = f"### Doing run {start_run}"
elif num_args == 2:
    start_run, stop_run = sorted(sys.argv[1:3])
    run_str = f"### Doing runs {start_run} to {stop_run}"
else:
    sys.exit(f"Usage: {sys.argv[0]} <start_run> [<stop_run>]")


cwd = os.getcwd()
input_dir = cwd[:cwd.rindex('/') + 1]
    
with open('framedist_vals', 'r') as f: 
    framedist = np.array(f.readline().split(), int) #these are the random values from perl
    #am using to make sure I get the same result across the two. 
    
date = datetime.datetime.now().astimezone()
date_string = date.strftime("%a %d %b %Y %I:%M:%S %p %Z").rstrip()
#print("############################################################")
#print(f"### Started {sys.argv[0]} on {date_string}")
#print(run_str)
    
#//line 99

# define defaults
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
sphere_radius = 70. # radius of boundary sphere in cm
num_protons_per_run = 1.e6 # number of proton primaries 
                                # in a simulatin run (from Jonathan)
                                # his email said 1e7, but looking at the 
                                # input files it's really 1e6!!!
detector_area = 4. * (130.e-4 * 512.)**2 # detector area in cm2, 
                                              # assuming 130 um pixels
texp_run = num_protons_per_run/3.14159/proton_flux/(sphere_radius**2)
                         # total exposure time in sec for this run
texp_frame = .005   # frame exposure time in sec
mu_per_frame = num_protons_per_run * texp_frame / texp_run
     # minimum ionizing pa

    
#print(f"### There are {num_protons_per_run} primaries in this run.")
#print(f"### The run exposure time is {texp_run} sec.")
#print(f"### The frame exposure time is {texp_frame} sec")
#print(f"### for an expected mean of {mu_per_frame} primaries per frame.")

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

# hash of codes indexed by particle type indexed
ptypes = {
     'proton': 0,
     'gamma': 1,
     'electron': 2,
     'neutron': 3,
     'pi+': 4,
     'e+': 5,
     'pi-': 6,
     'nu_mu': 7,
     'anti_nu_mu': 8,
     'nu_e': 9,
     'kaon+': 10,
     'mu+': 11,
     'deuteron': 12,
     'kaon0L': 13,
     'lambda': 14,
     'kaon-': 15,
     'mu-': 16,
     'kaon0S': 17,
     'alpha': 18,
     'anti_proton': 19,
     'triton': 20,
     'anti_neutron': 21,
     'sigma-': 22,
     'sigma+': 23,
     'He3': 24,
     'anti_lambda': 25,
     'anti_nu_e': 26,
     'anti_sigma-': 27,
     'xi0': 28,
     'anti_sigma+': 29,
     'xi-': 30,
     'anti_xi0': 31,
     'C12': 32,
     'anti_xi-': 33,
     'Li6': 34,
     'Al27': 35,
     'O16': 36,
     'Ne19': 37,
     'Mg24': 38,
     'Li7': 39,
     'He6': 40,
     'Be8': 41,
     'Be10': 42,
     'unknown': 99
}


# initialize rng
rng = np.random.RandomState(1234)


# temporary variables
x, y = 0, 0


"""
Comments...
line 246
"""

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
    if isinstance(data, tuple):
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

#######################################
# Main loop.
# Step through Geant4 output data files.
# For each one, create a frame per primary,
# find blobs and MIPs, and then find events in that frame.
# Change start_run to parallelize things (i.e. do chunks of 10 runs
# in parallel).
# delte variables later

for this_run in range(start_run, stop_run + 1):
    # see if there are files for this run
            
    infiles = glob.glob(f"{input_dir}input/{this_run}_detector?")
    if len(infiles) != 4:
        print (f"### Found something other than 4 datafiles for {this_run}, skipping.")
        continue

    # initialize event piddles, which will be written out or used later
    runid = np.zeros(0, dtype=int)
    detectorid = np.zeros(0, dtype=int)
    primid = np.zeros(0, dtype=int)
    actx = np.zeros(0, dtype=int)
    acty = np.zeros(0, dtype=int)
    phas = np.zeros((0,25), dtype=float)
    pha = np.zeros(0, dtype=float)
    ptype = np.zeros(0, dtype=int)
    energy = np.zeros(0, dtype = float)     # in keV
    evttype = np.zeros(0, dtype=int)
    blobdist = np.zeros(0, dtype=float)
    mipdist = np.zeros(0, dtype=float)
    pattern = np.zeros(0, dtype=int)
    vfaint = np.zeros(0, dtype=int)
    # assign frames and times to each primary
    # to start, assume mean of one primary per second
    evt_time = np.zeros(0, dtype=float)
    evt_frame = np.zeros(0, dtype=int)
    pix_time = np.zeros(0, dtype=float)
    pix_frame = np.zeros(0, dtype=int)
    
    # initialize structures to hold the secondary particle columns
    # piddles for the numeric columns; these are enough for now
    run = np.zeros(0, dtype=int)        # Geant4 run (* in *_detector[0123])
    detector = np.zeros(0, dtype=int)   # WFI detector (? in *_detector?)
    eid = np.zeros(0, dtype=int)        # primary ID
    particleid = np.zeros(0, dtype=int) # interacting particle ID
    parentid = np.zeros(0, dtype=int)   # don't really need probably

    # initialize piddles to hold the energy deposition (per pixel) columns
    # some of this will be written out the pixel list
    xdep = np.zeros(0, dtype=int)
    ydep = np.zeros(0, dtype=int)
    endep = np.zeros(0, dtype=float)
    rundep = np.zeros(0, dtype=int)
    detectordep = np.zeros(0, dtype=int)
    eiddep = np.zeros(0, dtype=int)
    framedep = np.zeros(0, dtype=int)
    piddep = np.zeros(0, dtype=int)
    ptypedep = np.zeros(0, dtype=int)
    cprocdep = np.zeros(0, dtype=int)
    blobid = np.zeros(0, dtype=int)

    # initialize piddles to hold frame-specific things to go in FITS table
    frame_frame = np.zeros(0, dtype=int)
    frame_time = np.zeros(0, dtype=float)
    frame_runid = np.zeros(0, dtype=int)
    frame_npix = np.zeros(0, dtype=int)
    frame_npixmip = np.zeros(0, dtype=int)
    frame_nevt = np.zeros(0, dtype=int)
    frame_nevtgood = np.zeros(0, dtype=int)
    frame_nevt27 = np.zeros(0, dtype=int)
    frame_nblob = np.zeros(0, dtype=int)
    frame_nprim = np.zeros(0, dtype=int)

    # initialize piddles to hold blob-specific things to go in FITS table
    blob_frame = np.zeros(0, dtype=int)
    blob_blobid = np.zeros(0, dtype=int)
    blob_cenx = np.zeros(0, dtype=float)
    blob_ceny = np.zeros(0, dtype=float)
    blob_cenxcl = np.zeros(0, dtype=float)
    blob_cenycl = np.zeros(0, dtype=float)
    blob_npix = np.zeros(0, dtype=int)
    blob_energy = np.zeros(0, dtype=float)
    blob_energycl = np.zeros(0, dtype=float)
    
    # initialize things for the running frames which we will
    # randomly populate
    # frame settings
    # we know there are $num_protons_per_run, so generate enough
    # random frames to hold them
    ##framedist = rng.poisson(mu_per_frame, int(2*num_protons_per_run/mu_per_frame))
    
    cumframedist = framedist.cumsum()
    # get the total number of frames needed to capture all the primaries; 
    # will write this to FITS header so we can combine runs
    numtotframes = which(cumframedist >= num_protons_per_run)[0] + 1
    # this is wrong, because it will remove the last bin which we need
    # it's also unnecessary
    #cumframedist = cumframedist[cumframedist <= num_protons_per_run]
    
    # running variables
    numevents = 0
    numtotblobs = 0

    
    
    # loop through the four quadrant data files for this run
    # now combine all four, since single primary can produce signal
    # in multiple quadrants
    for infile in infiles:

        #print(f"### Reading {infile}")
        
        rc = match('[0-9]+_detector([0-9]+)', infile) #extracts the detector name
        this_detector = int(rc.group(1))
        
        ptype = {}
        cproc = {}
        
        with open(infile, 'r') as IN:
            # step through the input file and accumulate primaries
            for line in IN: #switched to a for loop because of built-in __iter__ method
                if match('^\s*#', line): #skip comments #added ability to have arbritrary whitespace before '#'
                    continue
                if match('^\s*$', line): #skip blank lines:
                    continue
                if not match(',', line): #could be if ',' not in line
                    continue

                fields = line.rstrip().split(',')

                if match('[a-zA-Z]', fields[0]): # if the first column is a string, then this is a particle line
                    # retain the primary for this interaction
                    this_eid = int(float(fields[1]))
                    eid = np.append(eid, this_eid)
                    # particle type and interaction type are hashes so
                    # that the pixel-specific read can pick them up
                    # doesn't matter if the particle ID is re-used from
                    # primary to primary, since this will reset it
                    ptype[int(fields[2])] = fields[0]
                    cproc[int(fields[2])] = fields[4]
                else: # if the first column is a number, then this is a pixel hit line
                    #print(fields)
                    if float(fields[2]) <= splitth: # skip it if less than split threshold is deposited,
                        continue             #since that is ~ the lower threshold of pixels we'll get

                    tmp_x, tmp_y = int(fields[0]), int(fields[1])
                    if tmp_x<xydep_min or tmp_y<xydep_min or tmp_x>xydep_max or tmp_y>xydep_max:
                        continue # skip it if it's outside the 512x512 region of a quad

                    xdep = np.append(xdep, tmp_x)
                    ydep = np.append(ydep, tmp_y)
                    endep = np.append(endep, float(fields[2]))
                    rundep = np.append(rundep, this_run)
                    detectordep = np.append(detectordep, this_detector)
                    eiddep = np.append(eiddep, this_eid)
                    framedep = np.append(framedep, 0)
                    piddep = np.append(piddep, int(fields[3]))
                    # %ptype is hash of particle type strings indexed by the id
                    # %ptypes is (constant) hash of my own particle type IDs indexed
                    # by the string (confused yet?)
                    ptypedep = np.append(ptypedep, ptypes[ptype.get(int(fields[3]), 'unknown')])
                    blobid = np.append(blobid, 0)#np.zeros(0, dtype = int));
                    #print f"### piddep = {fields[3]}"
                    #print f"### ptype[piddep] = {ptype[fields[3]]}"
                    #print "### ptypes{ptype{piddep}} = {ptypes[ptype[fields[3]]]}"         
    # done loop through quadrant data files for this run
    uniq_eid = np.unique(eid)
    numprimaries = uniq_eid.size
    primary_flux = numprimaries / texp_run / detector_area
    numpixels = endep.size
    #print(f"### Run {this_run}: found {numprimaries} primaries that interacted.")
    #print(f"### Run {this_run}: that's {primary_flux} protons per sec per cm2.")
    #print(f"### Run {this_run}: found {numpixels} pixels with deposited energy")
    
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
        #print(f"Doing primary {primary})"
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
    #print(f"### Run {this_run}: number of primaries in 1 2 3 4 quads: {num_in_different_quadrants}")
    # min and max X,Y values

    # figure out the unique frames that are populated by 
    # primaries that interacted
    frame_frame = np.append(frame_frame, np.unique(framedep))
    numframes = frame_frame.size
    # now we can append to frame piddles since we know how many there are
    frame_runid = np.append(frame_runid, np.zeros(numframes) + this_run)
    frame_time = np.append(frame_time, np.zeros(numframes, dtype=float))
    frame_npix = np.append(frame_npix, np.zeros(numframes))
    frame_npixmip = np.append(frame_npixmip, np.zeros(numframes))
    frame_nevt = np.append(frame_nevt, np.zeros(numframes))
    frame_nevtgood = np.append(frame_nevtgood, np.zeros(numframes))
    frame_nevt27 = np.append(frame_nevt27, np.zeros(numframes))
    frame_nblob = np.append(frame_nblob, np.zeros(numframes))
    frame_nprim = np.append(frame_nprim, np.zeros(numframes))

    pct_interact = 100. * numprimaries / num_protons_per_run
    #print(f"### Run {this_run}: generated {numtotframes} total frames,")
    #print(f"### of which {numframes} frames with the {numprimaries}")
    #print("### interacting primaries will be written.")
    #print(f"### {num_protons_per_run} total primaries were simulated.")
    #print(f"### {numprimaries} or {pct_interact}% of these produced a WFI interaction.")
    
    # loop through frames and make a raw frame for each

    ptype = np.zeros(0, dtype=int)
    for i in range(numframes):
        # set the frame ID
        frame = frame_frame[i]

        # set the frame time
        frame_time[i] = frame * texp_frame

        # keep us updated
        if i % 100 == 0:
            #print(f"### Run {this_run}: done {i} of {numframes} frames")
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
        pixel_indx = which(framedep==frame)
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
            fitsfile = f"output/run{this_run}_frame{frame}_img.fits"
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
        # blobadjust decrements blob IDs when we chuck one, so the IDs are continuous from 1
        blobadjust = 0
        for j in range(1, np.amax(blobimg) + 1):
            indx = which(blobimg == j)
            #indx2d = whichND(blobimg == j) dsn
            # set the blob to zeros and skip it if there are too few elements
            # if it contains a MIP, it's a good blob irregardless
            if (indx.size < npixthresh and np.amax(img.flat[indx]) < mipthresh):
                blobimg.flat[indx] = 0
                blobadjust += 1
                continue
            blobimg.flat[indx] -= blobadjust

            # this is the running blobid which we need to add to blob piddles
            blob_blobid = np.append(blob_blobid, numtotblobs +j- blobadjust)
            blob_frame = np.append(blob_frame, frame)
            blob_npix = np.append(blob_npix, indx.size)
            # calculate unclipped blob centroid and summed energy
            tmp_en = img.flat[indx]
            tmp_toten = np.sum(tmp_en)
            tmp_wtd_x = img_xvals.flat[indx] * tmp_en
            tmp_wtd_y = img_yvals.flat[indx] * tmp_en
            blob_cenx = np.append(blob_cenx, xoff + tmp_wtd_x.sum() / tmp_toten)
            blob_ceny = np.append(blob_ceny, yoff + tmp_wtd_y.sum() / tmp_toten)
            blob_energy = np.append(blob_energy, tmp_toten)
            # calculate clipped blob centroid and summed energy
            tmp_en = np.clip(tmp_en, None, clip_energy)
            tmp_toten = tmp_en.sum()
            tmp_wtd_x = img_xvals.flat[indx] * tmp_en
            tmp_wtd_y = img_yvals.flat[indx] * tmp_en
            blob_cenxcl = np.append(blob_cenxcl, xoff + tmp_wtd_x.sum() / tmp_toten)
            blob_cenycl = np.append(blob_cenycl, yoff + tmp_wtd_y.sum() / tmp_toten)
            blob_energycl = np.append(blob_energycl, tmp_toten)
            
        # record number of blobs in this frame
        frame_nblob[i] = np.amax(blobimg)

        # if we found some blobs, change their IDs so they reflect the running total
        # for this run, increase that running blob number 
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
                    evt_phas = np.zeros((num_localmax, 25), float)
                    
                    phas_order = iter(range(1, 25))
                    evt_phas[:,0] = img[evty, evtx]  
                    for yi in range(-1,2):
                        for xi in range(-1,2):
                            if (xi, yi) != (0,0):
                                evt_phas[:,next(phas_order)] = img[evty + yi,evtx + xi]              
                    
                    # get outer 5x5 PHAS piddle
                    # first pad the image with zeros all around
                    tmp_img = np.pad(img, 1)
                    # now evtx and evty are too low by one, so correct for that in index
                    for yi in range(-1,4):
                        for xi in range(-1,4):
                            if not(-1 < xi < 3 and -1 < yi < 3):
                                evt_phas[:, next(phas_order)] = tmp_img[evty + yi,evtx + xi]
                    
                    # set the particle type for all these events
                    # based on whatever ptype produced the local max
                    evt_ptype = np.copy(ptypeimg[evty,evtx])
                    evt_primid = np.copy(primidimg[evty,evtx])
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

        runid = np.append(runid, np.zeros(numevents_thisframe, float) + this_run)
        detectorid = np.append(detectorid, np.zeros(numevents_thisframe, float) + this_run)
        primid = np.append(primid, evt_primid)
        evt_frame = np.append(evt_frame, np.zeros(numevents_thisframe, float) + frame)
        
        #actx = np.append(actx, np.zeros(numevents_thisframe, float))
        #acty = acty->append(zeros(long, $numevents_thisframe))
        #phas = phas->glue(1, zeros(long, 25, $numevents_thisframe))
        actx = np.append(actx, evtx)
        acty = np.append(acty, evty)
        phas = np.concatenate([phas, evt_phas])
        pha = np.append(pha, np.zeros(numevents_thisframe))
        ptype = np.append(ptype, evt_ptype)
        evttype = np.append(evttype, np.zeros(numevents_thisframe))
        energy = np.append(energy, np.zeros(numevents_thisframe))
        blobdist = np.append(blobdist, np.zeros(numevents_thisframe))
        mipdist = np.append(mipdist, np.zeros(numevents_thisframe))
        pattern = np.append(pattern, np.zeros(numevents_thisframe))
        vfaint = np.append(vfaint, np.zeros(numevents_thisframe))
        
        # step through all events to determine
        # EPIC-pn pattern, summed PHA, VFAINT flag, minimum distance to a
        # blob and mip.
        for j in range(numevents, numevents+numevents_thisframe):
            # below is already done in event finding and pasted above
            # get X,Y of center pixel
            x = actx[j]
            y = acty[j]
            # below is deprecated, since we've assumed this in event finding
            ## eliminate events on edges so we can use 5x5
            #next if ($x<2 or $x>1021 or $y<2 or $y>1021);

            # get ACIS flt grade; "which" returns the indices of
            # non-central pixels greater than or equal to event threshold,
            # these are used as bits to raise 2 to the power, and summed
            # (this can probably be removed from the loop if I'm smart)
            indx = which(phas[j,1:9] >= splitth)
            fltgrade = pow(2,indx).sum()
            # convert to EPIC-pn pattern (from look-up table)
            pattern[j] = epicpn_pattern_table[fltgrade]

            # sum 3x3 pixels over split threshold to get PHA
            # this is an ACIS PHA, not Suzaku
            # (this can probably be removed from the loop if I'm smart)
            pha[j] = phas[j][phas[j] >= splitth].sum()
            # apply gain correction for this node
            # get the gain parameters from the node (0-3)
            #print(f"{phas[j]} {gain_intercept} {gain_slope} {pha[j]}")
            pha[j] = (pha[j] - gain_intercept) * gain_slope / evperchan

            # convert pha to energy
            # (this can probably be removed from the loop if I'm smart)
            energy[j] = pha[j] * evperchan / 1000
            #print(f"{(phas[j])} {gain_intercept} {gain_slope} {pha[j]} {energy[j]}")

            # perform simple VFAINT filtering; also update the EPIC-pn pattern
            # of doubles, triples, and quads based on it
            # get outer 16 pixels and if any are above split flag it

            noutermost = which(phas[j, 9:25] > splitth).size #maybe don't need which
            if noutermost > 0:
                vfaint[j] = 1
                # EDM Fri Jan  4 11:24:20 EST 2019 
                # Reinstated the 5x5 filtering on PATTERN.
                # EDM Thu May 31 13:46:28 EDT 2018 
                # change to remove 5x5 criterion on EPIC-pn patterns
                # for doubles, triples, quads
                if pattern[j] > 0: 
                    pattern[j] = 13

            # get minimum distance to a blob
            # first find delta X,Y from this event to list of blob pixels
            delta_blobs = indx_blobs - np.array([x,y])
            #print(delta_blobs)
            # square that, sum it, square root it to get the distance
            if delta_blobs.shape[0] > 0:
                blobdist[j] = np.amin(np.sqrt(np.square(delta_blobs).sum(axis = 1)))
            # unless there aren't any blobs, in which case set blobdist to -1
            else:
                blobdist[j] = -1
            #print(f"{blobdist[j]}")

            # get minimum distance to a mip
            # first find delta X,Y from this event to list of mip pixels
            delta_mips = indx_mips - np.array([x,y])
            #print $delta_mips;
            # square that, sum it, square root it to get the distance
            if delta_mips.shape[0] > 0:
                
                mipdist[j] = np.amin(np.sqrt(np.square(delta_mips).sum(axis = 1)))
            # unless there aren't any mips, in which case set mipdist to -1
            else:
                mipdist[j] = -1

            #print(f"{mipdist[j]}")

            # we really want ACTX,Y in real WFI coords in the event list, 
            # so fix that here; $x and $y are children of $actx and $acty
            # for this event, which is why this works
            actx[j] += xoff
            acty[j] += yoff

            # add info to frame piddles
            frame_nevt[i] += 1
            if energy[j] < mipthresh:
                frame_nevtgood[i] += 1
            if energy[j] >= 2 and energy[j] < 7:
                frame_nevt27[i] += 1
                

        # increment the number of events and hit the next frame
        numevents += numevents_thisframe
    # done loop through frames

    # segregate the reds and greens as defined now by EPIC-pn pattern
    # reds are bad patterns only
    indx_goodpatterns, indx_badpatterns = which_both(pattern < 13)
    # cyans are singles and doubles, which are "best" of the good
    indx_goodbest, indx_badbest = which_both(pattern < 5);
    # combine indices for the filters we want via intersection
    indx_reds = indx_badpatterns
    indx_greens = np.intersect1d(indx_goodpatterns, indx_badbest)
    indx_cyans = np.intersect1d(indx_goodpatterns, indx_goodbest)
    # determine (event) count rate in cts/cm2/s/keV in the important band
    indx_goodenergy = which((energy > 2) & (energy <= 7))
    flux_goodenergy = indx_goodenergy.size / texp_run / detector_area / 5.

    evttype[indx_reds] = 3
    evttype[indx_greens] = 4
    evttype[indx_cyans] = 6
    
    
    hdr = {'NFRAMES': numtotframes, 'NBLOBS': numtotblobs}

    outevtfile = f"output/geant_events_{this_run}_evt.fits"
    wfits(
        (['ACTX','ACTY','BLOBDIST','DETID','ENERGY','EVTTYPE','FRAME','MIPDIST','PATTERN','PHA','PHAS','PRIMID','PTYPE','RUNID','VFAINT'],
           [actx,  acty, blobdist, detectorid, energy, evttype, evt_frame, mipdist, pattern, pha, phas, primid, ptype, runid, vfaint]),
           outevtfile, hdr=hdr)
    
    # write out pixel list
    #print(f"### Writing pixel list for run {this_run}.")
    outpixfile = f"output/geant_events_{this_run}_pix.fits"
    xdep += x_offset
    ydep += y_offset
    indx_sorted = np.argsort(framedep)
    wfits( 
    (['ACTX', 'ACTY', 'BLOBID', 'ENERGY', 'FRAME', 'PRIMID', 'PTYPE', 'RUNID'],
     [xdep[indx_sorted], ydep[indx_sorted], blobid[indx_sorted], endep[indx_sorted], framedep[indx_sorted], eiddep[indx_sorted], ptypedep[indx_sorted], rundep[indx_sorted]]), 
    outpixfile, hdr=hdr)

    # write out frame list
    #print(f"### Writing frame list for run {this_run}.")
    outframefile = f"output/pgeant_events_{this_run}_frames.fits"
    indx_sorted = np.argsort(frame_frame)
    wfits(
    (['FRAME', 'NBLOB', 'NEVT', 'NEVT27', 'NEVTGOOD', 'NPIX', 'NPIXMIP', 'NPRIM', 'RUNID', 'TIME'], 
     [frame_frame[indx_sorted], frame_nblob[indx_sorted], frame_nevt[indx_sorted], frame_nevt27[indx_sorted], frame_nevtgood[indx_sorted], frame_npix[indx_sorted], frame_npixmip[indx_sorted], frame_nprim[indx_sorted], frame_runid[indx_sorted], frame_time[indx_sorted]]),
    outframefile, hdr = hdr)  

    # write out blob list
    #print(f"### Writing blob list for run {this_run}.")
    outblobfile = f"output/pgeant_events_{this_run}_blobs.fits"
    indx_sorted = np.argsort(blob_blobid)
    wfits(
    (['BLOBID', 'CENX', 'CENXCL', 'CENY', 'CENYCL', 'ENERGY', 'ENERGYCL', 'FRAME', 'NPIX'],
    [blob_blobid[indx_sorted], blob_cenx[indx_sorted], blob_cenxcl[indx_sorted], blob_ceny[indx_sorted], blob_cenycl[indx_sorted], blob_energy[indx_sorted], blob_energycl[indx_sorted], blob_frame[indx_sorted], blob_npix[indx_sorted]]),
    outblobfile, hdr = hdr)

    #print(f"### Finished run {this_run}.")
    

# done loop through runs


rf_name = 'frame_full'
sense = 6
sense1 = 17

def split(line, typ=int):
    return [typ(i) for i in line.split(',')[:-1]]

def parse_la(line, typ = float):
    splited = [[typ(x) for x in y.split(',')[:-1]] for y in line.split('|')[:-1]]
    shape, inds, vals = splited
    arr = np.zeros([int(x) for x in shape[::-1]])
    arr.flat[inds] = vals
    return arr


to_test = {}
tests = ['pattern', 'indx_goodpatterns', 'indx_badpatterns', 'indx_goodbest', 'indx_badbest', 'indx_reds', 
         'indx_greens', 'indx_cyans', 'energy', 'indx_goodenergy', 'evttype', 'flux_goodenergy', 'numtotframes', 'numtotblobs']

fits_tests = [f"output/geant_events_{this_run}_evt.fits", f"output/geant_events_{this_run}_pix.fits"]


with open(f"{rf_name}.txt", 'r') as FH:
    state = 0
    for line in FH:
        if state <= -1:
            arr = parse_la(line)
            exec(f"perl_{tests[state]} = arr")
        elif state <= 10:
            exec(f"perl_{tests[state]} = np.array(line.split(',')[:-1], float)")
        else:
            exec(f"perl_{tests[state]} = float(line.rstrip())")
        state += 1

    for test in tests:
        to_test[test] = [eval('perl_'+ test), eval(test)]

        
for test, items in to_test.items():
    if isinstance(items[0], np.ndarray):
        if items[0].shape != items[1].shape:
            raise ValueError(f"Different shapes of {test}: {items[0].shape} vs {items[1].shape}")

        if not (np.around(items[0], sense) == np.around(items[1], sense)).all():
            i = np.where(np.around(items[0], sense) != np.around(items[1], sense))
            if len(i) > 1:
                j = [(i[0][x], i[0][x]) for x in range(len(i[0]))]
            else:
                j = list(i[0])
            #print(test, i)
            raise ValueError(f"Different values of {test} found at index {j}: {items[0][i]} vs {items[1][i]}")
    else:
        if round(items[0], sense1) != round(items[1], sense1):
            raise ValueError(f"Different values of {test}: {items[0]} vs {items[1]}")
                        

fits_files = [x[7:] for x in glob.glob("output/geant_events*")]

def evt_blobs_test(file):
    python_evt = Table.read('output/' + file)
    perl_evt = Table.read('output/p' + file)
    pylist = [[[round(z,1) for z in y] if isinstance(y, np.ndarray) else round(y,1) for y in x] for x in list(perl_evt)]
    pllist = [[[round(z,1) for z in y] if isinstance(y, np.ndarray) else round(y,1) for y in x] for x in list(python_evt)]
    for i,j in zip(pylist, pllist):
        if i != j:
            raise ValueError(f"{i} \n {j}")
              
def frames_test(file):
    python_frames = Table.read('output/' + file)
    perl_frames = Table.read('output/p' + file)
    for i, j in zip(list(perl_frames), list(python_frames)):
        if not list(i) == list(j):
            raise ValueError(i, j)
              
def pix_tests(file):
    python_pix = Table.read('output/' + file)
    perl_pix = Table.read('output/p' + file)
    plist = [[round(y, 2) for y in x] for x in list(perl_pix)]
    for i in python_pix:
        row = [round(x, 2) for x in i]
        if not row in plist:
            raise ValueError(row)
        else:
            plist.remove(row)

fits_tests = {'evt': evt_blobs_test, 'blobs': evt_blobs_test, 'frames': frames_test, 'pix': pix_tests}

for file in fits_files:
    for test in fits_tests:
        if test in file:
            fits_tests[test](file)
                
            
print("All good python")
