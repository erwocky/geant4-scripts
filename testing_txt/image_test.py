import numpy as np
from skimage import measure
from astropy.io import fits
import sys, glob, os
import datetime
import re

num_args = len(sys.argv) - 1
if num_args == 1:
    start_run = int(sys.argv[1])
    stop_run = start_run
else:
    sys.exit(f"Usage: {sys.argv[0]} [<run>]")

perl_files = glob.glob("output/run*")

cwd = os.getcwd()
input_dir = cwd[:cwd.rindex('/') + 1]

with open('framedist_vals', 'r') as f: 
    framedist = np.array(f.readline().split(), int) #these are the random values from perl
    #am using to make sure I get the same result across the two. 

date = datetime.datetime.now().astimezone()
date_string = date.strftime("%a %d %b %Y %I:%M:%S %p %Z").rstrip()
    
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

"""
Comments...
"""

epicpn_pattern_table=np.array([
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

"""#my $xdim; my $ydim;
#my $blobimg; my $mipimg;
#my $writeit; my $fitsfile;"""


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
def wfits(img, fitsfile, ow=True):
    hdu = fits.PrimaryHDU()
    hdu.data= img
    hdu.writeto(fitsfile, overwrite = ow)
def which(condition):
    if len(condition.shape) != 1:
        condition = condition.flat
    return np.where(condition)[0]

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
        #print (f"### Found something other than 4 datafiles for {this_run}, skipping.")
        continue

    # initialize event piddles, which will be written out or used later
    runid = np.zeros(0, dtype=int)
    detectorid = np.zeros(0, dtype=int)
    primid = np.zeros(0, dtype=int)
    actx = np.zeros(0, dtype=int)
    acty = np.zeros(0, dtype=int)
    phas = np.zeros((25,0), dtype=float)
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
    #$cumframedist = $cumframedist->where($cumframedist<=$num_protons_per_run);
    
    # running variables
    numevents = 0
    numtotblobs = 0

    
    
    # loop through the four quadrant data files for this run
    # now combine all four, since single primary can produce signal
    # in multiple quadrants
    for infile in infiles:
        
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
    
    # loop through frames and make a raw frame for each
    
    for i in range(numframes):
        # set the frame ID
        frame = frame_frame[i]

        # set the frame time
        frame_time[i] = frame * texp_frame

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
        this_blobid = blobid[pixel_indx]
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
        fitsfile = f"output/run{this_run}_frame{frame}_img.fits"
        if fitsfile in perl_files:
            fitsfile = f"output/prun{this_run}_frame{frame}_img.fits"
            wfits(img, fitsfile)

for fitsfile in perl_files:
    perl_image = fits.getdata(fitsfile, ext=0)
    pythonfile = "{}/p{}".format(*fitsfile.split('/'))
    python_image = fits.getdata(pythonfile, ext=0)
    file = match('output/(.+)_img', fitsfile).group(1)
    
    if python_image.shape != perl_image.shape:
        raise ValueError(f"Different shapes from {file}: {python_image.shape} vs {perl_image.shape}")
    i = 0    
    for j, k in zip(python_image.flat, perl_image.flat):
        if j != round(k, 6):
            raise ValueError(f"Different values in {file}: {j} vs {k}")
            
print("All good python")
