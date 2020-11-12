if ($#ARGV < 0) {
    die "Usage: $0 [<level (1-6)>] <run>\n";
} elsif ($#ARGV == 0) {
    $level = $ARGV[0];
    $start_run = 1;
    $stop_run = $start_run;
} elsif ($#ARGV >= 1) {
    $level = $ARGV[0];
    $start_run = $ARGV[1];
    $stop_run = $start_run;
    
}
    
use PDL;
use PDL::NiceSlice;
use PDL::GSL::INTEG;
use PDL::Image2D;
use PDL::ImageND;
use PDL::Ufunc;
use PDL::Fit::Gaussian;
use PDL::GSL::RNG;
use PDL::NiceSlice;
use warnings;
#use EDM;
#use Astro::FITS::CFITSIO qw( :constants );
#Astro::FITS::CFITSIO::PerlyUnpacking(0);

#$PDL::debug = 0;

sub print_all {
    for (my $i = 0; $i < $_[0] -> nelem; $i++) {
        my $g = $_[0] -> flat-> index($i);
        print FH "$g,";
    }
    print FH "\n";
    
};

sub print_nz {
    my $shape = $_[0] -> shape;
    my $ind = which($_[0] != 0);
    my $vals = $_[0]->where($_[0] != 0);
    
    foreach $piddle (($shape, $ind, $vals)) {
        for (my $i = 0; $i < $piddle -> nelem; $i++) {
            my $g = $piddle -> flat-> index($i);
            print FH "$g,";
        }
        print FH "|";
    }
    print FH "\n";
};

$pwd = `pwd`;
$input_dir = substr($pwd, 0, -12);
$wf_name = 'blobfind';

my $date = `date`;
chomp $date;
#print "############################################################\n";
#print "### Started $0 on $date\n";
print "### $0: Level $level: Run $start_run\n";

# define defaults
my $evtth = .1;          # 100 eV for WFI Geant4 simulations
my $splitth = $evtth;    # same as evtth for WFI Geant4 simulations
my $npixthresh = 5;      # minimum number of pixels in a blob
my $mipthresh = 15.;     # minimum ionizing particle threshold in keV
my $clip_energy = 22.;   # maximum pixel value reported by WFI, in keV
my $skip_writes = -1;    # writes FITS images for every skip_writes primary; set 
                         # to -1 to turn off writes
my $evperchan = 1000.;   # why not? PHA here is really PI
my $mipthresh_chan = $mipthresh * 1000. / $evperchan; 
                         # minimum ionizing particle threshold in PHA units
my $spec_maxkev = 100.;
my $numchans = int(($spec_maxkev*1000.) / $evperchan);
my $gain_intercept = 0.; # use this for Geant4 data
my $gain_slope = 1000.;  # use this for Geant4 data (pixel PH units are keV)
#my @gain_intercepts = (0., 0., 0., 0.);         # in ADU
#my @gain_slopes = ( 1000., 1000., 1000., 1000. ); # in eV/ADU
# rate and frame defaults
my $proton_flux = 4.1 * 7./5.;   # protons/s/cm2; 7/5 accounts for alphas, etc.
my $sphere_radius = 70.; # radius of boundary sphere in cm
my $num_protons_per_run = 1.e6; # number of proton primaries 
                                # in a simulatin run (from Jonathan)
                                # his email said 1e7, but looking at the 
                                # input files it's really 1e6!!!
my $detector_area = 4. * (130.e-4 * 512.)**2; # detector area in cm2, 
                                              # assuming 130 um pixels
my $texp_run = $num_protons_per_run/3.14159/$proton_flux/($sphere_radius**2);
                         # total exposure time in sec for this run
my $texp_frame = .005;   # frame exposure time in sec
my $mu_per_frame = $num_protons_per_run * $texp_frame / $texp_run; 
                         # mean number of primaries per frame

#print "### There are $num_protons_per_run primaries in this run.\n";
#print "### The run exposure time is $texp_run sec.\n";
#print "### The frame exposure time is $texp_frame sec,\n";
#print "### for an expected mean of $mu_per_frame primaries per frame.\n";

# conversion from ACIS 256 flt grades to EPIC-pn pattern
# see 'pnpatterns.png' for the patterns
# they include 0 (singles), 1-4 (doubles), 5-8 (triples), 9132 (triples)
#
# EDM Fri Jan  4 11:25:55 EST 2019 
# Re-change to below: once again I am enforcing the 5x5 filtering,
# thus using the real EPIC-pn PATTERN, as Esra is.
# EDM Thu May 31 13:43:32 EDT 2018 
# Change to below: no longer enforcing the 5x5 outer pixel criterion
# on doubles, triples, quads.  This is to more closely resemble what
# Tanja is doing (i.e. this is a "pre-filter" that naturally excludes
# events withing 2 pixels of MIP, which we want to do specifically later.
# Note that the VFAINT flag is still set, so this can be used to 
# remove those events if wanted.  But they won't be set to 13.
# old bits:
# but this isn't all, the doubles, triples, and quads
# also can't have an outer 5x5 pixel above split; so VFAINT is run
# on doubles, triples, quads, but not on singles
# I set everything else to 13 although I don't think that's
# the official grade
my $epicpn_pattern_table=pdl( 
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
);

# hash of codes indexed by particle type indexed
my %ptypes = (
     'proton' => 0, 'gamma' => 1, 'electron' => 2, 
     'neutron' => 3, 'pi+' => 4, 'e+' => 5, 'pi-' => 6,
     'nu_mu' => 7, 'anti_nu_mu' => 8, 'nu_e' => 9,
     'kaon+' => 10, 'mu+' => 11, 'deuteron' => 12,
     'kaon0L' => 13, 'lambda' => 14, 'kaon-' => 15,
     'mu-' => 16, 'kaon0S' => 17, 'alpha' => 18,
     'anti_proton' => 19, 'triton' => 20,
     'anti_neutron' => 21, 'sigma-' => 22, 'sigma+' => 23,
     'He3' => 24, 'anti_lambda' => 25, 'anti_nu_e' => 26,
     'anti_sigma-' => 27, 'xi0' => 28, 'anti_sigma+' => 29,
     'xi-' => 30, 'anti_xi0' => 31, 'C12' => 32,
     'anti_xi-' => 33, 'Li6' => 34, 'Al27' => 35,
     'O16' => 36, 'Ne19' => 37, 'Mg24' => 38, 
     'Li7' => 39, 'He6' => 40, 'Be8' => 41, 'Be10' => 42
);

# initialize rng
my $rng = PDL::GSL::RNG->new('taus');
$rng->set_seed(1234);

# temporary variables
my $x = 0.; my $y = 0.;
my $xdim; my $ydim;
my $blobimg; my $mipimg;
my $writeit; my $fitsfile;


# Offsets to add to Jonathan's coords, one per detector, to get
# them into 512 pixel quadrants.  He has a three-pixel gap between
# detectors in X and Y, with the center at 0,0, so the detectors
# run in (his) X and Y from -513->-2 and 2->513.  He also has an
# extra 77 pixels around the edges that I need to chop off, but
# (according to a 2018-10-31 email exchange with him) I will not
# eliminate secondaries that interact in those regions if they also
# produce charge in the good regions.
# So ACTX,ACTY needs to be 1024+3 = 1027x1027, 0-1026
# add 513 to Jonathan's coords, and lop off anything with ACTX or
# ACTY < 0 or > 1026.
# I think we also want an RAWX,RAWY and quadrant number

# 0-511; indexed by detector number 0-3
my $x_offset = 513;
my $y_offset = 513;
my $imgsize = 1027;
my $actmin = 0;
my $actmax = 1026;
my $xydep_min = -513;
my $xydep_max = 513;

# initialize the blobdist running histogram
# this histogram is the number of pixels available in each radial bin
# linear bins, this should be the same as plot_things.pl!:
my $numbins_blobhist = 300;
my $minrad_blobhist = 0.;
my $maxrad_blobhist = 600.;
my $binsize_blobhist = ($maxrad_blobhist-$minrad_blobhist)/$numbins_blobhist;
# needs some value to bin, so just make it a zero
my ($blobhist_bins,$blobhist_vals) = hist(zeros(long, 1),
     $minrad_blobhist,$maxrad_blobhist,$binsize_blobhist);
# really only wanted the bins, so reset the histogram
$blobhist_vals .= 0;

# initialize the mipdist running histogram
# this histogram is the number of pixels available in each radial bin
# linear bins, this should be the same as plot_things.pl!:
my $numbins_miphist = 300;
my $minrad_miphist = 0.;
my $maxrad_miphist = 600.;
my $binsize_miphist = ($maxrad_miphist-$minrad_miphist)/$numbins_miphist;
# needs some value to bin, so just make it a zero
my ($miphist_bins,$miphist_vals) = hist(zeros(long, 1),
     $minrad_miphist,$maxrad_miphist,$binsize_miphist);
# really only wanted the bins, so reset the histogram
$miphist_vals .= 0;

#######################################
# Main loop.
# Step through Geant4 output data files.
# For each one, create a frame per primary,
# find blobs and MIPs, and then find events in that frame.
# Change $start_run to parallelize things (i.e. do chunks of 10 runs
# in parallel).
for (my $this_run=$start_run;$this_run<=$stop_run;$this_run++) {
    # see if there are files for this run
    my @infiles = `/bin/ls ${input_dir}input/${this_run}_detector[0123]`;
    if ($#infiles != 3) {
        #print "### Found something other than 3 datafiles for $this_run, skipping.\n";
        next;
    }

    # initialize event piddles, which will be written out or used later
    my $runid = zeros(long, 0);
    my $detectorid = zeros(long, 0);
    my $primid = zeros(long, 0);
    my $actx = zeros(long, 0);
    my $acty = zeros(long, 0);
    my $phas = zeros(double, 25, 0);
    my $pha = zeros(double, 0);
    my $ptype = zeros(long, 0);
    my $energy = zeros(double, 0);     # in keV
    my $evttype = zeros(long, 0);
    my $blobdist = zeros(double, 0);
    my $mipdist = zeros(double, 0);
    my $pattern = zeros(ushort, 0);
    my $vfaint = zeros(ushort, 0);
    # assign frames and times to each primary
    # to start, assume mean of one primary per second
    my $evt_time = zeros(double, 0);
    my $evt_frame = zeros(long, 0);
    my $pix_time = zeros(double, 0);
    my $pix_frame = zeros(long, 0);

    # initialize structures to hold the secondary particle columns
    # piddles for the numeric columns; these are enough for now
    my $run = zeros(long, 0);        # Geant4 run (* in *_detector[0123])
    my $detector = zeros(long, 0);   # WFI detector (? in *_detector?)
    my $eid = zeros(long, 0);        # primary ID
    my $particleid = zeros(long, 0); # interacting particle ID
    my $parentid = zeros(long, 0);   # don't really need probably

    # initialize piddles to hold the energy deposition (per pixel) columns
    # some of this will be written out the pixel list
    my $xdep = zeros(long, 0);
    my $ydep = zeros(long, 0);
    my $endep = zeros(double, 0);
    my $rundep = zeros(long, 0);
    my $detectordep = zeros(long, 0);
    my $eiddep = zeros(long, 0);
    my $framedep = zeros(long, 0);
    my $piddep = zeros(long, 0);
    my $ptypedep = zeros(long, 0);
    my $cprocdep = zeros(long, 0);
    my $blobid = zeros(long, 0);

    # initialize piddles to hold frame-specific things to go in FITS table
    my $frame_frame = zeros(long, 0);
    my $frame_time = zeros(double, 0);
    my $frame_runid = zeros(long, 0);
    $frame_npix = zeros(long, 0);
    $frame_npixmip = zeros(long, 0);
    my $frame_nevt = zeros(long, 0);
    my $frame_nevtgood = zeros(long, 0);
    my $frame_nevt27 = zeros(long, 0);
    my $frame_nblob = zeros(long, 0);
    $frame_nprim = zeros(long, 0);

    # initialize piddles to hold blob-specific things to go in FITS table
    my $blob_frame = zeros(long, 0);
    my $blob_blobid = zeros(long, 0);
    my $blob_cenx = zeros(double, 0);
    my $blob_ceny = zeros(double, 0);
    my $blob_cenxcl = zeros(double, 0);
    my $blob_cenycl = zeros(double, 0);
    my $blob_npix = zeros(long, 0);
    my $blob_energy = zeros(double, 0);
    my $blob_energycl = zeros(double, 0);

    # initialize things for the running frames which we will
    # randomly populate
    # frame settings
    # we know there are $num_protons_per_run, so generate enough
    # random frames to hold them
    my $framedist = $rng->ran_poisson($mu_per_frame, int(2*$num_protons_per_run/$mu_per_frame))->long;
    my $cumframedist = $framedist->dcumusumover;
    # get the total number of frames needed to capture all the primaries; 
    # will write this to FITS header so we can combine runs
    my $numtotframes = which($cumframedist>=$num_protons_per_run)->index(0) + 1;
    # this is wrong, because it will remove the last bin which we need
    # it's also unnecessary
    #$cumframedist = $cumframedist->where($cumframedist<=$num_protons_per_run);

    # running variables
    my $numevents = 0;
    my $numtotblobs = 0;

    # loop through the four quadrant data files for this run
    # now combine all four, since single primary can produce signal
    # in multiple quadrants
    foreach my $infile (@infiles) {

        chomp $infile;
        #print "### Reading $infile\n";
        $infile =~ /([0-9]+)_detector([0-9]+)/;
        my $this_detector = $2;
        my %ptype;
        my %cproc;
        my $this_eid;
        open (IN, $infile) or die "### ERROR: Cannot open $infile: $!\n";

        # step through the input file and accumulate primaries
        while (<IN>) {
            next if (/^#/);          # skip comments
            next if (/^\s*$/);       # skip blank lines
            next if (not /,/);       # skip lines with just line counts
            chomp;
            my @fields = split(/,/);
            if ($fields[0] =~ /[a-zA-Z]/) {
                # if the first column is a string, then this is a particle line
                # retain the primary for this interaction
                $this_eid = $fields[1];
                $run = $run->append($this_run);
                $detector = $detector->append($this_detector);
                $eid = $eid->append($this_eid);
                $particleid = $particleid->append($fields[2]);
                $parentid = $parentid->append($fields[3]);
                # particle type and interaction type are hashes so
                # that the pixel-specific read can pick them up
                # doesn't matter if the particle ID is re-used from
                # primary to primary, since this will reset it
                $ptype{$fields[2]} = $fields[0];
                $cproc{$fields[2]} = $fields[4];
            }else{
                # if the first column is a number, then this is a pixel hit line
                # skip it if less than split threshold is deposited, since that is ~ the
                # lower threshold of pixels we'll get
                next if ($fields[2] <= $splitth);
                # skip it if it's outside the 512x512 region of a quad
                my $tmp_x = $fields[0];
                my $tmp_y = $fields[1];
                next if ($tmp_x<$xydep_min or $tmp_y<$xydep_min or $tmp_x>$xydep_max or $tmp_y>$xydep_max);
                $xdep = $xdep->append($tmp_x);
                $ydep = $ydep->append($tmp_y);
                $endep = $endep->append($fields[2]);
                $rundep = $rundep->append($this_run);
                $detectordep = $detectordep->append($this_detector);
                $eiddep = $eiddep->append($this_eid);
                $framedep = $framedep->append(0);
                $piddep = $piddep->append($fields[3]);
                # %ptype is hash of particle type strings indexed by the id
                # %ptypes is (constant) hash of my own particle type IDs indexed
                # by the string (confused yet?)
                $ptypedep = $ptypedep->append($ptypes{$ptype{$fields[3]}});
                $blobid = $blobid->append(pdl(long, 0));
                #print "### piddep = $fields[3]\n";
                #print "### ptype{piddep} = ".$ptype{$fields[3]}."\n";
                #print "### ptypes{ptype{piddep}} = ".$ptypes{$ptype{$fields[3]}}."\n";
            }
        }
        close IN;
    }
     
    my $uniq_eid = $eid->uniq;
    my $numprimaries = $uniq_eid->nelem;
    my $primary_flux = $numprimaries / $texp_run / $detector_area;
    #print "### Run $this_run: found $numprimaries primaries that interacted.\n";
    #print "### Run $this_run: that's $primary_flux protons per sec per cm2.\n";
    my $numpixels = $endep->nelem;
    #print "### Run $this_run: found $numpixels pixels with deposited energy\n";

    # loop through unique primaries to sort them into frames
    # also figure out the number of primaries with signal in 
    # multiple quadrants, just as a diagnostic
    @num_in_different_quadrants = ( 0, 0, 0, 0 );
    for (my $i=0;$i<$numprimaries;$i++) {
        my $primary = $uniq_eid -> index($i); #(($i));
        my $indx = which($eiddep==$primary);
        #print "#####################\n";
        #print "Doing primary $primary.\n";
        #print "eiddep: ".$eiddep->index($indx)."\n";
        #print "OLD framedep ".$framedep->index($indx)."\n";
        # assign each primary to a frame
        # first get the frame ID (indexed starts at 0)
        #print "".$cumframedist->((-1))."\n";
        my $frame = which($cumframedist>=$primary)-> index(0);
        #print "THIS IS FRAME $frame\n";
        # then set the primary and pixel piddles
        $framedep->index($indx) .= $frame;
        #print "NEW framedep ".$framedep->index($indx)."\n";
        my $num_quadrants = $detectordep->index($indx)->uniq->nelem;
        $num_in_different_quadrants[$num_quadrants-1]++;
    }
    
    $frame_frame = $frame_frame->append($framedep->uniq);
    my $numframes = $frame_frame->nelem;
    # now we can append to frame piddles since we know how many there are
    $frame_runid = $frame_runid->append(zeros($numframes) + $this_run);
    $frame_time = $frame_time->append(zeros(double, $numframes));
    $frame_npix = $frame_npix->append(zeros($numframes));
    $frame_npixmip = $frame_npixmip->append(zeros($numframes));
    $frame_nevt = $frame_nevt->append(zeros($numframes));
    $frame_nevtgood = $frame_nevtgood->append(zeros($numframes));
    $frame_nevt27 = $frame_nevt27->append(zeros($numframes));
    $frame_nblob = $frame_nblob->append(zeros($numframes));
    $frame_nprim = $frame_nprim->append(zeros($numframes));

    my $pct_interact = 100. * $numprimaries / $num_protons_per_run;

    #$numframes = ($level <= 5) ? 1: $numframes;
    $numframes = 1;
    for (my $i=0;$i<$numframes;$i++) {

        # set the frame ID
        my $frame = $frame_frame->index($i);

        # set the frame time
        $frame_time->index($i) .= $frame * $texp_frame;

        # keep us updated
        if ($i % 5000 == 0) {
            #print "### Run $this_run: done $i of $numframes frames\n";
        }

        # are we writing out?
        if ($skip_writes > 0 and $i % $skip_writes == 0) {
            $writeit = 1;
        } else {
            $writeit = 0;
        }

        my $pixel_indx = which($framedep==$frame);
        my $x = $xdep->index($pixel_indx)->sever;
        my $y = $ydep->index($pixel_indx)->sever;
        my $en = $endep->index($pixel_indx)->sever;
        my $pixptype = $ptypedep->index($pixel_indx)->sever;
        my $pixprimid = $eiddep->index($pixel_indx)->sever;
        # we want to populate this so don't sever it
        my $this_blobid = $blobid->index($pixel_indx);
        my $xoff = $x->min + $x_offset - 2;
        my $yoff = $y->min + $y_offset - 2;
        $x -= $x->min;
        $x += 2;
        $y -= $y->min;
        $y += 2;
        $xdim = $x->max+3;
        $ydim = $y->max+3;
        
        my $img = zeros(double, $xdim, $ydim);
        my $ptypeimg = zeros(long, $xdim, $ydim);
        my $primidimg = zeros(long, $xdim, $ydim);
        my $img_xvals = xvals($xdim, $ydim);
        my $img_yvals = yvals($xdim, $ydim);
        #for ($j=0;$j<$en->nelem;$j++) {
        #     $img($x(($j)),$y(($j))) .= $en(($j));
        #     $ptypeimg($x(($j)),$y(($j))) .= $pixptype(($j));
        #}
        # better way to do the above mapping to an image without a loop
        # indexND wants a 2xN piddle, where N is the number of pixels
        my $coos = $x->cat($y)->xchg(0,1);
        $img->indexND($coos) .= $en;
        $ptypeimg->indexND($coos) .= $pixptype;
        $primidimg->indexND($coos) .= $pixprimid;

        #print "$pixprimid\n$x\n$y\n$en\n";
        #print "$img\n";

        # add to some frame piddles
        $frame_npix->index($i) .= $pixel_indx->nelem;
        $frame_npixmip->index($i) .= which($en >= $mipthresh)->nelem;
        $frame_nprim->index($i) .= $pixprimid->uniq->nelem;
        
        $writeit = 0;
        if ($writeit) {
               $fitsfile = "output/run${this_run}_frame${frame}_img.fits";
               print "### Writing raw image $fitsfile.\n";
               wfits $img, $fitsfile;
        }
        
        #############################################
        # find blobs

        # segment image into blobs
        # original IDL code
        # blobRegions=label_region(phimg * (phimg GE evtthresh), /ulong)
        # PDL 'cc8compt' is equivalentish to IDL 'label_regions'
        $blobimg = cc8compt($img>$evtth);

        # the first blob is 1, not 0 (which indicates no blob)
        # blobadjust decrements blob IDs when we chuck one, so the IDs are continuous from 1
        my $blobadjust = 0;
        for (my $j=1;$j<=$blobimg->max;$j++) { 
            my $indx = which($blobimg==$j);
            my $indx2d = whichND($blobimg==$j);
            # set the blob to zeros and skip it if there are too few elements
            # if it contains a MIP, it's a good blob irregardless
            if ($indx->nelem < $npixthresh and $img->flat->index($indx)->max < $mipthresh) {
                $blobimg->flat->index($indx) .= 0;
                $blobadjust++;
                next;
            }
            $blobimg->flat->index($indx) -= $blobadjust;

            # this is the running blobid which we need to add to blob piddles

            $blob_blobid = $blob_blobid->append($numtotblobs + $j - $blobadjust);
            $blob_frame = $blob_frame->append($frame);
            $blob_npix = $blob_npix->append($indx->nelem);
            # calculate unclipped blob centroid and summed energy
            my $tmp_en = $img->flat->index($indx);
            my $tmp_toten = $tmp_en->sum;
            my $tmp_wtd_x = $img_xvals->flat->index($indx) * $tmp_en;
            my $tmp_wtd_y = $img_yvals->flat->index($indx) * $tmp_en;
            $blob_cenx = $blob_cenx->append($xoff + $tmp_wtd_x->sum / $tmp_toten);
            $blob_ceny = $blob_ceny->append($yoff + $tmp_wtd_y->sum / $tmp_toten);
            $blob_energy = $blob_energy->append($tmp_toten);
            # calculate clipped blob centroid and summed energy
            $tmp_en = $tmp_en->clip(undef, $clip_energy);
            $tmp_toten = $tmp_en->sum;
            $tmp_wtd_x = $img_xvals->flat->index($indx) * $tmp_en;
            $tmp_wtd_y = $img_yvals->flat->index($indx) * $tmp_en;
            $blob_cenxcl = $blob_cenxcl->append($xoff + $tmp_wtd_x->sum / $tmp_toten);
            $blob_cenycl = $blob_cenycl->append($yoff + $tmp_wtd_y->sum / $tmp_toten);
            $blob_energycl = $blob_energycl->append($tmp_toten);
        }
        
        open(FH, '>', $wf_name."1".".txt") or die "can't write";
        
        foreach (($blob_blobid, $blob_frame, $blob_npix, $blob_cenx, $blob_ceny, $blob_energy, $blob_cenxcl, $blob_cenycl, $blob_energycl)) {
            print_all $_;
        }
        close FH;
        
        if ($level >= 2) {
            $frame_nblob->index($i) .= $blobimg->max;
            $blobimg->where($blobimg>0) += $numtotblobs;
            $numtotblobs = ($blobimg->max > 0) ? $blobimg->max : $numtotblobs;

            # mark each pixel in a blob
            $this_blobid .= $blobimg->indexND($coos);
            $blobimg->where($blobimg>0) .= 1;
            $mipimg = $img->copy;
            $mipimg->where($mipimg<$mipthresh) .= 0.;
            $mipimg->where($mipimg>0) .= 1;
            my $indx_thcross = which($img > $evtth)->long;
            my $indx2d_thcross = whichND($img > $evtth)->long;
            my $num_thcross = $indx_thcross->nelem;
            
            my $evtx; my $evty; my $evt_phas; my $evt_ptype; my $evt_primid;
            my $localmax;
            my $num_localmax;
            if ($num_thcross > 0) {
                # get piddles containing X,Y coords of threshold crossings
                $evtx = $indx2d_thcross->index(0);
                $evty = $indx2d_thcross->index(1);
                my $tmp_evtx = $evtx + $xoff;
                my $tmp_evty = $evty + $yoff;
                no warnings;
                my ($indx_border,$indx_notborder) = which_both(
                    $tmp_evtx<=1 | ($tmp_evtx>=510 & $tmp_evtx<=516) | $tmp_evtx>=1025 | 
                    $tmp_evty<=1 | ($tmp_evty>=510 & $tmp_evty<=516) | $tmp_evty>=1025 | 
                    $img->index2d($evtx,$evty)<0
                );
                use warnings;
                my $num_notborder = $indx_notborder->nelem;
                
                open(FH, '>', $wf_name."2.txt") or die "can't write";

                foreach (($frame_nblob, $blobimg, $mipimg)) {
                    print_nz $_;
                }
                foreach (($tmp_evtx, $tmp_evty, $indx_border, $indx_notborder)) {
                    print_all $_;
                }
                
                print FH "$numtotblobs\n";
                close FH;
                
                if ($level >= 3) {
                    if ($num_notborder > 0) {
                        $evtx = $evtx->index($indx_notborder);
                        $evty = $evty->index($indx_notborder);
                        $indx_thcross = $indx_thcross->index($indx_notborder);
                        # find local maxima
                        # make a copy of the threshold crossing piddle, which we will
                        # subsitute with the largest neighbor value and then compare to the original
                        $localmax = $img->copy;
                        # compare each neighbor in the original with the central value in the copy
                        # and replace the central value if the neighbor is larger
                        # indx is index (into just thcross piddle) of pixels which have a larger neighbor
                        # NB: formerly:
                        #  left, upper left, upper, and upper right neighbors are
                        #  incremented before comparison to avoid equal split; these
                        #  pixels take precedence for local max
                        # but this is wrong, because PHAS are floating point, and there is
                        # very little chance two neighboring pixels will have the same PHA value
                        # see primary 1900 in run 1 for evidence of that
                        #
                        # left neighbor
                        my $indx = which($img->index2d($evtx-1,$evty) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx-1,$evty)->flat->index($indx);
                        # upper left neighbor
                        $indx = which($img->index2d($evtx-1,$evty+1) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx-1,$evty+1)->flat->index($indx);
                        # upper neighbor
                        $indx = which($img->index2d($evtx,$evty+1) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx,$evty+1)->flat->index($indx);
                        # upper right neighbor
                        $indx = which($img->index2d($evtx+1,$evty+1) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx+1,$evty+1)->flat->index($indx);
                        # right neighbor
                        $indx = which($img->index2d($evtx+1,$evty) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx+1,$evty)->flat->index($indx);
                        # lower right neighbor
                        $indx = which($img->index2d($evtx+1,$evty-1) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx+1,$evty-1)->flat->index($indx);
                        # lower neighbor
                        $indx = which($img->index2d($evtx,$evty-1) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx,$evty-1)->flat->index($indx);
                        # lower left neighbor
                        $indx = which($img->index2d($evtx-1,$evty-1) > $localmax->index2d($evtx,$evty));
                        $localmax->flat->index($indx_thcross->index($indx)) .= $img->index2d($evtx-1,$evty-1)->flat->index($indx);
                        # finally compare the original central pixel pulseheight
                        # with the maximum neighbor; if the former is greater than
                        # or equal to the latter, it's a local maximum
                        our $indx_localmax = which($img->flat->index($indx_thcross) >= $localmax->flat->index($indx_thcross));
                        $num_localmax = $indx_localmax->nelem;
                        
                    }
                    open(FH, '>', $wf_name."3.txt") or die "can't write";
                    foreach (($localmax,)) {
                        print_nz $_;
                    }
                    foreach (($evtx, $evty, $indx_thcross, $indx_localmax)) {
                        print_all $_;
                    }
                    print FH "$num_localmax\n";
                    close FH;
                    
                    if ($level >= 4) {
                        if ($num_localmax > 0) {
                             $evtx = $evtx->index($indx_localmax);
                             $evty = $evty->index($indx_localmax);
                             $indx_thcross = $indx_thcross->index($indx_localmax);
                             # get 3x3 PHAS piddle in correct
                             # 1D order (central pixel first)
                             $evt_phas = zeros(double, 25, $num_localmax);
                             $evt_phas((0),:) .= $img->index2d($evtx,$evty);
                             $evt_phas((1),:) .= $img->index2d($evtx-1,$evty-1);
                             $evt_phas((2),:) .= $img->index2d($evtx,$evty-1);
                             $evt_phas((3),:) .= $img->index2d($evtx+1,$evty-1);
                             $evt_phas((4),:) .= $img->index2d($evtx-1,$evty);
                             $evt_phas((5),:) .= $img->index2d($evtx+1,$evty);
                             $evt_phas((6),:) .= $img->index2d($evtx-1,$evty+1);
                             $evt_phas((7),:) .= $img->index2d($evtx,$evty+1);
                             $evt_phas((8),:) .= $img->index2d($evtx+1,$evty+1);
                             # get outer 5x5 PHAS piddle
                             # first pad the image with zeros all around
                             our $tmp_img = zeros(double, $img->shape+2);
                             $tmp_img(1:$img->shape->((0)),1:$img->shape->((1))) .= $img;
                             # now evtx and evty are too low by one, so correct for that in index
                             $evt_phas((9),:) .= $tmp_img->index2d($evtx-1,$evty-1);
                             $evt_phas((10),:) .= $tmp_img->index2d($evtx,$evty-1);
                             $evt_phas((11),:) .= $tmp_img->index2d($evtx+1,$evty-1);
                             $evt_phas((12),:) .= $tmp_img->index2d($evtx+2,$evty-1);
                             $evt_phas((13),:) .= $tmp_img->index2d($evtx+3,$evty-1);
                             $evt_phas((14),:) .= $tmp_img->index2d($evtx-1,$evty);
                             $evt_phas((15),:) .= $tmp_img->index2d($evtx+3,$evty);
                             $evt_phas((16),:) .= $tmp_img->index2d($evtx-1,$evty+1);
                             $evt_phas((17),:) .= $tmp_img->index2d($evtx+3,$evty+1);
                             $evt_phas((18),:) .= $tmp_img->index2d($evtx-1,$evty+2);
                             $evt_phas((19),:) .= $tmp_img->index2d($evtx+3,$evty+2);
                             $evt_phas((20),:) .= $tmp_img->index2d($evtx-1,$evty+3);
                             $evt_phas((21),:) .= $tmp_img->index2d($evtx,$evty+3);
                             $evt_phas((22),:) .= $tmp_img->index2d($evtx+1,$evty+3);
                             $evt_phas((23),:) .= $tmp_img->index2d($evtx+2,$evty+3);
                             $evt_phas((24),:) .= $tmp_img->index2d($evtx+3,$evty+3);

                             # set the particle type for all these events
                             # based on whatever ptype produced the local max
                             $evt_ptype = $ptypeimg->index2d($evtx,$evty)->sever->long;
                             $evt_primid = $primidimg->index2d($evtx,$evty)->sever->long;
                        } else {
                             #print "### Found no events, skipping processing.\n";
                             next;
                        } # end if there are localmaxes
                        open(FH, '>', $wf_name."4.txt") or die "can't write";
                        
                        foreach (($img, $tmp_img, $evt_phas)) {
                            print_nz $_;
                        }
                        
                        foreach (($evtx, $evty, $indx_thcross, $evt_ptype, $evt_primid)) {
                            print_all $_;
                        }
                        close FH;
                    } else {
                        #print "### Found no threshold crossings off border, skipping this frame.\n";
                        next;
                    } # end if there are threshold crossing not on borders
                } else {
                    #print "### Found no threshold crossings, skipping this frame.\n";
                    next;
                } # end if there are any threshold crossings
                
                if ($level >= 5) {
                    my $indx_blobs = whichND($blobimg>0);
                    ## get X,Y coords of any pixels in a mip (don't matter which mip)
                    my $indx_mips = whichND($mipimg>0);

                    my $numevents_thisframe = $num_localmax;
                    #print "### Found $numevents_thisframe events, processing.\n";

                    # process the detected events

                    # append events to running piddles
                    $runid = $runid->append(zeros(long, $numevents_thisframe)+pdl(long, $this_run));
                    $detectorid = $detectorid->append(zeros(long, $numevents_thisframe)+pdl(long, $this_run));
                    $primid = $primid->append($evt_primid);
                    $evt_frame = $evt_frame->append(zeros(long, $numevents_thisframe)+pdl(long, $frame));
                    #$actx = $actx->append(zeros(long, $numevents_thisframe));
                    #$acty = $acty->append(zeros(long, $numevents_thisframe));
                    #$phas = $phas->glue(1, zeros(long, 25, $numevents_thisframe));
                    $actx = $actx->append($evtx);
                    $acty = $acty->append($evty);
                    $phas = $phas->glue(1, $evt_phas);
                    $pha = $pha->append(zeros(double, $numevents_thisframe));
                    $ptype = $ptype->append($evt_ptype);
                    $evttype = $evttype->append(zeros(long, $numevents_thisframe));
                    $energy = $energy->append(zeros(double, $numevents_thisframe));
                    $blobdist = $blobdist->append(zeros(double, $numevents_thisframe));
                    $mipdist = $mipdist->append(zeros(double, $numevents_thisframe));
                    $pattern = $pattern->append(zeros(ushort, $numevents_thisframe));
                    $vfaint = $vfaint->append(zeros(ushort, $numevents_thisframe));
                    
                    open(FH, '>', $wf_name."5.txt") or die "can't write";
                        
                    #foreach (($img, $tmp_img, $evt_phas)) {
                    #    print_nz $_;
                    #}
                    foreach (($indx_blobs, $indx_mips, $phas)) {
                        print_nz $_;
                    }
                         
                    foreach (($runid, $detectorid, $evt_frame, $pha)) {
                        print_all $_;
                    }
                    print FH "$numevents_thisframe\n";
                    close FH;
                    
                    if ($level >= 6) {
                        for (my $j=$numevents; $j<($numevents+$numevents_thisframe); $j++) {
                            # below is already done in event finding and pasted above
                            # get X,Y of center pixel
                            $x = $actx(($j));
                            $y = $acty(($j));
                            # below is deprecated, since we've assumed this in event finding
                            ## eliminate events on edges so we can use 5x5
                            #next if ($x<2 or $x>1021 or $y<2 or $y>1021);

                            # get ACIS flt grade; "which" returns the indices of
                            # non-central pixels greater than or equal to event threshold,
                            # these are used as bits to raise 2 to the power, and summed
                            # (this can probably be removed from the loop if I'm smart)
                            my $indx = which($phas(1:8,($j))>=$splitth);
                            my $fltgrade = (pow(2,$indx))->sum;
                            # convert to EPIC-pn pattern (from look-up table)
                            $pattern(($j)) .= $epicpn_pattern_table(($fltgrade));

                            # sum 3x3 pixels over split threshold to get PHA
                            # this is an ACIS PHA, not Suzaku
                            # (this can probably be removed from the loop if I'm smart)
                            $pha(($j)) .= ($phas(:,($j))->where($phas(:,($j))>=$splitth))->sum;
                            # apply gain correction for this node
                            # get the gain parameters from the node (0-3)
                            #print "".($phas(:,($j)))." $gain_intercept $gain_slope ".$pha($j)."\n";
                            $pha(($j)) .= ($pha(($j)) - $gain_intercept) * $gain_slope / $evperchan;

                            # convert pha to energy
                            # (this can probably be removed from the loop if I'm smart)
                            $energy(($j)) .= $pha(($j)) * $evperchan / 1000.;
                            #print "".($phas(:,($j)))." $gain_intercept $gain_slope ".$pha($j),$energy($j)."\n";

                            # perform simple VFAINT filtering; also update the EPIC-pn pattern
                            # of doubles, triples, and quads based on it
                            # get outer 16 pixels and if any are above split flag it
                            my $noutermost = which($phas(9:24,($j)) > $splitth)->nelem;
                            if ($noutermost > 0) {
                                $vfaint(($j)) .= 1;
                                # EDM Fri Jan  4 11:24:20 EST 2019 
                                # Reinstated the 5x5 filtering on PATTERN.
                                # EDM Thu May 31 13:46:28 EDT 2018 
                                # change to remove 5x5 criterion on EPIC-pn patterns
                                # for doubles, triples, quads
                                if ($pattern(($j)) > 0) { 
                                    $pattern(($j)) .= 13;
                                }
                            }

                            # get minimum distance to a blob
                            # first find delta X,Y from this event to list of blob pixels
                            my $delta_blobs = $indx_blobs - pdl($x,$y);
                            #print $delta_blobs;
                            # square that, sum it, square root it to get the distance
                            if ($delta_blobs->shape->((1)) > 0) {
                                $blobdist(($j)) .= $delta_blobs->pow(2)->sumover->sqrt->min;
                                # unless there aren't any blobs, in which case set blobdist to -1
                            } else {
                                $blobdist(($j)) .= -1;
                            }
                            #print "".$blobdist(($j))."\n";

                            # get minimum distance to a mip
                            # first find delta X,Y from this event to list of mip pixels
                            our $delta_mips = $indx_mips - pdl($x,$y);
                            
                            #print $delta_mips;
                            # square that, sum it, square root it to get the distance
                            if ($delta_mips->shape->((1)) > 0) {
                                $mipdist(($j)) .= $delta_mips->pow(2)->sumover->sqrt->min;
                                # unless there aren't any mips, in which case set mipdist to -1
                            } else {
                                $mipdist(($j)) .= -1;
                            }
                            #print "".$mipdist(($j))."\n";

                            # we really want ACTX,Y in real WFI coords in the event list, 
                            # so fix that here; $x and $y are children of $actx and $acty
                            # for this event, which is why this works
                            $x += $xoff;
                            $y += $yoff;

                            # add info to frame piddles
                            $frame_nevt(($i))++;
                            if ($energy(($j)) < $mipthresh) {
                                $frame_nevtgood(($i))++;
                                if ($energy(($j)) >= 2 and $energy(($j)) < 7) {
                                    $frame_nevt27(($i))++;
                                }
                            }
                        } # done loop through events

                        open(FH, '>', $wf_name."6.txt") or die "can't write";
                        
                        foreach (($frame_nevt, $frame_nevtgood, $frame_nevt27)) {
                            print_nz $_;
                        }

                        foreach (($pattern, $pha, $energy, $blobdist, $mipdist)) {
                            print_all $_;
                        }
                        print FH "$x\n$y\n";
                        close FH;

                    } 
                }
            }
        }
    }
}

print "All done perl\n";

my $p =  `python3 blob_test.py $level $start_run`;
print $p;
if (($p cmp "All good python\n") == 0) {
    print "Test passed\n";
} else {
    print "Test failed\n";
}
