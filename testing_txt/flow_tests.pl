use PDL;
use warnings;

if ($#ARGV < 0) {
    die "Usage: $0 <test_number> [<files>]\n";
} elsif ($#ARGV == 0) {
    $inp_file = '1_detector0';
} elsif ($#ARGV >= 1) {
    $inp_file = $ARGV[1];
}

if (($ARGV[0] cmp "all") == 0) {
    @do_tests = ('1','2', '3');
    @td = ('1','2', '3');
} else {
    @do_tests = ($ARGV[0],);
    @td = ($ARGV[0],);
}

@passed = (0,0,0);
$failed = 0;

$pwd = `pwd`;
$input_dir = substr($pwd, 0, -12);
$infile =  $input_dir.'input/'.$inp_file; 
$infile =~ /([0-9]+)_detector([0-9]+)/;
$this_run = $1;
$this_detector = $2;

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

$evtth = .1;          # 100 eV for WFI Geant4 simulations
$splitth = $evtth;    # same as evtth for WFI Geant4 simulations
$xydep_min = -513;
$xydep_max = 513;


sub test_1 {
    my (%arg) = (
        'wf_name' => 'flow_test1',
    );
    
    open(IN, $infile) or die "### ERROR: Cannot open $infile: $!\n";
    open(FH, '>', $arg{'wf_name'}.'.txt') or die "heck 2";
    while (<IN>) {
       next if (/^#/);          # skip comments
       next if (/^\s*$/);       # skip blank lines
       next if (not /,/);       # skip lines with just line counts
       print FH $_;
       chomp;
       my @fields = split(/,/);
       if ($fields[0] =~ /[a-zA-Z]/) {
           print FH "if\n";
       } else {
           print FH "else\n";
       }
    }
    close(IN);
    close(FH);
    
    print "All done perl\n";
    my $p =  `python3 flow_tests.py 1 $inp_file`;
    print $p;
    if (($p cmp "All good python\n") == 0) {
        $passed[0] = 1;
    }
};

sub test_2 {
    my (%arg) = (
        'wf_name' => 'flow_test2',
    );
    
    open(IN, $infile) or die "### ERROR: Cannot open $infile: $!\n";
    open(FH, '>', $arg{'wf_name'}.'.txt') or die "can't write";
    while (<IN>) {
       next if (/^#/);          # skip comments
       next if (/^\s*$/);       # skip blank lines
       next if (not /,/);       # skip lines with just line counts
       print FH $_;
       chomp;
       my @fields = split(/,/);
       if ($fields[0] =~ /[a-zA-Z]/) {
           print FH "if\n";
       } else {
           next if ($fields[2] <= $splitth);
           my $tmp_x = $fields[0];
           my $tmp_y = $fields[1];
           next if ($tmp_x<$xydep_min or $tmp_y<$xydep_min or $tmp_x>$xydep_max or $tmp_y>$xydep_max);
           print FH "else\n";
       }
    }
    close(IN);
    close(FH);
    
    print "All done perl\n";
    my $p =  `python3 flow_tests.py 2 $inp_file`;
    print $p;
    if (($p cmp "All good python\n") == 0) {
        $passed[1] = 1;
    }
};


sub test_3 {
    my (%arg) = (
        'wf_name' => 'flow_test3',
    );
    
    my $eid = zeros(long, 0);
    my $ptype = zeros(long, 0);
    my $xdep = zeros(long, 0);
    my $ydep = zeros(long, 0);
    my $endep = zeros(double, 0);
    my $rundep = zeros(long, 0);
    my $detectordep = zeros(long, 0);
    my $eiddep = zeros(long, 0);
    my $framedep = zeros(long, 0);
    my $piddep = zeros(long, 0);
    my $ptypedep = zeros(long, 0);
    my $blobid = zeros(long, 0);
    # step through the input file and accumulate primaries
    my $this_eid;
    my %ptype;
    my %cproc;

    open(IN, $infile) or die "### ERROR: Cannot open $infile: $!\n";
    open(FH, '>', $arg{'wf_name'}.'.txt') or die "can't write";
    while (<IN>) {
       next if (/^#/);          # skip comments
       next if (/^\s*$/);       # skip blank lines
       next if (not /,/);       # skip lines with just line counts
       print FH $_;
       chomp;
       my @fields = split(/,/);
       if ($fields[0] =~ /[a-zA-Z]/) {
            $this_eid = $fields[1];
            $eid = $eid->append($this_eid);
            $ptype{$fields[2]} = $fields[0];
            $cproc{$fields[2]} = $fields[4];
            print FH "if\n";
       } else {
            next if ($fields[2] <= $splitth);
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
            $ptypedep = $ptypedep->append($ptypes{$ptype{$fields[3]}});
            $blobid = $blobid->append(pdl(long, 0));
            print FH "else\n";
       } 

    };
    @deps = ($xdep, $ydep, $endep, $rundep, $detectordep, 
             $eiddep, $framedep, $piddep, $ptypedep, $blobid);
    
    print FH "next\n";
    print FH $eid;
    print FH "\nnext\n";
    foreach (keys %ptype) {
        print FH "$_,$ptype{$_}\n";
    }
    print FH "next\n";
    foreach (keys %cproc) {
        print FH "$_,$cproc{$_}\n";
    }
    print FH "next\n";
    
    foreach my $dep (@deps) {
        for ($i = 0; $i < $dep->getdim(0); $i++){
            print FH $dep -> index($i), ","
        }
        print FH "\nnext\n";
    }
    
    
    close IN;
    close FH;
    
    print "All done perl\n";
    my $p =  `python3 flow_tests.py 3 $inp_file`;
    print $p;
    if (($p cmp "All good python\n") == 0) {
        $passed[2] = 1;
    }
    
};



print "Testing $inp_file at $input_dir"."input\n";


foreach (@td) {
    print "Doing test $_....\n";
    eval "test_$_";
    #eval "test_$_";
};

foreach (@do_tests) {
    #print "Evaluating test: $_\n";
    #print eval "$passed[$_ - 1]";
    #print "\n";
    if (eval "$passed[$_ - 1] == 0") {
        print "Test $_ failed\n";
        $failed = 1;
    }
};

if (not $failed) {
    print "All tests passed\n";
}
