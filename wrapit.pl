#!/usr/bin/env perl

# EDM Thu Nov  5 12:16:05 EST 2020
# Generic, pseudo-parallelization wrapper script.

use File::Basename;
use strict;

our $usage =<<EOUsage;
Usage:

   $0 <executable> <arg1> [ <arg2> <arg3> ... ]

where:
     <executable> is the script to run, with correct path from the
          current directory
     <arg1> etc. is a list of arguments that <executable> is expecting.
          These will be fed one-at-a-time to <executable> in parallel forks
          every WAITSEC seconds until MAXPROCS simultaneous forks are
          running. It will then check every WAITSEC seconds to see if any
          have finished, at which point it will fork a new one.
EOUsage

unless ($#ARGV>0) { die "$usage\n"; }

our $EXEC = shift @ARGV;
our @ARGS = @ARGV;

our $MAXPROCS = 20;
our $WAITSEC = 2;
our $PS = "/bin/ps auxww";

## -- uniq @ARGS, sorted in numerical order
#@ARGS = sort {$a <=> $b} keys %{{ map { $_ => 1 } @ARGS }}; 
#print "### Doing ".($#ARGS+1)." run(s): @ARGS\n";

# -- initialize some things
our @pids = ();
our $current = 0;
our $forkit = 0;
our $pid;
# -- step through runs until we're done
while ($#ARGS>=0) {
     # -- check to see if any children have terminated; if so, set forkit
     my @tmppids = ();
     foreach $pid (@pids) {
          # -- check if PID has exited
          my $sig = waitpid $pid, 1;
          if ($sig == 0) { 
               #print "### Process $pid returned $sig, exists, keeping.\n";
               push @tmppids, $pid;
          # -- otherwise don't keep it and force a fork
          } elsif ($sig == $pid) {
               print "### Process $pid returned $sig, ended, harvesting.\n";
               $forkit = 1; 
          } else {
               die "### ERROR: Process $pid returned $sig, has an issue: $?\n";
          }
     }
     @pids = @tmppids;
     # -- check to see how many forks we've done and compare to max allowed
     if ($#pids < ($MAXPROCS-1))  { $forkit = 1; }
     print "### Running ".($#pids+1)." processes, forkit = $forkit\n";

     # -- fork it?
     if ($forkit) {
          if ($pid = fork) {
               # -- need to pop off the run here too
               my $run = shift(@ARGS);
               # -- add current PID to list
               push @pids, $pid;
               print "### PARENT: Fork is started, child has PID $pid.\n";
          } elsif (defined $pid) {
               # -- pop off the next run
               my $run = shift(@ARGS);
               my $command = "${EXEC} ${run}";
               #my $command = "sleep 10";
               print "### CHILD STARTED: $command\n";
               system("$command");
               print "### CHILD ENDING:  $command\n";
               exit;
          } else {
               die "### ERROR: Cannot fork: $!\n";
          }
     }

     # -- reset forkit
     $forkit = 0;

     # -- nap for a bit
     sleep $WAITSEC;
}

# -- wait for all children to terminate
print "### Waiting for children to terminate.\n";
foreach $pid (@pids) {
     waitpid $pid, 0;
}
exit;
