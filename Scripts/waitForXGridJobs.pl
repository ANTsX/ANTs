#!/usr/bin/perl
#
# Usage: waitForXGridJobs.pl [-verbose] [-delay N] [-xgridflags string] job_IDs
#
# Minimal call would be (e.g., to wait for a single job ID #200):
#    waitForXGridJobs.pl 200
#
# delay: Poll interval in seconds (default is 30)
# verbose: Be verbose
# xgridflags: Flags to pass to XGrid.  Odds are you are going to want to put this
#             string in quotes as you'll be passing in several flags as this one string.
#
# More complex call:
#   waitForXGridJobs.pl -delay 60 -xgridflags "-p password -h hostname.local" 100 101 102 103 104
#
#
# Craig Stark, June 2010
#
# Code based loosly on waitForSGEJobs.pl

use Getopt::Long;
if ($ARGV[0] =~ /help/) {
   die "syntax: waitForXGridJobs.pl [-verbose] [-delay N] [-xgridflags string] job_IDs\n"; }
$result = GetOptions("verbose" => \$verbose, "delay:n" => \$delay, "xgridflags=s" => \$xgridflags);
if (!$result) { die "Invalid command line options.  Try waitForXGridJobs.pl -help\n"; }
if (!$delay) { $delay = 30; }
if (!$xgridflags) { $xgridflags = ''; }

$njobs = scalar(@ARGV);
if ($njobs == 0) {
	die "No jobs specified.  Try waitForXGridJobs.pl -help\n";
}
@jobIDs = @ARGV;

if ($verbose) {
	print "  Poll interval is $delay, and xgrid flags are: $xgridflags\n";
	print "  Waiting for " . scalar(@jobIDs) . " jobs: @jobIDs\n";
}

# Check for user stupidity
if ($delay < 10) {
    print STDERR "Sleep period is too short, will poll queue once every 10 seconds\n";
    $delay = 10;
}
elsif ($delay > 3600) {
    print STDERR "Sleep period is too long, will poll queue once every 60 minutes\n";
    $delay = 3600;
}

# Now check on all of our jobs
my $jobsIncomplete = 1;

# Set to 1 for any job in an error state
my $haveErrors = 0;

# Find the line number in the output we'll use
#$statline = 0;
#@result = `xgrid $xgridflags -job attributes -id $jobIDs[i]`;
#for ($i=0; $i<scalar(@result); $i++) {

#foreach $i (0..$#result) {
#	if ($result[$i] =~ /jobStatus/) {
#		$statline = $i;
#		if ($verbose) {
#			print "Status indication found on line $statline\n";
#		}
#	}
#}
#if ($statline == 0) {
#	die "Could not find jobStatus in xgrid -job attributes output";
#}


while ($jobsIncomplete) {

    # Jobs that are still showing up in qstat
    $jobsIncomplete = 0;

    foreach my $job (@jobIDs) {
	# Get the status of this ID
	@result = `xgrid $xgridflags -job attributes -id $job`;
	$statline=0;
		foreach $i (0..$#result) {
			if ($result[$i] =~ /jobStatus/) {
				$statline = $i;
				#if ($verbose) {
				#	print "Status indication found on line $statline\n";
				#}
			}
		}
		if ($statline == 0) {
			die "Could not find jobStatus in xgrid -job attributes output";
		}

	if ($result[1] =~ /InvalidJobIdentifier/) {
		die "Invalid job number ($job) specified";
	}
	if ($result[$statline]=~/Finished/) {
		if ($verbose) {
			print "Job $job finished\n";
		}
	}
	elsif ($result[$statline]=~/Failed/) {
		if ($verbose) { print "Job $job failed ****\n"; }
		$haveErrors++;
	}
		else {  # Should be "Running" or "Pending" by this point
		if ($verbose) { print "Job $job still going\n"; }
		$jobsIncomplete++;
	}

    }

    if ($jobsIncomplete && $verbose) {
	    my $timestamp = `date`;
	    chomp $timestamp;
	    print "  ($timestamp) Still waiting for $jobsIncomplete jobs - sleeping $delay\n\n";
		`sleep $delay`;
	}


}
