#!/usr/bin/perl -w

use strict;

# Usage: waitForSGEQJobs.pl <verbose [1 or 0]> <delay in seconds in range 10-600> [job IDs]
#
#
# Takes as args a string of qsub job IDs and periodically monitors them. Once they all finish, it returns 0
#
# If any of the jobs go into error state, an error is printed to stderr and the program waits for the non-error
# jobs to finish, then returns 1
#

# Usual qstat format - check this at run time
# job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID


# First thing to do is parse our input

my ( $verbose, $delay, @jobIDs ) = @ARGV;

# Check for user stupidity
if( $delay < 10 )
  {
  print STDERR "Sleep period is too short, will poll queue once every 10 seconds\n";
  $delay = 10;
  }
elsif( $delay > 3600 )
  {
  print STDERR "Sleep period is too long, will poll queue once every 60 minutes\n";
  $delay = 3600;
  }

print "  Waiting for " . scalar( @jobIDs ) . " jobs: @jobIDs\n";

my $user=`whoami`;

my $qstatOutput = `qstat -u $user`;

if( !scalar(@jobIDs) || !$qstatOutput )
  {
  # Nothing to do
  exit 0;
  }

my @qstatLines = split("\n", $qstatOutput);
# my @header = split('\s+', trim($qstatLines[0]));

# Position in qstat output of tokens we want
# Here we hardcode the values that work at UVa
my $jobID_Pos = 0;
my $statePos = 9;

# foreach my $i (0..$#header)
#   {
#   if ( $header[$i] eq "Job ID" )
#     {
#     $jobID_Pos = $i;
#     }
#   elsif ($header[$i] eq "state")
#     {
#     $statePos = $i;
#     }
#   }


# If we can't parse the job IDs, something is very wrong
# if ($jobID_Pos < 0 || $statePos < 0)
#   {
#   die "Cannot find job-ID and state field in qstat output, cannot monitor jobs\n";
#   }

# Now check on all of our jobs
my $jobsIncomplete = 1;

# Set to 1 for any job in an error state
my $haveErrors = 0;

while( $jobsIncomplete )
  {

  # Jobs that are still showing up in qstat
  $jobsIncomplete = 0;

  foreach my $job (@jobIDs)
    {
    # iterate over all user jobs in the queue
    qstatLine: foreach my $line ( @qstatLines )
      {
      # trim string for trailing white space so that the tokens are in the correct sequence
      # We are being paranoid by matching tokens to job-IDs this way. Less elegant than a
      # match but also less chance of a false-positive match

      my @tokens = split( '\s+', trim( $line ) );

      # The qstat command only prints the first 15 characters of the job
      # so we only compare the first 15 characters

      my $job_short = substr( $job, 0, 15 );

      if( @tokens > 0 && ( $tokens[$jobID_Pos] =~ m/$job_short/ ) )
        {
	# Check status - there's no error state in PBS
        # so we simply skip over this check
#        if( $tokens[$statePos] =~ m/E/ )
#         {
#         $haveErrors = 1;
#         }
#       else
#         {
          $jobsIncomplete = $jobsIncomplete + 1;
#	  }
        if( $verbose )
          {
	  print "    Job $job is in state $tokens[$statePos]\n";
	  }
	}
      last qstatLine if ( @tokens > 0 && ( $tokens[$jobID_Pos] =~ m/$job_short/ ) );
      }
    }


  if( $jobsIncomplete )
    {
    if( $verbose )
      {
      my $timestamp = `date`;
      chomp $timestamp;
      print "  ($timestamp) Still waiting for $jobsIncomplete jobs\n\n";
      }

    # Use of backticks rather than system permits a ctrl+c to work
    `sleep $delay`;
    $qstatOutput = `qstat -u $user`;
    @qstatLines = split("\n", $qstatOutput);
    }
  }

if ($haveErrors) {
    print "  No more jobs to run - some jobs had errors\n\n";
    exit 1;
}
else {
    print "  No more jobs in queue\n\n";
    exit 0;
}




sub trim {

    my ($string) = @_;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
