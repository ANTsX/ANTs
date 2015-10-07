#!/usr/bin/perl -w

use strict;

# Usage: waitForSlurmJobs.pl <verbose [1 or 0]> <delay in seconds in range 10-600> [job IDs]
#
#
# Takes as args a string of sbatch job IDs and periodically monitors them. Once they all finish, it returns 0
#
# If any of the jobs go into error state, an error is printed to stderr and the program waits for the non-error
# jobs to finish, then returns 1

my ( $verbose, $delay, @jobIDs ) = @ARGV;

my %COMPLETION_STATES = map {
    $_ => 1
} qw( COMPLETED );

my %FAILURE_STATES = map {
    $_ => 1
} qw(
    CANCELLED FAILED NODE_FAIL PREEMPTED TIMEOUT
    );

# Validate that the delay is within the acceptable range
if ($delay < 10) {
    print STDERR "Sleep period is too short, will poll queue once every 10 seconds\n";
    $delay = 10;
} elsif ($delay > 3600) {
    print STDERR "Sleep period is too long, will poll queue once every 60 minutes\n";
    $delay = 3600;
}


print "  Waiting for " . scalar( @jobIDs ) . " jobs: @jobIDs\n";

my $errorsEncountered = 0;

wait_for_all_jobs_to_complete(@jobIDs);

if ($errorsEncountered) {
    print "  No more jobs to run - some jobs had errors\n\n";
    exit 1;
}
else {
    print "  No more jobs in queue\n\n";
    exit 0;
}



sub wait_for_all_jobs_to_complete {
    my @pendingJobs = update_pending_jobs(@_);

    while (@pendingJobs) {
        if ($verbose) {
            my $timestamp = `date`;
            chomp $timestamp;
            printf "  ($timestamp) Still waiting for %d jobs\n\n", scalar(@pendingJobs);
        }

        # Use of backticks rather than system permits a ctrl+c to work
        `sleep $delay`;

        @pendingJobs = update_pending_jobs(@pendingJobs);
    };
}

sub update_pending_jobs {
    my (@jobsToQuery) = @_;
    my %jobStatuses = query_job_statuses(@jobsToQuery);

    if (!scalar(@jobsToQuery) || !%jobStatuses)
    {
        # No more jobs remain
        return ();
    }

    if ($verbose) {
        while (my ($job, $status) = each(%jobStatuses)) {
            print("    Job $job is in state $status\n");
        }
    }

    my @terminatedJobs = grep {
        !exists($jobStatuses{$jobsToQuery[$_]}) || exists($COMPLETION_STATES{$jobStatuses{$jobsToQuery[$_]}})
    } 0..$#jobsToQuery;

    my @failedJobs = grep {
        exists($jobStatuses{$jobsToQuery[$_]}) && exists($FAILURE_STATES{$jobStatuses{$jobsToQuery[$_]}})
    } 0..$#jobsToQuery;

    if (@failedJobs) {
        $errorsEncountered = 1;
    }

    push @terminatedJobs, @failedJobs;
    foreach my $index (reverse(@terminatedJobs)) {
        splice @jobsToQuery, $index, 1;
    }

    return @jobsToQuery;
}

sub query_job_statuses {
    my (@jobsToQuery) = @_;
    my $user = trim(`whoami`);
    my $squeueOutput = qx/squeue --noheader --user="$user" --format="%i,%T" --jobs=${\join(',', @jobsToQuery)}/;
    my $exitcode = $? >> 8;
    my %jobStatuses = ();

    if ($exitcode == 0) {
        %jobStatuses = map {
            my @statusParts = split(",", $_);
            $statusParts[0] => $statusParts[1];
        } split("\n", trim($squeueOutput));
    }

    return %jobStatuses;
}

sub trim {
    my ($string) = @_;

    $string =~ s/^\s+//;
    $string =~ s/\s+$//;

    return $string;
}
