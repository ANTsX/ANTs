#!/usr/bin/env perl

use strict;

my $masterBranchLabel = "master";

if ($#ARGV < 0) {
    print qq{
$0 <tag>

Run this from the ANTs/ directory to tag a release.

Given a tag "vX.Y.Z", the script will

  1. Check user is on the ${masterBranchLabel} branch, with no uncommitted changes
  2. Update Version.cmake to set the version major, minor, patch
  3. Commit Version.cmake and push to origin ${masterBranchLabel}
  4. git tag with the release tag and push the tags

};

    exit 1;
}

my $tag = $ARGV[0];

if (!($tag =~ m/^v[0-9]+\.[0-9]+\.[0-9]+/) ) {
    print "Tags for release should be in the format vX.Y.Z where X,Y,Z are integers\n";
    exit(1);
}

# Check that we are ready to update version

my $cleanTree = qq{On branch ${masterBranchLabel}
Your branch is up to date with 'origin/${masterBranchLabel}'.

nothing to commit, working tree clean
};

my $status = `git status -b`;

if (!($status eq $cleanTree)) {
    print "Must be up to date and on branch $masterBranchLabel to tag a release.\ngit status output:\n";
    # run it again so output is colorized
    system("git status -b");
    exit(1);
}

# Multiple tags get weird
chomp(my $existingTag = `git tag --points-at HEAD`);
if ($existingTag) {
    print "There is already a tag $existingTag on this commit. Exiting\n";
    exit(1);
}

# Also don't allow a duplicate tag. git will stop this later, but less messy to check here
my @allTags = `git tag`;

chomp(@allTags);

foreach my $repoTag (@allTags) {
    if ($repoTag eq ${tag}) {
        print "The tag $tag already exists. Exiting \n";
        exit(1);
    }
}

# Check tag matches Version.cmake
open(my $inFH, "<", "Version.cmake");
my $versionDotCmake = do { local $/; <$inFH> };
close($inFH);

my ($tagVersionMajor,$tagVersionMinor,$tagVersionPatch) = ($tag =~ m/^v([0-9]+)\.([0-9]+)\.([0-9]+)/);

$versionDotCmake =~ s/set\(\$\{PROJECT_NAME\}_VERSION_MAJOR "[0-9]+"\)/set\(\$\{PROJECT_NAME\}_VERSION_MAJOR "${tagVersionMajor}"\)/
    or die("Cannot find version information in Version.cmake");

$versionDotCmake =~ s/set\(\$\{PROJECT_NAME\}_VERSION_MINOR "[0-9]+"\)/set\(\$\{PROJECT_NAME\}_VERSION_MINOR "${tagVersionMinor}"\)/
    or die("Cannot find version information in Version.cmake");

$versionDotCmake =~ s/set\(\$\{PROJECT_NAME\}_VERSION_PATCH "[0-9]+"\)/set\(\$\{PROJECT_NAME\}_VERSION_PATCH "${tagVersionPatch}"\)/
    or die("Cannot find version information in Version.cmake");

$versionDotCmake =~ s/set\(\$\{PROJECT_NAME\}_RELEASE_VERSION 0\)/set\(\$\{PROJECT_NAME\}_RELEASE_VERSION 1\)/;

open(my $outFH, ">", "Version.cmake");
print $outFH $versionDotCmake;
close($outFH);

print "Modified version information in Version.cmake. Review git diff:\n";
system("git diff");
print "\nProceed with tag? Please enter (yY/nN):";
chomp(my $proceed = <STDIN>);

if (!(lc($proceed) eq "y")) {
    print "Reverting changes\n";
    system("git checkout Version.cmake");
    exit(1);
}

print("\nUpdating Version.cmake\n");
system("git add Version.cmake");
system("git commit -m \"Updating version for release $tag\"");
print("\nPushing changed Version.cmake\n");
system("git push origin $masterBranchLabel") == 0
    or die("Could not push updated Version.cmake");
print("\nApplying tag\n");
system("git tag -a $tag -m \"Tagging release $tag\"");
system("git push --tags");

print("\nUpdating Version.cmake for development");

$versionDotCmake =~ s/set\(\$\{PROJECT_NAME\}_RELEASE_VERSION 1\)/set\(\$\{PROJECT_NAME\}_RELEASE_VERSION 0\)/;

open(my $outFH, ">", "Version.cmake");
print $outFH $versionDotCmake;
close($outFH);

system("git add Version.cmake");
system("git commit -m \"[skip ci] Updating version for development post $tag\"");
print("\nPushing changed Version.cmake\n");
system("git push origin $masterBranchLabel") == 0 
    or die("Could not update Version.cmake post release");
