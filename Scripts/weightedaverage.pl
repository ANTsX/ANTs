#!/usr/bin/perl -w
    if ($#ARGV >= 0) { $who = join(' ', $ARGV[0]); }
    if ($#ARGV >= 1) { $imagedimension = join(' ', $ARGV[1]); }
    if ($#ARGV >= 2) { $outputname = join(' ', $ARGV[2]); }
 if ($#ARGV < 2) {
print " usage: \n   perl weightedaverage.pl  *files.jpg  imagedimension  weight1  ... weightN  \n ";
exit(1);
}

    opendir(DIR, $who);
    @filelist = glob("$who");
    $basect=3;
    $ct=$basect;

$dir = `pwd`;
$dir=substr($dir,0,length($dir)-1)."/";

    foreach $name (@filelist)
    {
    $weight=0;
     if ($#ARGV >= $ct) { $weight = join(' ', $ARGV[$ct]); }
     print " count ".$ct." & w=  ".$weight." \n ";
     $tempname=" temp.nii ";
     $exe = " ImageMath ".$imagedimension." ".$tempname." m ".$name." ".$weight." \n ";
     system $exe;
     $exe = " ImageMath ".$imagedimension."  ".$outputname." +  ".$tempname."  ".$outputname." \n ";
    if ($ct > $basect ) { system $exe; }
     else {
       $exe=" ImageMath ".$imagedimension." ".$outputname." m ".$name." ".$weight." \n ";
       print " YEE \n";
       system $exe;
         }
      $ct=$ct+1;
    }
    closedir(DIR);
