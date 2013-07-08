#!/usr/bin/perl -w
    if ($#ARGV >= 0) { $who = join(' ', $ARGV[0]); }

    opendir(DIR, $who);
    @filelist = glob("$who");
#    @filelist2 = glob("$who2");
    $ct=0;

$dir = `pwd`;
$dir=substr($dir,0,length($dir)-1)."/";

$pathpre="/mnt/aibs1/avants/bin/ants/";
    foreach $name (@filelist)
    {
        $exe = $pathpre."ANTS 3 -m PR[templatecontrol.nii,".$name.",1,2] -o ".$name." -t SyN[3] -r Gauss[2,0] -i 100x100x30 ";
        $exe2 = "/mnt/pkg/sge-root/bin/lx24-x86/qsub -q mac /mnt/data1/avants/Data/code/qsub3.sh  ".$dir." ".$exe." ";
        print "$exe\n";
      system $exe2;
# USE qstat TO VERIFY THE JOB IS SUBMITTED
        $ct=$ct+1;
    }
    closedir(DIR);
