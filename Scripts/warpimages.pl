#!/usr/bin/perl -w
    if ($#ARGV >= 0) { $who = join(' ', $ARGV[0]); }
    if ($#ARGV >= 1) { $template = join(' ', $ARGV[1]); }
    if ($#ARGV >= 2) { $dim = join(' ', $ARGV[2]); } else { $dim=0; }
    if ($#ARGV >= 3) { $outnaming = join(' ', $ARGV[3]); } else { $outnaming=""; }
    opendir(DIR, $who);
    @filelist = glob("$who");
#    @filelist2 = glob("$who2");
    $ct=0;
print " EXPECTING IMAGE DIMENSION ".$dim." \n \n \n ";

$dir = `pwd`;
$dir=substr($dir,0,length($dir)-1)."/";
$pathpre="/mnt/aibs1/avants/bin/ants/";

$qsname="superfresh.sh";

  foreach $name (@filelist)
    {
      $pre=$outnaming.$name;
      $pre=~s/\.[^.]*$//;
	print $pre." \n ";
	$prog=$pathpre."WarpImageMultiTransform ".$dim;
	$exe=$prog."  ".$name." ".$pre."deformed.nii ".$pre."Warp.nii ".$pre."Affine.txt -R ".$template;
        $exe2 = "/mnt/pkg/sge-root/bin/lx24-x86/qsub -q mac ".$qsname."  ".$dir." ".$exe." ";
        print $exe."\n";
      system $exe2;
# USE qstat TO VERIFY THE JOB IS SUBMITTED
	$ct=$ct+1;
    }



$user=`whoami`;
chomp($user);
print " you are ".$user." = a jingle-jangle-jungle ";
$atcfn=0;
if ($atcfn)  { $waitcmd="vq -s | grep ".$user; }
else { $waitcmd=" qstat -u ".$user." | grep ".substr($qsname,0,4); }
      $continuewaiting=1;
if ($continuewaiting) { system "sleep 40";  }
      while($continuewaiting)
      {
	system $waitcmd;
	$output = `$waitcmd`;

	if ( $output eq "" )
	{
	    $continuewaiting=0;
	    print "Jobs all finished!\n";
	}
	else { print " wait ";system "sleep 30"; }
       }

    closedir(DIR);
